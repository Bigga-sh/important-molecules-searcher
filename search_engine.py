"""
Importable search engine for the molecule GUI.
Thin wrapper around the core substructure-matching logic.
"""
import csv
import os
import requests
from rdkit import Chem, RDLogger

RDLogger.DisableLog('rdApp.*')

DATABASES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'databases')


def get_smiles_from_cas(cas: str) -> str | None:
    """Convert a CAS number to SMILES via CIR, then PubChem as fallback."""
    cas = cas.strip()
    try:
        r = requests.get(
            f'https://cactus.nci.nih.gov/chemical/structure/{cas}/smiles',
            timeout=10,
        )
        if r.status_code == 200:
            s = r.text.strip()
            if s and not s.startswith('<'):
                return s
    except Exception:
        pass
    try:
        r = requests.get(
            f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{cas}'
            f'/property/CanonicalSMILES/json',
            timeout=10,
        )
        if r.status_code == 200:
            return r.json()['PropertyTable']['Properties'][0].get('CanonicalSMILES')
    except Exception:
        pass
    return None


def _category_from_filename(filename: str) -> str:
    name = os.path.splitext(os.path.basename(filename))[0]
    name = name.replace('_data', '')
    name = name.replace('cosing_ingredients_with_smiles', 'cosing_ingredients')
    return name.replace('_', ' ').title()


def load_csv_database(csv_file: str) -> list[dict]:
    compounds = []
    category = _category_from_filename(csv_file)
    try:
        with open(csv_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                name   = (row.get('Name') or row.get('INCI_Name') or row.get('name') or '').strip()
                cas    = (row.get('CAS') or row.get('CAS No') or row.get('cas') or '').strip()
                smiles = (row.get('SMILES') or row.get('smiles') or '').strip()
                if smiles and smiles not in ('Not found', 'N/A', ''):
                    compounds.append({
                        'Name': name,
                        'CAS': cas,
                        'SMILES': smiles,
                        'Category': category,
                    })
    except Exception as exc:
        print(f"Warning: could not load {csv_file}: {exc}")
    return compounds


def load_all_databases(
    databases_folder: str = None,
) -> tuple[list[dict], dict[str, int]]:
    """
    Load all CSVs from the databases folder.
    Returns (all_compounds, {filename: compound_count}).
    """
    if databases_folder is None:
        databases_folder = DATABASES_DIR

    all_compounds: list[dict] = []
    stats: dict[str, int] = {}

    if not os.path.isdir(databases_folder):
        return all_compounds, stats

    for filename in sorted(os.listdir(databases_folder)):
        if not filename.endswith('.csv'):
            continue
        path = os.path.join(databases_folder, filename)
        compounds = load_csv_database(path)
        all_compounds.extend(compounds)
        stats[filename] = len(compounds)

    return all_compounds, stats


def search_substructure(
    query_smiles: str,
    all_compounds: list[dict],
) -> list[dict]:
    """
    Return all compounds that contain query_smiles as a substructure.
    Raises ValueError for an invalid query SMILES.
    """
    query_mol = Chem.MolFromSmiles(query_smiles)
    if query_mol is None:
        raise ValueError(f"Invalid SMILES: {query_smiles!r}")

    results = []
    for compound in all_compounds:
        try:
            mol = Chem.MolFromSmiles(compound['SMILES'], sanitize=False)
            if mol:
                Chem.SanitizeMol(mol)
                if mol.HasSubstructMatch(query_mol):
                    results.append(compound)
        except Exception:
            continue
    return results


def get_aromatic_position_presets(mol) -> list[dict]:
    """
    For each 6-membered aromatic ring with at least one external substituent,
    return a dict with ortho/meta/para atom-index sets relative to that substituent.

    Each entry: {anchor, anchor_label, ipso, ortho, meta, para, ring}
    where all position values are frozenset[int] of atom indices in mol.
    """
    ri = mol.GetRingInfo()
    presets = []

    for ring in ri.AtomRings():
        if len(ring) != 6:
            continue
        if not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            continue

        ring_set = set(ring)
        ring_list = list(ring)
        n = len(ring_list)

        for pos, atom_idx in enumerate(ring_list):
            atom = mol.GetAtomWithIdx(atom_idx)
            external = [nb for nb in atom.GetNeighbors() if nb.GetIdx() not in ring_set]
            if not external:
                continue

            ortho = frozenset({ring_list[(pos + 1) % n], ring_list[(pos - 1) % n]})
            meta  = frozenset({ring_list[(pos + 2) % n], ring_list[(pos - 2) % n]})
            para  = frozenset({ring_list[(pos + n // 2) % n]})

            ext_syms = '/'.join(nb.GetSymbol() for nb in external)
            presets.append({
                'ring':         tuple(ring),
                'anchor':       atom_idx,
                'anchor_label': f"atom {atom_idx} ({ext_syms})",
                'ipso':         frozenset({atom_idx}),
                'ortho':        ortho,
                'meta':         meta,
                'para':         para,
            })

    return presets


def search_with_attachment_constraints(
    query_smiles: str,
    all_compounds: list[dict],
    allowed_query_atoms: set[int] | None = None,
) -> list[dict]:
    """
    Substructure search with positional constraints.

    allowed_query_atoms: atom indices (in the query molecule) where extra
        substituents in the parent compound are permitted.  Any atom NOT in
        this set must match exactly — no additional neighbours allowed.
        If None, falls back to standard substructure search.

    Example: for aniline with allowed_query_atoms = {2, 6} (ortho carbons),
        only ortho-functionalised aniline derivatives are returned.
    """
    query_mol = Chem.MolFromSmiles(query_smiles)
    if query_mol is None:
        raise ValueError(f"Invalid SMILES: {query_smiles!r}")

    results = []
    for compound in all_compounds:
        try:
            mol = Chem.MolFromSmiles(compound['SMILES'], sanitize=False)
            if not mol:
                continue
            Chem.SanitizeMol(mol)
            matches = mol.GetSubstructMatches(query_mol)
            if not matches:
                continue

            if allowed_query_atoms is None:
                # No constraint — equivalent to basic substructure search.
                results.append(compound)
                continue

            # At least one match must satisfy the constraint.
            for match in matches:
                match_set = set(match)
                ok = True
                for query_idx, mol_idx in enumerate(match):
                    extra = [
                        nb for nb in mol.GetAtomWithIdx(mol_idx).GetNeighbors()
                        if nb.GetIdx() not in match_set
                    ]
                    if extra and query_idx not in allowed_query_atoms:
                        ok = False
                        break
                if ok:
                    results.append(compound)
                    break

        except Exception:
            continue

    return results

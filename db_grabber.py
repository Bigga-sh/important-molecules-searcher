"""
Database grabber — fetch CAS / SMILES / IUPAC from PubChem for a list of names.
Importable by the GUI; no stdin prompts.
"""
import csv
import os
import re
import time
import requests

DATABASES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'databases')

_CAS_PAT = re.compile(r'^\d{2,7}-\d{2}-\d$')


def clean_name(name: str) -> str:
    name = re.sub(r'^E\s*\d+[a-z]*\s*', '', name, flags=re.IGNORECASE)
    if ',' in name:
        name = name.split(',')[0]
    if '/' in name and not any(x in name.lower() for x in ['1,2-diol', "5′-"]):
        name = name.split('/')[0]
    name = re.sub(r'\s*\([^)]*\)', '', name)
    return name.strip()


def pubchem_lookup(name: str) -> tuple[str | None, str | None, str | None]:
    """Return (cas, smiles, iupac) or (None, None, None) on failure."""
    search = clean_name(name)
    if not search or search.lower() == 'quantum satis':
        return None, None, None

    try:
        r = requests.get(
            'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'
            f'{requests.utils.quote(search)}/property/CanonicalSMILES,IUPACName/json',
            timeout=15,
        )
        if r.status_code != 200:
            return None, None, None
        props = r.json()['PropertyTable']['Properties'][0]
        smiles = props.get('CanonicalSMILES', '')
        iupac  = props.get('IUPACName', '')

        r2 = requests.get(
            'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'
            f'{requests.utils.quote(search)}/synonyms/json',
            timeout=15,
        )
        cas = ''
        if r2.status_code == 200:
            syns = r2.json()['InformationList']['Information'][0].get('Synonym', [])
            hits = [s for s in syns if _CAS_PAT.match(s)]
            if hits:
                cas = hits[0]

        return cas or None, smiles or None, iupac or None

    except Exception:
        return None, None, None


def get_smiles_cir(identifier: str) -> str | None:
    try:
        r = requests.get(
            f'https://cactus.nci.nih.gov/chemical/structure/'
            f'{requests.utils.quote(identifier)}/smiles',
            timeout=10,
        )
        if r.status_code == 200:
            s = r.text.strip()
            if s and not s.startswith('<') and len(s) > 2:
                return s
    except Exception:
        pass
    return None


def get_smiles_opsin(iupac_name: str) -> str | None:
    try:
        r = requests.get(
            f'https://opsin.ch.cam.ac.uk/opsin/{requests.utils.quote(iupac_name)}.smi',
            timeout=10,
        )
        if r.status_code == 200:
            s = r.text.strip()
            if s and not s.startswith('<') and len(s) > 2:
                return s
    except Exception:
        pass
    return None


def parse_compound_names(text: str) -> list[str]:
    """Parse unique compound names from a block of text (one per line)."""
    seen: set[str] = set()
    names: list[str] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if any(line.startswith(x) for x in ['E-number', 'Name', '▼', 'PART', 'DEFINITIONS']):
            continue
        if 'quantum satis' in line.lower():
            continue
        if re.match(r'^\(\d+\)', line) or line.startswith('('):
            continue
        if line.endswith(':') or re.match(r'^\d+\.\s+', line):
            continue
        if re.match(r'^E\s*\d+$', line, re.IGNORECASE):
            continue
        if len(line) > 2 and line not in seen:
            seen.add(line)
            names.append(line)
    return names


def fetch_compounds(
    names: list[str],
    progress_callback=None,
) -> list[dict]:
    """
    Fetch CAS, SMILES, IUPAC for each name.
    progress_callback(i, total, name) is called after each lookup.
    """
    results = []
    total = len(names)

    for i, name in enumerate(names):
        cas, smiles, iupac = pubchem_lookup(name)

        # Fallbacks for missing SMILES
        if not smiles:
            if iupac:
                smiles = get_smiles_opsin(iupac)
            if not smiles:
                smiles = get_smiles_cir(name)
            if not smiles and cas:
                smiles = get_smiles_cir(cas)

        results.append({
            'Name': name,
            'CAS': cas or 'Not found',
            'SMILES': smiles or 'Not found',
            'IUPAC_Name': iupac or '',
        })

        if progress_callback:
            progress_callback(i + 1, total, name)

        time.sleep(0.3)  # rate-limit PubChem

    return results


def save_database(results: list[dict], output_path: str) -> tuple[int, int]:
    """
    Write results to a CSV file.
    Returns (total_rows, rows_with_smiles).
    """
    parent = os.path.dirname(output_path)
    if parent:
        os.makedirs(parent, exist_ok=True)

    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['Name', 'CAS', 'SMILES', 'IUPAC_Name'])
        writer.writeheader()
        writer.writerows(results)

    has_smiles = sum(1 for r in results if r['SMILES'] not in ('Not found', ''))
    return len(results), has_smiles

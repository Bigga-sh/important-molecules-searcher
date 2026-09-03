"""
Molecule Substructure Search — main Streamlit page.
Run with:  streamlit run app.py
"""
import os
import sys
from io import BytesIO

import pandas as pd
import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from search_engine import (
    get_smiles_from_cas,
    load_all_databases,
    search_with_attachment_constraints,
    get_aromatic_position_presets,
)

try:
    from streamlit_ketcher import st_ketcher
    _KETCHER = True
except ImportError:
    _KETCHER = False

# ── Page config ───────────────────────────────────────────────────────
st.set_page_config(
    page_title="Molecule DB Searcher",
    page_icon="⚗️",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Session state defaults ────────────────────────────────────────────
for key, default in [
    ("search_results", None),
    ("drawn_smiles", ""),
    ("last_query_smiles", ""),
    ("cas_fetched_smiles", ""),
    ("page_num", 0),
    ("page_size", 24),
]:
    if key not in st.session_state:
        st.session_state[key] = default

# ── Database loading (cached across reruns) ───────────────────────────
@st.cache_resource(show_spinner="Loading chemical databases…")
def _load_dbs():
    return load_all_databases()

all_compounds, db_stats = _load_dbs()


def _draw_with_atom_indices(mol, width: int = 450, height: int = 320) -> Image.Image:
    """Draw molecule with atom-index labels; tries Cairo then atom-map-number fallback."""
    try:
        d = rdMolDraw2D.MolDraw2DCairo(width, height)
        d.drawOptions().addAtomIndices = True
        d.DrawMolecule(mol)
        d.FinishDrawing()
        return Image.open(BytesIO(d.GetDrawingText()))
    except Exception:
        pass
    # Fallback: encode indices as atom map numbers (always renders)
    mol_copy = Chem.RWMol(mol)
    for a in mol_copy.GetAtoms():
        a.SetAtomMapNum(a.GetIdx())
    return Draw.MolToImage(mol_copy.GetMol(), size=(width, height))


# ── Sidebar ───────────────────────────────────────────────────────────
with st.sidebar:
    st.markdown("## ⚗️ Molecule DB Searcher")
    st.metric("Compounds loaded", f"{len(all_compounds):,}")
    st.caption(f"From {len(db_stats)} database files")
    if db_stats:
        with st.expander("Database breakdown"):
            for fname, count in db_stats.items():
                label = fname.replace('_data.csv', '').replace('_', ' ').title()
                st.markdown(f"- **{label}**: {count:,}")

# ── Title ─────────────────────────────────────────────────────────────
st.title("Substructure Search")
st.markdown(
    "Draw or enter a **query molecule** to find all database compounds "
    "that **contain it as a substructure** (i.e., parent compounds you could "
    "synthesise starting from your fragment)."
)

if not all_compounds:
    st.error(
        "No databases found. Go to **Database Manager** to build them first, "
        "or place CSV files in the `databases/` folder."
    )
    st.stop()

# ── Input tabs ────────────────────────────────────────────────────────
tab_draw, tab_smiles, tab_cas = st.tabs(
    ["✏️  Draw molecule", "  SMILES input", "  CAS lookup"]
)

query_smiles: str | None = None

with tab_draw:
    if _KETCHER:
        st.markdown("Draw your query fragment in the editor below.")
        drawn = st_ketcher(
            value=st.session_state.drawn_smiles,
            height=520,
            key="ketcher_main",
        )
        if drawn:
            st.session_state.drawn_smiles = drawn
        query_smiles = st.session_state.drawn_smiles or None
    else:
        st.warning(
            "Ketcher editor not available. "
            "Install it with:  `pip install streamlit-ketcher`\n\n"
            "Use the **SMILES input** tab in the meantime."
        )

with tab_smiles:
    raw = st.text_input(
        "SMILES string",
        placeholder="e.g.  NCc1ccccc1   (benzylamine)",
        key="smiles_input",
    )
    if raw:
        query_smiles = raw.strip()

with tab_cas:
    cas_input = st.text_input(
        "CAS number",
        placeholder="e.g.  100-46-9",
        key="cas_input",
    )
    fetch_col, _ = st.columns([1, 4])
    with fetch_col:
        if st.button("Fetch SMILES", key="fetch_cas_btn"):
            if cas_input:
                with st.spinner("Querying CIR / PubChem…"):
                    fetched = get_smiles_from_cas(cas_input.strip())
                st.session_state.cas_fetched_smiles = fetched or ""
    if st.session_state.cas_fetched_smiles:
        st.success(f"SMILES: `{st.session_state.cas_fetched_smiles}`")
        query_smiles = st.session_state.cas_fetched_smiles
    elif cas_input and not st.session_state.cas_fetched_smiles:
        pass
    elif cas_input:
        st.error("Could not resolve a SMILES for this CAS number.")

# ── Query preview, constraints and search ────────────────────────────
st.divider()

allowed_atoms_for_search: set[int] | None = None
search_clicked = False

if query_smiles:
    mol_preview = Chem.MolFromSmiles(query_smiles)
    if mol_preview is None:
        st.warning("⚠️ The SMILES is not valid — please check the structure.")
    else:
        col_img, col_info, col_btn = st.columns([1, 2, 1])
        with col_img:
            preview_img = Draw.MolToImage(mol_preview, size=(220, 160))
            st.image(preview_img, caption="Query fragment")
        with col_info:
            st.markdown(f"**SMILES:** `{query_smiles}`")
            st.markdown(f"**Heavy atoms:** {mol_preview.GetNumAtoms()}")
            st.markdown(f"**Database size:** {len(all_compounds):,} compounds")
        with col_btn:
            st.markdown("<br><br>", unsafe_allow_html=True)
            search_clicked = st.button(
                "🔍  Search", type="primary", use_container_width=True
            )

        # ── Position constraints expander ──────────────────────────────
        with st.expander("🔬 Position constraints (optional — click to expand)"):
            n_atoms = mol_preview.GetNumAtoms()

            def _atom_lbl(i: int) -> str:
                return f"{i}: {mol_preview.GetAtomWithIdx(i).GetSymbol()}"

            st.markdown(
                "Choose **which atoms of your query can bear extra substituents** "
                "in the parent compounds.  Atoms you do **not** select must appear "
                "exactly as drawn — no additional groups at those sites.\n\n"
                "_Example for aniline (Nc₁ccccc₁): selecting only the two ortho "
                "ring carbons restricts results to ortho-functionalised derivatives._"
            )

            col_idx, col_sel = st.columns([1, 1.4])

            with col_idx:
                idx_img = _draw_with_atom_indices(mol_preview)
                st.image(idx_img, caption="Atom indices")

            with col_sel:
                # Preset buttons come first so their session_state write is
                # visible to the multiselect rendered just below.
                presets = get_aromatic_position_presets(mol_preview)
                if presets:
                    st.markdown("**Quick ring presets:**")
                    for p in presets[:2]:   # cap at 2 rings to keep UI tidy
                        st.caption(f"Relative to {p['anchor_label']}")
                        pc1, pc2, pc3, pc4 = st.columns(4)
                        with pc1:
                            if st.button("Ortho", key=f"btn_o_{p['anchor']}"):
                                st.session_state["atom_pos_multiselect"] = sorted(p["ortho"])
                                st.rerun()
                        with pc2:
                            if st.button("Meta", key=f"btn_m_{p['anchor']}"):
                                st.session_state["atom_pos_multiselect"] = sorted(p["meta"])
                                st.rerun()
                        with pc3:
                            if st.button("Para", key=f"btn_p_{p['anchor']}"):
                                st.session_state["atom_pos_multiselect"] = sorted(p["para"])
                                st.rerun()
                        with pc4:
                            if st.button("O + P", key=f"btn_op_{p['anchor']}"):
                                st.session_state["atom_pos_multiselect"] = sorted(
                                    p["ortho"] | p["para"]
                                )
                                st.rerun()

                selected_atoms: list[int] = st.multiselect(
                    "Allowed substitution atoms:",
                    options=list(range(n_atoms)),
                    format_func=_atom_lbl,
                    key="atom_pos_multiselect",
                    help=(
                        "Select atom indices where extra groups are permitted in "
                        "the parent compound. Leave empty for standard search."
                    ),
                )
                allowed_atoms_for_search = set(selected_atoms) if selected_atoms else None

                if st.button("Clear constraint", key="btn_clear_constraint"):
                    st.session_state["atom_pos_multiselect"] = []
                    st.rerun()

            if allowed_atoms_for_search:
                labels = ", ".join(_atom_lbl(i) for i in sorted(allowed_atoms_for_search))
                st.info(f"Constraint **active** — extra substituents allowed only at: {labels}")
            else:
                st.caption("No constraint set — standard substructure search.")

        # ── Execute search ─────────────────────────────────────────────
        if search_clicked:
            with st.spinner(f"Searching {len(all_compounds):,} compounds…"):
                try:
                    hits = search_with_attachment_constraints(
                        query_smiles, all_compounds, allowed_atoms_for_search
                    )
                    st.session_state.search_results = hits
                    st.session_state.last_query_smiles = query_smiles
                    st.session_state.page_num = 0   # reset to page 1 on new search
                except ValueError as exc:
                    st.error(str(exc))
                    st.session_state.search_results = None

# ── Results ───────────────────────────────────────────────────────────
if st.session_state.search_results is not None:
    results = st.session_state.search_results
    st.divider()

    if not results:
        st.info("No compounds in the database contain this substructure.")
    else:
        st.success(f"Found **{len(results):,}** matching compounds")
        df = pd.DataFrame(results)

        # Category breakdown
        with st.expander("Results by category", expanded=True):
            cat_df = (
                df["Category"]
                .value_counts()
                .rename_axis("Category")
                .reset_index(name="Count")
            )
            st.bar_chart(cat_df.set_index("Category"))

        # Full data table
        st.subheader("All matches")
        st.dataframe(
            df[["Name", "CAS", "Category", "SMILES"]],
            use_container_width=True,
            height=420,
        )
        st.download_button(
            "⬇️  Download CSV",
            data=df.to_csv(index=False).encode("utf-8"),
            file_name="substructure_results.csv",
            mime="text/csv",
        )

        # ── Paginated structure grid ───────────────────────────────────
        PAGE_SIZES = [24, 48, 96]
        total      = len(results)
        page_size  = st.session_state.page_size
        page_num   = st.session_state.page_num

        total_pages = max(1, (total + page_size - 1) // page_size)
        page_num    = min(page_num, total_pages - 1)   # clamp after result-set shrinks
        start_idx   = page_num * page_size
        end_idx     = min(start_idx + page_size, total)

        st.subheader(f"Structure viewer — {start_idx + 1}–{end_idx} of {total:,}")

        # Pagination controls row
        c_ps, c_prev, c_info, c_next = st.columns([1.2, 1, 2.5, 1])
        with c_ps:
            new_ps = st.selectbox(
                "Per page",
                PAGE_SIZES,
                index=PAGE_SIZES.index(page_size) if page_size in PAGE_SIZES else 0,
                key="page_size_select",
                label_visibility="collapsed",
            )
            if new_ps != page_size:
                st.session_state.page_size = new_ps
                st.session_state.page_num  = 0
                st.rerun()
        with c_prev:
            if st.button("◀ Prev", disabled=(page_num == 0), key="prev_page_btn"):
                st.session_state.page_num = page_num - 1
                st.rerun()
        with c_info:
            st.markdown(
                f"<div style='text-align:center;padding-top:6px'>"
                f"Page <b>{page_num + 1}</b> / {total_pages}</div>",
                unsafe_allow_html=True,
            )
        with c_next:
            if st.button("Next ▶", disabled=(page_num >= total_pages - 1), key="next_page_btn"):
                st.session_state.page_num = page_num + 1
                st.rerun()

        # Build and render the grid for the current page
        page_results = results[start_idx:end_idx]
        mols, legends = [], []
        for r in page_results:
            try:
                m = Chem.MolFromSmiles(r["SMILES"], sanitize=False)
                if m:
                    Chem.SanitizeMol(m)
                    mols.append(m)
                    name = r["Name"]
                    if len(name) > 35:
                        name = name[:33] + "…"
                    legends.append(f"{name}\n{r['CAS']}\n[{r['Category']}]")
            except Exception:
                continue

        if mols:
            per_row = min(4, len(mols))
            with st.spinner("Rendering structures…"):
                grid_img = Draw.MolsToGridImage(
                    mols,
                    molsPerRow=per_row,
                    subImgSize=(350, 300),
                    legends=legends,
                )
            st.image(grid_img, use_container_width=True)
        else:
            st.warning("No valid structures to display on this page.")

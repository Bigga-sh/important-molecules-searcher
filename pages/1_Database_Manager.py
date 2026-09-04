"""
Database Manager — Streamlit page.
Add new compound lists and fetch their data from PubChem.
"""
import os
import sys

import pandas as pd
import streamlit as st

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from db_grabber import (
    DATABASES_DIR,
    fetch_compounds,
    parse_compound_names,
    save_database,
)
from search_engine import load_all_databases

st.set_page_config(
    page_title="Database Manager",
    page_icon="🗄️",
    layout="wide",
)

st.title("🗄️ Database Manager")
st.markdown(
    "Manage the compound databases used by the search engine. "
    "Upload a `.txt` file with compound names to fetch their data from PubChem "
    "and add them to the searchable database."
)

# ── Current database overview ─────────────────────────────────────────
st.header("Current Databases")

_, db_stats = load_all_databases()

if db_stats:
    df_stats = pd.DataFrame(
        [{"File": k, "Compounds (with SMILES)": v} for k, v in db_stats.items()]
    )
    st.dataframe(df_stats, use_container_width=True, hide_index=True)
    st.metric("Total searchable compounds", f"{sum(db_stats.values()):,}")
else:
    st.warning(
        f"No CSV databases found in `{DATABASES_DIR}`. "
        "Upload a compound list below to create the first one."
    )

st.divider()

# ── Upload new compound list ──────────────────────────────────────────
st.header("Add / Refresh a Database")

upload_tab, existing_tab = st.tabs(
    ["📤  Upload a new .txt file", "🔄  Refresh an existing list"]
)

def _run_fetch(names: list[str], out_name: str, key_suffix: str = ""):
    """Shared fetch-and-save logic used by both tabs."""
    estimate_min = len(names) * 0.6 / 60
    st.caption(
        f"Estimated time: ~{estimate_min:.0f}–{estimate_min * 2:.0f} min "
        f"for {len(names):,} compounds"
    )

    if not st.button("🚀  Fetch from PubChem and save", type="primary", key=f"fetch_{out_name}{key_suffix}"):
        return

    out_path = os.path.join(DATABASES_DIR, out_name + ".csv")
    progress_bar = st.progress(0.0)
    status_text = st.empty()

    def on_progress(i, total, name):
        progress_bar.progress(i / total)
        status_text.markdown(f"`[{i}/{total}]` **{name}**")

    try:
        results = fetch_compounds(names, progress_callback=on_progress)
    except Exception as exc:
        st.error(f"Error during fetch: {exc}")
        return

    total_rows, with_smiles = save_database(results, out_path)
    progress_bar.progress(1.0)
    status_text.empty()

    st.success(
        f"✅ Saved `{out_name}.csv` — "
        f"{total_rows:,} compounds, {with_smiles:,} with SMILES"
    )

    col1, col2, col3 = st.columns(3)
    col1.metric("Total saved", f"{total_rows:,}")
    col2.metric("With SMILES", f"{with_smiles:,}")
    col3.metric("Missing SMILES", f"{total_rows - with_smiles:,}")

    df_result = pd.DataFrame(results)
    no_smiles = df_result[df_result["SMILES"] == "Not found"]
    if not no_smiles.empty:
        with st.expander(f"Compounds with missing SMILES ({len(no_smiles):,})"):
            st.dataframe(
                no_smiles[["Name", "CAS", "IUPAC_Name"]],
                use_container_width=True,
            )

    csv_bytes = df_result.to_csv(index=False).encode("utf-8")
    st.download_button(
        "⬇️  Download results CSV",
        data=csv_bytes,
        file_name=out_name + ".csv",
        mime="text/csv",
    )

    # Bust the cache so the search page sees the new data on next load
    st.cache_resource.clear()
    st.info("🔄 Database cache cleared. Reload the **Search** page to use the updated data.")


with upload_tab:
    uploaded = st.file_uploader(
        "Upload compound list (.txt)",
        type=["txt"],
        help="Plain text with one compound name per line. "
             "E-numbers and parenthetical notes are cleaned automatically.",
    )
    if uploaded:
        raw_text = uploaded.read().decode("utf-8", errors="replace")
        names = parse_compound_names(raw_text)

        st.info(f"Parsed **{len(names):,}** unique names from `{uploaded.name}`")
        with st.expander("Preview parsed names"):
            preview = names[:60]
            st.write(preview)
            if len(names) > 60:
                st.caption(f"… and {len(names) - 60} more")

        default_out = os.path.splitext(uploaded.name)[0] + "_data"
        out_name = st.text_input(
            "Output filename (without .csv)",
            value=default_out,
            key="upload_outname",
        )
        if os.path.exists(os.path.join(DATABASES_DIR, out_name + ".csv")):
            st.warning(f"⚠️  `{out_name}.csv` already exists — it will be overwritten.")

        if names and out_name:
            _run_fetch(names, out_name)


with existing_tab:
    txt_files = []
    if os.path.isdir(DATABASES_DIR):
        txt_files = sorted(f for f in os.listdir(DATABASES_DIR) if f.endswith('.txt'))

    if not txt_files:
        st.info("No `.txt` files found in the databases folder.")
    else:
        selected = st.selectbox("Select a compound list to refresh", txt_files)
        if selected:
            txt_path = os.path.join(DATABASES_DIR, selected)
            with open(txt_path, 'r', encoding='utf-8', errors='replace') as fh:
                raw = fh.read()
            names = parse_compound_names(raw)
            st.info(f"Found **{len(names):,}** compounds in `{selected}`")

            default_out = os.path.splitext(selected)[0] + "_data"
            out_name = st.text_input(
                "Output filename (without .csv)",
                value=default_out,
                key="refresh_outname",
            )
            if names and out_name:
                _run_fetch(names, out_name)

# ── Active Pharmaceutical Ingredients (API) database ─────────────────
st.divider()
st.header("💊 Active Pharmaceutical Ingredients (API) Database")

api_out_name = "active_pharmaceutical_ingredients_data"
api_out_path = os.path.join(DATABASES_DIR, api_out_name + ".csv")
api_txt_path = os.path.join(DATABASES_DIR, "active_pharmaceutical_ingredients.txt")

if os.path.exists(api_out_path):
    with open(api_out_path, "r", encoding="utf-8") as _f:
        existing_api_count = sum(1 for _ in _f) - 1
    st.info(f"Local API database: **{existing_api_count:,} compounds with SMILES**.")

st.markdown(
    "**Primary source — bundled list (recommended):** ~300 WHO Essential Medicines "
    "and common APIs. Fetches SMILES from PubChem (reliable, ~15 min)."
)

if st.button("🏥  Build API database from bundled WHO/Essential list", type="primary"):
    if not os.path.exists(api_txt_path):
        st.error(f"Bundled list not found at `{api_txt_path}`")
    else:
        with open(api_txt_path, "r", encoding="utf-8") as fh:
            raw = fh.read()
        names = parse_compound_names(raw)
        st.info(f"Parsed **{len(names)}** drug names from the bundled list.")
        _run_fetch(names, api_out_name, key_suffix="_who")

st.markdown("---")
st.markdown(
    "**Optional — ChEMBL live download (~4,200 approved drugs):** "
    "Requires ChEMBL's API to be available. Checks status before attempting."
)

if st.button("🌐  Try ChEMBL download (requires server to be up)"):
    import requests as _req
    import csv as _csv2
    import time as _time2

    CHEMBL_HOST = "https://www.ebi.ac.uk"
    CHEMBL_BASE = f"{CHEMBL_HOST}/chembl/api/data"
    HEADERS     = {"User-Agent": "MoleculeDBSearcher/1.0", "Accept": "application/json"}
    PAGE_SIZE   = 200

    session = _req.Session()
    session.headers.update(HEADERS)

    # Status check first
    status_text = st.empty()
    status_text.markdown("Checking ChEMBL availability…")
    try:
        probe = session.get(
            f"{CHEMBL_BASE}/molecule",
            params={"max_phase": 4, "format": "json", "limit": 1},
            timeout=20,
        )
        if probe.status_code != 200 or not probe.text.strip():
            st.error(
                f"ChEMBL is currently unavailable (HTTP {probe.status_code}). "
                "Use the bundled list above instead."
            )
            st.stop()
        total = probe.json()["page_meta"]["total_count"]
    except Exception as exc:
        st.error(f"ChEMBL unreachable: {exc}. Use the bundled list above instead.")
        st.stop()

    status_text.markdown(f"ChEMBL online — downloading **{total:,}** approved drugs…")
    progress_bar = st.progress(0.0)
    all_rows: list[dict] = []

    # Start from page 1 with the full page size (the probe used limit=1,
    # so following its 'next' URL would download one molecule at a time)
    next_url = f"{CHEMBL_BASE}/molecule"
    next_params = {"max_phase": 4, "format": "json", "limit": PAGE_SIZE, "offset": 0}

    while next_url:
        try:
            r = session.get(next_url, params=next_params, timeout=60)
            next_params = None  # only send params on the first request; use next URL after
            if r.status_code != 200:
                st.warning(f"ChEMBL returned {r.status_code} mid-download; saving partial results.")
                break
            data = r.json()
        except Exception as exc:
            st.warning(f"Request failed mid-download: {exc}. Saving partial results.")
            break

        for mol in data.get("molecules", []):
            structs = mol.get("molecule_structures") or {}
            smi = structs.get("canonical_smiles", "")
            if not smi:
                continue
            props = mol.get("molecule_properties") or {}
            all_rows.append({
                "Name":       mol.get("pref_name") or mol.get("molecule_chembl_id", ""),
                "CAS":        "",
                "SMILES":     smi,
                "IUPAC_Name": "",
                "ChEMBL_ID":  mol.get("molecule_chembl_id", ""),
                "MolFormula": props.get("full_molformula", ""),
                "MW":         str(props.get("full_mwt", "")),
            })

        raw_next = data.get("page_meta", {}).get("next")
        if not raw_next or not data.get("molecules"):
            break
        next_url = CHEMBL_HOST + raw_next if raw_next.startswith("/") else raw_next

        progress_bar.progress(min(len(all_rows) / total, 1.0))
        status_text.markdown(f"Downloaded **{len(all_rows):,}** / {total:,}…")
        _time2.sleep(1.0)

    if all_rows:
        with open(api_out_path, "w", newline="", encoding="utf-8") as fh:
            writer = _csv2.DictWriter(fh, fieldnames=list(all_rows[0].keys()))
            writer.writeheader()
            writer.writerows(all_rows)
        progress_bar.progress(1.0)
        status_text.empty()
        st.success(f"✅ Saved **{len(all_rows):,}** APIs from ChEMBL to `{api_out_name}.csv`")
        st.cache_resource.clear()
        st.info("Reload the Search page to include the new database.")
    else:
        st.error("No compounds retrieved.")

# ── Manual SMILES lookup helper ───────────────────────────────────────
st.divider()
st.header("Single Compound Lookup")
st.markdown("Quickly look up one compound by name to preview its data from PubChem.")

lookup_name = st.text_input("Compound name", placeholder="e.g. benzyl alcohol")
if st.button("Look up", key="single_lookup"):
    if lookup_name:
        from db_grabber import pubchem_lookup, get_smiles_cir, get_smiles_opsin
        with st.spinner("Querying PubChem…"):
            cas, smiles, iupac = pubchem_lookup(lookup_name)
        if smiles:
            st.success("Found!")
        else:
            with st.spinner("Trying CIR / OPSIN fallbacks…"):
                smiles = get_smiles_cir(lookup_name) or (
                    get_smiles_opsin(iupac) if iupac else None
                )
            if smiles:
                st.warning("Found via fallback (CIR/OPSIN)")
            else:
                st.error("Not found in any source.")

        col_a, col_b = st.columns(2)
        col_a.markdown(f"**Name:** {lookup_name}")
        col_a.markdown(f"**CAS:** {cas or 'Not found'}")
        col_a.markdown(f"**SMILES:** `{smiles or 'Not found'}`")
        col_a.markdown(f"**IUPAC:** {iupac or '—'}")

        if smiles:
            from rdkit import Chem
            from rdkit.Chem import Draw
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                with col_b:
                    st.image(Draw.MolToImage(mol, size=(280, 200)), caption=lookup_name)

# ── COSING Update ─────────────────────────────────────────────────────
st.divider()
st.header("🇪🇺 COSING Ingredients Update")
st.markdown(
    "The **COSING** database (EU Commission cosmetic ingredients) requires a manual "
    "download because the EU Commission does not expose a stable machine-readable URL. "
    "Steps:\n"
    "1. Go to the [COSING portal](https://ec.europa.eu/growth/tools-databases/cosing/) "
    "→ *Download* → save the CSV file.\n"
    "2. Upload it below — the app will process it and refresh `cosing_ingredients_with_smiles.csv`."
)

cosing_upload = st.file_uploader(
    "Upload raw COSING CSV",
    type=["csv"],
    key="cosing_upload",
    help="The official COSING export file (COSING_Ingredients.csv or similar).",
)

if cosing_upload:
    import csv as _csv
    import io

    raw_bytes = cosing_upload.read()
    raw_text  = raw_bytes.decode("utf-8", errors="replace")

    # Detect how many header rows to skip (COSING files have 9 meta-rows before the CSV header)
    lines = raw_text.splitlines()
    header_skip = 0
    for i, line in enumerate(lines):
        if "INCI name" in line or "CAS No" in line:
            header_skip = i
            break

    st.info(f"Detected data starts at line {header_skip + 1}. File has {len(lines):,} rows.")

    # Save raw file to databases folder
    raw_out = os.path.join(DATABASES_DIR, "COSING_Ingredients.csv")
    with open(raw_out, "wb") as fh:
        fh.write(raw_bytes)
    st.caption(f"Saved raw file to `{raw_out}`")

    # Parse entries
    sys.path.insert(0, os.path.join(ROOT, "databases"))
    from cosing_search import read_cosing_csv, get_smiles_from_pubchem, get_smiles_from_cir, get_smiles_from_opsin

    import time as _time

    entries = read_cosing_csv(raw_out)
    st.info(f"Parsed **{len(entries):,}** entries with valid CAS numbers.")

    # Only process entries not already in existing output (incremental update)
    existing_out = os.path.join(DATABASES_DIR, "cosing_ingredients_with_smiles.csv")
    already_done: set[str] = set()
    existing_rows: list[dict] = []

    if os.path.exists(existing_out):
        with open(existing_out, "r", encoding="utf-8") as fh:
            for row in _csv.DictReader(fh):
                cas = row.get("CAS", "").strip()
                if cas:
                    already_done.add(cas)
                existing_rows.append(row)
        st.caption(f"Existing database has {len(existing_rows):,} entries ({len(already_done):,} unique CAS). Will only process new ones.")

    new_entries = [e for e in entries if e["CAS"] not in already_done]
    st.info(f"**{len(new_entries):,}** new entries to process (skipping {len(entries) - len(new_entries):,} already known).")

    if new_entries and st.button("🚀 Process & update COSING database", type="primary"):
        progress_bar = st.progress(0.0)
        status_text  = st.empty()
        new_results: list[dict] = []

        for i, entry in enumerate(new_entries):
            inci   = entry["INCI_Name"]
            cas    = entry["CAS"]
            iupac  = entry.get("IUPAC_Name", "")

            smiles = get_smiles_from_pubchem(inci, cas)
            if not smiles and iupac:
                smiles = get_smiles_from_opsin(iupac)
            if not smiles:
                smiles = get_smiles_from_cir(inci)
            if not smiles and cas:
                smiles = get_smiles_from_cir(cas)

            new_results.append({
                "INCI_Name":  inci,
                "CAS":        cas,
                "IUPAC_Name": iupac,
                "SMILES":     smiles or "Not found",
            })

            progress_bar.progress((i + 1) / len(new_entries))
            status_text.markdown(f"`[{i+1}/{len(new_entries)}]` **{inci[:60]}**")
            _time.sleep(0.3)

        # Merge with existing rows and write
        all_rows = existing_rows + new_results
        with open(existing_out, "w", newline="", encoding="utf-8") as fh:
            writer = _csv.DictWriter(
                fh, fieldnames=["INCI_Name", "CAS", "IUPAC_Name", "SMILES"]
            )
            writer.writeheader()
            writer.writerows(all_rows)

        progress_bar.progress(1.0)
        status_text.empty()

        found  = sum(1 for r in new_results if r["SMILES"] != "Not found")
        st.success(
            f"✅ Updated COSING database — "
            f"added {len(new_results):,} new entries ({found:,} with SMILES). "
            f"Total: {len(all_rows):,} entries."
        )
        st.cache_resource.clear()
        st.info("Reload the Search page to use the updated COSING data.")

    elif not new_entries:
        st.success("✅ No new entries — the COSING database is already up to date.")

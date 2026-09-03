# Molecule Substructure Database Searcher

A Streamlit web application for searching curated chemical databases by molecular substructure.  
Draw (or enter) a query fragment and instantly find all parent compounds across multiple industry-relevant databases.

---

## Features

- **Interactive drawing** via the Ketcher editor, or enter SMILES / CAS number directly
- **Substructure search** across 7 000+ compounds from multiple databases
- **Position constraints** — restrict results to parent compounds functionalised only at specific sites (e.g. ortho-only aniline derivatives, or amine-modified compounds with intact ring)
- **Paginated structure viewer** (24 / 48 / 96 per page, full navigation)
- **Database Manager** page: build new databases from any compound name list via PubChem, or update COSING incrementally
- **Download results** as CSV

---

## Included databases

| Database | Source | Notes |
|---|---|---|
| Herbicides | Manual curated list + PubChem | ~500 compounds |
| Fungicides | Manual curated list + PubChem | ~300 compounds |
| Insecticides | Manual curated list + PubChem | ~300 compounds |
| Synergists & Safeners | Manual curated list + PubChem | ~40 compounds |
| Wood Preservatives | Manual curated list + PubChem | ~50 compounds |
| Food Additives | Manual curated list + PubChem | ~900 compounds |
| COSING Ingredients | EU Commission COSING database | ~5 000 cosmetic ingredients |
| Active Pharmaceutical Ingredients | WHO Essential Medicines 23rd Ed. | ~300 drugs |

---

## Requirements

- Python 3.10 or newer
- Internet connection (for PubChem lookups in the Database Manager; search itself is fully offline)

---

## Installation

```bash
pip install -r requirements.txt
```

> **Note:** RDKit may take a minute to install. If `pip install rdkit` fails, try `pip install rdkit-pypi`.

---

## Running the app

**Windows** — double-click `Launch App.bat`

**Any platform** — run from the terminal:

```bash
streamlit run app.py
```

The app opens automatically in your browser at `http://localhost:8501`.

---

## Usage

### Basic substructure search

1. Go to the **Draw molecule** tab and sketch your fragment in the Ketcher editor  
   (or use **SMILES input** / **CAS lookup** tabs)
2. A preview of your fragment and its atom count appears below
3. Click **Search**
4. Results show a category breakdown chart, a full data table (downloadable as CSV), and a paginated structure grid

### Position-constrained search

Use this to filter *where* on the parent compound the query fragment may differ.

1. Draw / enter your query molecule
2. Expand **Position constraints** below the preview
3. The left panel shows your molecule with **atom index labels**
4. For aromatic rings, use the **Ortho / Meta / Para / O+P** quick-preset buttons  
   — these auto-populate the selection based on ring topology
5. Or select any combination of atoms manually in the multiselect
6. Click **Search**

**Example:** search for aniline (`Nc1ccccc1`) with **Ortho** selected  
→ returns only ortho-substituted aniline derivatives; meta/para positions must remain as drawn.

**Example:** select only atom 0 (N) for aniline  
→ returns compounds where the amine is further functionalised (NHR, NR₂…) but the benzene ring is unchanged.

### Adding new databases

1. Navigate to **Database Manager** in the left sidebar
2. Upload a `.txt` file with one compound name per line  
   (or select an existing name list from the `databases/` folder)
3. Click **Fetch from PubChem** — the app retrieves CAS and SMILES for each compound
4. The new database CSV is saved and becomes immediately searchable

---

## Project structure

```
database_searcher/
├── app.py                       # Main Streamlit page (search UI)
├── search_engine.py             # Substructure search logic (importable)
├── db_grabber.py                # PubChem fetcher for database building
├── pages/
│   └── 1_Database_Manager.py   # Database build / update UI
├── databases/
│   ├── *.txt                    # Compound name source lists
│   ├── *_data.csv               # Built compound databases (Name, CAS, SMILES)
│   └── COSING_Ingredients.csv   # Raw COSING source (for incremental updates)
├── requirements.txt
├── Launch App.bat               # Windows one-click launcher
└── README.md
```

---

## Troubleshooting

| Problem | Fix |
|---|---|
| `ModuleNotFoundError: rdkit` | Run `pip install rdkit` (or `rdkit-pypi`) |
| `streamlit-ketcher` not found | Run `pip install streamlit-ketcher` |
| Ketcher editor shows warning | Use the **SMILES input** tab instead |
| No databases found on startup | Run Database Manager to build them, or ensure `databases/*_data.csv` files are present |
| PubChem lookup times out | PubChem enforces rate limits; the fetcher waits 300 ms between requests automatically |

# FUS Gene Expression: Mechanosensitive Ion Channels in the Human Brain

Analysis of mechanosensitive ion channel gene expression across human brain regions using data from the [Allen Human Brain Atlas](https://human.brain-map.org/). This project supports research into focused ultrasound (FUS) neuromodulation by characterizing the spatial expression patterns of ion channels thought to mediate ultrasound's effects on neural tissue.

## Background

Focused ultrasound is hypothesized to exert neuromodulatory effects via mechanical activation of ion channels. Key channel families include:

- **Piezo** (PIEZO1/FAM38A, PIEZO2/FAM38B) — primary mechanotransducers
- **K2P / Two-pore potassium** (TREK-1/KCNK2, TREK-2/KCNK10, TRAAK/KCNK4) — mechanically gated leak channels
- **TRP** (TRPV1-4, TRPA1, TRPC1/5/6, TRPM4/7) — polymodal sensory channels
- **TRPP / Polycystin** (PKD2/TRPP2) — flow/pressure-sensitive channels

See [Bhatt et al. 2024 (PMC10423872)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10423872/) for a review.

## Analyses

1. **Regional expression profiling** — Expression levels of mechanosensitive channel genes across brain structures
2. **Hierarchical clustering** — Grouping brain regions by similarity of mechanosensitive channel expression
3. **Cell-type enrichment** — Association of channel genes with specific cell types (PV, SST, VIP, pyramidal, astrocytes)

## Data

Uses the Allen Human Brain Atlas normalized microarray (6 donors) and RNA-seq (2 donors) datasets. Data files are not included in this repo due to size — download from [Allen Brain Atlas API](http://human.brain-map.org/static/download).

Place downloaded zip files in the `data/` directory.

## Setup

```bash
conda activate py311env
pip install -r requirements.txt
```

## Project Structure

```
├── src/
│   ├── gene_lists.py          # Mechanosensitive channel gene definitions
│   ├── data_loader.py         # Allen Brain Atlas data loading utilities
│   ├── clustering.py          # Hierarchical clustering analysis
│   └── cell_type_markers.py   # Cell-type marker gene definitions and enrichment
├── notebooks/
│   ├── 01_data_exploration.ipynb
│   ├── 02_regional_expression.ipynb
│   ├── 03_hierarchical_clustering.ipynb
│   └── 04_cell_type_enrichment.ipynb
├── figures/
├── data/                      # Place Allen atlas zip files here
└── requirements.txt
```

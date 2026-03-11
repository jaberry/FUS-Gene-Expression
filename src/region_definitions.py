"""
Brain region definitions for anxiety/depression circuit analysis.

Maps Allen Human Brain Atlas structure acronyms to higher-level regions
relevant to mood/anxiety circuits and FUS targeting.

Laterality is collapsed — left and right samples for the same structure
share the same acronym in SampleAnnot.csv, so no special handling needed.
"""
from __future__ import annotations

# ---------------------------------------------------------------------------
# Subregion → parent region mapping
# Each key is the Allen acronym (as it appears in SampleAnnot structure_acronym).
# The value is the parent region label for grouping.
# ---------------------------------------------------------------------------

SUBREGION_TO_PARENT: dict[str, str] = {
    # Amygdala
    "BLA":   "Amygdala",
    "BMA":   "Amygdala",
    "CeA":   "Amygdala",
    "ATZ":   "Amygdala",   # amygdalohippocampal transition zone

    # BNST
    "BST":   "BNST",

    # Hippocampus
    "CA1":   "Hippocampus",
    "CA2":   "Hippocampus",
    "CA3":   "Hippocampus",
    "CA4":   "Hippocampus",
    "DG":    "Hippocampus",
    "S":     "Hippocampus",    # subiculum
    "PHG-cos": "Hippocampus",  # parahippocampal gyrus (bank of cos)
    "PHG-l":   "Hippocampus",  # parahippocampal gyrus (lateral)

    # Striatum (dorsal)
    "HCd":   "Striatum",      # head of caudate
    "BCd":   "Striatum",      # body of caudate
    "TCd":   "Striatum",      # tail of caudate
    "Pu":    "Striatum",      # putamen

    # Nucleus Accumbens
    "Acb":   "Nucleus Accumbens",

    # DLPFC (approximated via middle frontal gyrus)
    "MFG-s": "DLPFC",         # middle frontal gyrus, superior bank
    "MFG-i": "DLPFC",         # middle frontal gyrus, inferior bank
    "SFG-l": "DLPFC",         # superior frontal gyrus, lateral bank

    # Cingulate cortex
    "CgGf-i": "Cingulate Cortex",   # cingulate, frontal, inferior
    "CgGf-s": "Cingulate Cortex",   # cingulate, frontal, superior
    "CgGp-i": "Cingulate Cortex",   # cingulate, parietal, inferior
    "CgGp-s": "Cingulate Cortex",   # cingulate, parietal, superior
    "CgGr-i": "Cingulate Cortex",   # cingulate, retrosplenial, inferior
    "CgGr-s": "Cingulate Cortex",   # cingulate, retrosplenial, superior
    "SCG":    "Cingulate Cortex",    # subcallosal cingulate (BA25)

    # Thalamus
    "DTA":   "Thalamus",      # anterior group
    "DTM":   "Thalamus",      # medial group (includes MD)
    "DTLd":  "Thalamus",      # lateral group, dorsal
    "DTLv":  "Thalamus",      # lateral group, ventral
    "DTP":   "Thalamus",      # posterior group
    "ILc":   "Thalamus",      # intralaminar nuclei, caudal
    "ILr":   "Thalamus",      # intralaminar nuclei, rostral
    "Pa":    "Thalamus",      # paraventricular nuclei
    "R":     "Thalamus",      # reticular nucleus
    "LGd":   "Thalamus",      # dorsal lateral geniculate
    "MG":    "Thalamus",      # medial geniculate

    # Habenula
    "Hl":    "Habenula",      # lateral habenular nucleus
    "Hm":    "Habenula",      # medial habenular nucleus

    # Insula
    "SIG":   "Insula",        # short insular gyri
    "LIG":   "Insula",        # long insular gyri

    # Hypothalamus
    "AHA":   "Hypothalamus",  # anterior hypothalamic area
    "ARH":   "Hypothalamus",  # arcuate nucleus
    "DMH":   "Hypothalamus",  # dorsomedial hypothalamic nucleus
    "LHA":   "Hypothalamus",  # lateral hypothalamic area, anterior
    "LHM":   "Hypothalamus",  # lateral hypothalamic area, mammillary
    "LHT":   "Hypothalamus",  # lateral hypothalamic area, tuberal
    "PHA":   "Hypothalamus",  # posterior hypothalamic area
    "PVH":   "Hypothalamus",  # paraventricular nucleus
    "VMH":   "Hypothalamus",  # ventromedial hypothalamic nucleus
}

# All subregion acronyms
ALL_SUBREGIONS = list(SUBREGION_TO_PARENT.keys())

# Parent regions in display order
PARENT_REGIONS = [
    "Amygdala",
    "BNST",
    "Hippocampus",
    "Striatum",
    "Nucleus Accumbens",
    "DLPFC",
    "Cingulate Cortex",
    "Thalamus",
    "Habenula",
    "Insula",
    "Hypothalamus",
]

# Human-readable labels for subregions (for plot annotations)
SUBREGION_LABELS: dict[str, str] = {
    "BLA": "Basolateral Amy",
    "BMA": "Basomedial Amy",
    "CeA": "Central Amy",
    "ATZ": "Amy-Hpc Transition",
    "BST": "BNST",
    "CA1": "CA1",
    "CA2": "CA2",
    "CA3": "CA3",
    "CA4": "CA4",
    "DG":  "Dentate Gyrus",
    "S":   "Subiculum",
    "PHG-cos": "Parahipp (cos)",
    "PHG-l":   "Parahipp (lat)",
    "HCd": "Caudate (head)",
    "BCd": "Caudate (body)",
    "TCd": "Caudate (tail)",
    "Pu":  "Putamen",
    "Acb": "Nucl. Accumbens",
    "MFG-s": "MFG (superior)",
    "MFG-i": "MFG (inferior)",
    "SFG-l": "SFG (lateral)",
    "CgGf-i": "Cing frontal (inf)",
    "CgGf-s": "Cing frontal (sup)",
    "CgGp-i": "Cing parietal (inf)",
    "CgGp-s": "Cing parietal (sup)",
    "CgGr-i": "Cing retrospl (inf)",
    "CgGr-s": "Cing retrospl (sup)",
    "SCG":    "Subcallosal Cing",
    "DTA":  "Thal anterior",
    "DTM":  "Thal medial (MD)",
    "DTLd": "Thal lat dorsal",
    "DTLv": "Thal lat ventral",
    "DTP":  "Thal posterior",
    "ILc":  "Thal intralam (caud)",
    "ILr":  "Thal intralam (rost)",
    "Pa":   "Thal paraventr",
    "R":    "Thal reticular",
    "LGd":  "LGN (dorsal)",
    "MG":   "MGN",
    "Hl":   "Lat Habenula",
    "Hm":   "Med Habenula",
    "SIG":  "Short Insular Gyri",
    "LIG":  "Long Insular Gyri",
    "AHA":  "Ant Hypothal",
    "ARH":  "Arcuate Nucl",
    "DMH":  "Dorsomed Hypothal",
    "LHA":  "Lat Hypothal (ant)",
    "LHM":  "Lat Hypothal (mamm)",
    "LHT":  "Lat Hypothal (tub)",
    "PHA":  "Post Hypothal",
    "PVH":  "Paraventric Hypothal",
    "VMH":  "Ventromed Hypothal",
}

# Color palette for parent regions
REGION_COLORS: dict[str, str] = {
    "Amygdala":          "#E74C3C",
    "BNST":              "#C0392B",
    "Hippocampus":       "#3498DB",
    "Striatum":          "#2ECC71",
    "Nucleus Accumbens": "#27AE60",
    "DLPFC":             "#F39C12",
    "Cingulate Cortex":  "#E67E22",
    "Thalamus":          "#9B59B6",
    "Habenula":          "#8E44AD",
    "Insula":            "#1ABC9C",
    "Hypothalamus":      "#E91E63",
}


def get_subregions_for_parent(parent: str) -> list[str]:
    """Return list of Allen acronyms belonging to a parent region."""
    return [k for k, v in SUBREGION_TO_PARENT.items() if v == parent]


def get_parent(acronym: str) -> str | None:
    """Return the parent region for an Allen acronym, or None if not mapped."""
    return SUBREGION_TO_PARENT.get(acronym)


def get_label(acronym: str) -> str:
    """Return a human-readable label for an Allen acronym."""
    return SUBREGION_LABELS.get(acronym, acronym)

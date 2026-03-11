"""
Mechanosensitive ion channel gene families relevant to focused ultrasound neuromodulation.

Gene families and specific members based on:
    Bhatt et al. 2024 - PMC10423872
    "Bhatt, D.L., et al. Focused ultrasound neuromodulation..."

Notes:
    - PIEZO1/PIEZO2 are listed under old symbols FAM38A/FAM38B in some Allen
      Human Brain Atlas datasets.
    - K2P channels use KCNK nomenclature in HUGO/Allen datasets.
    - The "mechanosensitive" subset of K2P is TREK-1 (KCNK2), TREK-2 (KCNK10),
      and TRAAK (KCNK4). Other K2P members included for comparison.
"""

# ---------------------------------------------------------------------------
# Piezo family
# ---------------------------------------------------------------------------
PIEZO_GENES = {
    "FAM38A": "PIEZO1",   # primary FUS mechanotransducer
    "FAM38B": "PIEZO2",   # somatosensory mechanotransducer
}

# ---------------------------------------------------------------------------
# K2P / two-pore-domain potassium channels (KCNK family)
# ---------------------------------------------------------------------------
K2P_MECHANOSENSITIVE = {
    "KCNK2":  "TREK-1",   # strongly mechanosensitive
    "KCNK10": "TREK-2",   # strongly mechanosensitive
    "KCNK4":  "TRAAK",    # strongly mechanosensitive
}

K2P_OTHER = {
    "KCNK1":  "TWIK-1",
    "KCNK3":  "TASK-1",
    "KCNK5":  "TASK-2",
    "KCNK6":  "TWIK-2",
    "KCNK7":  "KCNK7",
    "KCNK9":  "TASK-3",
    "KCNK12": "THIK-2",
    "KCNK13": "THIK-1",
    "KCNK15": "TASK-5",
    "KCNK16": "TALK-1",
    "KCNK17": "TALK-2",
    "KCNK18": "TRESK",
}

K2P_ALL = {**K2P_MECHANOSENSITIVE, **K2P_OTHER}

# ---------------------------------------------------------------------------
# TRP channels — subfamilies implicated in mechanosensation / FUS
# ---------------------------------------------------------------------------
TRP_GENES = {
    # TRPV subfamily
    "TRPV1": "TRPV1",     # heat + mechanical
    "TRPV2": "TRPV2",     # mechanical stretch
    "TRPV3": "TRPV3",
    "TRPV4": "TRPV4",     # osmotic / mechanical
    # TRPA subfamily
    "TRPA1": "TRPA1",     # activated by low-intensity US in astrocytes
    # TRPC subfamily
    "TRPC1": "TRPC1",     # stretch-activated, US calcium responses
    "TRPC3": "TRPC3",
    "TRPC4": "TRPC4",
    "TRPC5": "TRPC5",     # mechanosensitive
    "TRPC6": "TRPC6",     # mechanosensitive
    # TRPM subfamily
    "TRPM2": "TRPM2",
    "TRPM3": "TRPM3",
    "TRPM4": "TRPM4",     # mechanosensitive
    "TRPM7": "TRPM7",     # Mg2+-permeable mechanosensor
    "TRPM8": "TRPM8",
}

# ---------------------------------------------------------------------------
# TRPP / Polycystin channels (mechanosensitive to flow / pressure)
# ---------------------------------------------------------------------------
TRPP_GENES = {
    "PKD1":   "TRPP1/PC1",
    "PKD2":   "TRPP2/PC2",   # primary mechanosensitive polycystin
    "PKD2L1": "TRPP3",
}

# ---------------------------------------------------------------------------
# Convenience collections
# ---------------------------------------------------------------------------

# Core mechanosensitive genes most directly implicated in FUS
FUS_CORE_GENES = {
    **PIEZO_GENES,
    **K2P_MECHANOSENSITIVE,
    "TRPV1": "TRPV1",
    "TRPV4": "TRPV4",
    "TRPA1": "TRPA1",
    "TRPC1": "TRPC1",
    "TRPM4": "TRPM4",
    "TRPM7": "TRPM7",
    "PKD2":  "TRPP2/PC2",
}

# All genes across all families
ALL_MECHANOSENSITIVE = {
    **PIEZO_GENES,
    **K2P_ALL,
    **TRP_GENES,
    **TRPP_GENES,
}

# Family labels for each gene (for annotation in plots)
GENE_FAMILY = {}
for g in PIEZO_GENES:
    GENE_FAMILY[g] = "Piezo"
for g in K2P_MECHANOSENSITIVE:
    GENE_FAMILY[g] = "K2P (mechano)"
for g in K2P_OTHER:
    GENE_FAMILY[g] = "K2P (other)"
for g in TRP_GENES:
    if g.startswith("TRPV"):
        GENE_FAMILY[g] = "TRPV"
    elif g.startswith("TRPA"):
        GENE_FAMILY[g] = "TRPA"
    elif g.startswith("TRPC"):
        GENE_FAMILY[g] = "TRPC"
    elif g.startswith("TRPM"):
        GENE_FAMILY[g] = "TRPM"
for g in TRPP_GENES:
    GENE_FAMILY[g] = "TRPP"


def get_display_name(gene_symbol: str) -> str:
    """Return a human-friendly display name for a gene symbol."""
    if gene_symbol in ALL_MECHANOSENSITIVE:
        alias = ALL_MECHANOSENSITIVE[gene_symbol]
        if alias != gene_symbol:
            return f"{alias} ({gene_symbol})"
    return gene_symbol


def get_family(gene_symbol: str) -> str:
    """Return the channel family label for a gene symbol."""
    return GENE_FAMILY.get(gene_symbol, "Unknown")

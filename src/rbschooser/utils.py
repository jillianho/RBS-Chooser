"""Utility functions for sequence processing and data parsing."""

import math
from typing import Dict, Set, Tuple

from bioe234_mentor.homeworks.RBSChooser.tools import (
    translate,
    edit_distance,
    hairpin_counter,
)

from .data import RBSOption


def reverse_complement(seq: str) -> str:
    """
    Compute the reverse complement of a DNA sequence.
    
    Args:
        seq: DNA sequence string.
        
    Returns:
        The reverse complement of the input sequence. Ambiguous bases
        (like 'N') are preserved as 'N'.
        
    Example:
        >>> reverse_complement("ATGC")
        'GCAT'
    """
    complement_map = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(complement_map.get(base.upper(), "N") for base in reversed(seq))


def parse_genbank(filename: str) -> Dict[str, Dict[str, str]]:
    """
    Parse a GenBank file to extract genes with their sequences and locations.
    
    Extracts the genome sequence from the ORIGIN section and parses all CDS
    (coding sequence) features with their qualifiers (gene name, locus tag).
    For each gene, returns the UTR (50 bases upstream) and CDS sequences.
    
    Args:
        filename: Path to the GenBank format file.
        
    Returns:
        Dictionary mapping locus tags to gene information dicts containing:
        - 'gene': Gene name
        - 'strand': '+' or '-' (forward/reverse strand)
        - 'start': 1-based start position
        - 'end': 1-based end position
        - 'cds': Coding sequence
        - 'utr_50': 50 bp upstream sequence in coding orientation
        
    Example:
        >>> genes = parse_genbank("sequence.gb")
        >>> genes["b0001"]["gene"]
        'thrA'
    """
    with open(filename, "r") as f:
        lines = f.readlines()

    # Extract genome from ORIGIN section
    genome_chunks = []
    in_origin = False

    for line in lines:
        if line.startswith("ORIGIN"):
            in_origin = True
            continue
        if in_origin:
            if line.startswith("//"):
                break
            # ORIGIN lines format: "     1 atgcat... tgca"
            parts = line.split()
            if len(parts) >= 2:
                # First token is line number; join sequence fragments
                genome_chunks.append("".join(parts[1:]))

    genome = "".join(genome_chunks).upper()

    # Parse CDS features
    genes = {}
    i = 0

    while i < len(lines):
        line = lines[i]

        if line.startswith("     CDS"):
            location = line.strip().split()[-1]

            # Handle complement (reverse strand)
            strand = "+"
            if location.startswith("complement("):
                strand = "-"
                location = location.replace("complement(", "").replace(")", "")

            # Skip complex locations like join(...)
            if ".." not in location:
                i += 1
                continue

            parts = location.split("..")
            if len(parts) != 2 or (not parts[0].isdigit()) or (not parts[1].isdigit()):
                i += 1
                continue

            start = int(parts[0])  # 1-based inclusive
            end = int(parts[1])    # 1-based inclusive

            # Read qualifiers
            locus_tag = None
            gene_name = None

            i += 1
            while i < len(lines) and lines[i].startswith("                     "):
                qualifier = lines[i]
                if "/locus_tag=" in qualifier:
                    locus_tag = qualifier.split('"')[1]
                if "/gene=" in qualifier:
                    gene_name = qualifier.split('"')[1]
                i += 1

            if locus_tag is None:
                continue
            if gene_name is None:
                gene_name = locus_tag

            # Extract CDS sequence
            cds_seq = genome[start - 1:end]

            # Extract 50nt UTR in coding strand orientation
            if strand == "+":
                utr_start = max(0, (start - 1) - 50)
                utr_end = start - 1
                utr = genome[utr_start:utr_end]
            else:
                # For reverse strand, extract and reverse complement
                utr_raw = genome[end:min(len(genome), end + 50)]
                utr = reverse_complement(utr_raw)
                cds_seq = reverse_complement(cds_seq)

            genes[locus_tag] = {
                "gene": gene_name,
                "strand": strand,
                "start": start,
                "end": end,
                "cds": cds_seq,
                "utr_50": utr,
            }

        else:
            i += 1

    return genes


def get_top_expressed_genes(
    filename: str, top_percent: float = 0.05
) -> Set[str]:
    """
    Extract locus tags of the top expressed genes from a proteomics file.
    
    Reads a two-column file (locus_tag, abundance_value) and returns the
    locus tags of genes in the top percentile by expression level.
    
    Args:
        filename: Path to proteomics data file with columns:
            [locus_tag, abundance_value]
        top_percent: Fraction of genes to return by abundance (default 0.05 = top 5%).
        
    Returns:
        Set of locus tag strings for top-expressed genes.
        
    Example:
        >>> top_genes = get_top_expressed_genes("proteome.txt", top_percent=0.1)
        >>> len(top_genes)
        150  # If there are 1500 total genes
    """
    abundance = []

    with open(filename, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            
            cols = line.split()
            if len(cols) < 2:
                continue

            # Extract locus tag (last component of identifier)
            locus_tag = cols[0].split(".")[-1]
            
            try:
                value = float(cols[1])
            except ValueError:
                continue

            abundance.append((locus_tag, value))

    # Sort by abundance descending
    abundance.sort(key=lambda x: x[1], reverse=True)

    if len(abundance) == 0:
        return set()

    # Compute number to keep
    n_keep = max(1, math.ceil(len(abundance) * top_percent))
    return set(tag for tag, _ in abundance[:n_keep])


def first_six_amino_acids(cds: str) -> str:
    """
    Translate a CDS and return the first six amino acids.
    
    Args:
        cds: A coding sequence starting with ATG.
        
    Returns:
        String of amino acids. May be shorter than 6 if CDS is very short.
        
    Example:
        >>> first_six_amino_acids("ATGGACTATAAAGAC")
        'MDYK'
    """
    if not cds:
        return ""
    protein = translate(cds)
    return protein[:6]


def occlusion_score(utr: str, cds_start: str = "") -> int:
    """
    Predict secondary structure interference in the UTR region.
    
    Counts hairpin structures in the UTR and the beginning of the CDS.
    More hairpins = more structural occlusion = worse RBS accessibility.
    
    Args:
        utr: The UTR sequence (typically 50 bp).
        cds_start: The first ~30 bases of the CDS (optional, for context).
        
    Returns:
        Number of hairpin structures detected.
        
    Example:
        >>> score = occlusion_score("AGGAGGUAAA")
        >>> score >= 0
        True
    """
    # Analyze UTR + first 30 bp of CDS for secondary structures
    window = utr + cds_start[:30] if cds_start else utr
    return hairpin_counter(window)


def peptide_distance(opt_first_six: str, input_first_six: str) -> int:
    """
    Calculate edit distance between two 6-AA peptide sequences.
    
    Lower distance = more similar N-terminal sequences = lower risk
    of protein misfolding when replacing the RBS.
    
    Args:
        opt_first_six: First six amino acids of the source RBS option.
        input_first_six: First six amino acids of the target protein.
        
    Returns:
        Edit distance (number of single-character edits).
        
    Example:
        >>> peptide_distance("MATVGR", "MATVGS")
        1
    """
    return edit_distance(opt_first_six, input_first_six)


def option_rank(
    opt: RBSOption, input_first_six: str, cds: str
) -> Tuple[int, int, str]:
    """
    Compute a ranking tuple for an RBS option.
    
    Lower-scoring options are "better" candidates. The ranking considers:
    1. Occlusion score (hairpin count in UTR)
    2. Peptide distance (amino acid similarity)
    3. Gene name (for tie-breaking)
    
    Args:
        opt: The RBSOption to rank.
        input_first_six: First six amino acids of the target gene.
        cds: The CDS sequence (for context).
        
    Returns:
        A tuple (occlusion, peptide_dist, gene_name) that sorts well
        with Python's min() function for ascending order.
        
    Example:
        >>> # Lower tuple values = better rank
        >>> rank = option_rank(option, "MATVGR", "ATGGCC...")
    """
    occlusion = occlusion_score(opt.utr, cds[:30])
    peptide_dist = peptide_distance(opt.first_six_aas, input_first_six)
    return (occlusion, peptide_dist, opt.gene_name)

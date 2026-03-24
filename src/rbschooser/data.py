"""Data models for RBS Chooser."""

from dataclasses import dataclass


@dataclass(frozen=True)
class RBSOption:
    """
    Represents a candidate RBS (Ribosomal Binding Site) option.
    
    This dataclass stores information about a potential RBS sequence from
    a high-expression source gene that can be used to improve translation
    initiation of target genes.
    
    Attributes:
        utr: The untranslated region (UTR) sequence containing the RBS
            and spacer region (typically 50 nt upstream of start codon).
        cds: The coding sequence (CDS) starting with ATG start codon.
        gene_name: Name or locus tag identifier for the source gene.
        first_six_aas: The first six amino acids translated from the CDS,
            used for peptide similarity comparison.
    """
    utr: str
    cds: str
    gene_name: str
    first_six_aas: str

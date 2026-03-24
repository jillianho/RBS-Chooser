"""Validation functions for DNA sequences."""


def validate_cds(cds: str) -> None:
    """
    Validate a coding sequence (CDS) for correctness.
    
    Checks that the CDS:
    - Is a string containing only A, T, C, G bases
    - Has length that is a multiple of 3 (complete codons)
    - Starts with the ATG start codon
    - Contains no internal stop codons (TAA, TAG, TGA)
    
    Args:
        cds: The DNA coding sequence to validate.
        
    Raises:
        ValueError: If the CDS fails any validation check.
        
    Example:
        >>> validate_cds("ATGGCC")  # Valid
        >>> validate_cds("ATGC")    # Raises ValueError: not multiple of 3
    """
    if not isinstance(cds, str):
        raise ValueError("CDS must be a string.")
    
    cds = cds.upper()
    
    if len(cds) == 0:
        raise ValueError("CDS is empty.")
    
    # Check for valid bases
    for base in cds:
        if base not in {"A", "T", "C", "G"}:
            raise ValueError(
                f"Invalid base '{base}' found in CDS (only A/T/C/G allowed)."
            )
    
    # Check length is multiple of 3
    if len(cds) % 3 != 0:
        raise ValueError(
            f"CDS length ({len(cds)}) is not a multiple of 3."
        )
    
    # Check starts with ATG
    if not cds.startswith("ATG"):
        raise ValueError(
            f"CDS does not start with ATG start codon. Found '{cds[:3]}'."
        )
    
    # Check for internal stop codons
    stop_codons = {"TAA", "TAG", "TGA"}
    codons = [cds[i:i+3] for i in range(0, len(cds), 3)]
    
    if any(codon in stop_codons for codon in codons[:-1]):
        raise ValueError(
            "CDS contains an internal stop codon (TAA/TAG/TGA)."
        )

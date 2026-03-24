"""RBSChooser: Optimal RBS selection for synthetic biology."""

from .data import RBSOption
from .core import RBSChooser
from .validators import validate_cds
from .utils import (
    reverse_complement,
    parse_genbank,
    get_top_expressed_genes,
    first_six_amino_acids,
    occlusion_score,
    peptide_distance,
)

__version__ = "1.0.0"
__author__ = "Jillian Ho"

__all__ = [
    "RBSChooser",
    "RBSOption",
    "validate_cds",
    "reverse_complement",
    "parse_genbank",
    "get_top_expressed_genes",
    "first_six_amino_acids",
    "occlusion_score",
    "peptide_distance",
]

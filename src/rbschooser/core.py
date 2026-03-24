"""Core RBSChooser class for selecting optimal ribosomal binding sites."""

from typing import List, Optional, Set, Union

from .data import RBSOption
from .validators import validate_cds
from .utils import (
    parse_genbank,
    get_top_expressed_genes,
    first_six_amino_acids,
    option_rank,
)

from bioe234_mentor.homeworks.RBSChooser.tools import translate


class RBSChooser:
    """
    Selects optimal Ribosomal Binding Site (RBS) sequences for genes.
    
    This class implements a two-step workflow:
    1. initiate(): Load and preprocess genomic and proteomics data once
    2. run(): For each target gene, find the best matching RBS from high-expression
       source genes.
    
    The algorithm ranks RBS options based on:
    - Secondary structure occlusion (fewer hairpins = better)
    - N-terminal peptide similarity (lower edit distance = safer)
    
    Attributes:
        options: Cached list of RBSOption objects (initialized by initiate()).
        
    Example:
        >>> chooser = RBSChooser()
        >>> chooser.initiate(genbank_file="genome.gb",
        ...                  proteomics_file="proteome.txt")
        >>> result = chooser.run("ATGGCCATTGTAATGGGCC...")
        >>> print(result.gene_name)
        'b0001'
    """
    
    options: List[RBSOption] = []
    
    @classmethod
    def initiate(
        cls,
        genbank_file: str = "sequence.gb",
        proteomics_file: str = "511145-WHOLE_ORGANISM-integrated.txt",
        top_percent: float = 0.05,
    ) -> None:
        """
        Initialize the RBSChooser by loading and preprocessing genomic data.
        
        This expensive operation is done once and the results are cached in
        the class-level options list, allowing run() to be called multiple times
        without reprocessing.
        
        Args:
            genbank_file: Path to GenBank format file with genome annotation.
            proteomics_file: Path to proteomics abundance file (locus_tag, value).
            top_percent: Fraction of genes to consider (by expression level).
                         Default 0.05 = top 5% highly expressed genes.
                         
        Raises:
            FileNotFoundError: If either input file does not exist.
            ValueError: If files cannot be parsed correctly.
            
        Example:
            >>> RBSChooser.initiate(
            ...     genbank_file="E_coli_MG1655.gb",
            ...     proteomics_file="proteome_data.txt",
            ...     top_percent=0.1
            ... )
        """
        # Step 1: Parse genomic annotations
        genes = parse_genbank(genbank_file)
        
        # Step 2: Identify top-expressed genes
        top_tags = get_top_expressed_genes(proteomics_file, top_percent=top_percent)
        
        # Step 3: Merge: keep only genes in both datasets
        high_expression_genes = {
            tag: genes[tag] for tag in top_tags if tag in genes
        }
        
        # Step 4: Build RBSOption objects
        options_list = []
        for locus_tag, data in high_expression_genes.items():
            protein = translate(data["cds"])
            options_list.append(
                RBSOption(
                    utr=data["utr_50"],
                    cds=data["cds"],
                    gene_name=locus_tag,
                    first_six_aas=protein[:6],
                )
            )
        
        cls.options = options_list
    
    def run(
        self,
        cds: str,
        ignores: Optional[Union[str, Set[str], List[str], Set[RBSOption]]] = None,
    ) -> RBSOption:
        """
        Select the best matching RBS option for a target coding sequence.
        
        Validates the input CDS, compares its N-terminal peptide to all
        available RBS options, and returns the option with the lowest
        combined ranking score (best structured RBS + closest peptide match).
        
        Args:
            cds: The target coding sequence starting with ATG.
            ignores: Gene names, locus tags, or RBSOption objects to exclude
                     from the selection. Can be a string for a single gene,
                     or a set/list for multiple. Default is None (no exclusions).
                     
        Returns:
            The best-matching RBSOption for the input CDS.
            
        Raises:
            ValueError: If the input CDS is invalid (see validate_cds).
            RuntimeError: If no valid options remain after applying ignores.
            
        Example:
            >>> chooser = RBSChooser()
            >>> chooser.initiate()
            >>> result = chooser.run("ATGGCCATTGTAA...")
            >>> print(f"Use RBS from {result.gene_name}")
            Use RBS from b0001
            
            >>> # Exclude specific genes
            >>> result2 = chooser.run("ATGGCCATTGTAA...", ignores={"b0001", "b0002"})
        """
        # Validate the input CDS
        validate_cds(cds)
        
        # Compute first six amino acids of the target
        input_first_six = first_six_amino_acids(cds)
        
        # Normalize ignores parameter
        if ignores is None:
            ignores_set = set()
        elif isinstance(ignores, str):
            ignores_set = {ignores}
        elif isinstance(ignores, (set, list, tuple)):
            ignores_set = set(ignores)
        else:
            ignores_set = {ignores}
        
        # Check if ignores contains RBSOption objects or strings
        contains_options = any(isinstance(x, RBSOption) for x in ignores_set)
        
        # Filter available options
        if contains_options:
            # Exclude by object identity/equality
            available = [opt for opt in self.options if opt not in ignores_set]
        else:
            # Treat ignores as a set of gene names / locus tags
            available = [
                opt for opt in self.options if opt.gene_name not in ignores_set
            ]
        
        if len(available) == 0:
            raise RuntimeError(
                "No RBS options left after applying ignores. "
                "Reduce exclusions or reinitialize with different parameters."
            )
        
        # Find the best option (lowest ranking score)
        best_option = min(available, key=lambda opt: option_rank(opt, input_first_six, cds))
        
        return best_option

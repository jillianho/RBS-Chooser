# RBSChooser

**Optimal Ribosomal Binding Site (RBS) Selection for Synthetic Biology**

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

## Overview

RBSChooser is a Python tool for selecting optimal Ribosomal Binding Sites (RBS) to improve protein expression in engineered genetic systems. Given a target coding sequence, it recommends the best RBS from a database of high-expression source genes by analyzing:

1. **Secondary structure** – UTR hairpin structures that can block ribosome access
2. **Peptide similarity** – N-terminal amino acid composition that affects translation initiation

This tool is particularly useful for:
- Optimizing heterologous protein expression in *E. coli*
- Rational RBS engineering based on genomic and proteomics data
- Synthetic biology applications requiring precise translational control

## Features

- **Genomic data parsing** – Extracts genes, annotations, and sequences from GenBank files
- **Proteomics integration** – Identifies naturally high-expressing genes as RBS sources
- **Sequence validation** – Comprehensive CDS validation (start codon, stop codons, frame)
- **Intelligent ranking** – Scores RBS candidates on structural and sequence properties
- **Flexible filtering** – Exclude specific genes from candidate pool
- **Type-safe design** – Full type hints for better IDE support and reliability

## Installation

### From source (recommended for development)

```bash
git clone https://github.com/jillianho/RBSChooser
cd RBSChooser
pip install -e .
```

### Install with development tools

```bash
pip install -e ".[dev]"
```

This includes pytest, coverage, linting, and type checking tools.

## Quick Start

```python
from rbschooser import RBSChooser

# Initialize with genomic and proteomics data
chooser = RBSChooser()
chooser.initiate(
    genbank_file="E_coli_MG1655.gb",
    proteomics_file="proteome_data.txt",
    top_percent=0.05  # Use top 5% expressed genes
)

# Find best RBS for your target gene
target_cds = "ATGGCCATTGTAATGGGCCGCTGCAAGGGT"
result = chooser.run(target_cds)

print(f"Selected RBS from: {result.gene_name}")
print(f"RBS sequence: {result.utr}")
print(f"Recommended CDS: {result.cds}")
```

## Usage Guide

### 1. Initialize RBSChooser

The `initiate()` class method loads and preprocesses data once:

```python
RBSChooser.initiate(
    genbank_file="sequence.gb",           # GenBank file with genome annotation
    proteomics_file="proteome.txt",       # Expression levels (locus_tag, value)
    top_percent=0.05                      # Fraction of genes to consider
)
```

**Input file formats:**

- **GenBank**: Standard NCBI format with ORIGIN section and CDS features
- **Proteomics**: Two-column tab/space-separated format:
  ```
  locus_tag.1  abundance_value
  locus_tag.2  abundance_value
  ...
  ```

### 2. Run Selection

For each target gene, call `run()` with the CDS:

```python
chooser = RBSChooser()

# Simple case
best_rbs = chooser.run("ATGATGATG...")

# Exclude specific genes
best_rbs = chooser.run("ATGATGATG...", ignores={"b0001", "b0002"})

# Exclude multiple genes
best_rbs = chooser.run(
    "ATGATGATG...",
    ignores={"b0001", "b0002", "thrA"}  # By locus tag or gene name
)
```

### 3. Access Results

The returned `RBSOption` contains:

```python
result = chooser.run(target_cds)

# Attributes
result.gene_name        # Source gene identifier (locus tag)
result.utr              # 50 bp UTR sequence + RBS + spacer
result.cds              # Coding sequence from source gene
result.first_six_aas    # First 6 amino acids (for verification)
```

## Ranking Algorithm

RBSChooser ranks options using a multi-criterion scoring system:

```
Rank = (occlusion_score, peptide_distance, gene_name)
```

**Lower scores are better.** The ranking considers:

1. **Occlusion Score** – Number of hairpin structures in the UTR
   - More hairpins = worse RBS accessibility
   
2. **Peptide Distance** – Edit distance between first 6 amino acids
   - Lower distance = safer protein-level changes
   
3. **Gene Name** – Alphabetical tie-breaker for reproducibility

## API Reference

### Core Classes

#### `RBSChooser`

Main class for RBS selection.

**Methods:**
- `initiate(genbank_file, proteomics_file, top_percent)` – Load and preprocess data
- `run(cds, ignores=None)` – Select best RBS for a target CDS

#### `RBSOption` (dataclass)

```python
@dataclass(frozen=True)
class RBSOption:
    utr: str              # 50 bp UTR with RBS
    cds: str              # Coding sequence
    gene_name: str        # Gene identifier
    first_six_aas: str    # First 6 amino acids
```

### Utility Functions

```python
from rbschooser import (
    validate_cds,              # Validate CDS correctness
    reverse_complement,        # Compute reverse complement
    parse_genbank,            # Parse GenBank files
    get_top_expressed_genes,  # Extract top genes by expression
    first_six_amino_acids,    # Translate and get first 6 AAs
    occlusion_score,          # Compute UTR structure penalty
    peptide_distance,         # Calculate peptide similarity
)
```

## Testing

Run the test suite:

```bash
pytest                          # Run all tests
pytest -v                       # Verbose output
pytest --cov                    # With coverage report
pytest tests/test_validators.py # Specific test file
```

## Development

### Code style

This project uses [Black](https://github.com/psf/black) for code formatting and [isort](https://pycqa.github.io/isort/) for import organization.

```bash
black src/ tests/
isort src/ tests/
```

### Type checking

Run [mypy](http://mypy-lang.org/) for static type analysis:

```bash
mypy src/
```

### Linting

Check code quality with [pylint](https://pylint.pycqa.org/):

```bash
pylint src/rbschooser/
```

## Project Structure

```
RBSChooser/
├── README.md                 # This file
├── LICENSE                   # MIT License
├── pyproject.toml           # Package configuration
├── .gitignore               # Git ignore patterns
├── src/rbschooser/          # Main package
│   ├── __init__.py          # Package exports
│   ├── core.py              # RBSChooser class
│   ├── data.py              # RBSOption dataclass
│   ├── validators.py        # CDS validation
│   └── utils.py             # Helper functions
├── tests/                   # Unit tests
│   ├── __init__.py
│   ├── test_validators.py
│   └── test_utils.py
└── notebooks/               # Jupyter notebooks (examples)
    └── RBSChooser_Example.ipynb
```

## Dependencies

- **Python 3.8+**
- **bioe234-mentor** – Provides translate(), edit_distance(), hairpin_counter()

To install dependencies:

```bash
pip install -e .
```

## Background

**What is an RBS?**

The Ribosomal Binding Site (RBS) is a ~10 bp sequence in the mRNA 5' UTR that forms base pairs with the 16S rRNA of the ribosome. It determines the efficiency of translation initiation.

**Why optimize it?**

- **Problem**: Heterologous genes often express poorly in *E. coli* due to suboptimal RBS sequences
- **Solution**: Replace the native RBS with one from a naturally high-expressing gene
- **Benefit**: Often increases protein yield 5–100 fold without genetic engineering

**How does RBSChooser help?**

Instead of random selection or expensive experimental screening, RBSChooser:
1. Prefilters to naturally high-expressing genes (evolutionary selection)
2. Ranks candidates by structural and sequence properties (rational design)
3. Identifies safe candidates with similar N-terminal peptides (protein safety)

## Citation

If you use RBSChooser in your research, please cite:

```bibtex
@software{ho2026rbschooser,
  author = {Ho, Jillian},
  title = {RBSChooser: Optimal Ribosomal Binding Site Selection},
  year = {2026},
  url = {https://github.com/jillianho/RBSChooser}
}
```

## License

MIT License – See [LICENSE](LICENSE) for details.

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit changes (`git commit -m 'Add amazing feature'`)
4. Push to branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request
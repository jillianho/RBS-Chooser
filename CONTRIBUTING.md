# Contributing to RBSChooser

Thank you for your interest in contributing! This document provides guidelines and instructions.

## Code of Conduct

Be respectful, inclusive, and constructive in all interactions.

## Getting Started

### 1. Fork and Clone

```bash
git clone https://github.com/your-username/RBSChooser.git
cd RBSChooser
```

### 2. Set Up Development Environment

```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -e ".[dev]"
```

### 3. Create a Feature Branch

```bash
git checkout -b feature/your-feature-name
```

## Development Workflow

### Code Style

We follow [PEP 8](https://pep8.org/) with these tools:

- **Black** – Code formatting (line length 88)
- **isort** – Import organization
- **pylint** – Code quality
- **mypy** – Static type checking

Format your code before committing:

```bash
black src/ tests/
isort src/ tests/
pylint src/rbschooser/
mypy src/
```

### Testing

Write tests for new features. Run the test suite:

```bash
pytest -v                  # Verbose output
pytest --cov              # With coverage
pytest tests/test_*.py    # Specific tests
```

**Test coverage target:** ≥80% for new code

### Type Hints

All public functions should have type hints:

```python
def my_function(param1: str, param2: int) -> bool:
    """Docstring here."""
    return True
```

### Docstrings

Use Google-style docstrings:

```python
def process_sequence(seq: str, quality: float = 0.95) -> str:
    """
    Process a DNA sequence with quality filtering.
    
    Args:
        seq: DNA sequence string.
        quality: Quality threshold (0.0-1.0).
        
    Returns:
        Processed sequence.
        
    Raises:
        ValueError: If quality is out of range.
        
    Example:
        >>> process_sequence("ATGC", 0.9)
        'ATGC'
    """
    pass
```

## Commit Guidelines

- Use clear, descriptive commit messages
- Reference issues when relevant: `Fixes #42`
- Keep commits focused and atomic
- Example: `Add peptide_distance validation tests`

```bash
git add src/
git commit -m "Add peptide_distance validation tests

- Add unit tests for exact matches
- Add tests for single substitutions
- Fixes #42"
```

## Pull Request Process

1. Update documentation and tests
2. Ensure all tests pass: `pytest`
3. Check code style: `black`, `isort`, `pylint`
4. Create a Pull Request with:
   - Clear title and description
   - Reference to related issues
   - Summary of changes
   - Any breaking changes noted

## Reporting Issues

**Security vulnerability?** Please email privately instead of opening an issue.

**Bug report:** Include:
- Python version
- Operating system
- Minimal reproducible example
- Expected vs actual behavior

**Feature request:** Describe:
- Use case
- Why it's needed
- Potential implementation approach

## Documentation

- Update README.md for user-facing changes
- Update docstrings for API changes
- Add examples in docstrings
- Update type hints

## Areas for Contribution

- **Testing** – Increase coverage, edge case testing
- **Documentation** – Clarifications, examples, tutorials
- **Performance** – Optimize bottlenecks
- **Features** – New RBS ranking criteria, integration with tools
- **Packaging** – CI/CD, wheels, conda packages

## Questions?

Open an issue or discussion for questions about:
- Development setup
- Architecture decisions
- Design patterns

---

Thank you for making RBSChooser better! 🎉

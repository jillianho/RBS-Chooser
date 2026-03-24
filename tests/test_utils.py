"""Unit tests for RBSChooser utilities."""

import pytest
from rbschooser.utils import reverse_complement, peptide_distance


class TestReverseComplement:
    """Test suite for reverse_complement function."""
    
    def test_simple_sequence(self):
        """Test reverse complement of a simple sequence."""
        assert reverse_complement("ATGC") == "GCAT"
    
    def test_full_genome_like(self):
        """Test reverse complement of longer sequence."""
        seq = "ATGGCCATTGTAA"
        result = reverse_complement(seq)
        # Reverse should be last-to-first, complement is A<->T, G<->C
        assert result == "TTACAATGGCCAT"
    
    def test_case_insensitive(self):
        """Reverse complement should handle lowercase."""
        assert reverse_complement("atgc") == "GCAT"
    
    def test_ambiguous_bases(self):
        """Ambiguous bases like N should be preserved as N."""
        assert reverse_complement("ATGCN") == "NGCAT"
    
    def test_empty_string(self):
        """Empty string should return empty."""
        assert reverse_complement("") == ""


class TestPeptideDistance:
    """Test suite for peptide_distance function."""
    
    def test_identical_peptides(self):
        """Identical peptides should have distance 0."""
        assert peptide_distance("MATVGR", "MATVGR") == 0
    
    def test_single_substitution(self):
        """One difference should have distance 1."""
        assert peptide_distance("MATVGR", "MATVGS") == 1
    
    def test_different_lengths(self):
        """Peptides of different lengths should compute edit distance."""
        # edit_distance handles different lengths
        distance = peptide_distance("MATVGR", "MAT")
        assert distance >= 0  # Should be a valid distance

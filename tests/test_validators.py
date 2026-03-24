"""Unit tests for RBSChooser validators."""

import pytest
from rbschooser.validators import validate_cds


class TestValidateCDS:
    """Test suite for CDS validation."""
    
    def test_valid_cds(self):
        """Valid CDS should pass without raising."""
        validate_cds("ATGGCCATTGTAATGGGCCGCTGCAAGGGT")
    
    def test_empty_cds(self):
        """Empty CDS should raise ValueError."""
        with pytest.raises(ValueError, match="empty"):
            validate_cds("")
    
    def test_invalid_base(self):
        """CDS with invalid bases should raise ValueError."""
        with pytest.raises(ValueError, match="Invalid base"):
            validate_cds("ATGGCCN")
    
    def test_not_multiple_of_three(self):
        """CDS length not multiple of 3 should raise ValueError."""
        with pytest.raises(ValueError, match="not a multiple of 3"):
            validate_cds("ATGC")
    
    def test_missing_start_codon(self):
        """CDS without ATG start should raise ValueError."""
        with pytest.raises(ValueError, match="does not start with ATG"):
            validate_cds("GCCATTGTAATGGGCCGCTGCAAGGGT")
    
    def test_internal_stop_codon(self):
        """CDS with internal stop codon should raise ValueError."""
        with pytest.raises(ValueError, match="internal stop codon"):
            validate_cds("ATGTAAGCCATTGTAATGGGCCGCTGCAAGGGT")
    
    def test_not_string(self):
        """Non-string CDS should raise ValueError."""
        with pytest.raises(ValueError, match="must be a string"):
            validate_cds(12345)

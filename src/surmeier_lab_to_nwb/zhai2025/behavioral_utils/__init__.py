"""
Behavioral data utilities for Zhai et al. 2025.

This package provides utilities for parsing and processing behavioral data,
particularly AIM (Abnormal Involuntary Movement) scoring data from Excel files.
"""

from .aim_parser import parse_aim_excel_to_wide_format
from .test_suite import EXPECTED_VALUES, test_aim_data_accuracy

__all__ = ["parse_aim_excel_to_wide_format", "test_aim_data_accuracy", "EXPECTED_VALUES"]

"""
Test suite for AIM data parsing accuracy.

This module contains test cases with hard-coded expected values extracted
directly from the Excel files to verify parsing accuracy.
"""

from pathlib import Path
from typing import Any, Dict

import pandas as pd

# Hard-coded test cases extracted directly from Excel sheets
EXPECTED_VALUES = [
    # Sheet 11202017, Animal 4124, Session 2017-11-20
    {
        "sheet_name": 11202017,
        "session_date": "2017-11-20",
        "animal_id": 4124,
        "score_type": "axial",
        "expected": {"20min": 3.0, "40min": 3.0, "60min": 4.0, "80min": 4.0, "100min": 3.0},
    },
    {
        "sheet_name": 11202017,
        "session_date": "2017-11-20",
        "animal_id": 4124,
        "score_type": "limb",
        "expected": {"20min": 1.5, "40min": 2.5, "60min": 2.0, "80min": 2.0, "100min": 1.0},
    },
    {
        "sheet_name": 11202017,
        "session_date": "2017-11-20",
        "animal_id": 4124,
        "score_type": "orolingual",
        "expected": {"20min": 2.0, "40min": 2.0, "60min": 2.0, "80min": 1.5, "100min": 3.0},
    },
    # Sheet 11202017, Animal 4124, Session 2017-11-22
    {
        "sheet_name": 11202017,
        "session_date": "2017-11-22",
        "animal_id": 4124,
        "score_type": "axial",
        "expected": {"20min": 4.0, "40min": 3.0, "60min": 3.0, "80min": 3.5, "100min": 4.0, "120min": 3.0},
    },
    {
        "sheet_name": 11202017,
        "session_date": "2017-11-22",
        "animal_id": 4124,
        "score_type": "orolingual",
        "expected": {"20min": 1.0, "40min": 1.0, "60min": 2.0, "80min": 2.0, "100min": 2.0, "120min": 1.5},
    },
    # Sheet 02232018, Animal 1590, Session 2018-02-23
    {
        "sheet_name": 2232018,
        "session_date": "2018-02-23",
        "animal_id": 1590,
        "score_type": "axial",
        "expected": {"20min": 4.0, "40min": 4.0, "60min": 4.0, "80min": 4.0},
    },
    {
        "sheet_name": 2232018,
        "session_date": "2018-02-23",
        "animal_id": 1590,
        "score_type": "limb",
        "expected": {"20min": 3.0, "40min": 3.0, "60min": 3.5, "80min": 3.5},
    },
    # Sheet 02232018, Animal 1944, Session 2018-03-04
    {
        "sheet_name": 2232018,
        "session_date": "2018-03-04",
        "animal_id": 1944,
        "score_type": "limb",
        "expected": {"20min": 2.0, "40min": 1.5, "60min": 2.0, "80min": 2.0},
    },
    # Sheet 11042020, Animal 2519, Multiple sessions (this was the problematic one)
    {
        "sheet_name": 11042020,
        "session_date": "2020-11-04",
        "animal_id": 2519,
        "score_type": "orolingual",
        "expected": {"20min": 2.0, "40min": 2.0, "60min": 2.0, "80min": 2.0},
    },
    {
        "sheet_name": 11042020,
        "session_date": "2020-11-06",
        "animal_id": 2519,
        "score_type": "orolingual",
        "expected": {"20min": 2.0, "40min": 2.5, "60min": 2.0, "80min": 2.0},
    },
    {
        "sheet_name": 11042020,
        "session_date": "2020-11-09",
        "animal_id": 2519,
        "score_type": "orolingual",
        "expected": {"20min": 2.0, "40min": 2.5, "60min": 2.0, "80min": 2.0, "100min": 1.5},
    },
    {
        "sheet_name": 11042020,
        "session_date": "2020-11-11",
        "animal_id": 2519,
        "score_type": "orolingual",
        "expected": {"20min": 2.5, "40min": 2.0, "60min": 2.0, "80min": 2.0, "100min": 1.5, "120min": 1.0},
    },
    # Sheet 03202018, Animal 2587, Session 2018-03-20
    {
        "sheet_name": 3202018,
        "session_date": "2018-03-20",
        "animal_id": 2587,
        "score_type": "axial",
        "expected": {"20min": 4.0, "40min": 4.0, "60min": 4.0, "80min": 4.0},
    },
    {
        "sheet_name": 3202018,
        "session_date": "2018-03-20",
        "animal_id": 2587,
        "score_type": "orolingual",
        "expected": {"20min": 3.0, "40min": 2.5, "60min": 2.5, "80min": 2.5},
    },
    # Sheet 09052018, Animal 6571, Session 2018-09-05
    {
        "sheet_name": 9052018,
        "session_date": "2018-09-05",
        "animal_id": 6571,
        "score_type": "axial",
        "expected": {"20min": 3.5, "40min": 4.0, "60min": 4.0, "80min": 4.0},
    },
    # Sheet 04012018, Animal 2828, Session 2018-04-01
    {
        "sheet_name": 4012018,
        "session_date": "2018-04-01",
        "animal_id": 2828,
        "score_type": "limb",
        "expected": {"20min": 1.0, "40min": 2.0, "60min": 2.0, "80min": 2.0},
    },
    # Additional random test cases for more comprehensive testing
    # Sheet 02232018, Animal 1945, Session 2018-02-25
    {
        "sheet_name": 2232018,
        "session_date": "2018-02-25",
        "animal_id": 1945,
        "score_type": "axial",
        "expected": {"20min": 4.0, "40min": 4.0, "60min": 4.0, "80min": 4.0},
    },
    # Sheet 11202017, Animal 4123, Session 2017-11-24
    {
        "sheet_name": 11202017,
        "session_date": "2017-11-24",
        "animal_id": 4123,
        "score_type": "orolingual",
        "expected": {"20min": 3.0, "40min": 3.0, "60min": 3.0, "80min": 3.0, "100min": 2.0},
    },
    # Sheet 09052018, Animal 5940, Session 2018-09-08
    {
        "sheet_name": 9052018,
        "session_date": "2018-09-08",
        "animal_id": 5940,
        "score_type": "limb",
        "expected": {"20min": 2.5, "40min": 3.5, "60min": 2.5, "80min": 3.0, "100min": 3.5},
    },
    # Sheet 04012018, Animal 2636, Session 2018-04-03
    {
        "sheet_name": 4012018,
        "session_date": "2018-04-03",
        "animal_id": 2636,
        "score_type": "axial",
        "expected": {"20min": 4.0, "40min": 4.0, "60min": 4.0, "80min": 4.0, "100min": 4.0},
    },
    # Sheet 7182018, Animal 5871, Session 2018-07-20
    {
        "sheet_name": 7182018,
        "session_date": "2018-07-20",
        "animal_id": 5871,
        "score_type": "orolingual",
        "expected": {"20min": 4.0, "40min": 3.0, "60min": 2.0, "80min": 2.5, "100min": 2.0},
    },
    # Sheet 11042020, Animal 2516, Session 2020-11-09
    {
        "sheet_name": 11042020,
        "session_date": "2020-11-09",
        "animal_id": 2516,
        "score_type": "limb",
        "expected": {"20min": 2.0, "40min": 2.0, "60min": 2.5, "80min": 2.0, "100min": 2.5},
    },
    # Sheet 03202018, Animal 2586, Session 2018-03-25
    {
        "sheet_name": 3202018,
        "session_date": "2018-03-25",
        "animal_id": 2586,
        "score_type": "orolingual",
        "expected": {"20min": 2.0, "40min": 2.0, "60min": 2.5, "80min": 3.0, "100min": 2.5},
    },
]


def test_aim_data_accuracy(csv_path: Path) -> Dict[str, Any]:
    """
    Test the tidy CSV data against hard-coded expected values.

    Parameters
    ----------
    csv_path : Path
        Path to the tidy CSV file to test

    Returns
    -------
    Dict[str, Any]
        Test results with pass/fail status and details
    """
    print(f"Testing AIM Data Accuracy: {csv_path.name}")
    print("=" * 60)

    # Load the CSV data
    df = pd.read_csv(csv_path)

    # Test results
    results = {"total_tests": len(EXPECTED_VALUES), "passed": 0, "failed": 0, "failures": []}

    for i, test_case in enumerate(EXPECTED_VALUES):
        test_num = i + 1
        print(
            f"\nTest {test_num}: Animal {test_case['animal_id']}, {test_case['session_date']}, {test_case['score_type']}"
        )

        # Find the matching row in CSV
        mask = (
            (df["sheet_name"] == test_case["sheet_name"])
            & (df["session_date"] == test_case["session_date"])
            & (df["animal_id"] == test_case["animal_id"])
            & (df["score_type"] == test_case["score_type"])
        )

        matching_rows = df[mask]

        if len(matching_rows) == 0:
            print(f"  ❌ FAIL: No matching row found")
            results["failed"] += 1
            results["failures"].append({"test": test_num, "reason": "No matching row found", "expected": test_case})
            continue
        elif len(matching_rows) > 1:
            print(f"  ❌ FAIL: Multiple matching rows found ({len(matching_rows)})")
            results["failed"] += 1
            results["failures"].append(
                {"test": test_num, "reason": f"Multiple matching rows ({len(matching_rows)})", "expected": test_case}
            )
            continue

        # Compare values
        csv_row = matching_rows.iloc[0]
        test_passed = True
        mismatches = []

        for time_point, expected_val in test_case["expected"].items():
            csv_val = csv_row[time_point]

            # Handle NaN comparisons
            if pd.isna(expected_val) and pd.isna(csv_val):
                continue
            elif pd.isna(expected_val) or pd.isna(csv_val):
                test_passed = False
                mismatches.append(f"{time_point}: expected {expected_val}, got {csv_val}")
            elif abs(float(expected_val) - float(csv_val)) > 0.001:
                test_passed = False
                mismatches.append(f"{time_point}: expected {expected_val}, got {csv_val}")

        if test_passed:
            print(f"  ✅ PASS: All values match")
            results["passed"] += 1
        else:
            print(f"  ❌ FAIL: Value mismatches - {'; '.join(mismatches)}")
            results["failed"] += 1
            results["failures"].append(
                {"test": test_num, "reason": "Value mismatches", "mismatches": mismatches, "expected": test_case}
            )

    return results


def print_test_summary(results: Dict[str, Any]):
    """Print test summary."""
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    print(f"Total tests: {results['total_tests']}")
    print(f"Passed: {results['passed']}")
    print(f"Failed: {results['failed']}")
    print(f"Success rate: {100 * results['passed'] / results['total_tests']:.1f}%")

    if results["failures"]:
        print("\nFAILED TESTS:")
        for failure in results["failures"]:
            print(f"  Test {failure['test']}: {failure['reason']}")
    else:
        print("\n🎉 ALL TESTS PASSED! 🎉")

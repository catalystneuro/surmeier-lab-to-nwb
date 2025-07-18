#!/usr/bin/env python3
"""
Main script to generate and test AIM parsing into tidy format.

This script uses the behavioral_utils package to parse AIM scoring data
from Excel files and validate the results against known test cases.
"""

from pathlib import Path

try:
    # Relative imports when used as module
    from .aim_parser import parse_aim_excel_to_wide_format
    from .test_suite import print_test_summary, test_aim_data_accuracy
except ImportError:
    # Direct imports when run as script
    from aim_parser import parse_aim_excel_to_wide_format
    from test_suite import print_test_summary, test_aim_data_accuracy


def main():
    """Main function to process AIM data and run tests."""

    # Define paths (go up one level since we're now in behavioral_utils)
    base_dir = Path(__file__).parent.parent
    excel_path = base_dir / "assets" / "AIM testing_CDGI KO.xlsx"
    genotype_path = base_dir / "assets" / "data_connections_D2_figures_3_4_6_7_8.csv"
    output_path = base_dir / "assets" / "aim_data_tidy.csv"

    print("AIM Data Processing and Testing")
    print("=" * 50)
    print(f"Excel file: {excel_path}")
    print(f"Genotype file: {genotype_path}")
    print(f"Output CSV: {output_path}")
    print()

    # Check if input files exist
    if not excel_path.exists():
        print(f"❌ ERROR: Excel file not found: {excel_path}")
        return False

    if not genotype_path.exists():
        print(f"❌ ERROR: Genotype file not found: {genotype_path}")
        return False

    try:
        # Parse the Excel data
        print("📊 Parsing AIM data from Excel...")
        df = parse_aim_excel_to_wide_format(
            excel_path=excel_path, genotype_csv_path=genotype_path, output_path=output_path
        )

        print(f"✅ Successfully parsed {len(df)} rows of AIM data")
        print(f"💾 Saved tidy CSV to: {output_path}")
        print()

        # Run tests
        print("🧪 Running accuracy tests...")
        results = test_aim_data_accuracy(output_path)
        print_test_summary(results)

        # Return success status
        success = results["failed"] == 0
        if success:
            print("\n🎉 All tests passed! AIM data parsing is accurate.")
        else:
            print(f"\n❌ {results['failed']} tests failed. Please check the output above.")

        return success

    except Exception as e:
        print(f"❌ ERROR during processing: {e}")
        return False


if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)

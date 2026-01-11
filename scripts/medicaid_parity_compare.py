import pandas as pd
from pathlib import Path


def compare_results(py_path, r_path, out_path):
    py_df = pd.read_csv(py_path)
    r_df = pd.read_csv(r_path)
    merged = py_df.merge(r_df, on=["weighting", "method"], suffixes=("_py", "_r"))
    merged["att_diff"] = merged["att_py"] - merged["att_r"]
    merged["se_diff"] = merged["se_py"] - merged["se_r"]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_path, index=False)


def main():
    script_dir = Path(__file__).resolve().parent
    py_csv = script_dir / "medicaid_python_results.csv"
    r_csv = script_dir / "medicaid_r_results.csv"
    cmp_csv = script_dir / "medicaid_comparison.csv"

    compare_results(py_csv, r_csv, cmp_csv)
    print(f"Comparison: {cmp_csv}")


if __name__ == "__main__":
    main()

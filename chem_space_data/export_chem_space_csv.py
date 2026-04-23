from pathlib import Path
import csv
import json
import pickle
import sys


def load_chem_space(path):
    with Path(path).open("rb") as handle:
        return pickle.load(handle)


def to_json_string(value):
    if value is None:
        return ""
    return json.dumps(value, default=float)


def export_csv(chem_space, output_path):
    fieldnames = [
        "probe ID",
        "CD type",
        "primary",
        "secondary",
        "dG_md",
        "ddG_md",
        "Kd_md",
        "Kd_SDS/Kd_PFOS",
    ]

    with Path(output_path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()

        for probe_id, entry in chem_space.items():
            if entry.get("ddG_md") is None:
                continue

            writer.writerow(
                {
                    "probe ID": probe_id,
                    "CD type": entry.get("CD", ""),
                    "primary": to_json_string(entry.get("primary")),
                    "secondary": to_json_string(entry.get("secondary")),
                    "dG_md": to_json_string(entry.get("dG_md")),
                    "ddG_md": to_json_string(entry.get("ddG_md")),
                    "Kd_md": to_json_string(entry.get("Kd_md")),
                    "Kd_SDS/Kd_PFOS": to_json_string(entry.get("Kd_SDS/Kd_PFOS")),
                }
            )


def main():
    input_path = sys.argv[1] if len(sys.argv) > 1 else "chem_space.pkl"
    output_path = sys.argv[2] if len(sys.argv) > 2 else "chem_space_export.csv"

    chem_space = load_chem_space(input_path)
    export_csv(chem_space, output_path)
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()

from pathlib import Path
import pickle
import sys


def load_chem_space(path):
    """Load the chem_space pickle from disk."""
    with Path(path).open("rb") as handle:
        return pickle.load(handle)


def mean_of(pair):
    """Return the mean value from a [mean, error] style entry."""
    if pair is None:
        return None
    return float(pair[0])


def summarize_counts(chem_space):
    """Print high-level counts that are usually the first sanity checks."""
    total = len(chem_space)
    observed = sum(entry["ddG_md"] is not None for entry in chem_space.values())
    generated = total - observed

    print(f"Total systems: {total}")
    print(f"Observed systems: {observed}")
    print(f"Generated systems: {generated}")
    print()


def summarize_cd_types(chem_space):
    """Show how many entries belong to each cyclodextrin family."""
    counts = {}
    for entry in chem_space.values():
        counts[entry["CD"]] = counts.get(entry["CD"], 0) + 1

    print("Counts by CD type:")
    for cd_type, count in sorted(counts.items()):
        print(f"  {cd_type}: {count}")
    print()


def summarize_data_coverage(chem_space):
    """Check which MD-derived fields are present in the file."""
    fields = ["dG_md", "ddG_md", "Kb_md", "Kd_md", "Kd_SDS/Kd_PFOS", "delta_dG", "delta_ddG"]

    print("Data coverage:")
    for field in fields:
        count = sum(entry[field] is not None for entry in chem_space.values())
        print(f"  {field}: {count}")
    print()


def print_example_entries(chem_space, limit=3):
    """Print a few entries so the structure is obvious without opening the pickle."""
    print(f"Example entries (first {limit}):")
    for key in list(chem_space)[:limit]:
        entry = chem_space[key]
        print(f"  {key}:")
        print(f"    CD = {entry['CD']}")
        print(f"    primary = {entry['primary']}")
        print(f"    secondary = {entry['secondary']}")
        print(f"    dG_md = {entry['dG_md']}")
        print(f"    ddG_md = {entry['ddG_md']}")
        print(f"    Kb_md = {entry['Kb_md']}")
        print(f"    Kd_md = {entry['Kd_md']}")
        print(f"    Kd_SDS/Kd_PFOS = {entry['Kd_SDS/Kd_PFOS']}")
    print()


def print_top_entries(chem_space, field, limit=5, reverse=False):
    """
    Print top entries for a numeric [mean, error] field.

    For dG-like values, more negative is typically stronger, so the default sort
    order is ascending. Set reverse=True if you want larger values first.
    """
    ranked = []
    for key, entry in chem_space.items():
        value = mean_of(entry[field])
        if value is not None:
            ranked.append((key, value, entry["CD"]))

    ranked.sort(key=lambda item: item[1], reverse=reverse)

    print(f"Top {limit} entries by {field}:")
    for key, value, cd_type in ranked[:limit]:
        print(f"  {key}: {value:.6f} ({cd_type})")
    print()


def main():
    # Default to the local chem_space.pkl, but allow a custom path from the CLI.
    path = sys.argv[1] if len(sys.argv) > 1 else "chem_space.pkl"
    chem_space = load_chem_space(path)

    print(f"Loaded {path}")
    print()

    summarize_counts(chem_space)
    summarize_cd_types(chem_space)
    summarize_data_coverage(chem_space)
    print_example_entries(chem_space)

    # More negative dG / ddG values are usually the interesting end of the range.
    print_top_entries(chem_space, "dG_md")
    print_top_entries(chem_space, "ddG_md")

    # For Kb and the Kd ratio, larger values are usually the interesting end.
    print_top_entries(chem_space, "Kb_md", reverse=True)
    print_top_entries(chem_space, "Kd_SDS/Kd_PFOS", reverse=True)


if __name__ == "__main__":
    main()

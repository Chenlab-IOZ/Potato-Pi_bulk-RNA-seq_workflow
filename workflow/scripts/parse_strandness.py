#!/usr/bin/env python3

import argparse
import re
from pathlib import Path


def parse_fraction(lines, marker):
    for line in lines:
        if marker in line:
            match = re.search(r":\s*([0-9]+(?:\.[0-9]+)?)", line)
            if match:
                return float(match.group(1))
    raise ValueError(f"Cannot find fraction for: {marker}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert RSeQC infer_experiment output to a TSV decision."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--stranded-threshold", type=float, required=True)
    parser.add_argument("--unstranded-threshold", type=float, required=True)
    args = parser.parse_args()

    lines = Path(args.input).read_text().splitlines()
    if not any("PairEnd" in line for line in lines):
        raise ValueError("RSeQC did not identify the data as paired-end")

    failed = parse_fraction(lines, "Fraction of reads failed to determine")
    fr_raw = parse_fraction(lines, '1++,1--,2+-,2-+')
    rf_raw = parse_fraction(lines, '1+-,1-+,2++,2--')
    determined = fr_raw + rf_raw
    if determined <= 0:
        raise ValueError("No reads were usable for strandness inference")

    fr_fraction = fr_raw / determined
    rf_fraction = rf_raw / determined

    if fr_fraction >= args.stranded_threshold:
        strandness = "FR"
        classification = "stranded"
    elif rf_fraction >= args.stranded_threshold:
        strandness = "RF"
        classification = "stranded"
    elif max(fr_fraction, rf_fraction) <= args.unstranded_threshold:
        strandness = "unstranded"
        classification = "unstranded"
    else:
        strandness = "unstranded"
        classification = "ambiguous_as_unstranded"

    with Path(args.output).open("w") as handle:
        handle.write(
            "sample\tlayout\tfailed_fraction\tfr_fraction\t"
            "rf_fraction\tstrandness\tclassification\n"
        )
        handle.write(
            f"{args.sample}\tpaired-end\t{failed:.6f}\t"
            f"{fr_fraction:.6f}\t{rf_fraction:.6f}\t"
            f"{strandness}\t{classification}\n"
        )


if __name__ == "__main__":
    main()

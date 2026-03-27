#!/usr/bin/env python3
import sys
import os


def parse_depth_file(depth_path: str) -> None:
    with open(depth_path, "r") as fh:
        header = fh.readline().rstrip("\n").split("\t")

        # Columns layout (0-indexed):
        #   0         : contigName
        #   1         : contigLen
        #   2         : totalAvgDepth
        #   3,4       : <sample1>.bam, <sample1>-var
        #   5,6       : <sample2>.bam, <sample2>-var
        #   ...
        # Every even index >= 3 is an average-depth column.
        n_abund = (len(header) - 3) // 2
        if n_abund < 1:
            sys.exit(
                f"ERROR: Expected at least one abundance column pair, "
                f"but only found {len(header)} columns in header."
            )

        # Map each sample to its 0-based column index and output file handle.
        samples: list[tuple[str, int]] = []
        for i in range(n_abund):
            col_idx = 3 + i * 2          # depth column (0-indexed)
            raw_name = header[col_idx]
            sample_name = raw_name.removesuffix(".bam")
            samples.append((sample_name, col_idx))

        # Open one output file per sample.
        handles: dict[str, object] = {}
        try:
            for sample_name, _ in samples:
                out_path = f"{sample_name}.abund"
                handles[sample_name] = open(out_path, "w")

            for line in fh:
                if not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                contig = fields[0]
                for sample_name, col_idx in samples:
                    depth_val = fields[col_idx]
                    handles[sample_name].write(f"{contig}\t{depth_val}\n")
        finally:
            for fh_out in handles.values():
                fh_out.close()

    written = [f"{name}.abund" for name, _ in samples]
    print(f"Written {len(written)} abundance file(s): {', '.join(written)}", file=sys.stderr)


def main() -> None:
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {os.path.basename(sys.argv[0])} <depth_file>")

    depth_path = sys.argv[1]
    if not os.path.isfile(depth_path):
        sys.exit(f"ERROR: depth file not found: {depth_path}")

    parse_depth_file(depth_path)


if __name__ == "__main__":
    main()
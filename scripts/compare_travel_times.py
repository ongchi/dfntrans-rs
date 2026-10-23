"""Compares the travel times of DFNTrans (C) and dfntrans-rs.

Usage: python3 -I scripts/compare_travel_times.py <partime.dat> <travel_times.dat> [tolerance]

Both files have the particles of the output (with `allparticles_output`,
all of them) in the order of the seeding. A particle matches if its
travel time differs by less than the relative tolerance (default 1e-4) and
its final position by less than 1e-3 m. With matrix diffusion, DFNTrans
writes no final position; then the advective and the total time are
compared.
"""

import sys


def read_rows(path, skip_columns):
    header = ""
    rows = []
    with open(path) as file:
        for line in file:
            if line.startswith("#"):
                header = line
                continue
            values = line.split()
            if values:
                rows.append([float(value) for value in values[skip_columns:]])
    return header, rows


def relative_difference(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def main():
    if len(sys.argv) not in (3, 4):
        sys.exit(__doc__)
    tolerance = float(sys.argv[3]) if len(sys.argv) == 4 else 1e-4
    _, c_rows = read_rows(sys.argv[1], 0)
    # Without the particle number of dfntrans-rs
    rust_header, rust_rows = read_rows(sys.argv[2], 1)
    has_diffusion = "diffusion" in rust_header
    print(f"particles: C {len(c_rows)}, Rust {len(rust_rows)}; matrix diffusion: {has_diffusion}")

    matches = 0
    largest_time_difference = 0.0
    for index, (c_row, rust_row) in enumerate(zip(c_rows, rust_rows)):
        if has_diffusion:
            # C: steps, weight, advective, total, diffusion, beta, length
            # Rust: steps, weight, advective, total, diffusion, x, y, z, beta, length
            time_difference = max(relative_difference(c_row[2], rust_row[2]),
                                  relative_difference(c_row[3], rust_row[3]))
            position_difference = 0.0
            description = (f"C times {c_row[2]:.9e} {c_row[3]:.9e}, "
                           f"Rust times {rust_row[2]:.9e} {rust_row[3]:.9e}")
        else:
            # C: steps, weight, time, x, y, z, beta, length
            # Rust: steps, weight, time, x, y, z, beta, length
            time_difference = relative_difference(c_row[2], rust_row[2])
            position_difference = max(abs(a - b) for a, b in zip(c_row[3:6], rust_row[3:6]))
            description = (f"C time {c_row[2]:.9e} at {c_row[3:6]}, "
                           f"Rust time {rust_row[2]:.9e} at {rust_row[3:6]}")
        if time_difference < tolerance and position_difference < 1e-3:
            matches += 1
            largest_time_difference = max(largest_time_difference, time_difference)
        else:
            print(f"particle {index + 1}: {description}")
    print(f"{matches} of {min(len(c_rows), len(rust_rows))} particles match; "
          f"largest relative time difference of a match {largest_time_difference:.3e}")
    sys.exit(0 if matches == len(c_rows) == len(rust_rows) else 1)


if __name__ == "__main__":
    main()

"""Adds the point array boundary__<name>_inner to an ASCII .vtu mesh: the
points of boundary__<name> without those at the smallest and the largest z.

Usage: python3 -I scripts/add_inner_boundary.py <input.vtu> <name> <output.vtu>

The reference comparison with DFNTrans uses it: a particle which starts at
a corner of a fracture moves along the fracture edge, and there the
rounding of the two codes decides if it leaves the fracture.
"""

import re
import sys


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    input_path, name, output_path = sys.argv[1:]
    text = open(input_path).read()

    points_match = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", text, re.S)
    coordinates = [float(value) for value in points_match.group(1).split()]
    z_values = coordinates[2::3]

    array_pattern = rf'(<DataArray type="Int32" Name="boundary__{name}"[^>]*>)(.*?)(</DataArray>)'
    array_match = re.search(array_pattern, text, re.S)
    if array_match is None:
        sys.exit(f"no point array boundary__{name}")
    values = [int(value) for value in array_match.group(2).split()]
    on_boundary = [z for z, value in zip(z_values, values) if value > 0]
    lowest, highest = min(on_boundary), max(on_boundary)
    inner = [1 if value > 0 and lowest < z < highest else 0 for z, value in zip(z_values, values)]

    header = array_match.group(1).replace(f"boundary__{name}", f"boundary__{name}_inner")
    new_array = header + "\n" + " ".join(map(str, inner)) + "\n" + array_match.group(3)
    end = array_match.end()
    open(output_path, "w").write(text[:end] + "\n" + new_array + text[end:])


if __name__ == "__main__":
    main()

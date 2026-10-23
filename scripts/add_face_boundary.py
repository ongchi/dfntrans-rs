"""Adds the point array boundary__<name> to an ASCII .vtu mesh: the points
of the array boundary__<source> whose coordinate <axis> is <value>, within 1e-6
of the size of the mesh.

Usage: python3 -I scripts/add_face_boundary.py <input.vtu> <output.vtu> <source> <name> <x|y|z> <value>

For example the faces of the box domain of a dfngen-rs mesh, from its
array boundary__mesh_surface.
"""

import re
import sys


def main():
    if len(sys.argv) != 7:
        sys.exit(__doc__)
    input_path, output_path, source, name, axis, value = sys.argv[1:]
    text = open(input_path).read()

    points_match = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", text, re.S)
    coordinates = [float(number) for number in points_match.group(1).split()]
    axis_values = coordinates["xyz".index(axis)::3]
    size = max(coordinates) - min(coordinates)

    array_pattern = rf'(<DataArray type="(\w+)" Name="boundary__{source}"[^>]*>)(.*?)(</DataArray>)'
    array_match = re.search(array_pattern, text, re.S)
    if array_match is None:
        sys.exit(f"no point array boundary__{source}")
    values = [int(float(number)) for number in array_match.group(3).split()]
    target = float(value)
    tolerance = 1e-6 * size
    face = [1 if flag > 0 and abs(coordinate - target) <= tolerance else 0
            for coordinate, flag in zip(axis_values, values)]

    header = array_match.group(1).replace(f"boundary__{source}", f"boundary__{name}")
    new_array = header + "\n" + " ".join(map(str, face)) + "\n" + array_match.group(4)
    end = array_match.end()
    open(output_path, "w").write(text[:end] + "\n" + new_array + text[end:])
    print(f"boundary__{name}: {sum(face)} points at {axis} = {target}")


if __name__ == "__main__":
    main()

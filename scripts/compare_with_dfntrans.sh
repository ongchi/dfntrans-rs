#!/bin/sh
# Runs a case with dfntrans-rs and with the C DFNTrans of dfnWorks, and
# compares the travel times of the particles.
#
# Usage: scripts/compare_with_dfntrans.sh <case.toml> <work_folder> <DFNTrans binary>
#
# The case must use the seeding at_nodes or random. Build DFNTrans without
# fused multiply-add, as gcc on x86-64 does, for example with clang:
#   make CFLAGS="-lm -w -O3 -std=gnu89 -Wno-error=return-mismatch -ffp-contract=off"
# Write the output of all particles (output.all_particles = true), else a
# particle which exits in only one of the two codes moves the rows.
set -eu

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <case.toml> <work_folder> <DFNTrans binary>" >&2
    exit 1
fi
case_file=$1
work_folder=$2
dfntrans_binary=$(cd "$(dirname "$3")" && pwd)/$(basename "$3")
repository=$(cd "$(dirname "$0")/.." && pwd)

cargo run --release --quiet --manifest-path "$repository/Cargo.toml" --example export_dfntrans -- "$case_file" "$work_folder"
(cd "$work_folder/c" && "$dfntrans_binary" control.dat > dfntrans.log 2>&1 < /dev/null)
python3 -I "$repository/scripts/compare_travel_times.py" "$work_folder/c/output/partime.dat" "$work_folder/rust/travel_times.dat"

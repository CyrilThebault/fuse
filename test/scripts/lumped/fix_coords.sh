#!/bin/bash
#
# fix_coords.sh
#
# Purpose
# -------
# Remove the unnecessary time dimension from the latitude, longitude,
# and hruId variables in a CAMELS-SPAT meteorological forcing file.
#
# Usage
# -----
#
#   ./fix_coords.sh <input_file> <output_file>
#
# Example
# -------
#
#   ./fix_coords.sh \
#       ../CAN_05BB001/input/work/time/CAN_05BB001_daymet_lumped.nc \
#       ../CAN_05BB001/input/work/coords/CAN_05BB001_daymet_lumped.nc
#
# Processing
# ----------
#
# The input variables:
#
#   latitude(time,hru)
#   longitude(time,hru)
#   hruId(time,hru)
#
# are replaced with variables that do not contain the time dimension:
#
#   latitude(hru)
#   longitude(hru)
#   hruId(hru)
#
# The coordinate values are taken from the first time step.
#
# Requirements
# ------------
#
# NCO utilities:
#
#   ncwa
#   ncks
#

set -euo pipefail


# ======================================================================
# Configuration
# ======================================================================

coord_vars="latitude,longitude,hruId"


# ======================================================================
# Command-line arguments
# ======================================================================

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <input_file> <output_file>" >&2
    exit 1
fi

input_file=$1
output_file=$2


# ======================================================================
# Check requirements
# ======================================================================

for command_name in ncwa ncks; do
    if ! command -v "${command_name}" >/dev/null 2>&1; then
        echo "ERROR: Required command not found: ${command_name}" >&2
        echo "Install NCO before running this script." >&2
        exit 1
    fi
done


# ======================================================================
# Check input and output paths
# ======================================================================

if [[ ! -f "${input_file}" ]]; then
    echo "ERROR: Input file not found: ${input_file}" >&2
    exit 1
fi

if [[ "${input_file}" == "${output_file}" ]]; then
    echo "ERROR: Input and output files must be different." >&2
    exit 1
fi

mkdir -p "$(dirname "${output_file}")"


# ======================================================================
# Define intermediate file
# ======================================================================

coord_file="${output_file%.nc}_coords.nc"


# Remove an old intermediate file if one exists.
rm -f "${coord_file}"


# ======================================================================
# Report operation
# ======================================================================

echo
echo "Fixing forcing coordinate variables"
echo "  Input : ${input_file}"
echo "  Output: ${output_file}"


# ======================================================================
# Create output copy
# ======================================================================

cp "${input_file}" "${output_file}"


# ======================================================================
# Extract coordinates without the time dimension
# ======================================================================

# Select the first time step and average over the singleton time
# selection. This removes the time dimension from the coordinate
# variables.

ncwa -O \
    -a time \
    -v "${coord_vars}" \
    -d time,0,0 \
    "${input_file}" \
    "${coord_file}"


# ======================================================================
# Replace the original coordinate variables
# ======================================================================

# Remove the original time-dependent coordinate variables from the
# output file.

ncks -O \
    -x \
    -v "${coord_vars}" \
    "${output_file}" \
    "${output_file}"


# Append the coordinate variables without the time dimension.

ncks -A \
    -v "${coord_vars}" \
    "${coord_file}" \
    "${output_file}"


# ======================================================================
# Clean up
# ======================================================================

rm -f "${coord_file}"


# ======================================================================
# Complete
# ======================================================================

echo
echo "Coordinate preprocessing completed successfully."
echo "Created: ${output_file}"
#

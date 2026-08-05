#!/bin/bash
#
# ------------------------------------------------------------------------------
# Prepare a mizuRoute hydrofabric
#
# This script creates a compact mizuRoute-compatible hydrofabric from the
# complete hydrologic fabric produced by create_hydrofabric.R.
#
# Processing steps:
#
#   1. Extract the variables required by mizuRoute.
#   2. Convert river length from kilometres to metres.
#   3. Convert HRU area from square kilometres to square metres.
#   4. Write a compact hydrofabric containing only the variables required by
#      mizuRoute.
#
# Intermediate files are written beneath:
#
#   <output_directory>/work/
#
# Usage:
#
#   prepare_mizuroute_hydrofabric.sh \
#       hydrofabric.nc \
#       hydrofabric_mizuroute.nc
#
# Arguments:
#
#   hydrofabric.nc
#       Complete hydrofabric created by create_hydrofabric.R.
#
#   hydrofabric_mizuroute.nc
#       Output hydrofabric for mizuRoute.
#
#
# Example:
#
#   hydrofabric_merit=test/CAN_05BB001/distributed/input/hydrofabric_merit.nc
#   hydrofabric_mizuRoute=test/CAN_05BB001/distributed/input/hydrofabric_mizuRoute.nc
#
#   test/scripts/distributed/prepare_mizuroute_hydrofabric.sh \
#       "${hydrofabric_merit}" \
#       "${hydrofabric_mizuRoute}"
#
# Requirements:
#
#   NCO utilities:
#
#     ncks
#     ncap2
#
# ------------------------------------------------------------------------------

set -euo pipefail

# ------------------------------------------------------------------------------
# Check command-line arguments
# ------------------------------------------------------------------------------

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <input_hydrofabric.nc> <output_hydrofabric.nc>" >&2
    exit 1
fi

INPUT_FILE=$1
OUTPUT_FILE=$2

# ------------------------------------------------------------------------------
# Check requirements
# ------------------------------------------------------------------------------

for command_name in ncks ncap2; do
    if ! command -v "${command_name}" >/dev/null 2>&1; then
        echo "ERROR: Required command not found: ${command_name}" >&2
        echo "Install NCO before running this script." >&2
        exit 1
    fi
done

# ------------------------------------------------------------------------------
# Check input file
# ------------------------------------------------------------------------------

if [[ ! -f "${INPUT_FILE}" ]]; then
    echo "ERROR: Input file not found: ${INPUT_FILE}" >&2
    exit 1
fi

# ------------------------------------------------------------------------------
# Create output and working directories
# ------------------------------------------------------------------------------

OUTPUT_DIR=$(dirname "${OUTPUT_FILE}")

WORK_DIR="${OUTPUT_DIR}/work"
WORK_FILE="${WORK_DIR}/hydrofabric_work.nc"

mkdir -p \
    "${OUTPUT_DIR}" \
    "${WORK_DIR}"

# ------------------------------------------------------------------------------
# Display workflow
# ------------------------------------------------------------------------------

echo
echo "Preparing mizuRoute hydrofabric"
echo
echo "Input hydrofabric : ${INPUT_FILE}"
echo "Work directory    : ${WORK_DIR}"
echo "Output hydrofabric: ${OUTPUT_FILE}"
echo

# ------------------------------------------------------------------------------
# Step 1: Extract variables required by mizuRoute
# ------------------------------------------------------------------------------

echo "Step 1: Extract required variables"

ncks -O -h \
    -v segId,NextDownID,hruId,hruSegId,lengthkm,slope,unitarea \
    "${INPUT_FILE}" \
    "${WORK_FILE}"

echo

# ------------------------------------------------------------------------------
# Step 2: Convert units
# ------------------------------------------------------------------------------

echo "Step 2: Convert units"

ncap2 -O -h -s '
    length = lengthkm * 1000.0;
    length@long_name = "River reach length";
    length@units = "m";

    area = unitarea * 1000000.0;
    area@long_name = "Unit catchment area";
    area@units = "m2";
' \
    "${WORK_FILE}" \
    "${WORK_FILE}"

echo

# ------------------------------------------------------------------------------
# Step 3: Create compact hydrofabric
# ------------------------------------------------------------------------------

echo "Step 3: Create compact hydrofabric"

ncks -O -h \
    -v segId,NextDownID,hruId,hruSegId,length,slope,area \
    "${WORK_FILE}" \
    "${OUTPUT_FILE}"

echo

# ------------------------------------------------------------------------------
# Complete
# ------------------------------------------------------------------------------

echo "mizuRoute hydrofabric preparation completed successfully."
echo "Created: ${OUTPUT_FILE}"

#!/bin/bash
#
# ------------------------------------------------------------------------------
# Create a lumped-to-distributed runoff mapping
#
# This script creates a mizuRoute runoff-remapping file for the case in which
# the hydrologic model is spatially lumped while the routing network is
# spatially distributed.
#
# In this configuration, the hydrologic model produces one runoff value for
# the entire basin at each time step. The routing model, however, contains
# multiple river reaches and associated river-network HRUs.
#
# mizuRoute normally remaps runoff from hydrologic-model grid cells or HRUs
# onto river-network HRUs using spatial-overlap weights. Here there is only
# one runoff cell, so every routing HRU receives runoff from that same cell.
#
#                     lumped runoff cell (1,1)
#                              |
#              +---------------+---------------+
#              |               |               |
#              v               v               v
#           routing HRU 1   routing HRU 2   routing HRU N
#
# Consequently, every routing HRU is assigned
#
#   polyid      = hruId
#   nOverlaps   = 1
#   weight      = 1.0
#   i_index     = 1
#   j_index     = 1
#
# The remapping variables use a contiguous ragged-array representation:
#
#   polyid
#       One entry for each routing HRU.
#
#   nOverlaps
#       Number of runoff-grid cells overlapping each routing HRU.
#
#   data
#       One entry for every routing-HRU/runoff-cell overlap.
#
# Therefore,
#
#   size(data) = sum(nOverlaps)
#
# For this lumped configuration, every routing HRU overlaps exactly one runoff
# cell. Therefore,
#
#   nOverlaps(:) = 1
#   size(data)    = number of routing HRUs
#
# The mapping variables are first created using the existing "hru" dimension
# because hruId is defined on that dimension. NCO requires the source and
# destination variables to share a dimension when values are assigned.
#
# After the mapping variables are created:
#
#   1. ncks extracts only the variables required by the remapping file; and
#   2. ncrename renames the "hru" dimension to "polyid", matching the default
#      mizuRoute remapping-file convention.
#
#
# Usage:
#
#   create_lumped_to_hru_mapping.sh \
#       hydrofabric_mizuRoute.nc \
#       lumped_to_hru.nc
#
# Arguments:
#
#   hydrofabric_mizuRoute.nc
#       Distributed mizuRoute hydrofabric containing hruId(hru).
#
#   lumped_to_hru.nc
#       Output runoff-remapping file.
#
# Example:
#
#   hydrofabric_mizuRoute=test/CAN_05BB001/distributed/input/hydrofabric_mizuRoute.nc
#   mapping_file=test/CAN_05BB001/distributed/input/lumped_to_hru.nc
#
#   test/scripts/distributed/create_lumped_to_hru_mapping.sh \
#       "${hydrofabric_mizuRoute}" \
#       "${mapping_file}"
#
# Dependencies:
#
#   NCO
#       ncap2
#       ncks
#       ncrename
#
# ------------------------------------------------------------------------------

set -euo pipefail

# ------------------------------------------------------------------------------
# Read command-line arguments
# ------------------------------------------------------------------------------

if [[ $# -ne 2 ]]; then
    echo "Usage:" >&2
    echo "  $0 hydrofabric_mizuRoute.nc lumped_to_hru.nc" >&2
    exit 1
fi

INPUT_FILE=$1
OUTPUT_FILE=$2

# ------------------------------------------------------------------------------
# Check inputs
# ------------------------------------------------------------------------------

if [[ ! -f "${INPUT_FILE}" ]]; then
    echo "ERROR: Input file not found:" >&2
    echo "  ${INPUT_FILE}" >&2
    exit 1
fi

for command_name in ncap2 ncks ncrename; do
    if ! command -v "${command_name}" >/dev/null 2>&1; then
        echo "ERROR: ${command_name} is not available." >&2
        echo "Install the NCO utilities before running this script." >&2
        exit 1
    fi
done

OUTPUT_DIR=$(dirname "${OUTPUT_FILE}")
OUTPUT_NAME=$(basename "${OUTPUT_FILE}")
WORK_DIR="${OUTPUT_DIR}/work"
WORK_FILE="${WORK_DIR}/${OUTPUT_NAME%.nc}_work.nc"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${WORK_DIR}"

rm -f "${WORK_FILE}"
rm -f "${OUTPUT_FILE}"

# ------------------------------------------------------------------------------
# Display workflow
# ------------------------------------------------------------------------------

echo
echo "Creating lumped-to-HRU runoff mapping"
echo
echo "Input hydrofabric : ${INPUT_FILE}"
echo "Output mapping    : ${OUTPUT_FILE}"
echo "Work file         : ${WORK_FILE}"
echo

# ------------------------------------------------------------------------------
# Create mapping variables
#
# polyid and nOverlaps are initially defined on the existing hru dimension
# because hruId is defined on that dimension. NCO requires the source and
# destination variables to share a dimension during assignment.
#
# The data dimension contains one record for each runoff-cell/HRU overlap.
# Consequently,
#
#   size(data) = sum(nOverlaps)
#
# In this lumped configuration, every routing HRU overlaps the single runoff
# cell exactly once. Therefore,
#
#   nOverlaps(:) = 1
#   size(data)    = size(hru)
#
# All overlap records point to runoff-grid cell (1,1).
# ------------------------------------------------------------------------------

ncap2 -O -h -s '
    defdim("data",$hru.size);

    polyid[hru]    = int(hruId);
    nOverlaps[hru] = 1;

    weight[data]   = 1.0;
    i_index[data]  = 1;
    j_index[data]  = 1;

    polyid@long_name = "River-network HRU identifier";
    polyid@units     = "1";

    nOverlaps@long_name = "Number of overlapping runoff cells";
    nOverlaps@units     = "1";

    weight@long_name = "Areal remapping weight";
    weight@units     = "1";

    i_index@long_name = "Runoff-grid x index";
    i_index@units     = "1";

    j_index@long_name = "Runoff-grid y index";
    j_index@units     = "1";
' \
"${INPUT_FILE}" \
"${WORK_FILE}"

# ------------------------------------------------------------------------------
# Extract only the variables required by the mizuRoute remapping file
#
# ncap2 preserves the original hydrofabric variables in the work file. These
# variables are not needed in the runoff-remapping file, so ncks extracts only
#
#   polyid
#   nOverlaps
#   weight
#   i_index
#   j_index
#
# into the final output file.
# ------------------------------------------------------------------------------

ncks -O -h \
    -v polyid,nOverlaps,weight,i_index,j_index \
    "${WORK_FILE}" \
    "${OUTPUT_FILE}"

# ------------------------------------------------------------------------------
# Rename the routing-HRU dimension
#
# The mapping variables were created on the original hru dimension so that
# polyid could be copied directly from hruId. Rename that dimension to polyid,
# which is the default mizuRoute dimension name for routing HRUs in a runoff-
# remapping file.
#
# The final dimensions are therefore
#
#   polyid = number of routing HRUs
#   data   = sum(nOverlaps)
# ------------------------------------------------------------------------------

ncrename -O -d hru,polyid "${OUTPUT_FILE}"

# ------------------------------------------------------------------------------
# Remove intermediate work file
# ------------------------------------------------------------------------------

rm -f "${WORK_FILE}"

# ------------------------------------------------------------------------------
# Complete
# ------------------------------------------------------------------------------

echo
echo "Lumped-to-HRU mapping creation completed successfully."
echo "Created: ${OUTPUT_FILE}"

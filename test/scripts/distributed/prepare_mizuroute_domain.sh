#!/bin/bash
#
# ------------------------------------------------------------------------------
# Prepare a mizuRoute domain
#
# This script prepares the routing-network input files required by mizuRoute
# from MERIT-Basins river and catchment shapefiles.
#
# The workflow supports the test configuration in which:
#
#   - hydrologic simulations are spatially lumped; and
#   - routing is performed on a distributed river network.
#
# Processing steps:
#
#   1. Create a complete hydrologic fabric from the river and catchment
#      shapefiles.
#   2. Prepare a compact mizuRoute hydrofabric by converting variables to the
#      required units and extracting the variables needed by mizuRoute.
#   3. Create a runoff-remapping file that maps the single lumped runoff cell
#      to every distributed river-network HRU.
#
# Intermediate files created while preparing the compact hydrofabric are
# written beneath:
#
#   <output_directory>/work/
#
# Usage:
#
#   test/scripts/distributed/prepare_mizuroute_domain.sh \
#       rivers.shp \
#       catchments.shp \
#       metadata.csv \
#       hydrofabric_merit.nc \
#       hydrofabric_mizuRoute.nc \
#       lumped_to_hru.nc
#
# Arguments:
#
#   rivers.shp
#       River network shapefile.
#
#   catchments.shp
#       Catchment polygon shapefile.
#
#   metadata.csv
#       CSV file describing the variables and metadata to include in the
#       complete hydrologic fabric.
#
#   hydrofabric_merit.nc
#       Complete hydrologic fabric created from the MERIT-Basins shapefiles.
#
#   hydrofabric_mizuRoute.nc
#       Compact hydrofabric containing the variables required by mizuRoute.
#
#   lumped_to_hru.nc
#       Runoff-remapping file that maps the single lumped runoff cell to all
#       distributed river-network HRUs.
#
# Example:
#
#   catchment_shp=test/CAN_05BB001/distributed/input/geospatial/shp/CAN_05BB001_distributed_basin.shp
#   river_shp=test/CAN_05BB001/distributed/input/geospatial/shp/CAN_05BB001_distributed_river.shp
#   metadata=test/metadata/merit_basins_shapefile_metadata.csv
#   hydrofabric_merit=test/CAN_05BB001/distributed/input/hydrofabric_merit.nc
#   hydrofabric_mizuRoute=test/CAN_05BB001/distributed/input/hydrofabric_mizuRoute.nc
#   mapping_file=test/CAN_05BB001/distributed/input/lumped_to_hru.nc
#
#   prepare_mizuroute_domain.sh \
#       "${river_shp}" \
#       "${catchment_shp}" \
#       "${metadata}" \
#       "${hydrofabric_merit}" \
#       "${hydrofabric_mizuRoute}" \
#       "${mapping_file}"
#
# ------------------------------------------------------------------------------

set -euo pipefail

# ------------------------------------------------------------------------------
# Check command-line arguments
# ------------------------------------------------------------------------------

if [[ $# -ne 6 ]]; then
    echo "Usage:" >&2
    echo "  $0 rivers.shp catchments.shp metadata.csv \\" >&2
    echo "     hydrofabric_merit.nc hydrofabric_mizuRoute.nc lumped_to_hru.nc" >&2
    exit 1
fi

RIVER_FILE=$1
CATCHMENT_FILE=$2
METADATA_FILE=$3
HYDROFABRIC_FILE=$4
ROUTING_FILE=$5
MAPPING_FILE=$6

# ------------------------------------------------------------------------------
# Locate preprocessing scripts
# ------------------------------------------------------------------------------

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

CREATE_SCRIPT="${SCRIPT_DIR}/create_hydrofabric.R"
PREPARE_SCRIPT="${SCRIPT_DIR}/prepare_mizuroute_hydrofabric.sh"
MAPPING_SCRIPT="${SCRIPT_DIR}/create_lumped_to_hru_mapping.sh"

# ------------------------------------------------------------------------------
# Check input files
# ------------------------------------------------------------------------------

for input_file in \
    "${RIVER_FILE}" \
    "${CATCHMENT_FILE}" \
    "${METADATA_FILE}"
do
    if [[ ! -f "${input_file}" ]]; then
        echo "ERROR: Input file does not exist: ${input_file}" >&2
        exit 1
    fi
done

# ------------------------------------------------------------------------------
# Check preprocessing scripts
# ------------------------------------------------------------------------------

if [[ ! -f "${CREATE_SCRIPT}" ]]; then
    echo "ERROR: Missing script: ${CREATE_SCRIPT}" >&2
    exit 1
fi

if [[ ! -x "${PREPARE_SCRIPT}" ]]; then
    echo "ERROR: Missing or non-executable script: ${PREPARE_SCRIPT}" >&2
    exit 1
fi

if [[ ! -x "${MAPPING_SCRIPT}" ]]; then
    echo "ERROR: Missing or non-executable script: ${MAPPING_SCRIPT}" >&2
    exit 1
fi

# ------------------------------------------------------------------------------
# Display workflow
# ------------------------------------------------------------------------------

echo
echo "Preparing mizuRoute domain"
echo
echo "River network        : ${RIVER_FILE}"
echo "Catchments           : ${CATCHMENT_FILE}"
echo "Metadata             : ${METADATA_FILE}"
echo "Complete hydrofabric : ${HYDROFABRIC_FILE}"
echo "mizuRoute hydrofabric: ${ROUTING_FILE}"
echo "Runoff mapping       : ${MAPPING_FILE}"
echo

# ------------------------------------------------------------------------------
# Step 1: Create complete hydrologic fabric
# ------------------------------------------------------------------------------

echo "Step 1: Create hydrologic fabric"

Rscript \
    "${CREATE_SCRIPT}" \
    "${CATCHMENT_FILE}" \
    "${RIVER_FILE}" \
    "${METADATA_FILE}" \
    "${HYDROFABRIC_FILE}"

echo

# ------------------------------------------------------------------------------
# Step 2: Prepare compact mizuRoute hydrofabric
# ------------------------------------------------------------------------------

echo "Step 2: Prepare mizuRoute hydrofabric"

"${PREPARE_SCRIPT}" \
    "${HYDROFABRIC_FILE}" \
    "${ROUTING_FILE}"

echo

# ------------------------------------------------------------------------------
# Step 3: Create lumped-to-HRU runoff mapping
# ------------------------------------------------------------------------------

echo "Step 3: Create lumped-to-HRU runoff mapping"

"${MAPPING_SCRIPT}" \
    "${ROUTING_FILE}" \
    "${MAPPING_FILE}"

echo

# ------------------------------------------------------------------------------
# Complete
# ------------------------------------------------------------------------------

echo "mizuRoute domain preparation completed successfully."
echo
echo "Created:"
echo "  ${HYDROFABRIC_FILE}"
echo "  ${ROUTING_FILE}"
echo "  ${MAPPING_FILE}"

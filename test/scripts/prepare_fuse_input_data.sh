#!/bin/bash
#
# ------------------------------------------------------------------------------
# Prepare FUSE input data
#
# This script prepares meteorological forcing and streamflow observations from
# CAMELS-SPAT for use with FUSE.
#
# Processing steps:
#
#   1. Standardize the time coordinates in the forcing and streamflow files.
#   2. Remove the unnecessary time dimension from forcing coordinate variables.
#   3. Compute daily mean air temperature from daily minimum and maximum
#      temperature.
#   4. Merge meteorological forcing and streamflow observations.
#   5. Compute Oudin potential evapotranspiration.
#   6. Add basin area and observed runoff depth.
#   7. Simplify global metadata.
#   8. Create a legacy-format FUSE input file.
#
# Intermediate files are written beneath:
#
#   <input_directory>/work/
#
# Usage:
#
#   prepare_fuse_input_data.sh forcing.nc streamflow.nc merged.nc
#
# Arguments:
#   forcing.nc     Original meteorological forcing file.
#   streamflow.nc  Original streamflow observation file.
#   merged.nc      Final merged FUSE input file.
#
# ------------------------------------------------------------------------------

set -euo pipefail

METADATA_FILE="../metadata/camels-spat-metadata.csv"

# ------------------------------------------------------------------------------
# Check command-line arguments
# ------------------------------------------------------------------------------

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 forcing_file streamflow_file merged_file" >&2
    exit 1
fi

FORCING_ORIG=$1
QOBS_ORIG=$2
MERGED_FILE=$3

# ------------------------------------------------------------------------------
# Locate the preprocessing scripts
#
# This assumes that the master script and the three sub-scripts are stored in
# the same directory.
# ------------------------------------------------------------------------------

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

TIME_SCRIPT="${SCRIPT_DIR}/standardize_time_coordinates.sh"
COORD_SCRIPT="${SCRIPT_DIR}/fix_coords.sh"
TEMP_SCRIPT="${SCRIPT_DIR}/compute_mean_temperature.sh"
OUDIN_SCRIPT="${SCRIPT_DIR}/compute_oudin_pet.R"
MERGE_SCRIPT="${SCRIPT_DIR}/merge_netcdf_files.sh"

# ------------------------------------------------------------------------------
# Check input files
# ------------------------------------------------------------------------------

for input_file in "${FORCING_ORIG}" "${QOBS_ORIG}" "${METADATA_FILE}"; do
    if [[ ! -f "${input_file}" ]]; then
        echo "ERROR: Input file does not exist: ${input_file}" >&2
        exit 1
    fi
done

# ------------------------------------------------------------------------------
# Check preprocessing scripts
# ------------------------------------------------------------------------------

for script_file in \
    "${TIME_SCRIPT}"  \
    "${COORD_SCRIPT}" \
    "${TEMP_SCRIPT}"  \
    "${OUDIN_SCRIPT}" \
    "${MERGE_SCRIPT}"
do
    if [[ ! -x "${script_file}" ]]; then
        echo "ERROR: Script is missing or not executable: ${script_file}" >&2
        exit 1
    fi
done

# ------------------------------------------------------------------------------
# Create output and working directories
#
# The work directory is placed inside the directory containing the final
# merged file.
# ------------------------------------------------------------------------------

INPUT_DIR=$(dirname "${MERGED_FILE}")

WORK_DIR="${INPUT_DIR}/work"
TIME_WORK_DIR="${WORK_DIR}/time"
COORD_WORK_DIR="${WORK_DIR}/coords"
TEMP_WORK_DIR="${WORK_DIR}/temperature"

mkdir -p \
    "${INPUT_DIR}" \
    "${TIME_WORK_DIR}" \
    "${COORD_WORK_DIR}" \
    "${TEMP_WORK_DIR}"

# ------------------------------------------------------------------------------
# Define intermediate filenames
# ------------------------------------------------------------------------------

FORCING_BASENAME=$(basename "${FORCING_ORIG}")
QOBS_BASENAME=$(basename "${QOBS_ORIG}")

FORCING_TIME="${TIME_WORK_DIR}/${FORCING_BASENAME}"
QOBS_TIME="${TIME_WORK_DIR}/${QOBS_BASENAME}"

FORCING_COORD="${COORD_WORK_DIR}/${FORCING_BASENAME}"

FORCING_TEMP="${TEMP_WORK_DIR}/${FORCING_BASENAME}"

# ------------------------------------------------------------------------------
# Display workflow
# ------------------------------------------------------------------------------

echo
echo "Preparing FUSE input data"
echo
echo "Forcing file:      ${FORCING_ORIG}"
echo "Streamflow file:   ${QOBS_ORIG}"
echo "Work directory:    ${WORK_DIR}"
echo "Merged output:     ${MERGED_FILE}"
echo

# ------------------------------------------------------------------------------
# Step 1: Standardize forcing time coordinates
# ------------------------------------------------------------------------------

echo "Step 1: Standardizing forcing and streamflow time coordinates"

"${TIME_SCRIPT}" \
    forcing \
    "${FORCING_ORIG}" \
    "${FORCING_TIME}"

"${TIME_SCRIPT}" \
    qobs \
    "${QOBS_ORIG}" \
    "${QOBS_TIME}"

echo

# ------------------------------------------------------------------------------
# Step 2: Remove the time dimension from forcing coordinate variables
#
# fix_coords.sh writes the processed file into COORD_WORK_DIR using the same
# basename as the input file.
# ------------------------------------------------------------------------------

echo "Step 2: Fixing forcing coordinate variables"

"${COORD_SCRIPT}" \
    "${FORCING_TIME}" \
    "${FORCING_COORD}"

echo


# ------------------------------------------------------------------------------
# Step 3: Compute daily mean temperature by averaging tmin and tmax
# ------------------------------------------------------------------------------

echo "Step 3: Computing daily mean temperature"

"${TEMP_SCRIPT}" \
    "${FORCING_COORD}" \
    "${FORCING_TEMP}"

echo


# ------------------------------------------------------------------------------
# Step 4: Merge forcing and streamflow observations
# ------------------------------------------------------------------------------

echo "Step 4: Merging forcing and streamflow observations"

"${MERGE_SCRIPT}" \
    "${FORCING_TEMP}" \
    "${QOBS_TIME}" \
    "${MERGED_FILE}"

echo

# ------------------------------------------------------------------------------
# Step 5: Compute Oudin PET
# ------------------------------------------------------------------------------

echo "Step 5: Compute Oudin PET"

Rscript "${OUDIN_SCRIPT}" "${MERGED_FILE}"

echo

# ------------------------------------------------------------------------------
# Step 6: Add variables runoff and basin area
# ------------------------------------------------------------------------------

echo "Step 6: Add runoff and basin area variables"

id_col=$(head -1 ${METADATA_FILE} | tr ',' '\n' | grep -nx Station_id | cut -d: -f1)
area_col=$(head -1 ${METADATA_FILE} | tr ',' '\n' | grep -nx Basin_area_km2 | cut -d: -f1)
echo id column: $id_col
echo area column: $area_col

station_id=$(basename "$(dirname "${INPUT_DIR}")")
station_id=${station_id#*_} # remove XXX_ (USA_, CAN_)
echo station id: $station_id

basin_area_km2=$(
    awk -F',' -v id_col="$id_col" -v station_id="$station_id" \
        '$id_col == station_id' \
        "${METADATA_FILE}" |
    cut -d',' -f"$area_col"
)
echo basin area: $basin_area_km2

# create new variables basin area and runoff depth
ncap2 -O -s "
    basin_area = ${basin_area_km2}f;
    runoff_obs[time] =
        q_obs * 86400.0f * 1000.0f /
        (basin_area * 1000000.0f);
" \
    "${MERGED_FILE}" "${MERGED_FILE}"

# Remove the original discharge variable
ncks -O -x -v q_obs "${MERGED_FILE}" "${MERGED_FILE}"

# add attrinutes
ncatted -O -h \
    -a long_name,basin_area,o,c,"basin area" \
    -a units,basin_area,o,c,"km2" \
    -a long_name,runoff_obs,o,c,"observed runoff depth" \
    -a units,runoff_obs,o,c,"mm day-1" \
    "${MERGED_FILE}"

echo

# ----------------------------------------------------------------------
# Step 7: Define concise global metadata for the final FUSE input file
# ----------------------------------------------------------------------

echo "Step 7: Clean up metadata"

creation_date=$(date -u +"%Y-%m-%dT%H:%M:%SZ")

ncatted -O -h \
    -a History,global,d,, \
    -a history,global,d,, \
    -a history_of_appended_files,global,d,, \
    -a Author,global,d,, \
    -a License,global,d,, \
    -a Source,global,d,, \
    -a easymore_hash,global,d,, \
    -a NCO,global,d,, \
    -a conventions,global,d,, \
    -a title,global,o,c,"FUSE meteorological forcing and streamflow observations" \
    -a Conventions,global,o,c,"CF-1.10" \
    -a institution,global,o,c,"University of Calgary" \
    -a source,global,o,c,"Meteorological forcing and streamflow observations derived from CAMELS-SPAT" \
    -a forcing_source,global,o,c,"Daymet meteorological forcing distributed through CAMELS-SPAT" \
    -a streamflow_source,global,o,c,"Water Survey of Canada streamflow observations distributed through CAMELS-SPAT" \
    -a references,global,o,c,"CAMELS-SPAT dataset: doi:10.20383/103.01306; Knoben et al. (2025): doi:10.5194/hess-29-5791-2025" \
    -a processing,global,o,c,"Time coordinates standardized; coordinate variables corrected; daily mean air temperature calculated from minimum and maximum temperature; meteorological forcing merged with streamflow observations." \
    -a history,global,o,c,"${creation_date}: Created by prepare_fuse_input_data.sh" \
    "${MERGED_FILE}"

echo

# ------------------------------------------------------------------------------
# Step 8: Convert files to legacy format to support old code
# ------------------------------------------------------------------------------

echo "Step 8: Convert files to legacy format"

base=$(basename "$MERGED_FILE" .nc)
LEGACY_FILE="$(dirname "$MERGED_FILE")/${base}_legacy.nc"

# remove the singleton hru dimension
ncwa -O -h -a hru "${MERGED_FILE}" "${LEGACY_FILE}"

# add latitude and longitude variables
ncecat -O -h -u latitude -v prcp,temp,pet,pet_oudin,runoff_obs "${LEGACY_FILE}" "${LEGACY_FILE}"
ncecat -O -h -u longitude -v prcp,temp,pet,pet_oudin,runoff_obs "${LEGACY_FILE}" "${LEGACY_FILE}"

# put the dimensions in the correct order
ncpdq -O -h -a time,latitude,longitude "${LEGACY_FILE}" "${LEGACY_FILE}"

# Create the coordinate variables directly from the merged file values
lat=$(ncks -H -C -v latitude -s '%g\n' "${MERGED_FILE}")
lon=$(ncks -H -C -v longitude -s '%g\n' "${MERGED_FILE}")
echo ${lat} ${lon}

ncap2 -O -h -s "
    latitude[latitude]=${lat};
    longitude[longitude]=${lon};
" "${LEGACY_FILE}" "${LEGACY_FILE}"

echo

# ------------------------------------------------------------------------------
# Complete
# ------------------------------------------------------------------------------

echo "FUSE input preparation completed successfully."
echo "Created: ${MERGED_FILE}"

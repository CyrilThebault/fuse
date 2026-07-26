#!/bin/bash
#
# ------------------------------------------------------------------------------
# Merge meteorological forcing and streamflow observations
#
# This script creates a FUSE input file by combining meteorological forcing
# variables with streamflow observations using the forcing file as the master
# time axis.
#
# The script:
#   1. Determines the period common to the forcing and observation files.
#   2. Extracts the required forcing variables over the complete forcing period.
#   3. Creates a streamflow variable spanning the complete forcing period and
#      initializes it to NaN.
#   4. Inserts observed streamflow into the corresponding locations of the
#      forcing time series.
#   5. Writes a merged NetCDF file suitable for FUSE simulations.
#
# Assumptions:
#   - Both files use the same calendar and time units.
#   - Time coordinates are daily and contain no missing records.
#   - Streamflow observations form a contiguous time series.
#
# Usage:
#   merge_netcdf_files.sh forcing.nc streamflow.nc merged.nc
#
# Arguments:
#   forcing.nc     Meteorological forcing file.
#   streamflow.nc  Streamflow observation file.
#   merged.nc      Output merged NetCDF file.
#
# ------------------------------------------------------------------------------

set -euo pipefail

met_file=$1
qobs_file=$2
merged_file=$3

# Variable names used in the forcing and streamflow files
met_vars="prcp,temp,pet"
qobs_var="q_obs"

script_dir=$(cd "$(dirname "$0")" && pwd)
input_dir=$(dirname "$merged_file")

work_dir="${input_dir}/work"
mkdir -p "$work_dir"

merged_work="${work_dir}/merged_work.nc"
qobs_insert_file="${work_dir}/qobs_insert.nc"

# ----------------------------------------------------------------------
# Get the start and end times
# ----------------------------------------------------------------------

start_met=$(ncks -H -C -s '%.17g\n' \
    -v time_bnds -d time,0,0 -d nbnds,0,0 "$met_file")

end_met=$(ncks -H -C -s '%.17g\n' \
    -v time_bnds -d time,-1,-1 -d nbnds,1,1 "$met_file")

echo "Start/end time for meteorological forcing: $start_met $end_met"

start_qobs=$(ncks -H -C -s '%.17g\n' \
    -v time_bnds -d time,0,0 -d nbnds,0,0 "$qobs_file")

end_qobs=$(ncks -H -C -s '%.17g\n' \
    -v time_bnds -d time,-1,-1 -d nbnds,1,1 "$qobs_file")

echo "Start/end time for streamflow observations: $start_qobs $end_qobs"

# ----------------------------------------------------------------------
# Get the common period
# ----------------------------------------------------------------------

start_time=$(
    printf "%s\n%s\n" "$start_met" "$start_qobs" |
        sort -g |
        tail -1
)

end_time=$(
    printf "%s\n%s\n" "$end_met" "$end_qobs" |
        sort -g |
        head -1
)

echo "Start/end time for common period: $start_time $end_time"

# ----------------------------------------------------------------------
# Check that the files overlap
# ----------------------------------------------------------------------

minimum_common_time=$(
    printf "%s\n%s\n" "$start_time" "$end_time" |
        sort -g |
        head -1
)

if [[ "$minimum_common_time" != "$start_time" ]]; then
    echo "ERROR: Forcing and streamflow observations do not overlap." >&2
    exit 1
fi

# ----------------------------------------------------------------------
# Convert common-period times to forcing-array indices
#
# This assumes:
#
#   - daily data;
#   - time is measured in days;
#   - no missing records in either time axis.
# ----------------------------------------------------------------------

start_index=$(
    awk -v start="$start_time" -v origin="$start_met" \
        'BEGIN {printf "%.0f\n", start-origin}'
)

end_index=$(
    awk -v end="$end_time" -v origin="$start_met" \
        'BEGIN {printf "%.0f\n", end-origin-1}'
)

echo "Streamflow insertion indices: $start_index to $end_index"

# ----------------------------------------------------------------------
# Extract the complete forcing time series
# ----------------------------------------------------------------------

ncks -O -C \
    -v "time,time_bnds,latitude,longitude,${met_vars}" \
    "$met_file" \
    "$merged_work"

# ----------------------------------------------------------------------
# Create a full-length streamflow vector initialized to missing
# ----------------------------------------------------------------------

missing_value=-9999.0

ncap2 -O -s "
    ${qobs_var}[time] = float(${missing_value});
    ${qobs_var}.set_miss(${missing_value});
" \
    "$merged_work" \
    "$merged_work"

# ----------------------------------------------------------------------
# Extract streamflow over the period that overlaps the forcing data
# ----------------------------------------------------------------------

qobs_start_index=$(
    awk -v start="$start_time" -v origin="$start_qobs" \
        'BEGIN {printf "%.0f\n", start-origin}'
)

qobs_end_index=$(
    awk -v end="$end_time" -v origin="$start_qobs" \
        'BEGIN {printf "%.0f\n", end-origin-1}'
)

echo "Observation extraction indices: ${qobs_start_index} to ${qobs_end_index}"

forcing_records=$(( end_index - start_index + 1 ))
qobs_records=$(( qobs_end_index - qobs_start_index + 1 ))

printf "Forcing insertion:      %6d to %-6d = %6d records\n" \
    "$start_index" \
    "$end_index" \
    "$forcing_records"

printf "Observation extraction: %6d to %-6d = %6d records\n" \
    "$qobs_start_index" \
    "$qobs_end_index" \
    "$qobs_records"

ncks -O -C \
    -v "time,${qobs_var}" \
    -d time,"${qobs_start_index}","${qobs_end_index}" \
    "$qobs_file" \
    "$qobs_insert_file"

echo Created $qobs_insert_file

# ----------------------------------------------------------------------
# Insert the observations into the full-length missing vector
# ----------------------------------------------------------------------

# Give the observation vector a temporary variable and dimension name
#   - this allows it to coexist with the full forcing time dimension.
ncrename -O \
    -d time,obs_time \
    -v time,obs_time \
    -v "${qobs_var},${qobs_var}_insert" \
    "$qobs_insert_file"


# Append the temporary observation vector to the working file
ncks -A \
    "$qobs_insert_file" \
    "$merged_work"

# Insert the observations into the full-length missing vector
ncap2 -O -s "
    ${qobs_var}(${start_index}:${end_index}) =
        ${qobs_var}_insert;
" \
    "$merged_work" \
    "$merged_work"

echo Inserted observations into $merged_work

# ----------------------------------------------------------------------
# Remove the temporary observation variable and coordinate
# ----------------------------------------------------------------------

ncks -O -x \
    -v "${qobs_var}_insert,obs_time" \
    "$merged_work" \
    "$merged_file"

echo "Cleaned up merged file: $merged_file"

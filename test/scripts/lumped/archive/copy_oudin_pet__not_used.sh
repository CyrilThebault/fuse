#!/bin/bash

set -euo pipefail

oudin_file=$1
output_file=$2

work_dir="$(dirname "$output_file")/work"
mkdir -p "$work_dir"

oudin_insert_file="${work_dir}/oudin_pet_insert.nc"

echo "  Oudin file : ${oudin_file}"
echo "  Output file: ${output_file}"
echo "  Work dir   : ${work_dir}"
echo "  Temp file  : ${oudin_insert_file}"
echo

# ------------------------------------------------------------
# Standardize time in the Oudin file
# ------------------------------------------------------------

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
TIME_SCRIPT="${SCRIPT_DIR}/standardize_time_coordinates.sh"

STANDARD_OUDIN="${work_dir}/$(basename "${oudin_file}")"

"${TIME_SCRIPT}" \
    legacy \
    "${oudin_file}" \
    "${STANDARD_OUDIN}"

# ------------------------------------------------------------
# Determine the time ranges
# ------------------------------------------------------------

start_met=$(ncks -H -C -s '%.17g\n' \
    -v time_bnds -d time,0,0 -d nbnds,0,0 "$output_file")

end_met=$(ncks -H -C -s '%.17g\n' \
    -v time_bnds -d time,-1,-1 -d nbnds,1,1 "$output_file")

start_oudin=$(ncks -H -C -s '%.17g\n' \
    -v time_bnds -d time,0,0 -d nbnds,0,0 "$STANDARD_OUDIN")

end_oudin=$(ncks -H -C -s '%.17g\n' \
    -v time_bnds -d time,-1,-1 -d nbnds,1,1 "$STANDARD_OUDIN")

start_time=$(
    printf "%s\n%s\n" "$start_met" "$start_oudin" |
        sort -g |
        tail -1
)

end_time=$(
    printf "%s\n%s\n" "$end_met" "$end_oudin" |
        sort -g |
        head -1
)

echo
echo "Time ranges"
printf "  Merged file : start = %.17g, end = %.17g\n" \
    "$start_met" "$end_met"
printf "  Oudin file  : start = %.17g, end = %.17g\n" \
    "$start_oudin" "$end_oudin"
printf "  Overlap     : start = %.17g, end = %.17g\n" \
    "$start_time" "$end_time"

# ------------------------------------------------------------
# Convert overlap to array indices
# ------------------------------------------------------------

met_start_index=$(
    awk -v start="$start_time" -v origin="$start_met" \
        'BEGIN {printf "%.0f\n", start-origin}'
)

met_end_index=$(
    awk -v end="$end_time" -v origin="$start_met" \
        'BEGIN {printf "%.0f\n", end-origin-1}'
)

oudin_start_index=$(
    awk -v start="$start_time" -v origin="$start_oudin" \
        'BEGIN {printf "%.0f\n", start-origin}'
)

oudin_end_index=$(
    awk -v end="$end_time" -v origin="$start_oudin" \
        'BEGIN {printf "%.0f\n", end-origin-1}'
)

echo
echo "Array indices"
printf "  Merged file : %6d to %-6d (%6d records)\n" \
    "$met_start_index" \
    "$met_end_index" \
    "$((met_end_index - met_start_index + 1))"

printf "  Oudin file  : %6d to %-6d (%6d records)\n" \
    "$oudin_start_index" \
    "$oudin_end_index" \
    "$((oudin_end_index - oudin_start_index + 1))"
echo

# ------------------------------------------------------------
# Extract Oudin PET over the common period
# ------------------------------------------------------------

ncks -O -h -C \
    -v time,pet \
    -d time,"${oudin_start_index}","${oudin_end_index}" \
    "$STANDARD_OUDIN" \
    "$oudin_insert_file"
echo Created $oudin_insert_file

# Give temporary names so it can coexist with the destination time axis
ncrename -O -h \
    -d time,oudin_time \
    -v time,oudin_time \
    -v pet,pet_oudin_insert \
    "$oudin_insert_file"
echo Renamed variables in $oudin_insert_file

# Remove singleton latitude and longitude dimensions
ncwa -O \
    -a latitude,longitude \
    -v pet_oudin_insert \
    "$oudin_insert_file" \
    "$oudin_insert_file"
echo "Removed singleton spatial dimensions in ${oudin_insert_file}"

# fix the record dimension (oudin time) to enable merging time with other files
ncks -O --fix_rec_dmn oudin_time \
    "$oudin_insert_file" \
    "$oudin_insert_file"
echo Fixed time dimension in $oudin_insert_file

ncks -A -h \
    "$oudin_insert_file" \
    "$output_file"
echo Added $oudin_insert_file to $output_file

# ------------------------------------------------------------
# Create full-length Oudin PET and insert available values
# ------------------------------------------------------------

ncap2 -O -h -s "
    pet_oudin[time,hru] = float(-9999.0);
    pet_oudin.set_miss(-9999.0);

    pet_oudin(${met_start_index}:${met_end_index},0) =
        pet_oudin_insert(:);
" \
    "$output_file" \
    "$output_file"

# Remove temporary variables
ncks -O -h -x \
    -v pet_oudin_insert,oudin_time \
    "$output_file" \
    "$output_file"

ncatted -O -h \
    -a long_name,pet_oudin,o,c,"potential evapotranspiration estimated using the Oudin method" \
    -a units,pet_oudin,o,c,"mm/day" \
    "$output_file"

rm -f "$oudin_insert_file"

#!/bin/bash
#
# compute_mean_temperature
#
# Compute daily mean air temperature from daily minimum and maximum
# temperature:
#
#   temp = 0.5 * (tmin + tmax)
#
# Usage:
#
#   ./compute_mean_temperature.sh <input_file> <output_file>
#

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <input_file> <output_file>" >&2
    exit 1
fi

input_file=$1
output_file=$2

if ! command -v ncap2 >/dev/null 2>&1; then
    echo "ERROR: Required command not found: ncap2" >&2
    exit 1
fi

if [[ ! -f "${input_file}" ]]; then
    echo "ERROR: Input file not found: ${input_file}" >&2
    exit 1
fi

if [[ "${input_file}" == "${output_file}" ]]; then
    echo "ERROR: Input and output files must be different." >&2
    exit 1
fi

mkdir -p "$(dirname "${output_file}")"

echo
echo "Creating daily mean temperature"
echo "  Input : ${input_file}"
echo "  Output: ${output_file}"

ncap2 -O -s '
    temp[time,hru] = 0.5f*(tmin+tmax);
    temp@long_name = "daily mean air temperature";
    temp@units = "degrees C";
    temp@_FillValue = -9999.0f;
' \
    "${input_file}" \
    "${output_file}"

echo
echo "Mean-temperature preprocessing completed successfully."
echo "Created: ${output_file}"

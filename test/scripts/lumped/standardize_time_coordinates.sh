#!/bin/bash
#
# standardize_time_coordinates.sh
#
# Purpose
# -------
# Standardize the time coordinates in CAMELS-SPAT forcing,
# streamflow-observation, or legacy FUSE input files.
#
# Usage
# -----
#
#   ./standardize_time_coordinates.sh \
#       <forcing|qobs|legacy> \
#       <input_file> \
#       <output_file>
#
# Examples
# --------
#
#   ./standardize_time_coordinates.sh \
#       forcing \
#       ../CAN_05BB001/input/forcing/CAN_05BB001_daymet_lumped.nc \
#       ../CAN_05BB001/input/work/time/CAN_05BB001_daymet_lumped.nc
#
#   ./standardize_time_coordinates.sh \
#       qobs \
#       ../CAN_05BB001/input/q_obs/CAN_05BB001_daily_flow_observations.nc \
#       ../CAN_05BB001/input/work/time/CAN_05BB001_daily_flow_observations.nc
#
# Output convention
# -----------------
#
# The processed file contains:
#
#   double time(time)
#   double time_bnds(time,nbnds)
#
# with:
#
#   units    = "days since 1950-01-01 00:00:00"
#   calendar = "proleptic_gregorian"
#
# The time coordinate is reconstructed as the midpoint of time_bnds.
#
# Requirements
# ------------
#
# NCO utilities:
#
#   ncap2
#   ncatted
#   ncpdq
#   ncrename
#   ncks
#

set -euo pipefail


# ======================================================================
# Load shared function
# ======================================================================

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${script_dir}/common_time_coordinates.sh"


# ======================================================================
# Command-line arguments
# ======================================================================

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <forcing|qobs|legacy> <input_file> <output_file>" >&2
    exit 1
fi

file_type=$1
input_file=$2
output_file=$3

# ======================================================================
# Validate file type
# ======================================================================

case "${file_type}" in
    forcing|qobs|legacy)
        ;;
    *)
        echo "ERROR: File type must be 'forcing', 'qobs', or 'legacy'." >&2
        exit 1
        ;;
esac

# ======================================================================
# Check requirements
# ======================================================================

for command_name in ncap2 ncatted ncpdq ncrename ncks; do
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
# Prepare input-specific coordinate structure
# ======================================================================

case "${file_type}" in

    qobs)

        echo
        echo "Preparing streamflow observations"
        echo "  Input : ${input_file}"
        echo "  Output: ${output_file}"

        # Original observation structure:
        #
        #   time_bnds(nbnds,time)
        #
        # Required structure:
        #
        #   time_bnds(time,nbnds)
        #
        ncpdq -O \
            -a time,nbnds \
            "${input_file}" \
            "${output_file}"

        description="streamflow observations"
        units_per_day=1440
        ;;

    forcing)

        echo
        echo "Preparing meteorological forcing"
        echo "  Input : ${input_file}"
        echo "  Output: ${output_file}"

        # Work on a copy so the original forcing file is unchanged.
        cp "${input_file}" "${output_file}"

        # Original forcing structure:
        #
        #   time_bnds(time,nv)
        #
        # Required structure:
        #
        #   time_bnds(time,nbnds)
        #
        ncrename \
            -d nv,nbnds \
            "${output_file}"

        description="meteorological forcing"
        units_per_day=1
        ;;

    legacy)

        echo
        echo "Preparing legacy forcing"
        echo "  Input : ${input_file}"
        echo "  Output: ${output_file}"
        
        cp "${input_file}" "${output_file}"
        
        # Original convention:
        #
        #   time = 0, 1, 2, ...
        #   units = "days since 1987-01-01"
        #
        # Convert the numerical values to days since 1950-01-01.
        #
        # Number of days from 1950-01-01 to 1987-01-01:
        time_offset=13514
        
        echo "  Time offset: ${time_offset} days"
        
        ncap2 -O -s "
            defdim(\"nbnds\",2);
        
            time[time] = double(time) + ${time_offset}.0;
        
            time_bnds[time,nbnds] = 0.0;
            time_bnds(:,0) = time;
            time_bnds(:,1) = time + 1.0;
        " \
            "${output_file}" \
            "${output_file}"
        
        description="legacy forcing"
        units_per_day=1
        ;;

esac


# ======================================================================
# Apply common time-coordinate processing
# ======================================================================

common_time_coordinates \
    "${description}" \
    "${output_file}" \
    "${units_per_day}"


# ======================================================================
# Report result
# ======================================================================

echo
echo "Time coordinates in:"
echo "  ${output_file}"

ncks -m \
    -v time,time_bnds \
    "${output_file}"

echo
echo "Time-coordinate preprocessing completed successfully."

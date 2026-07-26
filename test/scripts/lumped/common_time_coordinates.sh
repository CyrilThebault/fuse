#!/bin/bash
#
# ----------------------------------------------------------------------
# Common time-coordinate processing
# ----------------------------------------------------------------------
#
# Arguments:
#
#   $1  Description used in status messages
#   $2  Processed working file
#   $3  Number of original time units per day
#
# The input file must already contain:
#
#   time(time)
#   time_bnds(time,nbnds)
#
# The routine:
#
#   1. Promotes time and time_bnds to double precision.
#   2. Removes missing-value attributes from the coordinates.
#   3. Converts time_bnds to days since 1950-01-01.
#   4. Reconstructs time as the midpoint of time_bnds.
#   5. Applies common coordinate metadata.
#
# Typical values:
#
#   observations: units_per_day = 1440
#                 because the original bounds are in minutes
#
#   forcing:      units_per_day = 1
#                 because the original bounds are already in days
#

common_time_coordinates()
{
    description=$1
    work_file=$2
    units_per_day=$3

    echo
    echo "Standardizing ${description}"
    echo "  File: ${work_file}"

    # Promote both coordinate variables to double precision.
    #
    # Left-hand casting changes the stored NetCDF variable type rather
    # than merely casting values during evaluation of an expression.
    ncap2 -O -s "
        time[time]            = double(time);
        time_bnds[time,nbnds] = double(time_bnds);
    " \
        "${work_file}" \
        "${work_file}"

    # Coordinate values are expected to be complete. Remove _FillValue
    # and missing_value attributes before performing arithmetic.
    #
    # This also avoids warnings caused by comparisons between NaN-valued
    # missing-value attributes.
    ncatted -O \
        -a _FillValue,time,d,, \
        -a missing_value,time,d,, \
        -a _FillValue,time_bnds,d,, \
        -a missing_value,time_bnds,d,, \
        "${work_file}"

    # Convert the interval bounds to days.
    #
    # For observations:
    #
    #     units_per_day = 1440 minutes/day
    #
    # For forcing:
    #
    #     units_per_day = 1 day/day
    #
    # The time coordinate is then reconstructed from the authoritative
    # interval bounds:
    #
    #     time(t) =
    #         0.5 * [time_bnds(t,0) + time_bnds(t,1)]
    #
    ncap2 -O -s "
        time_bnds = time_bnds/${units_per_day}.0;
        time      = 0.5*(time_bnds(:,0)+time_bnds(:,1));
    " \
        "${work_file}" \
        "${work_file}"

    # Apply common units, calendar, and descriptive metadata.
    ncatted -O \
        -a units,time,o,c,"days since 1950-01-01 00:00:00" \
        -a calendar,time,o,c,"proleptic_gregorian" \
        -a axis,time,o,c,"T" \
        -a bounds,time,o,c,"time_bnds" \
        -a long_name,time,o,c,"midpoint of the daily interval in local standard time" \
        -a time_basis,time,o,c,"local standard time" \
        -a units,time_bnds,o,c,"days since 1950-01-01 00:00:00" \
        -a calendar,time_bnds,o,c,"proleptic_gregorian" \
        -a long_name,time_bnds,o,c,"daily interval bounds in local standard time" \
        -a time_basis,time_bnds,o,c,"local standard time" \
        "${work_file}"

}

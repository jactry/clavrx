#!/bin/ksh
#   
# shift_date: Shift specified date (YYYYMMDD) by number of days
# e.g. shift_date 20000411 -42
#      = 20000229
# Uses internal Julian date based on nearest convenient year
# (calculated each time based on year and number of days specified).
# Will therefore accommodate any date within reason.
# No it doesn't fully handle the Gregorian adjustment of 1752
# (Britain and its then colonies adopted the Gregorian calendar in 1752 when September 2nd was followed by
# September 14th. Catholic countries had already changed over in 1582. For more information, see
# http://aa.usno.navy.mil/data/docs/JulianDate.html)
#
# Parameters:
# $1   shift_date   Date (YYYYMMDD)
# $2   shift_days   Days to add (integer, can be negative)
#
# Set environment variable DEBUG=TRUE to see diagnostic mesages during execution.
# William Robertson www.williamrobertson.net
# 21 Feb 2004 AD

# 19 Jan 2005 AD: Defining variables and assigning them values from funtion expressions in the same line
# e.g. typeset -Z2 var=$(func)
# does not seem to work on some platforms including Solaris and Sequent Dynix.

# 14 Apr 2007 AD: Integer division fails in (new?) version of ksh93 on CentOS5 when divisor is negative
# e.g. (( -4 / -2 )) gives 2 in other versions, fails with zero divide error on CentOS 5.
# Fixed by reversing both signs (thanks Loren Siebert)

# Define current script name, stripping any leading "*/" or "-" (in case run in
# current shell for some reason and $0 = "-sh")
program=${0##*(*/|-)}

# Allow references to unset variables:
set +u
alias integer="typeset -i"

# Define 'usage' text for error messages:
usage="Usage: ${program} [date YYYYMMDD] [days to add]
\n\tparameter 1 = date in YYYYMMDD format
\n\tparameter 2 = days to add (negative to subtract), integer
\ne.g:\n\t$program $(date +'%Y%m%d') -42"

# Array of month lengths for use by conversion functions
# (Feb 28 is default value only - gets reset to 29 when necessary)
set -A monthlengths dummy 31 28 31 30 31 30 31 31 30 31 30 31

# Set DEBUG=TRUE to enable diagnostic output.
typeset -ux DEBUG=${DEBUG:-FALSE}

function debug
{
	# Optional diagnostic message, only displayed if environment variable DEBUG=TRUE

	message="$*"
	[[ ${DEBUG} = TRUE ]] && print -u2 "${message#-}"
}

debug "${program}: Debugging messages enabled. Set environment variable DEBUG=FALSE to disable."

function bomb
{
	# Display error message and exit with non-zero exit status.

	# Pick up exit status from environment (reset to 2 if zero):
	exitstatus=$?
	(( exitstatus == 0 )) && exitstatus=2

	# Display message only if one was passed as parameter:
	if (( $# > 0 )); then
		print "${program}: $*" >&2
	fi

	debug Exiting with exit status ${exitstatus}
	exit ${exitstatus}
}

# Ensure two parameters were supplied:
(( $# == 2 )) || bomb ${usage}

# Ensure parameter 1 is a number and is 8 digits long:
shift_date=$1 && [[ ${#1} = 8 ]] && print ${shift_date} | grep >/dev/null "[0-9]\{8\}" || bomb "Invalid date '$1'\n${usage}"
[[ ${shift_date} != $1 ]] && bomb "Invalid date '$1'\n${usage}"

# Ensure parameter 2 is a number (note shell rounds decimals down):
# Some shells give fatal error of integer assignment fails, others silently set value to 0,
# so to validate $2 we assign to an integer variable and then do a text comparison (stripping
# any leading "+" off start of $2).
integer shift_days=$2 || bomb "Invalid number of days $2\n${usage}"
[[ ${shift_days} != ${2#+} ]] && bomb "Invalid number '$2'\n${usage}"

# Rather than hardcode arbitrary julian base year e.g. 1970, set each time
# using specified year and number of shift_days:

typeset -L4 base_year=${shift_date}  # Extract as leftmost 4 characters of $shiftdate

if (( shift_days != 1 )); then
	plural="s"
fi

if (( shift_days >= 0 ))
then
	# Positive shift: don't need base date in the past; set to year specified:
	# (redefine existing character variable as integer)
	sign="+"
	integer max_years=1
	integer base_year

	debug "Using ${base_year} as base year for julian date conversion."

	debug Adding ${shift_days} day${plural} to ${shift_date}
else
	# Negative shift: set julian base date as far in the past as it needs to be:
	sign=""
	# integer max_years=$(( (shift_days/-365) +1 ))
	# Latest version of ksh93 in CentOS 5 gives zero divide error when multiplying by a negative number,
	# so reverse both signs to give same result:
	integer max_years=$(( ( (shift_days * -1) / 365) +1 ))
	integer base_year=$(( base_year - max_years ))  # Redefine as integer

	debug "Negative shift ${shift_days} days = ${max_years} year(s) or less (${shift_days}/-365 +1):"
	debug "Using ${base_year} as base year for julian date conversion."

	debug Subtracting $(( shift_days * -1 )) day${plural} from ${shift_date}
fi

# Make read-only (some shells don't allow typeset -r var=value):
typeset -r max_years
typeset -r base_year

# Sanity check values derived above:
if (( base_year < 1 ))  # There was no year zero
then
	bomb "Cannot calculate BC dates."
fi

function leapyear
{
	# Test for a leap year.
	# Parameters:
	# $1 testyear YYYY
	# Returns actual boolean (not text) via return exit status.
	# Leap Year occurs every four years, except for years ending in 00, in which case
	# only if the year is divisible by 400.
	# ...however, this system was only standardised in 1752 (1582 in Catholic countries), so the 400 rule
	# does not apply to years before then.

	# The actual length of a year is 365 days, five hours, 48 minutes, 46 seconds.

	integer testyear=$1 TRUE=0 FALSE=1

	# Initially assume FALSE, only reset to TRUE if passes tests:
	integer leapyearanswer=${FALSE}

	if (( testyear % 4 == 0 ))
	then
		if (( testyear % 100 == 0 ))
		then
			if (( testyear < 1752 ))
			then
				leapyearanswer=${TRUE}

			elif (( testyear % 400 == 0 ))
			then
				leapyearanswer=${TRUE}
			fi
		else
			leapyearanswer=${TRUE}
		fi
	fi

	return ${leapyearanswer}
}

function yearlength
{
	# Return number of days in specified year:
	# Parameters:
	# $1 testyear YYYY

	integer result=365

	integer testyear=$1 && (( testyear > 0 )) && (( testyear < 10000 )) \
		|| bomb "$0: Invalid year '$1': year must be between 1 and 9999."

	if leapyear ${testyear}
	then
		(( result = 366 ))

	elif (( testyear == 1582 )) || (( testyear == 1752 ))
	then
		(( result = 355 ))
	fi

	print ${result}
}

function split_yyyy_mm_dd
{
	# Accept YYYYMMDD
	# Return YYYY MM DD, validated as integers in vaguely appropriate ranges.
	date_yyyymmdd=$1

	# Confirm length of date parameter is 8:
	[[ ${#date_yyyymmdd} = 8 ]] || bomb "$0: Invalid date '$1': must be 8 digits"

	# Break date string into YYYY MM DD elements:
	typeset -L4 yyyy=${date_yyyymmdd}
	typeset -R2 dd=${date_yyyymmdd}
	typeset -R4 mm=${date_yyyymmdd}
	typeset -L2 mm=${mm}

	debug $0: ${date_yyyymmdd} = ${yyyy} ${mm} ${dd}

	integer yyyy || bomb "$0: year '${year}' is not numeric."

	# Confirm yyyy is a number between 1 and 9999:
	if (( yyyy < 1 )) || (( yyyy > 9999 ))
	then
		bomb "$0: Invalid year '${yyyy}': must be between 1 and 9999"
	fi

	print ${yyyy} ${mm} ${dd}
}

function to_j
{
	# Parameters:
	# $1: date_yyyymmdd: date in format YYYYMMDD
	# Return date converted to julian format (days since 01/01/${base_year})

	[[ $1 = [0-9][0-9][0-9][0-9] ]] || bomb $1 is not 4 digits
	[[ $2 = [0-1][0-9] ]]           || bomb $2 is not 2 digits
	[[ $3 = [0-3][0-9] ]]           || bomb $3 is not 2 digits [0-3][0-9]

	integer yyyy=$1 || bomb "$0: Year $1 is not numeric"
	integer mm=$2 || bomb "$0: Month $2 is not numeric"
	integer dd=$3 || bomb "$0: Day $1 is not numeric"

	debug yyyy = $1 mm = $2 dd = $3

	# Confirm yyyy is a number between 1 and 9999:
	if (( yyyy < 1 )) || (( yyyy > 9999 ))
	then
		bomb "$0: Invalid year '${yyyy}': must be between 1 and 9999"
	fi

	# If yyyy is a leap year, set Feb length to 29:
	if leapyear ${yyyy}
	then
		monthlengths[2]=29
	else
		monthlengths[2]=28
	fi

	# Confirm mm is a number between 1 and 12:
	integer mm && (( 1 <= mm )) && (( mm <= 12 )) \
		|| bomb "$0: Invalid month '${mm}'"

	# Confirm dd is a number between 1 and length of month 'mm':
	integer dd && (( 1 <= dd )) && (( dd <= ${monthlengths[mm]} )) \
		|| bomb "$0: Invalid day '${dd}' for month ${mm}/${yyyy}"

	integer yearlength y=${base_year} j=0 m=1
	yearlength=$(yearlength ${y}) || bomb  # Cascade on failure

	debug "$0: Convert ${yyyy} ${mm} ${dd} to julian date, base = 01/01/${base_year}:"

	# Add whole years to total 'j', stopping at year before yyyy:
	# (y is counter starting at ${base_year}; yyyy is the year specified in date parameter)

	while (( y < yyyy ))
	do
		yearlength=$(yearlength ${y}) || bomb  # Cascade on failure

		debug "\tYear ${y}: incrementing j by ${yearlength}"

		(( j += yearlength ))
		(( y += 1 ))
	done

	debug "\tFinished years at ${y}: j = ${j}."

	# Set length of February (reset yearlength as loop may have overshot in last iteration):
	yearlength=$(yearlength ${y}) || bomb  # Cascade on failure
	if (( yearlength == 366 ))
	then
		(( monthlengths[2] = 29 ))
	else
		(( monthlengths[2] = 28 ))
	fi

	# Add whole months to total 'j', stopping at month before mm:
	# (m is counter starting at 1; mm is the month specified in date parameter)
	while (( m < mm ))
	do
		debug "\tMonth ${m}: incrementing j by ${monthlengths[m]}"

		(( j += ${monthlengths[m]} ))
		(( m += 1 ))
	done

	debug "\tFinished months at ${m}: j = ${j}."

	# Add days to total 'j':
	debug "\tAdding ${dd} days"
	(( j += dd ))

	debug "$0 $* = ${j}"

	if (( j < 1 ))
	then
		bomb "$0: Date conversion ${yyyy} ${mm} ${dd} to julian with" \
		"base year ${base_year} gave invalid result '${j}'."
	else
		print ${j}
	fi
}

function j_to_yyyymmdd
{
	# Parameters:
	# $1: date_j: date in julian format (days since 01/01/${base_year})
	# Return date converted to format YYYYMMDD

	integer date_j=$1 || bomb "$0: Invalid julian date '$1'"
	(( date_j > 0 )) || bomb "$0: julian date must be greater than 0 (01/01/${base_year})"

	integer yyyy=${base_year} mm=1 dd=0 yearlength monthlength

	debug "$0: Convert ${date_j} to YYYYMMDD using base date 01/01/${base_year}:"

	# Add whole years to total 'j', stopping when running total goes < 0:
	while (( date_j > 0 ))
	do
		if (( yyyy < 1 ))
		then
			bomb "Cannot calculate BC dates."
		elif (( yyyy > 9999 ))
		then
			bomb "Invalid year ${yyyy}: cannot be greater than 9999."
		fi

		(( yearlength = $(yearlength ${yyyy}) )) || bomb  # Cascade on failure

		if (( date_j <= yearlength ))
		then
			# Less than a whole year to add:
			# use yearlength to set length of Feb for use in Month stage, and break:
			debug "\t${date_j} <= ${yearlength}: Finished years: yyyy = ${yyyy}."
			break
		else
			(( yyyy += 1 ))
			(( date_j -= yearlength ))

			debug "\tYear $(( yyyy -1 )):" \
				"Reduced date_j by ${yearlength} to ${date_j}"
		fi
	done

	# Re-use final yearlength from loop to set length of February:
	if (( yearlength == 366 ))
	then
		(( monthlengths[2] = 29 ))
	else
		(( monthlengths[2] = 28 ))
	fi

	# Add whole months to total 'j', stopping at month before mm:
	while (( mm <= 12 )) && (( date_j > 0 ))
	do
		(( monthlength = ${monthlengths[mm]} ))

		if (( date_j <= monthlength ))
		then
			debug "\tMonth ${mm}:" \
				"${date_j} <= ${monthlength}: Finished months: mm = ${mm}."
			break
		else
			(( mm += 1 ))
			(( date_j -= monthlength ))

			debug "\tMonth $(( mm -1 )) length ${monthlength}:" \
				"Reduced date_j by ${monthlength} to ${date_j}"
		fi
	done

	# Assign remaining days to DD:
	dd=${date_j}
	debug dd = $dd

	# Finished integer arithmetic: change to zero-padded strings:
	typeset -Z4 yyyy
	typeset -Z2 mm dd  # left-padded with zeroes, 2 digits

	debug "$0 $* = ${yyyy} ${mm} ${dd}"

	# Sanity-check result before returning:
	if (( yyyy < 1 ))
	then
		bomb "$0: Date conversion ${date_j} to YYYYMMDD base year ${base_year}" \
			"gave invalid result '${yyyy}'."
	else
		print ${yyyy}${mm}${dd}
	fi
}

# Convert to julian date and confirm valid result:
# integer jdate=$(to_j ${shift_date}) || bomb
set -A date_elements $(split_yyyy_mm_dd ${shift_date}) || bomb "Failed to split ${shift_date}) into its components"

# Inelegant optimization:
# Short-circuit the whole procedure in simplest case where dd + shift_days > 0 < month_length
# e.g. if adding 2 days to 2004 03 02, answer will be 2004 03 04 (month unchanged)
integer test_result_dd
typeset -Z2 mm monthlength
mm=${date_elements[1]}
monthlength=${monthlengths[${mm}]}

(( test_result_dd = date_elements[2] + shift_days ))

debug Day ${date_elements[2]} + ${shift_days} = ${test_result_dd}
debug Date elements are: ${date_elements[*]}
debug Month lengths are: ${monthlengths[*]}
debug Length of month ${mm} is ${monthlengths[${mm}]} which should be the same as ${monthlength} 
debug Testing whether ${test_result_dd} is less than month length ${monthlength} \
	and required day ${test_result_dd} is in future

if (( test_result_dd < monthlength )) && (( test_result_dd > 0 )); then
	debug "${test_result_dd} is less than ${monthlength} and ${test_result_dd} is greater than zero"
	debug "Short circuit: day ${date_elements[2]} + ${shift_days} = ${test_result_dd}" \
		"is in range 1-${monthlength} for month ${mm}"
	debug Quick result: ${date_elements[0]} ${mm} ${test_result_dd}

	# Declared test_result_dd as an integer above for efficiency,
	# but make it 2-char zero-padded here immediately prior to display:
	typeset -Z2 test_result_dd
	print ${date_elements[0]}${mm}${test_result_dd}
	exit 0
else
	debug Short circuit not possible: applying full procedure.
fi

debug Defining jdate from by calling function to_j "${date_elements[*]}"
# Solaris reports success indiscriminately if var defined and set on same line e.g. integer jdate=$(func)
integer jdate
jdate=$(to_j ${date_elements[*]}) || bomb
if [[ $? = 0 ]]; then
	debug Call successful
else
	debug Call failed
fi

# Apply shift_days (+ or -) to jdate, and confirm valid result:
(( jdate += shift_days )) || bomb "Date arithmetic failed (${jdate} - ${shift_days})"

debug "$((jdate - shift_days)) ${sign}${shift_days} = ${jdate}"

# Convert back:
new_date_yyyymmdd=$(j_to_yyyymmdd ${jdate}) || bomb

# Confirm we have a result before returning:
if [[ -z ${new_date_yyyymmdd} ]]
then
	print "$0: Failed: $(date)" >&2
	exit 1
else
	print ${new_date_yyyymmdd}
fi

 

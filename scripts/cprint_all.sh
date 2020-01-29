#!/bin/bash
#
# Usage:
# ./script.sh [J_min] [J_step] [J_max] [arrang] [shift] [scale] [adiabatic] [cprint.out]
#

set -e
set -u

input_filename="input.d"
cmatrix_datafile="cmatrix_arrang=%c_n=%d_J=%d.bin"

assert_file ()
{
	if [ ! -e $1 ]
	then
		echo
		echo "$0, error: $1 not found"
		echo
		exit 666
	fi
}

assert_file $8
start_dir=$PWD

for J in $(seq $1 $2 $3)
do
	work_dir="$start_dir/J=$J"
	assert_file $work_dir

	for parity in $(seq -1 2 +1)
	do
		if [ "$parity" == "-1" ]
		then
			bin_dir="$work_dir/parity=-1/bin"
			assert_file $bin_dir

			channels_dir="$work_dir/parity=-1/channels"
			assert_file $channels_dir
		else
			bin_dir="$work_dir/parity=+1/bin"
			assert_file $bin_dir

			channels_dir="$work_dir/parity=+1/channels"
			assert_file $channels_dir
		fi

		echo ""
		cd $channels_dir

		ln -s "$bin_dir/$cmatrix_datafile" .

		echo "J = $J"                > $input_filename
		echo "arrang = $4"          >> $input_filename

		echo ""
		echo "energy_shift = $5"    >> $input_filename
		echo "energy_scale = $6"    >> $input_filename
		echo "print_adiabatic = $7" >> $input_filename

		$8 $input_filename
		rm -rf $cmatrix_datafile

		cd $start_dir
	done
done

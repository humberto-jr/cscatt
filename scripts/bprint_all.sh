#!/bin/bash
#
# Usage:
# ./script.sh [J_min] [J_step] [J_max] [bprint.out] [arrang]
#

set -e
set -u

input_filename="input.d"
bprint_datafile="basis_arrang=*_ch=*_J=*.bin"

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

assert_file $4
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

			basis_wavef_dir="$work_dir/parity=-1/basis_wavef"
			assert_file $basis_wavef_dir
		else
			bin_dir="$work_dir/parity=+1/bin"
			assert_file $bin_dir

			basis_wavef_dir="$work_dir/parity=+1/basis_wavef"
			assert_file $basis_wavef_dir
		fi

		echo ""
		cd $basis_wavef_dir

		ln -s "$bin_dir/$bprint_datafile" .

		echo "J = $J"       > $input_filename
		echo "arrang = $5" >> $input_filename

		$4 $input_filename
		rm -rf $bprint_datafile

		cd $start_dir
	done
done

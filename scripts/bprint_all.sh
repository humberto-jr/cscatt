#!/bin/bash
#
# Usage:
# ./script.sh [J_min] [J_step] [J_max] [bprint.out]
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
			if [ "$J" == "0" ]
			then
				continue
			fi

			parity_dir="$work_dir/parity=-1"
		else
			parity_dir="$work_dir/parity=+1"
		fi

		assert_file $parity_dir

		bin_dir="$parity_dir/bin"
		assert_file $bin_dir

		basis_wavef_dir="$parity_dir/basis_wavef"
		assert_file $basis_wavef_dir

		input="$parity_dir/$input_filename"
		assert_file $input

		echo "#"
		echo "# $basis_wavef_dir"

		cd $basis_wavef_dir

		rm -rf $bprint_datafile
		ln -s $bin_dir/$bprint_datafile .

		cp $input .

		$4 $input_filename
		rm -rf $bprint_datafile

		cd $start_dir
	done
done

#!/bin/bash
#
# Usage:
# ./script.sh [J_min] [J_step] [J_max] [shift] [scale] [adiabatic] [cprint.out]
#

set -e
set -u

input_filename="input.d"
cmatrix_datafile="cmatrix_arrang=*_n=*_J=*.bin"

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

assert_file $7
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

		channels_dir="$parity_dir/channels"
		assert_file $channels

		input="$parity_dir/$input_filename"
		assert_file $input

		echo "#"
		echo "# $channels_dir"

		cd $channels_dir

		rm -rf $cmatrix_datafile
		ln -s $bin_dir/$cmatrix_datafile .

		cp $input .

		shift_pattern=$(grep "energy_shift = " $input_filename)
		sed -i "s/$shift_pattern/energy_shift = $4/g" $input_filename

		scale_pattern=$(grep "energy_scale = " $input_filename)
		sed -i "s/$scale_pattern/energy_scale = $5/g" $input_filename

		print_pattern=$(grep "print_adiabatic = " $input_filename)
		sed -i "s/$print_pattern/print_adianatic = $6/g" $input_filename

		$7 $input_filename
		rm -rf $cmatrix_datafile

		cd $start_dir
	done
done

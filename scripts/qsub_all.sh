#!/bin/bash
set -e
set -u

start_dir=$PWD

for J in $(seq $1 $2 $3)
do
	for parity in $(seq -1 2 +1)
	do
		if [ "$parity" == "-1" ]
		then
			dir="J=$J/parity=-1"
			pbs_filename="$dir/job.pbs"
		else
			dir="J=$J/parity=+1"
			pbs_filename="$dir/job.pbs"
		fi

		if [ -e $pbs_filename ]
		then
			echo ""
			echo $pbs_filename

			cd $dir
			qsub job.pbs
			cd $start_dir
		fi
	done
done

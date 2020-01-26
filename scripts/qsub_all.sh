#!/bin/bash
set -e
set -u

for J in $(seq $1 $2 $3)
do
	for parity in $(seq -1 2 +1)
	do
		if [ "$parity" == "-1" ]
		then
			pbs_filename="J=$J/parity=-1/job.pbs"
		else
			pbs_filename="J=$J/parity=+1/job.pbs"
		fi

		if [ -e $pbs_filename ]
		then
			echo
			echo $pbs_filename
			qsub $pbs_filename
		fi
	done
done

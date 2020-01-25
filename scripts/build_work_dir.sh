#!/bin/bash
set -u
set -e

J_min=0
J_max=5
J_step=1

arrang=1

atom_a="35Cl"
atom_b="1H"
atom_c="1H"

J_min=0
J_max=5
J_step=1

v_min=0
v_max=2
v_step=1

j_min=0
j_max=2
j_step=2

lambda_max=4

r_min=0.5
r_max=5.0
rovib_grid_size=200

R_min=0.5
R_max=5.0
scatt_grid_size=200

use_omp=0

# executables
dbasis_exe="hahaha"
bprint_exe="cucuc"
cmatrix_exe="kakak"
cprint_exe="kikiikik"

# misc
input_filename="input.d"
slurm_filename="job.sh"
bprint_datafile="basis_arrang=*_ch=*_J=*.dat"
cmatrix_datafile="cmatrix_arrang=*_n=*_J=*.dat"

# slurm configuration
wall_time="5:00:00"
max_memory="2G"
nodes=6
mpi_cpus=12
omp_threads=8
modules=""

# OpenMP configuration
env_omp_threads='OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK'

# MPI configuration
mpi_proc_placement='I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0'
mpirun_call='mpirun -np $SLURM_NTASKS -genv OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK -genv I_MPI_PIN_DOMAIN=omp'

# end of inputs ###############################################################

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

build_dir ()
{
	if [ -d $1 ]
	then
		rm -rf "$1/*"
	else
		mkdir $1
	fi
}

assert_file $dbasis_exe

if [ "$bprint_exe"  != "" ]
then
	assert_file $bprint_exe
fi

if [ "$cmatrix_exe" != "" ]
then
	assert_file $cmatrix_exe
fi

if [ "$cprint" != "" ]
then
	assert_file $cprint_exe
fi

for J in $(seq $J_min $J_step $J_max)
do
	work_dir="J=$J"
	build_dir $work_dir

	for parity in $(seq -1 2 +1)
	do
		if [ "$J" == "0" ] && [ "$parity" == "-1" ]
		then
			continue
		else
			if [ "$parity" == "-1" ]
			then
				parity_dir="$work_dir/parity=-1"
				build_dir $parity_dir

				job_name="J=-$J"
			else
				parity_dir="$work_dir/parity=+1"
				build_dir $parity_dir

				job_name="J=+$J"
			fi
		fi

		scatt_wavef_dir="$parity_dir/scatt_wavef"
		build_dir $scatt_wavef_dir

		basis_wavef_dir="$parity_dir/basis_wavef"
		build_dir $basis_wavef_dir

		channels_dir="$parity_dir/channels"
		build_dir $channels_dir

		csection_dir="$parity_dir/csection"
		build_dir $csection_dir

		kmatrix_dir="$parity_dir/kmatrix"
		build_dir $kmatrix_dir

		smatrix_dir="$parity_dir/smatrix"
		build_dir $smatrix_dir

		bin_dir="$parity_dir/bin"
		build_dir $bin_dir

		filename="$parity_dir/$input_filename"
		input=$filename

		echo "arrang = $arrang"                    > $filename

		echo ""                                   >> $filename
		echo "mass_a = $atom_a"                   >> $filename
		echo "mass_b = $atom_b"                   >> $filename
		echo "mass_c = $atom_c"                   >> $filename

		echo ""                                   >> $filename
		echo "J = $J"                             >> $filename
		echo "J_parity = $parity"                 >> $filename

		echo ""                                   >> $filename
		echo "v_min  = $v_min"                    >> $filename
		echo "v_max  = $v_max"                    >> $filename
		echo "v_step = $v_step"                   >> $filename

		echo ""                                   >> $filename
		echo "j_min  = $j_min"                    >> $filename
		echo "j_max  = $j_max"                    >> $filename
		echo "j_step = $j_step"                   >> $filename

		echo ""                                   >> $filename
		echo "lambda_max = $lambda_max"           >> $filename

		echo ""                                   >> $filename
		echo "r_min = $r_min"                     >> $filename
		echo "r_max = $r_max"                     >> $filename
		echo "rovib_grid_size = $rovib_grid_size" >> $filename

		echo ""                                   >> $filename
		echo "R_min = $R_min"                     >> $filename
		echo "R_max = $R_max"                     >> $filename
		echo "scatt_grid_size = $scatt_grid_size" >> $filename

		echo ""                                   >> $filename
		echo "use_omp = $use_omp"                 >> $filename

		filename="$parity_dir/$slurm_filename"

		echo "#!/bin/sh"                                 > $filename
		echo '#SBATCH --job-name="'$job_name'"'         >> $filename

		echo ""                                         >> $filename
		echo "#SBATCH --mem=$max_memory"                >> $filename
		echo "#SBATCH --time=$wall_time"                >> $filename

		echo ""                                         >> $filename
		echo "#SBATCH --nodes=$nodes"                   >> $filename
		echo "#SBATCH --ntasks=$mpi_cpus"               >> $filename
		echo "#SBATCH --cpus-per-task=$omp_threads"     >> $filename

		echo ""                                         >> $filename
		echo "# Job script generated at $(date) by $0"  >> $filename

		echo ""                                         >> $filename
		echo "export $mpi_proc_placement"               >> $filename
		echo "export $env_omp_threads"                  >> $filename

		if [ "$modules" != "" ]
		then
			echo ""                                      >> $filename
			echo "module load $modules"                  >> $filename
		fi

		echo ""                                         >> $filename
		echo 'echo "Calculation starting at $(date)"'   >> $filename

		echo ""                                         >> $filename
		echo "cd $bin_dir/"                             >> $filename

		echo ""                                         >> $filename
		echo "$dbasis_exe $input"                       >> $filename

		if [ "$bprint_exe" != "" ]
		then
			echo ""                                      >> $filename
			echo "$bprint_exe $input"                    >> $filename
			echo "mv $bprint_datafile $basis_wavef_dir/" >> $filename
		fi

		if [ "$cmatrix_exe" != "" ]
		then
			echo ""                                      >> $filename
			echo "$mpirun_call $cmatrix_exe $input"      >> $filename
		fi

		if [ "$cprint_exe" != "" ]
		then
			echo ""                                      >> $filename
			echo "$cprint_exe $input"                    >> $filename
			echo "mv $cmatrix_datafile $channels_dir/"   >> $filename
		fi

		echo ""                                         >> $filename
		echo 'echo "Calculation ending at $(date)"'     >> $filename
	done

	echo $work_dir
done

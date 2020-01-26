#!/bin/bash
set -u
set -e

J_min=0
J_max=10
J_step=1

arrang=1

atom_a="35Cl"
atom_b="1H"
atom_c="1H"

v_min=0
v_max=2
v_step=1

j_min=0
j_max=8
j_step=2

lambda_max=4

r_min="0.5"
r_max="5.0"
rovib_grid_size=500

R_min="2.0"
R_max="500.0"
scatt_grid_size=2000

use_omp=1

energy_shift="0.0"
energy_scale="1.0"
print_adiabatic=0

# executables
dbasis_exe="/home/humberto/H2+Cl_minus/exec/dbasis.out"
bprint_exe=""
cmatrix_exe="/home/humberto/H2+Cl_minus/exec/cmatrix.out"
cprint_exe=""

# misc
input_filename="input.d"
slurm_filename="job.sh"
pbs_filename="job.pbs"
bprint_datafile="basis_arrang=*_ch=*_J=*.dat"
cmatrix_datafile="cmatrix_arrang=*_n=*_J=*.dat"

# batch job configuration
wall_time="12:00:00"
max_memory="2Gb"
nodes=10
mpi_cpus=10
omp_threads=12
queue_name="balalab"
modules="intelmpi/16.1.2"

# OpenMP configuration
env_omp_threads="OMP_NUM_THREADS=$omp_threads"

# MPI configuration
mpi_pin_domain_name="I_MPI_PIN_DOMAIN=omp"
mpi_proc_placement="I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0"

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

build_batch_job ()
{
	if [ -e $1 ]
	then
		echo ""                                         >> $1
		echo "export $env_omp_threads"                  >> $1

		if [ "$mpi_proc_placement" != "" ]
		then
			echo "export $mpi_proc_placement"            >> $1
		fi

		if [ "$modules" != "" ]
		then
			echo ""                                      >> $1
			echo "module load $modules"                  >> $1
		fi

		echo ""                                         >> $1
		echo 'bin_dir="'$bin_dir'"'                     >> $1
		echo 'channels_dir="'$channels_dir'"'           >> $1
		echo 'basis_wavef_dir="'$basis_wavef_dir'"'     >> $1

		echo ""                                         >> $1
		echo 'dbasis_exe="'$dbasis_exe'"'               >> $1
		echo 'bprint_exe="'$bprint_exe'"'               >> $1
		echo 'cmatrix_exe="'$cmatrix_exe'"'             >> $1
		echo 'cprint_exe="'$cprint_exe'"'               >> $1

		echo ""                                         >> $1
		echo 'mpirun_call="'$mpirun_call'"'             >> $1

		echo ""                                         >> $1
		echo 'input="'$input'"'                         >> $1

		echo ""                                         >> $1
		echo 'bprint_datafile="'$bprint_datafile'"'     >> $1
		echo 'cmatrix_datafile="'$cmatrix_datafile'"'   >> $1

		echo ""                                         >> $1
		echo 'echo "Calculation starting at $(date)"'   >> $1

		echo ""                                         >> $1
		echo 'cd $bin_dir/'                             >> $1

		echo ""                                         >> $1
		echo '$dbasis_exe $input'                       >> $1

		if [ "$bprint_exe" != "" ]
		then
			echo ""                                      >> $1
			echo '$bprint_exe $input'                    >> $1
			echo 'mv $bprint_datafile $basis_wavef_dir/' >> $1
		fi

		if [ "$cmatrix_exe" != "" ]
		then
			echo ""                                      >> $1
			echo '$mpirun_call $cmatrix_exe $input'      >> $1
		fi

		if [ "$cprint_exe" != "" ]
		then
			echo ""                                      >> $1
			echo '$cprint_exe $input'                    >> $1
			echo 'mv $cmatrix_datafile $channels_dir/'   >> $1
		fi

		echo ""                                         >> $1
		echo 'echo "Calculation ending at $(date)"'     >> $1
	fi
}

#assert_file $dbasis_exe

#if [ "$bprint_exe"  != "" ]
#then
#	assert_file $bprint_exe
#fi

#if [ "$cmatrix_exe" != "" ]
#then
#	assert_file $cmatrix_exe
#fi

#if [ "$cprint" != "" ]
#then
#	assert_file $cprint_exe
#fi

if [ "$mpi_cpus" == "1" ]
then
	mpirun_call=""
else
	if [ "$mpi_pin_domain_name" == "" ]
	then
		mpirun_call="mpirun -np $mpi_cpus -genv OMP_NUM_THREADS=$omp_threads"
	else
		mpirun_call="mpirun -np $mpi_cpus -genv OMP_NUM_THREADS=$omp_threads -genv $mpi_pin_domain_name"
	fi
fi

#pbs_ncpus=$(echo "($mpi_cpus + $omp_threads)/$nodes" | bc)
#pbs_mpiprocs=$(echo "$pbs_ncpus - $omp_threads" | bc)
#pbs_mpiprocs=$mpi_cpus

for J in $(seq $J_min $J_step $J_max)
do
	work_dir="$PWD/J=$J"
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
		echo "energy_shift = $energy_shift"       >> $filename
		echo "energy_scale = $energy_scale"       >> $filename
		echo "print_adiabatic = $print_adiabatic" >> $filename

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

		build_batch_job $filename

		filename="$parity_dir/$pbs_filename"

		echo "#!/bin/sh"                                                                                     > $filename
		echo "#PBS -l select=$nodes:ncpus=$omp_threads:mpiprocs=1:ompthreads=$omp_threads:mem=$max_memory"  >> $filename
		echo "#PBS -l walltime=$wall_time"                                                                  >> $filename

		if [ $queue_name != "" ]
		then
			echo "#PBS -q $queue_name"                                                                        >> $filename
		fi

		echo '#PBS -N "'$job_name'"'                                                                        >> $filename

		echo ""                                                                                             >> $filename
		echo "# Job script generated at $(date) by $0"                                                      >> $filename

		build_batch_job $filename
	done

	echo "J=$J"
done

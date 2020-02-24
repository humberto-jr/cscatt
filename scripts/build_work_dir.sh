#!/bin/bash
set -u
set -e

# Total angular momentum, J
J_min=0
J_max=10
J_step=1

# Arrangement (1 = a + bc, 2 = b + ac, 3 = c + ab)
arrang=1
atom_a="1H"
atom_b="16O"
atom_c="1H"

# Diatomic vibration quantum number, v
v_min=0
v_max=2
v_step=1

# Diatomic rotation quantum number, j
j_min=0
j_max=8
j_step=1

# Number of atom-diatom potential expansion coefficients
lambda_max=4

# Diatomic radial grid
r_min="0.6"
r_max="6.0"
rovib_grid_size=150

# Scattering radial grid
R_min="0.5"
R_max="20.0"
scatt_grid_size=250

# Use of OpenMP (0 = no, 1 = yes)
use_omp=0

# Printing properties for energies
energy_shift="0.0"
energy_scale="219474.63137054"
print_adiabatic=0

# Executables
dbasis_exe="/home/hsilva/workdir/H+OH/murrell1984/exe/dbasis.out"
bprint_exe="/home/hsilva/workdir/H+OH/murrell1984/exe/bprint.out"
cmatrix_exe="/home/hsilva/workdir/H+OH/murrell1984/exe/cmatrix.out"
cprint_exe="/home/hsilva/workdir/H+OH/murrell1984/exe/cprint.out"

# PES external datafiles and/or any other dependency needed (format: "file_a file_b etc")
pes_extern_data=""

# Batch job configuration (slurm and/or pbs)
wall_time="24:00:00"
max_memory="5Gb"
nodes=20
mpi_cpus=250
omp_threads=1
queue_name=""
modules=""

# Libraries to load by LD_LIBRARY_PATH (format: 'path_a:path_b:path_c:etc')
env_ld_path='$HOME/lib/gsl/lib'

# OpenMP configuration
env_omp_threads="OMP_NUM_THREADS=$omp_threads"

# MPI configuration
mpi_pin_domain="I_MPI_PIN_DOMAIN=omp"
mpi_proc_placement="I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0"

# Misc
input_filename="input.d"
slurm_filename="job.sh"
pbs_filename="job.pbs"
bprint_datafile="basis_arrang=*_ch=*_J=*.dat"
cmatrix_datafile="cmatrix_arrang=*_n=*_J=*.dat"

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
		rm -rf $1/*
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

		if [ "$env_ld_path" != "" ]
		then
			echo $env_ld_path                            >> $1
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
		echo "bprint_datafile=$bprint_datafile"         >> $1
		echo "cmatrix_datafile=$cmatrix_datafile"       >> $1

		echo ""                                         >> $1
		echo 'echo "# Calculation starting at $(date)"' >> $1

		echo 'echo ""'                                  >> $1

		echo ""                                         >> $1
		echo 'cd $bin_dir/'                             >> $1

		echo '$dbasis_exe $input'                       >> $1

		echo ""                                         >> $1
		echo 'echo ""'                                  >> $1

		if [ "$bprint_exe" != "" ]
		then
			echo ""                                      >> $1
			echo '$bprint_exe $input'                    >> $1
			echo 'mv $bprint_datafile $basis_wavef_dir/' >> $1
		fi

		echo ""                                         >> $1
		echo 'echo ""'                                  >> $1

		if [ "$cmatrix_exe" != "" ]
		then
			echo ""                                      >> $1
			echo '$mpirun_call $cmatrix_exe $input'      >> $1
		fi

		echo ""                                         >> $1
		echo 'echo ""'                                  >> $1

		if [ "$cprint_exe" != "" ]
		then
			echo ""                                      >> $1
			echo '$cprint_exe $input'                    >> $1
			echo 'mv $cmatrix_datafile $channels_dir/'   >> $1
		fi

		echo ""                                         >> $1
		echo 'echo ""'                                  >> $1

		echo 'echo "# Calculation ending at $(date)"'   >> $1
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

if [ "$cprint_exe" != "" ]
then
	assert_file $cprint_exe
fi

if [ "$mpi_cpus" == "0" ]
then
	mpirun_call=""
	mpi_cpus=$(($mpi_cpus + 1))
else
	if [ "$mpi_pin_domain" == "" ]
	then
		mpirun_call="mpirun -np $mpi_cpus -genv OMP_NUM_THREADS=$omp_threads"
	else
		mpirun_call="mpirun -np $mpi_cpus -genv OMP_NUM_THREADS=$omp_threads -genv $mpi_pin_domain"
	fi
fi

if [ "$env_ld_path" != "" ]
then
	env_ld_path='export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:'$env_ld_path'"'
fi

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

		if [ "$queue_name" != "" ]
		then
			echo "#PBS -q $queue_name"                                                                        >> $filename
		fi

		echo '#PBS -N "'$job_name'"'                                                                        >> $filename

		echo ""                                                                                             >> $filename
		echo "# Job script generated at $(date) by $0"                                                      >> $filename

		build_batch_job $filename

		if [ ! "$pes_extern_data" == "" ]
		then
			cp -f $pes_extern_data $bin_dir
		fi
	done

	echo "J=$J"
done

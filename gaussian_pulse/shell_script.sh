#!/bin/sh -l
# FILENAME: myjobsubmissionfile

#SBATCH -A physics		# Allocation name
##SBATCH -p queue-name		# Queue name
#SBATCH --nodes=1		# number of nodes
#SBATCH --ntasks=1		# number of MPI tasks
#SBATCH --time=00:10:00		# wall time
#SBATCH --job-name dir_script	# Job name
#SBATCH -o myjob.o		# Name of stdout output file
#SBATCH -e myjob.e		# Name of stderr error file


# module load mathematica/12.3
##cd $PBS_O_WORKDIR

# for i in -1 -2 -3 -4 -5
# do

#     cp -r ~/num_cal/band_modulation/kx_5_ky_5/. ~/num_cal/band_modulation/kx_0_ky_$i
#     cd ~/num_cal/band_modulation/kx_0_ky_$i
#     # parti = $1
#     # ntasks = $2
#     sed -i 's/-A .*/-A '$1'/' job_script
#     sed -i 's/--ntasks=.*/--ntasks='$2'/' job_script
#     sed -i 's/--job-name .*/--job-name kx_0_ky_'$i'/' job_script
#     sed -i 's/sk = .*/sk = 0;/' num_current_delay.nb
#     sed -i 's/ktx = .*/ktx = 0;/g' num_current_delay.nb
#     sed -i 's/kty = .*/kty = '$i'\/100;/g' num_current_delay.nb
#     sbatch job_script
# done

for i in 1 2 3 4
do

    cp -r ~/num_cal/gaussian_pulse/MgO/kx_0_ky_5/. ~/num_cal/gaussian_pulse/MgO/kx_0_ky_"$i"
    cd ~/num_cal/gaussian_pulse/MgO/kx_0_ky_"$i"
    # parti = $1
    # ntasks = $2
    sed -i 's/-A .*/-A '$1'/' job_script
    sed -i 's/--ntasks=.*/--ntasks='$2'/' job_script
    sed -i 's/--job-name .*/--job-name kx_0_ky_'$i'/' job_script
    # sed -i 's/sk = .*/sk = 0;/' num_current_delay.nb
    # sed -i 's/ktx = .*/ktx = '$i'\/100;/g' num_current_delay.nb
    sed -i 's/kty = .*/kty = '$i'\/100;/g' num_current_delay.nb
    sbatch job_script
done

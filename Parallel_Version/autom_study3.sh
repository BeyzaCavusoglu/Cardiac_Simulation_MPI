#!/bin/bash
#
# You should only work under the /scratch/users/<username> directory.
#
# Example job submission script
#
# -= Resources =-
#
#SBATCH --job-name=cardiacsim
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cores-per-socket=16
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --output=cardiacsim.out

module load gcc/9.3.0
module load openmpi/4.0.1
export PATH=/kuacc/apps/openmpi/4.0.1/bin/:$PATH
export LD_LIBRARY_PATH=/kuacc/apps/openmpi/4.0.1/lib/:$LD_LIBRARY_PATH
export PATH=/kuacc/apps/openmpi/4.0.1/include:$PATH
export OMPI_MCA_btl=^openib

################################################################################
##################### !!! DO NOT EDIT ABOVE THIS LINE !!! ######################
################################################################################

echo "Compiling..."
make clean
make cardiacsim
make cardiacsim_mpi_1d
make cardiacsim_mpi_openmp
make cardiacsim_mpi_2d

echo "Running Job...!"
echo "==============================================================================="
echo "Running compiled binary..."

# Redirect the output to a text file
output_file="cardiacsim_results_study3.txt"

{
  # Study 3
  echo ".....STUDY 3 - COMMUNICATION DISABLED....."

  # Serial version
  lscpu
  echo "Serial version..."
  ./cardiacsim -n 64 -t 100 -p 100
  echo ".................................................."

  # Parallel version openmp 1d
  echo "Parallel version with 1 process"
  mpirun -np 16 ./cardiacsim_mpi_1d -n 64 -t 100 -y 1 -x 16 -k
  echo ".................................................."

  echo "Parallel version with 2 processes"
  mpirun -np 16 ./cardiacsim_mpi_1d -n 64 -t 100 -y 16 -x 1 -k
  echo ".................................................."

  echo "Parallel version with 4 processes"
  mpirun -np 16 ./cardiacsim_mpi_1d -n 64 -t 100 -y 4 -x 4 -k
  echo ".................................................."

  echo "Parallel version with 8 processes"
  mpirun -np 16 ./cardiacsim_mpi_1d -n 64 -t 100 -y 2 -x 8 -k
  echo ".................................................."

  echo "Parallel version with 16 processes"
  mpirun -np 16 ./cardiacsim_mpi_1d -n 64 -t 100 -y 8 -x 2 -k
  echo ".................................................."
  echo ".................................................."
  echo ".................................................."
  echo ".................................................."
  echo ".................................................."
  echo "N SHRINKED TO 32..."


 # Parallel version openmp 1d
  echo "Parallel version with 1 process"
  mpirun -np 16 ./cardiacsim_mpi_1d -n 32 -t 100 -y 1 -x 16 -k
  echo ".................................................."

  echo "Parallel version with 2 processes"
  mpirun -np 16 ./cardiacsim_mpi_1d -n 32 -t 100 -y 16 -x 1 -k
  echo ".................................................."

  echo "Parallel version with 4 processes"
  mpirun -np 16 ./cardiacsim_mpi_1d -n 32 -t 100 -y 4 -x 4 -k
  echo ".................................................."

  echo "Parallel version with 8 processes"
  mpirun -np 16 ./cardiacsim_mpi_1d -n 32 -t 100 -y 2 -x 8 -k
  echo ".................................................."

  echo "Parallel version with 16 processes"
  mpirun -np 16 ./cardiacsim_mpi_1d -n 32 -t 100 -y 8 -x 2 -k


echo ".................................................."


} > "$output_file"

echo "Results saved to $output_file"

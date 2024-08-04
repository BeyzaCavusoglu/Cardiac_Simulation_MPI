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
output_file="cardiacsim_results.txt"

{
  # Study one
  echo ".....STUDY 1....."

  # Serial version
  lscpu
  echo "Serial version..."
  ./cardiacsim -n 256 -t 100 -p 100
  echo ".................................................."

  # Parallel version openmp 1d
  echo "Parallel version with 1 process"
  mpirun -np 1 ./cardiacsim_mpi_1d -n 64 -t 100 -y 1 -x 1
  echo ".................................................."

  echo "Parallel version with 2 processes"
  mpirun -np 2 ./cardiacsim_mpi_1d -n 64 -t 100 -y 1 -x 2
  echo ".................................................."

  echo "Parallel version with 4 processes"
  mpirun -np 4 ./cardiacsim_mpi_1d -n 64 -t 100 -y 1 -x 4
  echo ".................................................."

  echo "Parallel version with 8 processes"
  mpirun -np 8 ./cardiacsim_mpi_1d -n 64 -t 100 -y 1 -x 8
  echo ".................................................."

  echo "Parallel version with 16 processes"
  mpirun -np 16 ./cardiacsim_mpi_1d -n 64 -t 100 -y 1 -x 16

  # Study 2
  echo ".....STUDY 2....."


  # Parallel version openmp 2d
  echo "Parallel 2D version with 1 process"
  mpirun -np 1 ./cardiacsim_mpi_2d -n 64 -t 100 -y 1 -x 1
  echo ".................................................."

  echo "Parallel 2D version with 2 processes"
  mpirun -np 2 ./cardiacsim_mpi_2d -n 64 -t 100 -y 1 -x 2
  echo ".................................................."

  echo "Parallel 2D version with 4 processes"
  mpirun -np 4 ./cardiacsim_mpi_2d -n 64 -t 100 -y 1 -x 4
  echo ".................................................."

  echo "Parallel 2D version with 8 processes"
  mpirun -np 8 ./cardiacsim_mpi_2d -n 64 -t 100 -y 1 -x 8
  echo ".................................................."

  echo "Parallel 2D version with 16 processes"
  mpirun -np 16 ./cardiacsim_mpi_2d -n 64 -t 100 -y 1 -x 16


 

} > "$output_file"

echo "Results saved to $output_file"

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
  echo ".................................................."



  # Study 3
  echo ".....STUDY 3....."

#Different configuration of MPI+OpenMP
#[1 + 16] [2 + 8] [4 + 4] [8 + 2] [ 16 + 1]

export KMP_AFFINITY=verbose,compact

echo "MPI1 + OMP16"
export OMP_NUM_THREADS=16
mpirun -np 1 ./cardiacsim_mpi_openmp -n 64 -t 100 -y 1 -x 1

echo ".................................................."

echo "MPI2 + OMP8"
export OMP_NUM_THREADS=8
mpirun -np 2 ./cardiacsim_mpi_openmp -n 64 -t 100 -y 1 -x 2

echo ".................................................."

echo "MPI4 + OMP4"
export OMP_NUM_THREADS=4
mpirun -np 4 ./cardiacsim_mpi_openmp -n 64 -t 100 -y 1 -x 4


echo ".................................................."


echo "MPI8 + OMP2"
export OMP_NUM_THREADS=2
mpirun -np 8 ./cardiacsim_mpi_openmp -n 64 -t 100 -y 1 -x 8


echo ".................................................."


echo "MPI16 + OMP1"
export OMP_NUM_THREADS=1
mpirun -np 16 ./cardiacsim_mpi_openmp -n 64 -t 100 -y 1 -x 16

echo ".................................................."


} > "$output_file"

echo "Results saved to $output_file"

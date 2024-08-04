/*
 * Solves the Panfilov model using an explicit numerical scheme.
 * Based on code orginally provided by Xing Cai, Simula Research Laboratory
 * and reimplementation by Scott B. Baden, UCSD
 *
 * Modified and  restructured by Didem Unat, Koc University
 * MPI 2D VERSION ADDED BY BEYZA CAVUSOGLU, Parcorelab
 *
 */
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <math.h>
#include <sys/time.h>

#include <mpi.h>

using namespace std;

// 2D MPI IMPLEMENTATION FUNCTIONS:

/*
 -----------------INFORMATION BEGINS------------------

************** 1 *************
MPI_Scatter(
    void* send_data,
    int send_count,
    MPI_Datatype send_datatype,
    void* recv_data,
    int recv_count,
    MPI_Datatype recv_datatype,
    int root,
    MPI_Comm communicator)
******************************

**************** 2 ************

  MPI_Isend() AND MPI_Irecv() --  Begins a nonblocking send and receive

  The class of non-blocking protocols returns from the send or
  receive operations before it is semantically safe to do so.
  Thus copy operation may not be completed on return
  Non-blocking send and receive operations in MPI
  “I” stands for “Immediate”

  int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request)

  --> Processing continues immediately without waiting for the message to be
      copied out from the application buffer.


  int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source,
              int tag, MPI_Comm comm, MPI_Request *request)

  –-> Processing continues immediately without actually waiting for the message to
  be received and copied into the the application buffer.
***************************


  ------------------INFORMATION ENDS------------------
*/


// distribute_to_processes(dest, source, size_x, size_y);
void distribute_to_processes(double **source, double **dest, int size_x, int size_y, int rank, int m, int n )
{

  int area_size = size_x * size_y;
 // print2D(source, size_x, size_x);
  //printf("BEFORE SCATTER with size %d at %d\n", yarea_size, rank);
  MPI_Scatter(source[0], area_size, MPI_DOUBLE, dest[1], area_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //printf("after SCATTER at %d\n", rank);
  //if(rank == 0){
     // printf("process : %d", rank);
//  printf("sizex: %d, sizey: %d ",size_x, size_y);
  //print2D(dest, m+2, n+2);
  //printf("\n\n");
  //}
}

// reverse of the distribution function above
void gather_from_processes(double**dest, double **source, int size_x, int size_y, int m, int n)
{
  int area_size = size_x * size_y;
  MPI_Gather(dest[1], area_size, MPI_DOUBLE, source[0], area_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

// Utilities
//

// Timer
// Make successive calls and take a difference to get the elapsed time.
static const double kMicro = 1.0e-6;
double getTime()
{
  struct timeval TV;
  struct timezone TZ;

  const int RC = gettimeofday(&TV, &TZ);
  if (RC == -1)
  {
    cerr << "ERROR: Bad call to gettimeofday" << endl;
    return (-1);
  }

  return (((double)TV.tv_sec) + kMicro * ((double)TV.tv_usec));

} // end getTime()

// Allocate a 2D array
double **alloc2D(int m, int n)
{
  double **E;
  int nx = n, ny = m;
  E = (double **)malloc(sizeof(double *) * ny + sizeof(double) * nx * ny);
  assert(E);
  int j;
  for (j = 0; j < ny; j++)
    E[j] = (double *)(E + ny) + j * nx;
  return (E);
}

// Reports statistics about the computation
// These values should not vary (except to within roundoff)
// when we use different numbers of  processes to solve the problem
double stats(double **E, int m, int n, double *_mx)
{
  double mx = -1;
  double l2norm = 0;
  int i, j;
  for (j = 1; j <= m; j++)
    for (i = 1; i <= n; i++)
    {
      l2norm += E[j][i] * E[j][i];
      if (E[j][i] > mx)
        mx = E[j][i];
    }
  *_mx = mx;
  l2norm /= (double)((m) * (n));
  l2norm = sqrt(l2norm);
  return l2norm;
}

// External functions
extern "C"
{
  void splot(double **E, double T, int niter, int m, int n);
}
void cmdLine(int argc, char *argv[], double &T, int &n, int &px, int &py, int &plot_freq, int &no_comm, int &num_threads);

void simulate(double **E, double **E_prev, double **R,
              const double alpha, const int n, const int m, const double kk,
              const double dt, const double a, const double epsilon,
              const double M1, const double M2, const double b, int pos_X, int pos_Y, int px, int py)
{
  int i, j;
  /*
   * Copy data from boundary of the computational box
   * to the padding region, set up for differencing
   * on the boundary of the computational box
   * Using mirror boundaries
   */

  if (pos_X == 0)
    for (j = 1; j <= m; j++)
      E_prev[j][0] = E_prev[j][2];
  if (pos_X == px - 1)
    for (j = 1; j <= m; j++)
      E_prev[j][n + 1] = E_prev[j][n - 1];
  if (pos_Y == 0)
    for (i = 1; i <= n; i++)
      E_prev[0][i] = E_prev[2][i];
  if (pos_Y == py - 1)
    for (i = 1; i <= n; i++)
      E_prev[m + 1][i] = E_prev[m - 1][i];

  // Solve for the excitation, the PDE
  for (j = 1; j <= m; j++)
  {
    for (i = 1; i <= n; i++)
    {
      E[j][i] = E_prev[j][i] + alpha * (E_prev[j][i + 1] + E_prev[j][i - 1] - 4 * E_prev[j][i] + E_prev[j + 1][i] + E_prev[j - 1][i]);
    }
  }

  /*
   * Solve the ODE, advancing excitation and recovery to the
   *     next timtestep
   */
  for (j = 1; j <= m; j++)
  {
    for (i = 1; i <= n; i++)
      E[j][i] = E[j][i] - dt * (kk * E[j][i] * (E[j][i] - a) * (E[j][i] - 1) + E[j][i] * R[j][i]);
  }

  for (j = 1; j <= m; j++)
  {
    for (i = 1; i <= n; i++)
      R[j][i] = R[j][i] + dt * (epsilon + M1 * R[j][i] / (E[j][i] + M2)) * (-R[j][i] - kk * E[j][i] * (E[j][i] - b - 1));
  }
}

// Main program
int main(int argc, char **argv)
{
  /*
   *  Solution arrays
   *   E is the "Excitation" variable, a voltage
   *   R is the "Recovery" variable
   *   E_prev is the Excitation variable for the previous timestep,
   *      and is used in time integration
   */

  int num_processes, rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_processes); // Get the number of processes
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);          // Get process id

 // cout << "BURASI MPI 2D VERSIYON" << endl;

/*   if (rank == 0)
    printf("Hello from process %d of %d\n", rank, num_processes);
 */
  double **E, **R, **E_prev;

  // For MPI processes
  double **local_E, **local_E_prev, **local_R;

  // Various constants - these definitions shouldn't change
  const double a = 0.1, b = 0.1, kk = 8.0, M1 = 0.07, M2 = 0.3, epsilon = 0.01, d = 5e-5;

  double T = 1000.0;
  int m = 200, n = 200;
  int plot_freq = 0;
  int px = 1, py = 1;
  int no_comm = 0;
  int num_threads = 1;

  cmdLine(argc, argv, T, n, px, py, plot_freq, no_comm, num_threads);
  m = n;
  // Allocate contiguous memory for solution arrays
  // The computational box is defined on [1:m+1,1:n+1]
  // We pad the arrays in order to facilitate differencing on the
  // boundaries of the computation box
  E = alloc2D(m + 2, n + 2);
  E_prev = alloc2D(m + 2, n + 2);
  R = alloc2D(m + 2, n + 2);

  int i, j;
  // Initialization
  for (j = 1; j <= m; j++)
    for (i = 1; i <= n; i++)
      E_prev[j][i] = R[j][i] = 0;

  for (j = 1; j <= m; j++)
    for (i = n / 2 + 1; i <= n; i++)
      E_prev[j][i] = 1.0;

  for (j = m / 2 + 1; j <= m; j++)
    for (i = 1; i <= n; i++)
      R[j][i] = 1.0;

  double dx = 1.0 / n;

  // For time integration, these values shouldn't change
  double rp = kk * (b + 1) * (b + 1) / 4;
  double dte = (dx * dx) / (d * 4 + ((dx * dx)) * (rp + kk));
  double dtr = 1 / (epsilon + ((M1 / M2) * rp));
  double dt = (dte < dtr) ? 0.95 * dte : 0.95 * dtr;
  double alpha = d * dt / (dx * dx);

  cout << "Grid Size       : " << n << endl;
  cout << "Duration of Sim : " << T << endl;
  cout << "Time step dt    : " << dt << endl;
  cout << "Process geometry: " << px << " x " << py << endl;
  if (no_comm)
    cout << "Communication   : DISABLED" << endl;

  cout << endl;

  // Start the timer
  double t0 = getTime();

  // Simulated time is different from the integer timestep number
  // Simulated time
  double t = 0.0;
  // Integer timestep number
  int niter = 0;

  // 2D MPI IMPLEMENTATION:

  // size_X = GridSize_in_x_direction / process_geo_x
  int size_X = n / px;

  // pos_X = process_id / process_geo_x
  int pos_X = rank % px;

  // size_Y = GridSize_in_y_direction / process_geo_y
  int size_Y = m / py;

  // pos_Y = process_id / process_geo_x
  int pos_Y = rank / px;

  // Allocation for the local data of processes
  local_E = alloc2D(n + 2, m + 2);
  local_E_prev = alloc2D(n + 2, m + 2);
  local_R = alloc2D(n + 2, m + 2);

  // After the allocation for processes, data should be sent to processes memory using Scatter
   distribute_to_processes(E_prev,  local_E_prev,   size_X+1,   size_Y+2, rank, m, n);


  double *north, *south, *west, *east, *west_2, *east_2;

  int rank_north = rank - px;
  int rank_south = rank + px;
  int rank_west = rank + 1;
  int rank_east = rank - 1;

  /*
  Ghost cells are allocated for the information needed from
  other processes. Thus each processor allocates space to accommodate ghost region
  */

  west = (double *)malloc((size_Y+2)* sizeof(double));
  east = (double *)malloc((size_Y+2)* sizeof(double));
  west_2 = (double *)malloc((size_Y+2)* sizeof(double));
  east_2 = (double *)malloc((size_Y+2)* sizeof(double));

  // Necessary for MPI_Isend and MPI_Irecv and MPI_Waitall
  MPI_Request arr_request[4];

  // Necessary for MPI_Waitall
  MPI_Status arr_status[4];

  while (t < T)
  { // while time set by user is not exceeded

    // Ghost cell communications between neighbor processes should happen here after allocation.

    int request = 0;

    // North neighbors - receive and send
    if (pos_Y != 0)
    { // check if its the top-most one
    MPI_Irecv(local_E_prev[0], size_X+2, MPI_DOUBLE, rank_north, 0, MPI_COMM_WORLD, &arr_request[request++]);
    MPI_Isend(local_E_prev[1], size_X+2, MPI_DOUBLE, rank_north, 0, MPI_COMM_WORLD, &arr_request[request++]); 
    }

    // South neighbors - receive and send
    if (pos_Y != ((num_processes - 1) / px))
    { // check if its not exceeding the south-most one

    MPI_Irecv(local_E_prev[n-1], size_X+2, MPI_DOUBLE, rank_south, 0, MPI_COMM_WORLD, &arr_request[request++]);
    MPI_Isend(local_E_prev[n-2], size_X+2, MPI_DOUBLE, rank_south, 0, MPI_COMM_WORLD, &arr_request[request++]); 
    }

    // West neighbors - receive and send
    if (pos_X != 0)
    { // check if not the west most one

      for (int i=0; i < size_Y; i++){
        west_2[i] = local_E_prev[i+1][1];
      }

      MPI_Irecv(&west[0], size_Y, MPI_DOUBLE, rank_east, 0, MPI_COMM_WORLD, &arr_request[request++]);
      MPI_Isend(&west_2[0], size_Y, MPI_DOUBLE, rank_east, 0, MPI_COMM_WORLD, &arr_request[request++]);
    }
    // east neighbors - receive and send
    if (pos_X != ((num_processes - 1) % px))
    {
    for (int i=0; i < size_Y; i++){
      east_2[i] = local_E_prev[i+1][size_X];
      }
      MPI_Irecv(&east[0], size_Y, MPI_DOUBLE, rank_west, 0, MPI_COMM_WORLD, &arr_request[request++]);
      MPI_Isend(&east_2[0], size_Y, MPI_DOUBLE, rank_west, 0, MPI_COMM_WORLD, &arr_request[request++]);
    }


    // For the west part
    if (pos_X != 0)
    {
      for (int i = 0; i < size_Y; i++)
      {
        local_E_prev[i + 1][0] = west[i];
      }
    }
    // For the east part
    if (pos_X != ((num_processes - 1) % px))
    {
      for (int i = 0; i < size_Y; i++)
      {
        local_E_prev[i + 1][size_X + 1] = east[i];
      }
    }
  
    MPI_Waitall(request, arr_request, arr_status);
    //cout << "MPI communication completed in process " << rank << endl;

  t += dt;
  niter++;

  // simulate(E, E_prev, R, alpha, n, m, kk, dt, a, epsilon, M1, M2, b);
  simulate(local_E, local_E_prev, local_R, alpha, size_X, size_Y, kk, dt, a, epsilon, M1, M2, b, pos_X, pos_Y, px, py);

  // swap current E with previous E
  double **tmp = local_E;
  local_E = local_E_prev;
  local_E_prev = tmp;

  if (plot_freq)
  {
    int k = (int)(t / plot_freq);
    if ((t - k * plot_freq) < dt)
    {
      //splot(E, t, niter, m + 2, n + 2); //closed for performance tests
    }

  } 
}// end of while loop
  
  // Now after computation, the data needed to be gathered from processes to master
  gather_from_processes(E_prev,local_E_prev ,size_X+1, size_Y+2, m, n);
  
  // Master process only! Just do once
  if (rank == 0)
  {
    double time_elapsed = getTime() - t0;

    double Gflops = (double)(niter * (1E-9 * n * n) * 28.0) / time_elapsed;
    double BW = (double)(niter * 1E-9 * (n * n * sizeof(double) * 4.0)) / time_elapsed;

    cout << "Number of Iterations        : " << niter << endl;
    cout << "Elapsed Time (sec)          : " << time_elapsed << endl;
    cout << "Sustained Gflops Rate       : " << Gflops << endl;
    cout << "Sustained Bandwidth (GB/sec): " << BW << endl
         << endl;
    double mx;
    double l2norm = stats(E_prev, m, n, &mx);
    cout << "Max: " << mx << " L2norm: " << l2norm << endl;
  
  }

  MPI_Finalize();

  return 0;
}

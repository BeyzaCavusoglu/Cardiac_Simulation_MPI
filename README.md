# Cardiac_Simulation_MPI
The Aliev-Panfilov heart electrophysiology simulator using MPI and OpenMP.

## Description:
A cardiac electrophysiology simulator can be used for clinical diagnostic and therapeutic purposes. Cell simulators entail solving a coupled set of a equations: a system of Ordinary Differential Equations (ODEs) together with Partial Differential Equations (PDEs). 
The simulator models the propagation of electrical signals in the heart, and it incorporates a cell model describing the kinetics of a the membrane of a single cell. The PDE couples multiple cells into a system. There can be different cell models, with varying degrees of complexity. 
The simulator keeps track of two state variables that characterize the electrophysiology that are being simulated. Both are represented as 2D arrays. The first variable, called the excitation, is stored in the E[ ][ ] array. The second variable, called the recovery variable, is stored in the R[ ][ ] array. Lastly, Eprev is stored, the voltage at the previous timestep.This is necessary to advance the voltage over time. Since the method of finite differences is used to solve the problem, the variables E and R are discretized by considering the values only at a regularly spaced set of discrete points. Here is the formula for solving the PDE, where the E and Eprev refer to the voltage at current and previous timestep, respectively, and the constant α is defined in the simulator:

<img width="481" alt="a1" src="https://github.com/user-attachments/assets/a52f7654-9af1-4db4-8765-13453dc9c285">

Here is the formula for solving the ODE, where references to E and R correspond to whole arrays. Expressions involving whole arrays are pointwise operations (the value on i,j depends only on the value of i,j), and the constants kk, a, b, ϵ, M1 and M2 are defined in the simulator and dt is the time step size:

<img width="481" alt="a2" src="https://github.com/user-attachments/assets/57b6fe5f-6770-4cfa-ab4f-72e5ca0e4e0b">

# Serial Code:
A working serial simulator was already provided that uses the Aliev-Panfilov cell model. The simulator includes a plotting capability (using gnuplot) which can be used to debug the code, and also to observe the simulation dynamics. The plot frequency can be adjusted from command line. The timing results for performance tests will be taken when the plotting is disabled.

<img width="483" alt="a3" src="https://github.com/user-attachments/assets/d612d755-75c6-48af-9a59-7c3439519c7c">

*Note: This project is done for a course project, the serialized version of this code was already designed and given by the instructor. My implementation was to parallelize them using OPENMP and conducting the performance tests.*

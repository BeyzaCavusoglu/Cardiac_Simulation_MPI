include ./arch.gnu
# OPTIMIZATION = -fast
# OPTIMIZATION = -O3
# DEBUG += -g


app:		cardiacsim 

OBJECTS = cardiacsim.o splot.o cmdLine.o
MPI_OBJECTS = cardiacsim_mpi_1d.o splot.o cmdLine.o
OPENMP_OBJECTS = cardiacsim_mpi_openmp.o splot.o cmdLine.o
MPI_OBJECTS_2D = cardiacsim_mpi_2d.o splot.o cmdLine.o


cardiacsim:	        $(OBJECTS) 
		$(C++LINK) $(LDFLAGS) -o $@ $(OBJECTS)  $(LDLIBS)

cardiacsim_mpi_1d: 		    $(MPI_OBJECTS) 
		$(C++LINK) -g $(LDFLAGS) -o $@ $(MPI_OBJECTS) $(LDLIBS)

cardiacsim_mpi_openmp: 		    $(OPENMP_OBJECTS) 
		$(C++LINK) $(LDFLAGS) $(OPENMP) -o $@ $(OPENMP_OBJECTS) $(LDLIBS) -lgomp

cardiacsim_mpi_2d: 		    $(MPI_OBJECTS_2D) 
		$(C++LINK) $(LDFLAGS) -o $@ $(MPI_OBJECTS_2D) $(LDLIBS)

clean:	
	$(RM) *.o cardiacsim *~;
	$(RM) cardiacsim_mpi_1d *~;
	$(RM) cardiacsim_mpi_openmp *~;
	$(RM) cardiacsim_mpi_2d *~;
	$(RM) core;




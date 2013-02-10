FC = mpif90
FFLAGS = -O3
#FFLAGS = -O3 -g -traceback -check all

LIBNAME = libdstcltmtx.a

OBJ_DISTRIBUTE = $(LIBNAME)(complex_dstmtx.o) \
		$(LIBNAME)(double_dstmtx.o) \
		$(LIBNAME)(real_dstmtx.o) 

OBJ_COLLECT = $(LIBNAME)(complex_cltmtx.o) \
		$(LIBNAME)(double_cltmtx.o) \
		$(LIBNAME)(real_cltmtx.o) 

OBJ_EXTRA = $(LIBNAME)(findnode.o)

all: $(LIBNAME) 

$(LIBNAME): $(OBJ_DISTRIBUTE) $(OBJ_COLLECT) $(OBJ_EXTRA)

complex_dstmtx.o: dstmtx.F90
	$(FC) $(FFLAGS) -c -o $@ -DDISTRIBUTE=complex_dstmtx \
	  -DTIPO=complex*16 -DTYPE_MPI=MPI_DOUBLE_COMPLEX $<

double_dstmtx.o: dstmtx.F90
	$(FC) $(FFLAGS) -c -o $@ -DDISTRIBUTE=double_dstmtx \
	  -DTIPO=real*8 -DTYPE_MPI=MPI_DOUBLE_PRECISION $<

real_dstmtx.o: dstmtx.F90
	$(FC) $(FFLAGS) -c -o $@ -DDISTRIBUTE=real_dstmtx -DTIPO=real*4 \
	  -DTYPE_MPI=MPI_REAL $<

complex_cltmtx.o: cltmtx.F90
	$(FC) $(FFLAGS) -c -o $@ -DCOLLECT=complex_cltmtx \
	  -DTIPO=complex*16 -DTYPE_MPI=MPI_DOUBLE_COMPLEX $<

double_cltmtx.o: cltmtx.F90
	$(FC) $(FFLAGS) -c -o $@ -DCOLLECT=double_cltmtx -DTIPO=real*8 \
	  -DTYPE_MPI=MPI_DOUBLE_PRECISION $<

real_cltmtx.o: cltmtx.F90
	$(FC) $(FFLAGS) -c -o $@ -DCOLLECT=real_cltmtx \
	  -DTIPO=real*4 -DTYPE_MPI=MPI_REAL $<

findnode.o: findnode.F90
	$(FC) $(FFLAGS) -c -o $@  $<

clean:
	rm -f $(LIBNAME) *.o core*

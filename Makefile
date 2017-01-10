CC= gcc

WOPT= -Wall -fopenmp
#WOPT+= -D_WEIGH_L2
LIB_GSL= -L/opt/local/lib
INC_GSL= -I/opt/local/include/gsl
LIB_HP= -L/opt/local/lib
INC_HP= -I/opt/local/include
LIB_FITS= -L/opt/local/lib
INC_FITS= -I/opt/local/include
LIB_SHARP= -L/Users/allisonradmin/Documents/Cambridge/computing/libsharp-code/auto/lib
INC_SHARP= -I/Users/allisonradmin/Documents/Cambridge/computing/libsharp-code/auto/include

CFLAGS= $(WOPT) -I./src $(INC_GSL) $(INC_HP) $(INC_FITS) $(INC_SHARP)
LIBS= $(LIB_GSL) $(LIB_HP) $(LIB_FITS) $(LIB_SHARP)
LIBS+= -lgsl -lgslcblas -lchealpix -lsharp -lfftpack -lc_utils -lcfitsio -lm

COMMONO= src/common.o
HEO= src/healpix_extra.o
MASTERO= src/master.o
MAINO= src/main.o
OBJ= $(COMMONO) $(HEO) $(MASTERO) $(MAINO)

COMMVQEO= src/common_mvqe.o
CGO= src/cg.o
RNGO= src/rng.o
MAINMVQEO= src/main_mvqe.o
OBJMVQE= $(COMMONO) $(HEO) $(COMMVQEO) $(RNGO) $(CGO) $(MAINMVQEO)

EXEC= NaMaster MVQE
all: $(EXEC)

NaMaster : $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LIBS) -o $@

MVQE : $(OBJMVQE)
	$(CC) $(CFLAGS) $(OBJMVQE) $(LIBS) -o $@

clean :
	rm -f $(EXEC) $(OBJ) $(OBJMVQE)
cleaner :
	rm -f $(EXEC) $(OBJ) $(OBJMVQE) *~ src/*~

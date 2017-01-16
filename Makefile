CC= gcc
WOPT= -Wall -fopenmp
#WOPT+= -D_WEIGH_L2
LIB_GSL= -L/home/damonge/lib
INC_GSL= -I/home/damonge/include
LIB_HP=
INC_HP=
LIB_FITS=
INC_FITS= 
LIB_SHARP=
INC_SHARP=

CFLAGS= $(WOPT) -I./src $(INC_GSL) $(INC_HP) $(INC_FITS) $(INC_SHARP)
LIBS= $(LIB_GSL) $(LIB_HP) $(LIB_FITS) $(LIB_SHARP)
LIBS+= -lgsl -lgslcblas -lchealpix -lsharp -lfftpack -lc_utils -lcfitsio -lm

COMMONO= src/common.o
HEO= src/healpix_extra.o
MASTERO= src/master.o
MAINO= src/main.o
OBJ= $(COMMONO) $(HEO) $(MASTERO) $(MAINO)

EXEC= NaMaster
all: $(EXEC)

NaMaster : $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LIBS) -o $@

clean :
	rm -f $(EXEC) $(OBJ)
cleaner :
	rm -f $(EXEC) $(OBJ) *~ src/*~

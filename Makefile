CC= icc
WOPT= -Wall -fopenmp
LIB_GSL= -L/users/damonge/lib
INC_GSL= -I/users/damonge/include
LIB_HP=
INC_HP=
LIB_FITS=
INC_FITS=
LIB_SHARP=
INC_SHARP=

CFLAGS= $(WOPT) -I./src $(INC_GSL) $(INC_HP) $(INC_FITS) $(INC_SHARP)
LIBS= $(LIB_GSL) $(LIB_HP) $(LIB_FITS) $(LIB_SHARP)
LIBS+= -lUtilsDAM -lgsl -lgslcblas -lchealpix -lsharp -lfftpack -lc_utils -lcfitsio -lm

HEO= src/healpix_extra.o
MASTERO= src/master.o
MAINO= src/main.o
OBJ= $(HEO) $(MASTERO) $(MAINO)

EXEC= namaster
all: $(EXEC)

namaster : $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LIBS) -o $@

clean :
	rm -f $(EXEC) $(OBJ)
cleaner :
	rm -f $(EXEC) $(OBJ) *~ src/*~

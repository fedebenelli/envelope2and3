FC=gfortran

CFLAGS=-llapack -lblas -fno-align-commons

SRCS=Flash MulticompSRK_PR ThermoRoutines_RKPR main
OBJS=${SRCS:=.o}

MAIN=main

BIN=bin

all: ${MAIN}
	@echo Compilation ended!

${MAIN}: ${OBJS}
	-mkdir ${BIN}
	${FC} ${CFLAGS} -o ${BIN}/envelope2and3 ${OBJS}

.SUFFIXES: .f90 .f90.o

.f90.o:
	${FC} ${CFLAGS} -c $< 

clean:
	@rm *.o

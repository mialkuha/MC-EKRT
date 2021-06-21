IDIR =./include
SDIR =./src
ODIR =./obj
LDIR =./lib
BDIR =./bin
DBDIR =./bin/Debug
CC=/usr/bin/clang++

CFLAGS=-I$(IDIR) -I$(SDIR) -I. \
	   -I/usr/local/include -I/usr/local/include/boost \
	   -I/home/mikko/vcpkg/installed/x64-linux/include \
	   -march=corei7 -O1 -std=c++2a -m64 -pthread -g
WFLAGS=-Wall -Wextra -Wshadow -Wnon-virtual-dtor \
	   -Wold-style-cast -Wcast-align -Wunused \
	   -Woverloaded-virtual -Wpedantic -Wconversion \
	   -Wsign-conversion -Wdouble-promotion -Wformat=2 \
	   -pedantic -Weffc++ #-fsanitize=address,undefined
FLAGS=$(CFLAGS) $(WFLAGS)

LIBS=-L/usr/local/lib -lgsl -lpthread\
	   -lgslcblas -lm -lLHAPDF

_SRCS_CPP = ars.cpp hcubature.cpp linear_interpolator.cpp nn_coll.cpp nucleon.cpp nucleus_generator.cpp pqcd.cpp
SRCS_CPP = $(patsubst %,$(SDIR)/%,$(_SRCS_CPP))

_SRCS_H = ars.hpp converged.h cubature.h linear_interpolator.hpp linterp.h nn_coll.hpp nucleon.hpp nucleus_generator.hpp pqcd.hpp typedefs.hpp vwrapper.h
SRCS_H = $(patsubst %,$(IDIR)/%,$(_SRCS_H))
SRCS = $(SRCS_CPP) $(SRCS_H)

_OBJS = $(_SRCS_CPP:.cpp=.o)
OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

build : $(ODIR)/main.o $(OBJS)
	$(CC) $(OBJS) $(ODIR)/main.o $(FLAGS) $(LIBS) -o $(BDIR)/mcaa

debug : $(ODIR)/main.o $(OBJS)
	$(CC) -g $(OBJS) $(ODIR)/main.o $(FLAGS) $(LIBS) -o $(DBDIR)/mcaa

clean :
	rm $(OBJS) $(ODIR)/main.o

$(ODIR)/main.o : main.cpp $(OBJS)
	$(CC) -c $(FLAGS) $< -o $@

$(ODIR)/hcubature.o: $(SDIR)/hcubature.cpp
	$(CC) -c $(CFLAGS) $< -o $@

$(ODIR)/%.o: $(SDIR)/%.cpp
	$(CC) -c $(FLAGS) $< -o $@
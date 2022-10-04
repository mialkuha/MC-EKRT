IDIR =./include
SDIR =./src
ODIR =./obj
LDIR =./lib
BDIR =./bin
DBDIR =./bin/Debug
CC=clang++

CFLAGS=-I$(IDIR) -I$(SDIR) -I. \
	   -I/usr/local/include -I/usr/local/include/boost \
	   -I/home/mikko/vcpkg/installed/x64-linux/include \
	   -I/home/mikko/src/lhapdf/include \
	   -O1 -std=c++20 -m64 -pthread
WFLAGS=-Wall -Wextra -Wshadow -Wnon-virtual-dtor \
	   -Wold-style-cast -Wcast-align -Wunused \
	   -Woverloaded-virtual -Wpedantic -Wconversion \
	   -Wsign-conversion -Wdouble-promotion -Wformat=2 \
	   -pedantic -Weffc++
FLAGS=$(CFLAGS) $(WFLAGS)

LIBS=-L/usr/local/lib -lgsl -lpthread\
	   -lgslcblas -lm -ltbb -lLHAPDF\
	   -L/home/mikko/src/lhapdf/lib

_SRCS_CPP = ars.cpp hcubature.cpp histo.cpp linear_interpolator.cpp nn_coll.cpp nucleus_generator.cpp pqcd.cpp
SRCS_CPP = $(patsubst %,$(SDIR)/%,$(_SRCS_CPP))

_SRCS_H = ars.hpp converged.h cubature.h generic_helpers.hpp high_level_calcs.hpp histo.hpp io_helpers.hpp linear_interpolator.hpp linterp.h nn_coll.hpp nucleon.hpp nucleus_generator.hpp pqcd.hpp typedefs.hpp vwrapper.h
SRCS_H = $(patsubst %,$(IDIR)/%,$(_SRCS_H))
SRCS = $(SRCS_CPP) $(SRCS_H)

_OBJS = $(_SRCS_CPP:.cpp=.o)
OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

dir_guard=@mkdir -p $(@D)

build : $(ODIR)/main.o $(OBJS)
	$(dir_guard)
	@mkdir -p $(BDIR)
	$(CC) $(OBJS) $(ODIR)/main.o $(FLAGS) $(LIBS) -o $(BDIR)/mcaa

debug : $(ODIR)/main.o $(OBJS)
	$(dir_guard)
	$(CC) -g $(OBJS) $(ODIR)/main.o $(FLAGS) $(LIBS) -o $(DBDIR)/mcaa

clean :
	rm $(OBJS) $(ODIR)/main.o

$(ODIR)/main.o : main.cpp $(OBJS)
	$(dir_guard)
	$(CC) -c $(FLAGS) $< -o $@

$(ODIR)/hcubature.o: $(SDIR)/hcubature.cpp
	$(dir_guard)
	$(CC) -c $(CFLAGS) $< -o $@

$(ODIR)/%.o: $(SDIR)/%.cpp
	$(dir_guard)
	$(CC) -c $(FLAGS) $< -o $@
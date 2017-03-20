#
# Isotropic Distortion of Gaussian field
#

EXEC      = discrete_power
DIRS      = $(GSL_DIR)

all: $(EXEC)

INCLDIRS  = $(foreach dir, $(DIRS), -I$(dir)/include)
LIBDIRS   = $(foreach dir, $(DIRS), -L$(dir)/lib) 
CXXFLAGS := $(INCLDIRS) -O2

OBJS := discrete_power.o power.o

LIBS := $(LIBDIRS) -lm -lboost_program_options
LIBS += -lgsl -lgslcblas

discrete_power: $(OBJS)
	$(CXX) $(OBJS) $(LIBS) -o $@


corr.o: corr.cpp corr.h power.h
growth.o: growth.cpp
power.o: power.cpp power.h

.PHONY: clean
clean :
	rm -f $(EXEC) $(OBJS)


#
# Makefile
#

########################################################################### 
#
# Compiling and Linking Support
#
########################################################################### 

#
# Note: the variables CC, CFLAGS, LDFLAGS, etc., are usually left to be set
# outside of the makefile. Most system administrators set these variables so
# that compilations work "outside the box" with no developer intervention.
# Because we're compiling MPI programs using specific compilers and specific
# MPI libraries, they're redefined here for each system.
#

# --------------------------------
# Gordon Configuration
# --------------------------------

#
# Compilers
#

# Intel Compilers (Default)
CC = mpicc
CXX = mpicxx

CFLAGS = -g -Wall -shared-intel -DMPICH_IGNORE_CXX_SEEK -mpe=mpilog
LDFLAGS= -shared-intel -mpe=mpilog
#
# MPI Library
#
#CFLAGS = -I/opt/mpich-xl64/include -q64
#CXXFLAGS = -I/opt/mpich-xl64/include -q64
#LDFLAGS = -L/opt/mpich-xl64/lib -lmpich -q64

########################################################################### 
#
# File Declarations
#
########################################################################### 

#
# C and Object Files
#
C_FILES = conway.c
O_FILES = conway.o

#
# Main Targets
#
all: cgl

cgl: $(O_FILES) 
	$(CC) -o cgl $(O_FILES) $(LDFLAGS) 
#
# Housekeeping Targets
#
.PHONY: clean
clean:		
	/bin/rm -f core $(O_FILES) cgl

#
# Dependencies
#

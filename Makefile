# Makefile for CLASSIC.

# Usage:

# make mode=<mode>

#   or, for the default mode, just:

# make

# where <mode> is:
# 
# serial        - compiles code for serial runs using the GNU compiler on the local platform (default)
# parallel      - compiles code for parallel runs using the GNU compiler on the local platform
# ppp           - compiles code for parallel runs using the Intel compiler on the front-ends
# supercomputer - compiles code for parallel runs using the Intel compiler on the supercomputers

# The "cray" mode uses the native Cray compiler on the supercomputers. It still requires testing.

# ==================================================================================================================

# Modifications: Ed Chan, Sep 2018.
#   - Updated various compiler options.
#   - Added "cray" and "ppp" modes.
#   - All objects go into directories specific to the mode used.
#   - All executable names are labelled with the mode appended.

# Object files
OBJ = fileIOModule.o classic_params.o ctem_statevars.o class_statevars.o peatlands_mod.o \
	generalUtils.o APREP.o GRALB.o mvidx.o SNOW_ALBVAL.o SNOW_TRANVAL.o SNOALBA.o \
	TNPREP.o WFILL.o CANADD.o GRDRAN.o SNOALBW.o  TPREP.o WFLOW.o CANALB.o CLASSG.o \
	GRINFL.o SNOVAP.o TSOLVC.o WPREP.o wetland_methane.o ctemg1.o bio2str.o ctems1.o ctemg2.o \
	PHTSYN3.o CANVAP.o CLASSI.o CWCALC.o ICEBAL.o SUBCAN.o TSOLVE.o XIT.o CGROW.o CLASSS.o \
	DIASURFZ.o SCREENRH.o TFREEZ.o TSPOST.o CHKWAT.o CLASST.o DRCOEF.o SLDIAG.o TMCALC.o TSPREP.o \
	CLASSA.o CLASSW.o FLXSURFZ.o SNINFL.o TMELT.o TWCALC.o CLASSB.o CLASSZ.o GATPREP.o SNOADD.o \
	TNPOST.o WEND.o balcar.o mainres.o allocate.o phenolgy.o turnover.o mortality.o \
	disturb.o ctems2.o landuse_change_mod.o soil_ch4uptake.o \
	competition_mod.o hetres_mod.o ctemUtilities.o ctem.o outputManager.o prepareOutputs.o  \
	model_state_drivers.o read_from_job_options.o metModule.o main.o xmlParser.o xmlManager.o CLASSIC.o

# Object files (.o and .mod) are placed into directories named according to the mode. In the serial
# compiler case, these are placed into the "objectFiles" directory.
ODIR = objectFiles_$(mode)

# Set variables based on the "mode" specified.

# *** NB: These settings are particular to the ECCC platforms. Other users will need to find their own recipes.
#         The use most general option is "serial".

ifeq ($(mode), supercomputer)
	# Wrapper to the default Fortran compiler loaded via a module. The following is specific to the Intel compiler.
	# Include/library flags should NOT be specified as the appropriate ones should be available via the loaded module.
	COMPILER = ftn
	# Fortran Flags.
 	FFLAGS = -DPARALLEL -r8 -align array64byte -O1 -g -init=arrays -init=zero -traceback -module $(ODIR)
else ifeq ($(mode), cray)
	# Wrapper to the default Fortran compiler loaded via a module. The following is specific to the Cray compiler.
	# Include/library flags should NOT be specified as the appropriate ones should be available via the loaded module.
	COMPILER = ftn
	# Fortran Flags. Note: -g prevents most optimizations, -rm gives listings.
	FFLAGS = -DPARALLEL -s real64 -e0 -ez -O2 -g
else ifeq ($(mode), ppp)
	# Location of parallel netCDF library.
	NCLIBS=/fs/ssm/comm/eccc/ccrn/nemo/nemo-1.0/openmpi-1.6.5/enable-shared/intelcomp-2016.1.156/ubuntu-14.04-amd64-64
	# Open MPI wrapper to the default Fortran compiler. The following is specific to the Intel compiler.
	# Note: -xCORE-AVX2 generates specialized (optimized) code which must be supported by the system's processors.
	# Include/library flags for Open MPI should NOT be specified as the appropriate ones should be available via the wrapper.
	COMPILER = mpif90
	# Fortran Flags. Switches are as close as possible to those used for Intel on the Cray supercomputer.
	FFLAGS = -DPARALLEL -xCORE-AVX2 -r8 -align array64byte -O1 -g -init=arrays -init=zero -traceback -module $(ODIR)
	# Include Flags for netCDF.
	IFLAGS = -I$(NCLIBS)/include
	# Library Flags for netCDF.
	LFLAGS = -L$(NCLIBS)/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5
else ifeq ($(mode), parallel)
        # Parallel compiler.
        COMPILER = mpif90
        # Fortran Flags. The following is specific to the gfortran compiler.
        FFLAGS = -DPARALLEL -O3 -g -fdefault-real-8 -ffree-line-length-none -fbacktrace -ffpe-trap=invalid,zero,overflow -fbounds-check -J$(ODIR)  
        # Include Flags.
        IFLAGS =  -I${HOME}/PnetCDF/include 
        # Library Flags.
        LFLAGS = -L${HOME}/PnetCDF/lib -lnetcdf -ldl -lz -lm 
else
	# Serial compiler.
	COMPILER = gfortran
	ODIR = objectFiles
	# Fortran Flags.
	FFLAGS = -O3 -g -fdefault-real-8 -ffree-line-length-none -fbacktrace -ffpe-trap=invalid,zero,overflow -fbounds-check -J$(ODIR) #-Wall -Wextra
	# Include Flags.
	#IFLAGS = -I/home/rjm/Public/NETCDF/include -J$(ODIR)
	IFLAGS = -I/usr/include
	# Library Flags.
	#LFLAGS = -L/home/rjm/Public/NETCDF/lib -lnetcdff -L/home/rjm/Public/NETCDF/lib -fPIC -lnetcdf -lnetcdf -L/home/rjm/Public/NETCDF/lib
	LFLAGS = -lnetcdff -ldl -lz -lm 
endif

# Create required directory/.gitignore file, if missing.
VOID := $(shell mkdir -p $(ODIR))
VOID := $(shell [ ! -f $(ODIR)/.gitignore ] && cp objectFiles/.gitignore $(ODIR))

# RECIPES
# Compile object files from .F90 sources
$(ODIR)/%.o: src/%.F90
	$(COMPILER) $(FFLAGS) $(IFLAGS) -c $< -o $@

# Compile object files from .f90 sources
$(ODIR)/%.o: src/%.f90
	$(COMPILER) $(FFLAGS) $(IFLAGS) -c $< -o $@

# Compile object files from .f (Fortran 77) sources
$(ODIR)/%.o: src/%.f
	$(COMPILER) $(FFLAGS) $(IFLAGS) -c $< -o $@
	
# Properly reference the ODIR for the linking
OBJD = $(patsubst %,$(ODIR)/%,$(OBJ))

# Link objects together and put executable in the bin/ directory
CLASSIC: $(OBJD)
	$(COMPILER) $(FFLAGS) $(IFLAGS) -o bin/CLASSIC_$(mode) $(OBJD) $(LFLAGS)

# "make clean mode=supercomputer" removes all object files in objectFiles_supercomputer
clean:
	rm -f $(ODIR)/*.o $(ODIR)/*.mod

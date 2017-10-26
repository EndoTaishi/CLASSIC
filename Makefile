# Object files
OBJ = dataTransferModule.o fileIOModule.o ctem_params.o ctem_statevars.o  peatlands_mod.o \
	class_statevars.o generalUtils.o APREP.o CLASSBD.o GRALB.o mvidx.o SNOW_ALBVAL.o SNOW_TRANVAL.o SNOALBA.o \
	TNPREP.o WFILL.o CANADD.o CLASSD.o GRDRAN.o SNOALBW.o  TPREP.o WFLOW.o CANALB.o CLASSG.o \
	GRINFL.o SNOVAP.o TSOLVC.o WPREP.o wetland_methane.o ctemg1.o bio2str.o ctems1.o ctemg2.o \
	PHTSYN3.o CANVAP.o CLASSI.o CWCALC.o ICEBAL.o SUBCAN.o TSOLVE.o XIT.o CGROW.o CLASSS.o \
	DIASURFZ.o SCREENRH.o TFREEZ.o TSPOST.o CHKWAT.o CLASST.o DRCOEF.o SLDIAG.o TMCALC.o TSPREP.o \
	CLASSA.o CLASSW.o FLXSURFZ.o SNINFL.o TMELT.o TWCALC.o CLASSB.o CLASSZ.o GATPREP.o SNOADD.o \
	TNPOST.o WEND.o balcar.o ORDLEG.o mainres.o allocate.o phenolgy.o turnover.o mortality.o \
	disturb.o ctems2.o landuse_change_mod.o soil_ch4uptake.o \
	competition_mod.o hetres_mod.o ctemUtilities.o ctem.o outputManager.o io_driver.o  \
	model_state_drivers.o read_from_job_options.o main.o xmlParser.o xmlManager.o CLASSIC.o

# This is the directory our .o and .mod files will be put into.
ODIR = objectFiles

# Serial compiler
GC = gfortran
# Serial Include Flags
SIFLAGS = -O3 -g -I/usr/include -J$(ODIR)
# Serial Library Flags
SLFLAGS = -lnetcdff

# Parallel compiler
MPIC = mpif90
# Parallel Include Flags
PIFLAGS = -O3 -g -I${HOME}/PnetCDF/include -L${HOME}/PnetCDF/lib -J$(ODIR)
# Parallel Library Flags
PLFLAGS = -DPARALLEL -lpnetcdf

# Check to see if the model should be compiled to run in serial or parallel
# For a parallel run, use the command "make PARALLEL=true"
# For a serial run, use the command "make" 
ifeq ($(PARALLEL), true)
	COMPILER=$(MPIC)
	IFLAGS = $(PIFLAGS)
	LFLAGS = $(PLFLAGS)
else
	COMPILER=$(GC)
	IFLAGS = $(SIFLAGS)
	LFLAGS = $(SLFLAGS)
endif

# Additional Libraries (serial/parallel)
LFLAGS += -ldl -lz -lm -g -fdefault-real-8 -ffree-line-length-none -fbacktrace -ffpe-trap=invalid,zero,overflow  -fbounds-check
#-Wall -Wextra

# RECIPES
# Compile object files from .F90 sources
$(ODIR)/%.o: src/%.F90
	$(COMPILER) $(IFLAGS) -c $< -o $@ $(LFLAGS)

# Compile object files from .f90 sources
$(ODIR)/%.o: src/%.f90
	$(COMPILER) $(IFLAGS) -c $< -o $@ $(LFLAGS)

# Compile object files from .f (Fortran 77) sources
$(ODIR)/%.o: src/%.f
	$(COMPILER) $(IFLAGS) -c $< -o $@ $(LFLAGS)

# Properly reference the ODIR for the linking
OBJD = $(patsubst %,$(ODIR)/%,$(OBJ))

# Link objects together and put executable in the bin/ directory
CLASSIC: $(OBJD)
	$(COMPILER) $(IFLAGS) -o bin/CLASSIC $(OBJD) $(LFLAGS)

# "make clean" removes all object files and executables
clean:
	rm -f $(ODIR)/*.o $(ODIR)/*.mod *~ core CLASSIC

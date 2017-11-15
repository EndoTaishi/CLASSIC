# Object files
OBJ = fileIOModule.o ctem_params.o ctem_statevars.o  peatlands_mod.o \
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

# This is the directory our .o and .mod files will be put into. (except for supercomputer runs where they are just dumped in
# same directory as the Makefile
ODIR = objectFiles

# Check to see if the model should be compiled to run in supercomputer, serial or parallel mode
# For a supercomputer run, use the command "make mode=supercomputer"
# For a parallel run, use either the command "make mode=parallel"
# For a serial run, use the command "make mode=serial" or just "make"

ifeq ($(mode), supercomputer)
	# Reset the object directory only for the supercomputer
	ODIR =.
	# Supercomputer compiler
	COMPILER = ftn
	# Supercomputer Include Flags
	IFLAGS =
	# Supercomputer Library Flags
	LFLAGS = -DPARALLEL -lpnetcdf -r8 -align array64byte -O1 -g -p -init=arrays -init=zero -traceback
else
ifeq ($(mode), parallel)
	# Parallel compiler
	COMPILER = mpif90
	# Parallel Include Flags
	IFLAGS =  -I${HOME}/PnetCDF/include -L${HOME}/PnetCDF/lib -J$(ODIR)
	# Parallel Library Flags
	LFLAGS = -DPARALLEL -lpnetcdf
	# Additional flags
	LFLAGS += -O3 -ldl -lz -lm -g -fdefault-real-8 -ffree-line-length-none -fbacktrace -ffpe-trap=invalid,zero,overflow  -fbounds-check #-Wall -Wextra
else
	# Serial compiler
	COMPILER = gfortran
	# Serial Include Flags
	IFLAGS = -I/usr/include -J$(ODIR)
	# Serial Library Flags
	LFLAGS = -lnetcdff
	# Additional flags
	LFLAGS += -O3 -g -ldl -lz -lm -fdefault-real-8 -ffree-line-length-none -fbacktrace -ffpe-trap=invalid,zero,overflow  -fbounds-check #-Wall -Wextra
endif
endif

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

# "make {mode={parallel}{supercomputer}} clean" removes all object files and executables
clean:
	rm -f $(ODIR)/*.o $(ODIR)/*.mod *~ bin/CLASSIC

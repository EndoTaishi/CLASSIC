# Make file for GNU fortran
#

OBJ = dataTransferModule.o serialFileIOModule.o ctem_params.o ctem_statevars.o  peatlands_mod.o class_statevars.o APREP.o CLASSBD.o GRALB.o mvidx.o SNOW_ALBVAL.o SNOW_TRANVAL.o SNOALBA.o TNPREP.o WFILL.o CANADD.o CLASSD.o GRDRAN.o SNOALBW.o  TPREP.o WFLOW.o CANALB.o CLASSG.o GRINFL.o SNOVAP.o TSOLVC.o WPREP.o wetland_methane.o ctemg1.o bio2str.o ctems1.o ctemg2.o PHTSYN3.o CANVAP.o CLASSI.o CWCALC.o ICEBAL.o SUBCAN.o TSOLVE.o XIT.o CGROW.o CLASSS.o DIASURFZ.o SCREENRH.o TFREEZ.o TSPOST.o CHKWAT.o CLASST.o DRCOEF.o SLDIAG.o TMCALC.o TSPREP.o CLASSA.o CLASSW.o FLXSURFZ.o SNINFL.o TMELT.o TWCALC.o CLASSB.o CLASSZ.o GATPREP.o SNOADD.o TNPOST.o WEND.o balcar.o ORDLEG.o mainres.o allocate.o phenolgy.o turnover.o mortality.o disturb.o ctems2.o competition_map.o competition_unmap.o landuse_change_mod.o soil_ch4uptake.o  competition_mod.o hetres_mod.o ctem.o io_driver.o netcdf_tools.o model_state_drivers.o read_from_job_options.o main.o supportFunctions.o xmlParser.o xmlManager.o CLASSIC.o
#outputManager.o 
# Binary dir
BDIR = bin

# GNU Fortran
FC=gfortran
MPIC=mpif90
PIFLAGS = -I${HOME}/PnetCDF/include -L${HOME}/PnetCDF/lib
PLFLAGS = -lpnetcdf

NETCDFDIR=${HOME}/Public/NETCDF_GF

export LD_LIBRARY_PATH=$HOME/Public/NETCDF_GF/lib/

#FC := $(shell $(NETCDFDIR)/bin/nf-config --fc)
FFLAGS := $(shell $(NETCDFDIR)/bin/nf-config --fflags)
LDLIBS := $(shell $(NETCDFDIR)/bin/nf-config --flibs)
LDLIBS += $(shell $(NETCDFDIR)/bin/nc-config --libs)

# Optimized running of model ----------------
#FFLAGS += -O3 -ffree-line-length-none

# Debugging of model ----------------

#GNU
FFLAGS += -g -ffree-line-length-none -fbacktrace -ffpe-trap=invalid,zero,overflow -Wall -Wextra

# Other option is -Og which allows optimization but still debugs well.

#-----------------------

LDLIBS += -ldl -lz -lm -lhdf5_hl -lhdf5

# These are the rules to make the targets
#

%.o: src/%.f
	$(MPIC) $(FFLAGS) -c $< -o $@

%.o: src/%.f90
	$(MPIC) $(FFLAGS) -c $< -o $@

%.o: src/%.F90
#	$(MPIC) $(FFLAGS) -c $< -o $@
	$(MPIC) $(FFLAGS) $(PIFLAGS) -c $< -o $@

%.mod : src/%.f90
#	$(MPIC) $(FFLAGS) -c $< -mod $@
	$(MPIC) $(PIFLAGS) -c $< -mod $@

CLASSIC: $(OBJ)
#	$(MPIC) $(FFLAGS) -o $(BDIR)/CLASSIC $(OBJ) $(LDLIBS)
	$(MPIC) $(PIFLAGS) -o $(BDIR)/CLASSIC $(OBJ) $(LDLIBS) $(PLFLAGS)

.PHONY: clean

all: CLASSIC

clean:
	rm -f *.o *.mod *~ core CLASSIC

$(ODIR):
	mkdir $(ODIR)

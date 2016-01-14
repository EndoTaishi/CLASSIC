# Make file for PGI or GNU fortran
#
# Objects that must be created, now they get created in ODIR
#
OBJ = ctem_params.o ctem_statevars.o APREP.o CLASSBD.o GRALB.o mvidx.o SNOW_ALBVAL.o SNOW_TRANVAL.o SNOALBA.o  TNPREP.o WFILL.o CANADD.o CLASSD.o GRDRAN.o SNOALBW.o TPREP.o WFLOW.o CANALB.o CLASSG.o GRINFL.o SNOVAP.o TSOLVC.o WPREP.o wetland_methane.o ctemg1.o bio2str.o ctems1.o ctemg2.o PHTSYN3.o CANVAP.o CLASSI.o CWCALC.o ICEBAL.o SUBCAN.o TSOLVE.o XIT.o CGROW.o CLASSS.o DIASURFZ.o SCREENRH.o TFREEZ.o TSPOST.o CHKWAT.o CLASST.o DRCOEF.o SLDIAG.o TMCALC.o TSPREP.o CLASSA.o CLASSW.o FLXSURFZ.o SNINFL.o TMELT.o TWCALC.o CLASSB.o CLASSZ.o GATPREP.o SNOADD.o TNPOST.o WEND.o balcar.o GAUSSG.o ORDLEG.o TRIGL.o mainres.o allocate.o phenolgy.o turnover.o mortality.o disturb.o ctems2.o competition_map.o competition_unmap.o read_from_job_options.o landuse_change_mod.o competition_mod.o hetres_mod.o ctem.o io_driver.o runclass36ctem.o

# Binary dir
#
BDIR = ../bin

# PGI
FC=pgf90

# GNU Fortran
#FC=gfortran

# General running of model (PGI)

FFLAGS = -Bstatic -r8 -O2 -gopt -Mbyteswapio -Mbackslash -Mpreprocess -Kieee -uname -Ktrap=fp

export PGIMACH=linux86-64

# Debugging of model ----------------

# PGI
#FFLAGS = -r8 -Minform,warn -g -Mbyteswapio -Mbackslash -Mpreprocess -Kieee -uname -Ktrap=fp,align,denorm,unf -traceback -Mbounds

#export PGIMACH=linux86-64

#GNU
#FFLAGS = -g -fdefault-real-8 -ffree-line-length-none -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow -Waliasing -Wampersand -Wconversion -Wsurprising -Wintrinsics-std -Wno-tabs -Wintrinsic-shadow

#-----------------------


# These are the rules to make the targets
#

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.mod : %.f90
	$(FC) $(FFLAGS) -c $< -mod $@

CLASS36CTEM: $(OBJ)
	 $(FC) $(FFLAGS) -o $(BDIR)/CLASS36CTEM $(OBJ)

.PHONY: clean

clean:
	rm -f *.o *.mod *~ core CLASS36CTEM

$(ODIR):
	mkdir $(ODIR)	

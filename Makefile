# Make file for intel fortran
#
# Objects that must be created, now they get created in ODIR
#
OBJ = APREP.o CLASSBD.o GRALB.o SNOALBA.o TNPREP.o WFILL.o CANADD.o CLASSD.o GRDRAN.o SNOALBW.o  TPREP.o WFLOW.o CANALB.o CLASSG.o GRINFL.o SNOVAP.o TSOLVC.o WPREP.o CTEMG1.o BIO2STR.o CTEMS1.o CTEMG2.o PHTSYN.o CANVAP.o CLASSI.o CWCALC.o ICEBAL.o SUBCAN.o TSOLVE.o XIT.o CGROW.o CLASSS.o DIASURFZ.o RUNCLASS35CTEM.o TFREEZ.o TSPOST.o CHKWAT.o CLASST.o DRCOEF.o SLDIAG.o TMCALC.o TSPREP.o CLASSA.o CLASSW.o FLXSURFZ.o SNINFL.o TMELT.o TWCALC.o CLASSB.o CLASSZ.o GATPREP.o SNOADD.o TNPOST.o WEND.o CTEM.o BALCAR.o GAUSSG.o ORDLEG.o TRIGL.o LUC.o MAINRES.o HETRESV.o HETRESG.o ALLOCATE.o PHENOLGY.o TURNOVER.o MORTALTY.o DISTURB.o CTEMS2.o BIOCLIM.o EXISTENCE.o COMPETITION.o read_from_job_options.o

# Binary dir
#
BDIR = ../bin

# This is our compiler
#
FC=pgf90

# General running of model
#FFLAGS = -r8 -O2 -gopt -Mbyteswapio -Mbackslash -Mpreprocess -Kieee -uname -Ktrap=fp

# Debugging of model
FFLAGS = -r8 -O0 -Minform,warn -g -Mbyteswapio -Mbackslash -Mpreprocess -Kieee -uname -Ktrap=fp

export PGIMACH=-64

# These are the rules to make the targets
#
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

CLASS35CTEM: $(OBJ)
	 $(FC) $(FFLAGS) -o $(BDIR)/CLASS35CTEM $(OBJ)

.PHONY: clean

clean:
	rm -f *.o *~ core CLASS35CTEM

$(ODIR):
	mkdir $(ODIR)	

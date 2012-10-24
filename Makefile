# Make file for PGI or GNU fortran
#
# Objects that must be created, now they get created in ODIR
#
OBJ = APREP.o CLASSBD.o GRALB.o SNOALBA.o TNPREP.o WFILL.o CANADD.o CLASSD.o GRDRAN.o SNOALBW.o  TPREP.o WFLOW.o CANALB.o CLASSG.o GRINFL.o SNOVAP.o TSOLVC.o WPREP.o CTEMG1.o BIO2STR.o CTEMS1.o CTEMG2.o PHTSYN.o CANVAP.o CLASSI.o CWCALC.o ICEBAL.o SUBCAN.o TSOLVE.o XIT.o CGROW.o CLASSS.o DIASURFZ.o SCREENRH.o RUNCLASS36CTEM.o TFREEZ.o TSPOST.o CHKWAT.o CLASST.o DRCOEF.o SLDIAG.o TMCALC.o TSPREP.o CLASSA.o CLASSW.o FLXSURFZ.o SNINFL.o TMELT.o TWCALC.o CLASSB.o CLASSZ.o GATPREP.o SNOADD.o TNPOST.o WEND.o CTEM.o BALCAR.o GAUSSG.o ORDLEG.o TRIGL.o LUC.o MAINRES.o HETRESV.o HETRESG.o ALLOCATE.o PHENOLGY.o TURNOVER.o MORTALTY.o DISTURB.o CTEMS2.o BIOCLIM.o EXISTENCE.o COMPETITION.o COMPETITION_MAP.o COMPETITION_UNMAP.o read_from_job_options.o

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

CLASS36CTEM: $(OBJ)
	 $(FC) $(FFLAGS) -o $(BDIR)/CLASS36CTEM $(OBJ)

.PHONY: clean

clean:
	rm -f *.o *~ core CLASS36CTEM

$(ODIR):
	mkdir $(ODIR)	

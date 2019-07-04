# How do I ... / Something has gone wrong! {#howDoI}

1. @ref chgSoil
2. @ref failStart
3. @ref makeClean
4. @ref soilColourIndex

----

# Change the number/depth/etc of soil layers? {#chgSoil}

Information about the soil layers is taken in from the initialization netcdf file. model_state_drivers::read_modelsetup reads the netcdf for the number of soil layers (*ignd*). The *DELZ* variable in the intialization file is the thickness of each layer. So change the initialization file for the number and thicknesses of the soil layers as required and CLASSIC will simply read them in.

# My run starts and I get an immediate fail! {#failStart}

If you start a new run and it immediately fails with an output something like:

        Singularity jormelton-containerCLASSIC-master-latest.simg:~/Documents/CLASSIC> bin/CLASSIC configurationFiles/template_job_options_file.txt 254.85/53.98
         in process           0           1           1   254.84999999999999        53.979999999999997                1           1           1           1
          CANOPY ENERGY BALANCE         1       8092.77620036       8077.49643236
                 0.000000       1.875648      -2.045030      -0.189964       0.000000    8073.385790
                 0.000000       1.339229     257.925000
         died on   254.84999999999999        53.979999999999997 
         
This is a general indication of a problem in your initialization file. Often it is fixed by setting snow/liquid in the canopy to zero (see @ref initPhysProgVar). You may also get snow energy/water balance failing, again look to @ref initPhysProgVar for advice with setup.

# My run won't compile but I didn't do anything I can think of. {#makeClean}

Sometimes you need to run a '`make mode=???? clean`' (where mode is the one you used to compile the code, for site-level it is typically `serial`) to clean out old .mod and .o files that might be causing issues. Do that and then recompile.

# initFileConverter tells me I need a soil colour index {#soilColourIndex}

Soil colour index is described in @ref soilData and used in soilProperties.f90. A global file of soil colour index produced by Peter Lawrence (NCAR) can be obtained from ftp://ftp.cccma.ec.gc.ca/pub/jmelton/mksrf_soilcol_global_c090324.nc

# Gotcha for Windows users {#windowsEndings}

Text files created on DOS/Windows machines have different line endings than files created on Unix/Linux. DOS uses carriage return and line feed ("\r\n") as a line ending, while Unix uses just line feed ("\n"). You need to be careful about transferring files between Windows machines and Unix machines to make sure the line endings are translated properly.(from http://www.cs.toronto.edu/~krueger/csc209h/tut/line-endings.html)

Any CLASSIC tool that reads in ASCII expects the Linux/Unix line endings.

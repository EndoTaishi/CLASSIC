Create netcdf model initialization files from depricated INI (and CTM) file formats or fortran namelists
========

The initFileConverter program takes in ASCII CLASS initialization files (.INI) and, if required,
CTEM initialization files (.CTM) and converts them to the newer NetCDF format. It is also possible
to take in a fortran namelist format. This is the preferred method for site-level users.

# Compile

To compile, use the `make` command from the console.
This will generate the binary executable in the *bin* folder.

# Execute

To run the program, use the following command:

`bin/initFileConverter [file.INI] [file.CTM]

The expected format of the INI file is:
ASCII to NetCDF met file loader
========

The metASCIILoader program takes in ASCII met files and converts them to the newer NetCDF format.

# Compile

To compile, use the `make` command from the console.
This will generate the binary executable in the *bin* folder.

# Execute

To run the program, use the following command:

`bin/metASCIILoader [file.MET] [longitude] [latitude]`

Where *longitude* and *latitude* are the desired real values for the local file and the *file.MET* is your local met file; it should look like this:

` 10  0    1  1901   195.90   259.21    0.0000E+00     6.64   5.426E-03    1.09    99244.75`

` 10 30    1  1901   205.17   259.21    0.0000E+00     7.01   5.543E-03    1.06    99247.07`

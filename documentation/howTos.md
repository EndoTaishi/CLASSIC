# How do I ... {#howDoI}

1. @ref chgSoil
2. @ref windowsEndings

----

# Change the number/depth/etc of soil layers? {#chgSoil}

Information about the soil layers is taken in from the initialization netcdf file. model_state_drivers::read_modelsetup reads the netcdf for the number of soil layers (ignd). The ZBOT variable in the intialization file is the depth to the bottom of each layer. So change the initialization file for the number and depths of the soil layers as required and CLASSIC will simply read them in.

# Gotcha for Windows users {#windowsEndings}

Text files created on DOS/Windows machines have different line endings than files created on Unix/Linux. DOS uses carriage return and line feed ("\r\n") as a line ending, which Unix uses just line feed ("\n"). You need to be careful about transferring files between Windows machines and Unix machines to make sure the line endings are translated properly.(from http://www.cs.toronto.edu/~krueger/csc209h/tut/line-endings.html)

Any CLASSIC tool that reads in ASCII expects the Linux/Unix line endings.

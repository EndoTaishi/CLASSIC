# How do I ... {#howDoI}

1. @ref chgSoil

----

# Change the number/depth/etc of soil layers? {#chgSoil}

Information about the soil layers is taken in from the initialization netcdf file. model_state_drivers::read_modelsetup reads the netcdf for the number of soil layers (ignd). The ZBOT variable in the intialization file is the depth to the bottom of each layer. So change the initialization file for the number and depths of the soil layers as required and CLASSIC will simply read them in.

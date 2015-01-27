#!/bin/sh

# NOTE: this likely can only work on your local machine. Mike informs me that
# this will not work on the lxsrv2 or 3 servers since they need some updates for
# this to work. So compile locally then run on the servers. If this doesn't 
# compile on your local, you might need to get Mike to do the updates on your
# local machine. JM 25.06.2012

#export PGIVER=11.6
#export PGIMACH=linux86-64

#creator
#pgf90 -Bstatic creator_module_bigfile.f90 netcdf_creator_bigfile.f90 -I /users/tor/acrn/rls/local/aix64/netcdf4.1.3/include -o netcdf_creator_bf  -L/users/tor/acrn/rls/local/aix64/netcdf4.1.3/lib -lnetcdff -lnetcdf -L/usr/local/hdf5-1.8.7/lib -lhdf5_hl -lhdf5 -L/usr/local/zlib-1.2.5 -lz -lm -lhdf5_hl -lhdf5 -lz -L/usr/local/hdf5-1.8.7/lib -lhdf5_hl -lhdf5 -L/usr/local/zlib-1.2.5/lib -lz

#xlf90_r creator_module_bigfile.f90 netcdf_creator_bigfile.f90 -I /users/tor/acrn/rls/local/aix64/netcdf4.1.3/include -o netcdf_creator_bf  -L/users/tor/acrn/rls/local/aix64/netcdf4.1.3/lib -q64 -O -qarch=auto -qtune=auto -qextname -qrealsize=8 -qfree=f90 -qsuffix=f=f90

xlf90_r -lnetcdff -lnetcdf creator_module_bigfile.f90 netcdf_creator_bigfile.f90 -I/fs/dev/crb/had01/data/sundry/include  -o netcdf_creator_bf  -L/fs/dev/crb/had01/data/sundry/lib -qmaxmem=-1 -q64 -O -qarch=auto -qtune=auto -qrealsize=8 -qfree=f90 -qsuffix=f=f90
#writer
#pgf90 -Bstatic creator_module_bigfile.f90 netcdf_writer_bigfile.f90 -I /users/tor/acrn/rls/local/aix64/netcdf4.1.3/include -o netcdf_writer_bf  -L/users/tor/acrn/rls/local/aix64/netcdf4.1.3/lib -lnetcdff -lnetcdf -L/usr/local/hdf5-1.8.7/lib -lhdf5_hl -lhdf5 -L/usr/local/zlib-1.2.5 -lz -lm -lhdf5_hl -lhdf5 -lz -L/usr/local/hdf5-1.8.7/lib -lhdf5_hl -lhdf5 -L/usr/local/zlib-1.2.5/lib -lz

#xlf90_r creator_module_bigfile.f90 netcdf_writer_bigfile.f90 -I /users/tor/acrn/rls/local/aix64/netcdf4.1.3/include -o netcdf_writer_bf  -L/users/tor/acrn/rls/local/aix64/netcdf4.1.3/lib -q64 -O -qarch=auto -qtune=auto -qextname -qrealsize=8 -qfree=f90 -qsuffix=f=f90

xlf90_r -lnetcdff -lnetcdf creator_module_bigfile.f90 netcdf_writer_bigfile.f90 -I/fs/dev/crb/had01/data/sundry/include  -o netcdf_writer_bf  -L/fs/dev/crb/had01/data/sundry/lib -q64 -O -qarch=auto -qtune=auto -qmaxmem=-1 -qrealsize=8 -qfree=f90 -qsuffix=f=f90

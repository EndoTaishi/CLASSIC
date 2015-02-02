#!/bin/sh

# NOTE: this likely can only work on your local machine. Mike informs me that
# this will not work on the lxsrv2 or 3 servers since they need some updates for
# this to work. So compile locally then run on the servers. If this doesn't 
# compile on your local, you might need to get Mike to do the updates on your
# local machine. JM 25.06.2012

export PGIVER=11.6
export PGIMACH=linux86-64

#creator
pgf90 -Bstatic creator_module_bigfile.f90 netcdf_creator_bigfile.f90 -I /usr/local/netcdf-4.1.3/include -o netcdf_creator_bf  -L/usr/local/netcdf-4.1.3/lib -lnetcdff -lnetcdf -L/usr/local/hdf5-1.8.7/lib -lhdf5_hl -lhdf5 -L/usr/local/zlib-1.2.5 -lz -lm -lhdf5_hl -lhdf5 -lz -lcurl -L/usr/local/hdf5-1.8.7/lib -lhdf5_hl -lhdf5 -L/usr/local/zlib-1.2.5/lib -lz -lcurl

#writer
pgf90 -Bstatic creator_module_bigfile.f90 netcdf_writer_fast.f90 -I /usr/local/netcdf-4.1.3/include -o netcdf_writer_fast  -L/usr/local/netcdf-4.1.3/lib -lnetcdff -lnetcdf -L/usr/local/hdf5-1.8.7/lib -lhdf5_hl -lhdf5 -L/usr/local/zlib-1.2.5 -lz -lm -lhdf5_hl -lhdf5 -lz -lcurl -L/usr/local/hdf5-1.8.7/lib -lhdf5_hl -lhdf5 -L/usr/local/zlib-1.2.5/lib -lz -lcurl


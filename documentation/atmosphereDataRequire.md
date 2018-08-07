# Atmospheric Forcing Data {#forcingData}

At each physics time step, for each grid cell or modelled area, the following atmospheric forcing data are,

**Required**:
- FDLROW Downwelling longwave sky radiation \f$[ W m^{-2} ]\f$
- PREROW Surface precipitation rate \f$[kg m^{-2} s^{-1} ]\f$
  - CLASSIC is able to run with total incoming precipitation, partitioning it into rainfall and snowfall on the basis of empirically-derived equations. If the rainfall rate (RPREROW) and snowfall rate (SPREROW) are available, they could be used instead. The getMet and updateMet subroutines of the @ref model_state_drivers should be modified accordingly, and the job options switch IPCP should be set to 4.
- PRESROW Surface air pressure \f$[P_a ]\f$
- QAROW Specific humidity at reference height \f$[kg kg^{-1} ]\f$
- TAROW Air temperature at reference height [K]
  - For atmospheric models, the air temperature supplied to CLASS should be the lowest level air temperature extrapolated using the dry adiabatic lapse rate to the bottom of the atmosphere, i.e. to where the wind speed is zero and the pressure is equal to the surface pressure, Pa . For field data, the actual measured air temperature at the reference height should be used, since in this case the adiabatic extrapolation is performed within CLASS.
- VMODROW Wind speed at reference height \f$[m s^-1 ]\f$
  - Atmospheric models provide the zonal and meridional components of the wind velocity, but CLASS does not actually require information on wind direction. Thus, if only the scalar wind speed is available, either ULROW or VLROW can be set to it, and the other to zero. (Both of these terms, plus the scalar wind speed VMODROW, must be supplied to CLASS.)
- FSSROW Downwelling shortwave radiation incident on a horizontal surface \f$[W m^{-2} ]\f$
  -CLASS ordinarily requires that the forcing incoming shortwave radiation be partitioned into the visible and near-infrared components. If these are not available, however, they can each be roughly estimated as approximately half of the total incoming solar radiation.  Note: if the ISNOALB switch (see below) is set to 1, the incoming shortwave radiation is required in four wavelength bands, and the direct and diffuse components are required as well.


**Optional**:
- FCLOROW Fractional cloud cover [ ]
  - The fractional cloud cover is used in the calculation of the transmissivity of the vegetation canopy. If it is not available it can be estimated on the basis of the solar zenith angle and the occurrence of precipitation (COMBAK see the section on the RUNCLASS driver).
- FSIHROW Near infrared shortwave radiation incident on a horizontal surface \f$[W m^{-2} ]\f$ (see above for FSSROW)
- FSVHROW Visible shortwave radiation incident on a horizontal surface \f$[W m^{-2} ]\f$ (see above for FSSROW)
- ULROW Zonal component of wind velocity \f$[m s^{-1} ]\f$ (see above for VMODROW)
- VLROW Meridional component of wind velocity \f$[m s^{-1} ]\f$ (see above for VMODROW)


**Specified**:
- ZBLDROW Atmospheric blending height for surface roughness length averaging [m]
  - If the surface being modelled is a very heterogeneous one, care must be taken to ensure that the reference heights are greater than the “blending height”, the distance above the surface at which the atmospheric variables are not dominated by any one surface type. In principle this height depends on the length scale of the roughness elements; it may be as large as 50-100 m. In CLASSIC the blending height is used in averaging the roughness lengths over the modelled area, and is read in separately from ZRFMROW and ZRFHROW as ZBLDROW.
- ZRFHROW Reference height associated with forcing air temperature and humidity [m]
- ZRFMROW Reference height associated with forcing wind speed [m]
  - In atmospheric models the forcing wind speed, air temperature and specific humidity are obtained from the lowest modelled atmospheric layer, and thus the reference height will be the height above the “surface” (i.e. the location where the wind speed is zero and the pressure is equal to the surface pressure Pa) corresponding to that lowest layer. Some atmospheric models use a vertical co-ordinate system in which the momentum and thermodynamic levels are staggered, and if so, ZFRMROW and ZRFHROW will have different values. If that is the case, the switch ISLFD in the job options file should be set to 2, so that the subroutines src/FLXSURFZ.f and src/DIASURFZ.f are called, since the other options do not support different reference heights. In the case of field data, the reference height is the height above the ground surface at which the variables are measured. If the measurement height for wind speed is different from that for the air temperature and specific humidity, again the ISLFD switch in the job options file should be set to 2.
  - **Note** that neither ZRFHROW nor ZRFMROW may be smaller than the vegetation canopy height, as this will cause the model run to crash.

## Advisement regarding the physics timestep

The length of the time step should be carefully considered in assembling the forcing data. CLASS has been designed to run at a time step of 30 minutes or less, and the explicit prognostic time stepping scheme used for the soil, snow and vegetation variables is based on this assumption. Longer time steps may lead to the emergence of numerical instabilities in the modelled prognostic variables. The physics timestep can be changed in COMBAK

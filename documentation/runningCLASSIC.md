# Running CLASSIC {#runStandAloneMode}

Some points:
1. runclass36ctem.f is the primary driver via which most input data is read (CTEM's initialization files are read in io_driver.f90) and subroutines are called.
2. The model initialization files include one for CLASS (SITENAME.INI) and one for CTEM (SITENAME.CTM). Examples provided in the Benchmarks folder.
3. The file containing meteorological data (SITENAME.MET). Also provided in the Benchmarks folder.
4. A Makefile is included which compiles all .f and .f90 files and generates an executable called CLASS36CTEM. The standard Makefile, currently set up for pgf90 (but can also be setup for xlf and gfortran) may not be suitable for your Unix station and might require modification. You may need your system administrator’s help to do this. Note that when the code is compiled it generates .mod files that corresponds to the Fortran modules.

Other than the .INI, .CTM  and .MET files which are mandatory, some other CTEM related files are also required depending on which CTEM functionality is switched on. These required files are discussed further down in this user guide.


## The Input Files {#inputs}

### The MET file (Meteorological forcings) {#MET}

The .MET files contain (usually) half-hourly values of all seven variables that are needed to drive the CLASS land surface scheme. Please see the Benchmarks folder for an example .MET file.  The format of the .MET file cannot be changed since CLASS reads in these data using a fixed FORTRAN format (1X,I2,I3,I5,I6,2F9.2,E14.4,F9.2,E12.3,F8.2,F12.2,3F9.2,F9.4).

### The INI file (CLASS's initialization file) {#INI}

The .INI file is the initialization file for CLASS, and includes initialization data for soil moisture, temperature, and other related fields. Like the .MET file, the .INI file also uses fixed FORTRAN formats to read data, so caution must be taken when changing the .INI file fields. In addition, note that when the code is run with dynamic vegetation (when CTEM is switched on) several of the vegetation- related fields specified in the .INI file are replaced by those estimated by CTEM. These fields include:

   - Leaf area index (LAI)
   - Log of roughness length
   - Visible and near-infrared albedos
   - Canopy mass
   - Rooting depth

In addition, CTEM calculates stomatal resistance values so parameters related to CLASS’ stomatal resistance formulation and the associated parameters in the .INI file (RSMN, QA50, VPDA, VPDB, PSGA and PSGB) do not matter.

If CLASS is being run on its own, without CTEM, proper specification of vegetation related attributes and parameters is crucial.

Version 2.0 of CTEM adds the capability to divide each cell into several mosaics or ‘tiles’ as mentioned above. When running for multiple tiles, Both the CLASS and CTEM calculations are performed over each tile. The prognostic variables are saved for each tile in addition to grid averaged values.

A sample .INI file is found in the Benchmarks folder. The INI file has very speficic formatting. You can see it in runclass36ctem.f where the INI file is read and written (search for file units 10 and 100, respectively). Below is an explanation of the variables found in the INI file:


* ZRFM and ZRFH (m):
The reference height at which the climate variables (wind speed, temperature and specific humidity) are provided. If the model is driven by the atmospheric model forcing data, these heights would vary by time step. If the model is driven by field data (such as the sample file here), these heights refer to the measurement height of these variables.

* ZBLD (m):
Atmospheric blending height. Here it is assigned a value of 50 m.

* GC:
GCM surface description variable. For land surface (including inland water) it has a value of -1.

* Plant growth index:
CLASS uses specified maximum and minimum leaf area indices (LAIs) for each of its four plant functional types (PFTs). The LAI goes from this minimum to the maximum specified LAI at a specified rate and dependent on the soil temperature. Plant growth index determines what the LAI would be at initialization. This is little tricky if you aren’t familiar with CLASS, so it’s best left at zero. In addition, when coupled to CTEM this doesn’t matter because then CLASS uses LAI simulated by CTEM. When CTEM is switched on then all vegetation-related fields are estimated by CTEM. Competition between PFTs is now simulated and so if a user switches on this option then the fractional coverage of PFTs that grow in a grid cell is also simulated dynamically by CTEM. However, please read the discussion related to each of CTEM’s capabilities associated with each switch further down this user guide.

* NLTEST and NMTEST:
Number of grid cells and the number of mosaic tiles. NLTEST can not be greater than NLAT. Also, NMTEST must not be greater than NMOS. Typically this standalone code is run for 1 grid cell (NLTEST=1) with 1 or more mosaic tiles (NMTEST=1 or more).

* RSMN (s m-1), QA50 (W/m2), VPDA, VPDB, PSGA and PSGB:
Used in stomatal resistance calculation in CLASS subroutine CANALB. Note that these values do not matter when CTEM is switched on. Typical values for the four CLASS vegetation categories are listed above in section **Vegetation Data**.


* DRN:
Set it to 1 to allow drainage. Set it to 0 if drainage is suppressed at the bottom of the soil layer. Coupling of CLASS 3.6 to CTEM required a lot of tuning because it seems CLASS 3.6 is drier than CLASS 2.7. As a result a value of 0.1 has been lately used the CCCma Earth system model.

* SDEP (m):
Depth to bedrock in the soil profile


* XSLP, GRKF, WFSF, WFCI:
Parameters used when running MESH code, not used in Version 3.6 of CLASS.

* TBAR, THLQ and THIC:
Thin soil layers near the surface equilibrate quickly, but thicker, deeper layers respond more slowly. Long-term biases can be introduced into the simulation if the soil layers temperatures and moisture contents are not initialized accurately. For the moisture contents, it is better to err on the low side, since soil moisture recharge typically takes place on shorter scales than soil moisture loss. Very deep soil temperatures do not have a large effect on surface fluxes, but errors in their initial values can affect hydrological simulations. For rock or ice layers, THLQ and THIC should both be set to zero.

* RCAN (kg/m2):
Intercepted liquid water stored on canopy. RCAN can be initialized to zero.

* SCAN (kg/m2):
Intercepted frozen water stored on canopy
The vegetation canopy has a relatively small heat capacity and water storage capacity relative to the soil, so its temperature and intercepted water stores equilibrate quickly. TCAN can be initialized to the air temperature. SCAN can be initialized to zero.

* SNO (kg/m2):
Mass of snow pack.

* ALBS:
Snow albedo.

* RHOS (kg/m3):
Snow density

The above three variables and TSNO can be all initialized to 0 in snow-free conditions.

* TPND and ZPND:
Surface ponded water is a small term and is ephemeral in nature, so ZPNDROW and TPNDROW can both be initialized to zero.

* GRO:
Vegetation growth index, should be initialized to 1 during the growing season and to 0 otherwise. If CTEM is switched on, these values do not matter.

* DELZ and ZBOT:
The standard operational configuration for CLASS consists of three soil layers, of thicknesses 0.10 m, 0.25 m and 3.75 m, and thus of bottom depths 0.10, 0.35 and 4.10 m respectively. Version 3.6 supports other options: the third soil layer may be replaced with a larger number of thinner layers, and/or the bottom of the soil profile may be extended below 4.10 m. However, care must be taken not to make the soil layers too thin since this may lead to numerical instability.


### Typical values of vegetation-related fields for CLASS-only simulations {#classvals}

If you run CLASS without CTEM you will need to specify all vegetation-related parameters. The table below shows typical values you may want to use for different vegetation types.

\f[
\begin{tabular}{ | l | c | c | c | c | c | c | c | c | }
Biome & CLASS PFT & Visible albedo  & Near-infrared Albedo & Roughness length (m) \\
Needleleaf evergreen forest                 & 1 & 0.03 & 0.19 & 1.5  \\
Needleleaf deciduous forest                 & 1 & 0.03 & 0.19 & 1.0  \\
Broadleaf evergreen forest                  & 2 & 0.03 & 0.23 & 3.5  \\
Broadleaf cold deciduous forest             & 2 & 0.05 & 0.29 & 2.0  \\
Broadleaf tropical evergreen forest         & 2 & 0.03 & 0.23 & 3.0  \\
Broadleaf tropical drought deciduous forest & 2 & 0.05 & 0.29 & 0.8  \\
Broadleaf evergreen shrub                   & 4 & 0.03 & 0.19 & 0.05  \\
Broadleaf deciduous shrub                   & 2  & 0.05 & 0.29 & 0.15  \\
Broadleaf thorn shrub                       & 2 & 0.06 & 0.32 & 0.15  \\
Short grass and forbs                       & 4 & 0.06 & 0.34 & 0.02  \\
Long grass                                  & 4 & 0.05 & 0.31 & 0.08  \\
Arable                                      & 3 & 0.06 & 0.34 & 0.08  \\
Rice                                        & 3 & 0.06 & 0.36 & 0.08  \\
Sugarcane                                   & 3 & 0.05 & 0.31 & 0.35  \\
Maize                                       & 3 & 0.05 & 0.33 & 0.25  \\
Cotton                                      & 3 & 0.07 & 0.43 & 0.10  \\
Irrigated crop                              & 3 & 0.07 & 0.36 & 0.08  \\
Tundra                                      & 4 & 0.05 & 0.29 & 0.01  \\
Swamp                                       & 4 & 0.03 & 0.25 & 0.05  \\
Urban                                       & 5 & 0.09 & 0.15 & 1.35  \\
\end{tabular}
\f]


\f[
\begin{tabular}{ | l | c | c | c | c | c | c | c | c | }
Biome & CLASS PFT & Max. LAI (m2/m2) & Min. LAI (m2/m2) & Aboveground biomass (kg/m2) & Rooting depth (m) \\
Needleleaf evergreen forest                 & 1 & 2.0 & 1.6 & 25.0 & 1.0 \\
Needleleaf deciduous forest                 & 1 & 2.0 & 0.5 & 15.0 & 1.0 \\
Broadleaf evergreen forest                  & 2 & 8.0 & 8.0 & 50.0 & 5.0 \\
Broadleaf cold deciduous forest             & 6.0 & 0.5 & 20.0 & 2.0 \\
Broadleaf tropical evergreen forest         & 8.0 & 8.0 & 40.0 & 5.0 \\
Broadleaf tropical drought deciduous forest & 4.0 & 4.0 & 15.0 & 5.0 \\
Broadleaf evergreen shrub                   & 2.0 & 2.0 & 2.0 & 0.2 \\
Broadleaf deciduous shrub                   & 4.0 & 0.5 & 8.0 & 1.0 \\
Broadleaf thorn shrub                       & 3.0 & 3.0 & 8.0 & 5.0 \\
Short grass and forbs                       & 3.0 & 3.0 & 1.5 & 1.2 \\
Long grass                                  & 4.0 & 4.0 & 3.0 & 1.2 \\
Arable                                      & 4.0 & 0.0 & 2.0 & 1.2 \\
Rice                                        & 6.5 & 0.0 & 2.0 & 1.2 \\
Sugarcane                                   & 5.0 & 0.0 & 5.0 & 1.0 \\
Maize                                       & 4.0 & 0.0 & 5.0 & 1.5 \\
Cotton                                      & 5.0 & 0.0 & 2.0 & 2.0 \\
Irrigated crop                              & 4.0 & 0.0 & 2.0 & 5.0 \\
Tundra                                      & 1.5 & 1.5 & 0.2 & 0.1 \\
Swamp                                       & 1.5 & 1.5 & 1.0 & 5.0 \\
Urban                                       & - & - & - & - \\
\end{tabular}
\f]

### CLASS and CTEM PFTs {#classtoctem}

While CLASS models all physical processes related to energy and water balance for four PFTs, CTEM models its terrestrial ecosystem processes for nine PFTs. The nine PFTs of CTEM, however, line up with CLASS’ four PFTs perfectly as shown in Table 4. Needleleaf trees are divided into their evergreen and deciduous versions, broadleaf trees are divided into evergreen and cold and drought/dry deciduous versions, and crops and grasses are divided into their C3 and C4 versions based on their photosynthetic pathway.

\f[
\begin{tabular}{ | l | c | c | }
 & CTEM  PFTs & CLASS PFTs \\
 1. & Needleleaf Evergreen & Needleleaf \\
 2. & Needleleaf Deciduous & \\
 3. & Broadleaf Evergreen & Broadleaf \\
 4. & Broadleaf Cold Deciduous & \\
 5. & Broadleaf Drought/Dry Deciduous & \\
 6. & C3 Crop & Crops \\
 7. & C4 Crop & \\
 8. & C3 Grass & Grasses \\
 9. & C4 Grass & \\
\end{tabular}
\f]

<!-- start ignore of simple table example
\f[
\begin{tabular} { l c r }
  1 & 2 & 3 \\
  4 & 5 & 6 \\
  7 & 8 & 9 \\
\end{tabular}
\f]

Or can be done via the markdown syntax (http://www.stack.nl/~dimitri/doxygen/manual/markdown.html#md_toc) which is a bit easier.

end ignore -->

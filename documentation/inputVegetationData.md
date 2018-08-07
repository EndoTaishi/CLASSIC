# Input Vegetation Data {#vegetationData}

CLASSIC can be run with either dynamic vegetation (CTEM is turned on) or a physics only simulation (CLASS). The model inputs differ between the two simulations with some input vegetation data for a CLASS only run ignored when the CTEM is turned on. As well CTEM requires some additional inputs as described below.

# Required vegetation data for a physics only simulation (CLASS alone) {#vegCLASSonly}

The four main vegetation categories for the physics (CLASS) are needleleaf trees, broadleaf trees, crops and grass. Urban areas are also treated as “vegetation” in the CLASS code, and have associated values for FCANROT, ALVCROT, ALICROT and LNZ0ROT (see below). Thus these arrays have a third dimension of 5 rather than 4. For each of those the following data are required for each mosaic tile over each grid cell or modelled area (**NOTE**: When CTEM is turned on, i.e. dynamic vegetation is desired, the variables indicated in bold font are overwritten by CTEM during model run. FCANROT may be overwritten if land use change or if competition between PFTs is turned on.)

 1. ALICROT Average near-IR albedo of vegetation category when fully-leafed [ ]
 2. ALVCROT Average visible albedo of vegetation category when fully-leafed [ ]
 3. **CMASROT** Annual maximum canopy mass for vegetation category \f$[kg m^{-2} ]\f$
 4. **FCANROT** Annual maximum fractional coverage of modelled area [ ]
 5. **LNZ0ROT** Natural logarithm of maximum vegetation roughness length [ ]
 6. **PAMNROT** Annual minimum plant area index of vegetation category [ ]
 7. **PAMXROT** Annual maximum plant area index of vegetation category [ ]
 8. PSGAROT Soil moisture suction coefficient (used in stomatal resistance calculation) [ ]
 9. PSGBROT Soil moisture suction coefficient (used in stomatal resistance calculation) [ ]
 10. QA50ROT Reference value of incoming shortwave radiation (used in stomatal resistance calculation) \f$[W m^{-2} ]\f$
 11. **ROOTROT** Annual maximum rooting depth of vegetation category [m]
 12. RSMNROT Minimum stomatal resistance of vegetation category \f$[s m^{-1} ]\f$
 13. VPDAROT Vapour pressure deficit coefficient (used in stomatal resistance calculation) [ ]
 14. VPDBROT Vapour pressure deficit coefficient (used in stomatal resistance calculation) [ ]

In physics only runs (CLASS only), the vegetation is prescribed as follows (For full details of these calculations, see the documentation for subroutine src/APREP.f):
- CLASS models the physiological characteristics of trees as remaining constant throughout the year except for the leaf area index and plant area index, which vary seasonally between the limits defined by PAMXROT and PAMNROT.
- The areal coverage of crops varies from zero in the winter to FCANROT at the height of the growing season, and their physiological characteristics undergo a corresponding cycle.
- Grasses remain constant year-round.

Ideally the vegetation parameters should be measured at the modelled location. Of course this is not always possible, especially when running over a large modelling domain. As a guide, the table below provides generic values from the literature for the 20 categories of globally significant vegetation types. If more than one type of vegetation in a given category is present on the modelled area, the parameters for the category should be areally averaged over the vegetation types present.

\image html "landcovercat_table.png" ""
\image latex "landcovercat_table.png" ""

For the stomatal resistance parameters, typical values for the four principal vegetation types are given below:

\f[
\begin{tabular}{ | l | c | c | c | c | c | c | }
 & RSMN & QA50 & VPDA & VPDB & PSGA & PSGB \\
Needleleaf trees & 200.0 & 30.0 & 0.65 & 1.05 & 100.0 & 5.0 \\
Broadleaf trees & 125.0 & 40.0 & 0.50 & 0.60 & 100.0 & 5.0 \\
Crops & 85.0 & 30.0 & 0.50 & 1.00 & 100.0 & 5.0 \\
Grass & 100.0 & 30.0 & 0.50 & 1.00 & 100.0 & 5.0 \\
\end{tabular}
\f]

# Required vegetation data for a biogeochemical simulation (CTEM) {#vegCTEMtoo}

In addition to the CLASS variables described above, CTEM requires the following further variables:

COMBAK

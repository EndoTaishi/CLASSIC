# Overview of the Canadian Land Surface Scheme (CLASS) {#overviewCLASS}

The Canadian Land Surface Scheme, CLASS, was originally developed for use with the Canadian Global Climate Model or GCM (Verseghy, 1991 \cite Verseghy1991-635 ; Verseghy et al., 1993 \cite Verseghy1993-1ee ). The table at the end of this overview summarizes the development of CLASS from the late 1980’s onward.

The basic function of CLASS is to integrate the energy and water balances of the land surface forward in time from an initial starting point, making use of atmospheric forcing data to drive the simulation. When CLASS is run in coupled mode with a global or regional atmospheric model, the required forcing data are passed to it at each time step over each modeled grid cell from the atmospheric driver. CLASS then performs its internal calculations, evaluating a suite of prognostic and diagnostic variables such as albedo and surface radiative and turbulent fluxes, which are in turn passed back to the driver. CLASS can also be run in uncoupled or offline mode, using forcing data derived from field measurements, and the output values of its prognostic and diagnostic variables can then be validated against observations.

CLASS models separately the energy and water balances of the soil, snow, and vegetation canopy (see the diagram below). The basic prognostic variables consist of the temperatures and the liquid and frozen moisture contents of the soil layers; the mass, temperature, density, albedo and liquid water content of the snow pack; the temperature of the vegetation canopy and the mass of intercepted rain and snow present on it; the temperature and depth of ponded water on the soil surface; and an empirical vegetation growth index. These variables must be initialized, and a set of physical parameters describing the soil and vegetation existing on the modelled area must be assigned background values, at the beginning of the simulation (see @ref forcingData).

At each time step, CLASS calculates the bulk characteristics of the vegetation canopy on the basis of the vegetation types present over the modelled area. In a pre-processing step, each vegetation type is assigned representative values of parameters such as albedo, roughness length, annual maximum and minimum plant area index, rooting depth and so on (see @ref initProgVar). These values are then aggregated over four main vegetation categories identified by CLASS usually: needleleaf trees, broadleaf trees, crops, and grass (i.e. short vegetation). The physiological characteristics of the vegetation in each category are determined at the current time step using the aggregated background parameters and assumed annual or diurnal variation functions. These physiological characteristics are then aggregated to produce the bulk canopy characteristics for the current time step.

\image html "schematicDiagramOfClass.png" "Schematic Diagram Of CLASS"
\image latex "schematicDiagramOfClass.png" "Schematic Diagram Of CLASS"

In performing the surface flux calculations the modeled area is divided into up to four subareas: bare soil, vegetation over soil, snow over bare soil, and vegetation over snow. The fractional snow coverage is determined using the concept of a threshold snow depth. If the calculated snow depth is less than this value, the snow depth is set to the threshold value and the fractional snow cover is calculated on the basis of conservation of snow mass. The fluxes are calculated for each of the four subareas, and these and the prognostic variables are then areally averaged before being passed back to the atmospheric model.

Originally CLASS performed only one set of these calculations for each grid cell of the model domain. In more recent versions, a “mosaic” option has been added to handle sub-grid scale heterogeneity more effectively. When this option is utilized, each grid cell is divided into a user-specified number of mosaic “tiles”, and the CLASS calculations are performed in turn over each. The surface fluxes are averaged, but the prognostic variables are kept separate for each of the tiles of the mosaic between time steps (see [here for more](@ref compvsmosaic)).

In the CLASSIC offline driver, a gather-scatter operation is included in the driver, mimicking the practice in atmospheric models of “gathering” land surface points on latitude circles onto long vectors prior to the calculations (e.g. src/CLASSG.f or src/ctemg2.f), for improved computational efficiency on vector supercomputers. For CLASS, the mosaic tiles on each of the modelled grid cells are “gathered” onto long arrays prior to calling the CLASS subroutines (thus collapsing the first two dimensions of the arrays into one), and subsequently “scattered” back onto the grid cells before performing the diagnostic averaging calculations (e.g. src/CLASSS.f or src/ctems2.f).

## Development history of CLASS {#devHistory}

\f[
\begin{tabular}{ | l | l || l | }
1.0 & April 1989 & Basic thermal and hydrological model of snow and soil. \\
2.0 & August 1991 & Addition of vegetation thermal and hydrological model. \\
2.1 & May 1993 & Full vectorization of code to enable efficienr running on vector supercomputers. \\
2.2 & April 1994 & Augmentation of diagnostic calculations; incorporation of in-line comments throughout; development \\
    &            & of a parallel stand-alone version of the model for use with field data. \\
2.3 & December 1994 & Revisions to diagnostic calculations; new near-surface atmospheric stability functions. \\
2.4 & August 1995 & Complete set of water budget diagnostic calculations; parametrizations of organic soils and \\
    &            &rock soils; allowance for inhomegeneity between soil layers; incorporation of variable \\
    & & surface detention capacity. \\
2.5 & January 1996 & Completion of energy budget diagnostic calculations. \\
2.6 & August 1997 & Revisions to surface stability function calculations. \\
2.7 & December 1997 & Incorporation of variable soil permeable depth; calculation of soil thermal and hydraulic properties based on textural
composition; modified surface temperature iteration scheme. \\
3.0 & December 2002 & Improved treatment of soil evaporation; complete treatment of organic soils; new canopy conductance formulation; preliminary routines for lateral movement of soil water; enhanced snow density and snow interception; improved turbulent transfer from vegetation; mosaic formulation. \\
3.1 & April 2005 & Faster surface temperature iteration scheme; refinements to leaf boundary resistance formulation; improved treatment of snow sublimation and interception; transition to Fortran 90 and single precision variables. \\
3.2 & May 2006 & Option for multiple soil layers at depth; additional liquid water content of snow pack; revised radiation transmission in vegetation. \\
3.3 & December 2006 & Separate temperature profile curve fit for snow and soil; multiple-layer option for ice sheets; water and energy balance checks for each time step; modifications to soil hydraulic conductivity calculations. \\
3.4 & April 2008 & Streamline and clean up code; updated soil thermal conductivity calculations; revisions to handling of water stored on vegetation. \\
3.5 & December 2010 & Updated field capacity calculation; revised treatment of water on canopy; reworked calculation of baseflow. \\
3.6 & December 2011 & Revised ponding depth over organic soils; revised snow albedo refreshment threshold; new snow thermal conductivity algorithm; interface with Canadian Terrestrial Ecosystem Model (CTEM). \\
3.6.1 & December 2016 & New treatment of bare soil albedo; new optional four-band snow albedo formulation; fixes to guard against overshoots in water drawdown by evapotranspiration; upper limit on snow depth. \\
\end{tabular}
\f]

# Overview of the Canadian Terrestrial Ecosystem Model (CTEM) {#overviewCTEM}

Version 1 of the CTEM is the terrestrial carbon cycle component of the second generation Canadian Earth System Model (CanESM2) (Arora et al., 2011)\cite Arora2011-79f where it is coupled to version 2.7 of the Canadian Land Surface Scheme (CLASS). CTEM v. 2.0 (Melton and Arora, 2016) \cite Melton2016-zx was coupled to CLASS v. 3.6 (Verseghy, 2012) \cite Verseghy2012-c0e. The coupled CLASS--CTEM model is capable of being run online in the CanESM model or offline, driven by observation-based meteorological forcings. CTEM models terrestrial ecosystem processes for nine PFTs, two of which are crop PFTs (see table below), by tracking the flow of carbon through three living vegetation components (leaves, stem and roots) and two dead carbon pools (litter and soil).

| CLASS PFTs | CTEM PFTs ||| Peatland PFTs ||
|:---------------:|:---------:|----------------|----------------------|:----------------:|------------------|
| Needleleaf tree | Evergreen | Deciduous |  |  |  |
| Broadleaf tree | Evergreen | Cold Deciduous | Drought/Dry Decidous | Evergreen Shrubs | Deciduous Shrubs |
| Crop | \f$C_3\f$ | \f$C_4\f$ |  |  |  |
| Grass | \f$C_3\f$ | \f$C_4\f$ |  | Sedges |  |
<!-- \f[
\begin{table}[]
\caption{CTEM and peatland PFTs and their mapping to the CLASS PFTs}
\label{my-label}
\begin{tabular}{|l|l|l|l|l|l|}
\hline
\multicolumn{1}{|c|}{CLASS PFTs} & \multicolumn{3}{c|}{CTEM PFTs}                    & \multicolumn{2}{c|}{Peatland PFTs}  \\ \hline
Needleleaf tree                  & Evergreen & Deciduous      &                      &                  &                  \\ \hline
Broadleaf tree                   & Evergreen & Cold Deciduous & Drought/Dry Decidous & Evergreen Shrubs & Deciduous Shrubs \\ \hline
Crop                             & C$_3$     & C$_4$          &                      &                  &                  \\ \hline
Grass                            & C$_3$     & C$_4$          &                      & Sedges           &                  \\ \hline
\end{tabular}
\end{table}
\f] -->

The amount of carbon in these five carbon pools is simulated prognostically (see below). In the CLASSIC framework, CLASS uses structural vegetation attributes (including LAI, vegetation height, canopy mass and rooting depth) simulated by CTEM, and CTEM uses soil moisture, soil temperature and net radiation calculated by CLASS. Combined, CLASS and CTEM simulate the atmosphere--land fluxes of energy, water and \f$CO_2\f$.

Version 1.0 of CTEM is described in a collection of papers detailing parametrization of photosynthesis, autotrophic and heterotrophic respiration (Arora, 2003) \cite Arora2003-3b7; phenology, carbon allocation, biomass turnover and conversion of biomass to structural attributes (Arora and Boer, 2005) \cite Arora2005-6b1; dynamic root distribution (Arora and Boer, 2003) \cite Arora2003838; and disturbance (fire) (Arora and Boer, 2005) \cite Arora20052ac. These processes are modelled over prescribed fractional coverage of (typically) nine PFTs (Wang et al., 2006) \cite Wang2006-he and determine the structural vegetation dynamics including vegetation biomass, LAI, vegetation height, fraction of roots in each of the three soil layers, leaf onset and offset times and primary \f$CO_2\f$ fluxes of gross primary productivity (GPP) and NPP. CTEM v. 2.0 is described in Melton and Arora, 2016 \cite Melton2016-zx.

## Rate change equations for carbon pools {#CTEMRateChgEqns}

From the gross canopy photosynthesis rate (\f$G_{\text{canopy}}\f$, src/PHTSYN3.f), maintenance and growth respirations (\f$R_\mathrm{m}\f$ and \f$R_\mathrm{g}\f$, src/mainres.f), and
heterotrophic respiration components (\f$R_{\text{h,H}}\f$ and \f$R_{\text{h,D}}\f$, src/hetres_mod.f90), it is possible to estimate the change in carbon amount of the model's five pools.

When the daily NPP (\f$ G_{canopy} - R_\mathrm{m} - R_\mathrm{g}\f$) is positive, carbon is allocated to the plant's live carbon pools and the rate of change is given by

  \f[
  \frac{\mathrm{d}C_i}{\mathrm{d}t} = a_{fi} \left(G_{canopy}-R_\mathrm{m}-R_\mathrm{g} \right) - D_i - H_i - M_i \\ \quad i = {L, S, R}
   \f]
   <!-- {#rate_change_eqns_live_pools} -->

where \f$a_{fi}\f$ is the corresponding allocation fractions for each pool (stem, root and leaves) and \f$D_i\f$ is the litter produced from these components as explained in src/phenolgy.f90. \f$H_i\f$ is the loss associated with fire that releases \f\chem{CO_2}\f and other trace gases to the atmosphere and \f$M_i\f$ is the mortality associated with fire that contributes to the litter pool as explained in src/disturb.f90.

If the daily NPP is negative (\f$G_{canopy} < R_\mathrm{m}\f$, \f$R_\mathrm{g} = 0\f$), the rate of change is given by

\f[
 \frac{\mathrm{d}C_i}{\mathrm{d}t} = a_{fi}G_{canopy} - R_{m,i}  - D_i  - H_i - M_i, \\ \quad i = {L, S, R}
 \f]
<!-- \label{rate_change_eqns_live_pools2} -->

Negative NPP causes the plant to lose carbon from its live carbon pools due to respiratory costs in addition to the losses due to litter production (\f$D_i\f$) and disturbance (\f$H_i\f$, \f$M_i\f$).

The rate change equations for the litter and soil carbon pools are given by
\f{eqnarray*}{
\frac{\mathrm{d}C_\mathrm{D}}{\mathrm{d}t} &=& D_\mathrm{L} + D_\mathrm{S} +
D_\mathrm{R} + M_\mathrm{L} + M_\mathrm{R} + M_\mathrm{S} - H_\mathrm{D} -C_{\mathrm{D} \rightarrow \mathrm{H}} - R_{h,D} \\
\frac{\mathrm{d}C_\mathrm{H}}{\mathrm{d}t} &=& C_{\mathrm{D} \rightarrow
\mathrm{H}} - R_{h,H}
\f}
<!-- \label{rate_change_eqns_dead_pools}, -->

where \f$C_{\mathrm{D} \rightarrow \mathrm{H}}\f$ represents the transfer of
humified litter to the soil carbon pool and \f$H_\mathrm{D}\f$
is loss associated with burning of litter associated with fire that releases
\f\chem{CO_2}\f and other trace gases to the atmosphere.


A parametrization for competition between PFTs in an earlier version of CTEM is described by \cite Arora2006-pp \cite Arora2006-ax where it was evaluated at select locations. Here we present CTEM v. 2.0, which builds upon the model framework of CTEM v. 1.0 and can be run in two different modes, either (i) using specified fractional coverage of its nine PFTs, or (ii) allowing the fractional coverage of its seven non-crop PFTs to be dynamically determined based on competition between PFTs. The parametrization for simulating competition between PFTs is summarized in Sect. \ref{compmain}. The fire parametrization has also been refined in the new model version as described in Appendix \ref{fire}. The CLASS--CTEM modelling framework has the capability of representing the sub-grid scale variability of PFTs using either a composite or a mosaic configuration \cite Li2012-f7f \cite Melton2014-xk. In the composite (or single tile) configuration, the vegetation attributes for all PFTs present in a grid cell are averaged and used in energy and water balance calculations that determine the physical land surface conditions including soil moisture, soil temperature and thickness and fractional coverage of snow (if present). In the mosaic (or multi-tile) configuration each PFT is allocated its own tile for which separate energy and water balance calculations are performed. As a result, the simulated carbon balance evolves somewhat differently in the two configurations despite being driven with identical climate forcing (see \cite Melton2014-xk). The results presented in this paper are obtained using the composite configuration.

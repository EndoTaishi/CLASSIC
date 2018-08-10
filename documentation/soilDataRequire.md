# Input Soil Data {#soilData}

The following specifications are required for each modelled soil layer:

- DELZ Layer thickness [m]

The standard operational configuration for CLASS consists of three soil layers, of thicknesses 0.10 m, 0.25 m and 3.75 m, and thus of bottom depths 0.10, 0.35 and 4.10 m respectively. CLASSIC supports other options: the second and third soil layers may be replaced with a larger number of thinner layers, and/or the bottom of the soil profile may be extended below 4.10 m. However, because the temperature stepping scheme used in CLASS is of an explicit formulation, care must be taken not to make the layers too thin, since this may lead to numerical instability problems. As a rule of thumb, the thicknesses of layers should be limited to \f$\geq\f$ 0.10 m.

For each of the modelled soil layers on each of the mosaic tiles, the following texture data are required:

- CLAYROT Percentage clay content
- ORGMROT Percentage organic matter content
- SANDROT Percentage sand content

\image html "percentSand.png" "Percent Sand"
\image latex "percentSand.png" "Percent Sand"

- For mineral soils, the percentages of sand, clay and organic matter content need not add up to 100%, since the residual is assigned to silt content. If the exact sand, clay and organic matter contents are not known, estimates can be made for the general soil type on the basis of the standard USDA texture triangle shown above. Organic matter contents in mineral soils are typically not more than a few percent.
- If the layer consists of rock, SANDROT is assigned a flag value of -3. If it is part of a continental ice sheet, it is assigned a flag value of -4. In both cases, CLAYROT and ORGMROT are not used and are set to zero.
- Highly organic soils have different behaviour if the area is being modelled as a peatland:
  - If the peatland flag is 0, that is the tile is not being modelled as a peatland, and the soil layer is a fully organic one, SANDROT, CLAYROT and ORGMROT are used differently. The sand content is assigned a flag value of -2, and the organic matter content may be assigned a flag value of 1, 2 or 3 depending on whether the peat texture is fibric, hemic or sapric (see Letts et al. (2000) \cite Letts2000-pg). The current default is for the first layer to be assumed as fibric, the second as hemic and any lower layers as sapric. CLAYROT is not used and is set to zero.
  - If the tile is being treated as a peatland then the first soil layer is considered moss following Wu et al. (2016) \cite Wu2016-zt. The lower soil layers are treated such that the lower layers are assigned fibric, hemic or sapric characteristics (see @ref CLASSB.f)


SANDROT, CLAYROT and ORGMROT are utilized in the calculation of the soil layer thermal and hydraulic properties in @ref CLASSB.f. If measured values of these properties are available, they could be used instead (with modificiations to the code).

For each of the mosaic tiles over the modelled area, the following surface parameters must be specified:

- DRNROT Soil drainage index
  - The drainage index, DRNROT, is usually set to 0.005 except in cases of deep soils where it is desired to suppress drainage from the bottom of the soil profile (e.g. in bogs, or in deep soils with a high water table). In such cases it is set to 0.
- FAREROT Fractional coverage of mosaic tile on the modelled area (also discussed [here](@ref compvsmosaic))
- MIDROT Mosaic tile type identifier (1 for land surface, 0 for inland lake) (also discussed [here](@ref compvsmosaic))
- SDEPROT Soil permeable depth [m]
  - The soil permeable depth, i.e. the depth to bedrock, may be less than the modelled thermal depth of the soil profile. If the depth to bedrock occurs within a soil layer, CLASS assigns the specified mineral or organic soil characteristics to the part of the layer above bedrock, and values corresponding to rock to the portion below. All layers fully below SDEPROT are treated as rock.
- SOCIROT Soil colour index
  - The soil colour index is used to assign the soil albedo.  It ranges from 1 to 20; low values indicate bright soils and high values indicate dark (see Lawrence and Chase (2007) \cite Lawrence2007-bc).  The wet and dry visible and near-infrared albedos for the given index are obtained from lookup tables, which can be found in @ref CLASSB.f.

CLASS provides a means of accounting for the possibility of the depth to bedrock falling within a layer, and therefore of phase changes of water taking place in only the upper part of the layer, by introducing the variable TBASROT, which refers to the temperature of the lower part of the layer containing the bedrock. At the beginning of the time step the temperature of the upper part of the layer is disaggregated from the overall average layer temperature using the saved value of TBASROT. The heat flow between the upper part of the soil layer and the lower part is diagnosed from the heat flux at the top of the layer. The upper layer temperature and TBASROT are stepped ahead separately, and the net heat flux in the upper part of the layer is used in the phase change of water if appropriate. The upper layer temperature and TBASROT are re-aggregated at the end of the time step to yield once again the overall average layer temperature.

Two variables, assumed to be constant over the grid cell, are provided if required for atmospheric model runs:

- GGEOROW Geothermal heat flux \f$[W m^{-2} ]\f$
  - Unless the soil depth is very large and/or the run is very long, the geothermal heat flux can be set to zero. Since this is rarely used, this is not presently read in from the initialization file.
- Z0ORROW Orographic roughness length [m]
  - Z0ORROW is the surface roughness length representing the contribution of orography or other terrain effects to the overall roughness, which becomes important when the modelled grid cell is very large (e.g. in a GCM). For field studies it can be set to zero. It is presently not required as input for a model run and is set to zero in @ref model_state_drivers.read_initialstate

Four parameters are required for modelling lateral movement of soil water: GRKFROT, WFCIROT, WFSFROT and XSLPROT. However, the routines for interflow and streamflow modelling are not implemented in this version of CLASSIC, so unless the user is involved in this development, these parameters can be set to arbitrary values, since they will not be used.

- grclarea Area of grid cell \f$[km^{2} ]\f$

Area of the grid cell is required if CTEM is on. It is used for calculations relating to fire (@ref disturbance_scheme.disturb) and land use change (@ref landuse_change.luc). Both of those subroutines are not generally used (or meaningful) at the point scale so the grclarea can be set to an arbitrary number like 100 \f$km^{2}\f$.

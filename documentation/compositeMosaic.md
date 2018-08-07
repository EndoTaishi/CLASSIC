# Notes about composite versus mosaic running of the model {#compvsmosaic}

In the composite mode, the structural vegetation attributes (including leaf area index, vegetation height, rooting depth) of PFTs that exist in a grid cell are averaged in proportion to their fractional coverages and then used in the grid-averaged energy and water balance calculations. As a result the entire grid cell is characterized by land surface physical environment (including soil temperature, soil moisture, fractional snow cover, and net radiation) that is common to all PFTs. In contrast, in the mosaic mode a grid box is split into multiple tiles representing individual PFTs (or land unit, such as soil texture e.g. Melton et al. (2017) \cite Melton2017-gp) for each of which energy, water and carbon balance calculations are performed separately. The mosaic mode is, however, able to represent all forms of sub-grid scale variabilities. In principle, the mosaic mode may be used to represent tiles that are characterized by any chosen distinction such a lowlands vs. uplands, soil texture, vegetation, soil depth, etc.

Some variables that impact upon composite vs. mosaic model runs:

- FAREROT Fractional coverage of mosaic tile on the modelled area
- MIDROT Mosaic tile type identifier (1 for land surface, 0 for inland lake)

COMBAK


Although pretty clever and powerful, running the mosaic version and interpreting the model results can be a logical nightmare, especially with competition on where the fractions of different tiles/mosaic change with time. So if you are new at this, please consider running the model in the composite mode.

For more information see Melton and Arora (2014) \cite Melton2014-xk, Shrestha et al. (2016) \cite Shrestha2016-do, Melton et al. (2017) \cite Melton2017-gp.

\image html "compVsMosaic_MeltonArora_BG_2014.png" "Schematic representation of the composite and mosaic ap- proaches for the coupling of CLASS v 3.6 and CTEM v 1.2 models in a stand-alone mode. (From Melton and Arora, 2014)"
\image latex "compVsMosaic_MeltonArora_BG_2014.png" "Schematic representation of the composite and mosaic ap- proaches for the coupling of CLASS v 3.6 and CTEM v 1.2 models in a stand-alone mode. (From Melton and Arora, 2014)"

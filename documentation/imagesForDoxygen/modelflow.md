graph TD

  subgraph 
    contain(Singularity Container)
    style contain fill:#ffffdd,stroke:#ffffdd
    classic["<br/><br/><br/>   &ensp; &ensp;    &ensp; &ensp; CLASSIC  &ensp; &ensp; &ensp;      <br/> <br/> <br/> <br/>  "]
    style classic fill:#f9f,stroke:#333,stroke-width:6px
    moddep(Dependencies:<br/>make, libnetcdff-dev,<br/> git, gfortran,<br/> netcdf-bin)
  end

  subgraph Configuration Files
    nml[CLASSIC parameters<br/>namelist] --> classic
    jobo[Simulation joboptions] --> classic
    outxml[Output xml] --> classic
  end

  vared[OutputVariableEditor] -.-> outxml

  subgraph Supplementary info
    supl[bibliography,<br/>figures,<br/>documentation/*.md]
  end

  classic -.-> supl
  supl -.Doxygen.-> doc["Model<br/>Documentation<br/>(html,pdf)"]

  subgraph NetCDF Output
    classic -.-> half-hourly
    classic -.-> daily
    classic -.-> monthly
    classic --> annually
  end

  subgraph Model Geophysical Inputs
    met["Meteorological,<br/>{Land Cover Changes,<br/>GHGs,<br/>Population density,<br/>Lightning density,<br/>Peatland distribution,<br/>Wetland distribution}"] --> classic
  end

  init[NetCDF Initial<br/>Conditions File] --> classic

  classic --> rs[NetCDF Model<br/>restart file]

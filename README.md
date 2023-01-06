# A cross species map of neutrophil inflammatory responses

## Structure of this repository

    .
    ├── data                                 # Contains all data needed for the analysis.
    │   ├── metadata                            # metadata, part of this repo
    │   ├── processed                           # generated, content not part of vcs
    │   └── raw_data                            # downloaded, content not part of vcs
    ├── figures                              # Container for generated figures.
    ├── scripts                              # Contains all scripts and utils required.
    │   ├── experimental_validation             # code relating to the validation experiments
    │   ├── generate_figures.R                  # the script that can be run to generate the manuscript's figures
    │   ├── haemopedia                          # code pertaining the haemopedia analysis
    │   ├── inflammatory_meta                   # code to analyze the inflammatory dataset
    │   ├── utils                               # various utils for data retrieval, plotting, ...
    │   └── zymosan                             # zymosan analysis scripts
    ├── README.md
    ├── inflammatory_metaanlysis.Rproj
    └── .gitignore

## Figures

This is how figure generation works:

```mermaid
flowchart TD
    g1[Rmd files] -. calls build function .-> g2{figure_builder.R}
    g2 -- returns plot object --> g1
    g3[generate_figures.R] -. calls build function .-> g2
    g2 -- returns plot object --> g3
    g3 -- arranges subfigures --> g5(final PDF figure file)
    g4[config.R] -- provides palettes/... --> g2
    g1 -- generate --> g6[(Result data files)]
    g6 -- loaded by --> g2
```

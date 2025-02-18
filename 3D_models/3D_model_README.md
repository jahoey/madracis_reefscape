This directory contains data for colonies annotated in 3D photogrammetry models of reefs in Curaçao. Colony annotations and the 3D model are used to calculate various microhabitat measurements for each annotated colony. Many of these colonies were also genotyped and this process is described in the `/WGS` directory.

`/figures`: contains the microhabitat plots

`/input`: contains the colony-level annotation files for each plot, including non-genotyped colonies & some poorly-annotated colonies (filtered out later)

`/output`: contains the output after calculating the microhabitat measurements for each annotated colony

`/scripts`: contains the scripts for  visualizations
  * `filter_measures_for_quality.R`: filters measures output from unicorn package to only those whose annotations were actually on the coral
  * `plot_measures.R`: plots the quality-filtered measures output according to site, depth, and species & runs ANOVAs

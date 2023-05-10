# Phylogenetic_piepeline

## Table of contents
* [General info](#general-info)
* [Technologies](#technologies)
* [Quick Start](#quick-start)

## General info
This project allows you to automate steps of phylogentic analysis such as proteomes downloading from UniprotKB, Proteins Clustering, MSA, Maximum Likelihood Trees construction, Consensus Tree and Super Tree construction.  

	
## Technologies
Project is created with:
* [Clustal Omega](http://www.clustal.org/omega/)
* [IQTREE](http://www.iqtree.org/)
* [CD-HIT](https://github.com/weizhongli/cdhit)
* [Clann](http://chriscreevey.github.io/clann/)
* BioPython 1.81
	
## Quick Start
If all required execs are in PATH you can simply run:
```
$ python phylo_pipeline -t taxon_name --og out_group_organism_name
```

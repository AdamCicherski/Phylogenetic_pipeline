# Phylogenetic_piepeline

## Table of contents
* [General info](#general-info)
* [Technologies](#technologies)
* [Quick Start](#quick-start)

## General info
This script allows you to automate steps of phylogenetic analysis such as:
* proteome downloading from UniProtKB
*  Proteins Clustering 
*  MSA of each cluster in two variants
	* Only one protein from each organism
  	* Many protein from the same organism (paralogous proteins)   
*   Maximum Likelihood Tree construction for both variants of alignment
*   Consensus Tree construction in following variants:
  	* From all ML-Trees
  	* Only from trees with high bootstrap confidence level. 
*   Super Tree construction

	
## Technologies
The script requires:
* [Clustal Omega](http://www.clustal.org/omega/)
* [IQTREE](http://www.iqtree.org/)
* [CD-HIT](https://github.com/weizhongli/cdhit)
* [Clann](http://chriscreevey.github.io/clann/)
* [BioPython](https://biopython.org/)
	
## Quick Start
If all required execs are in PATH you can simply run:
```
$ python phylo_pipeline -t taxon_name -o output_dir --og out_group_organism_name
```

# MAGNET (MAny GeNE Trees) v0.1.4  :deciduous_tree::deciduous_tree::deciduous_tree:
Shell script pipeline for inferring ML gene trees for many loci (e.g. genomic SNP data)

## LICENSE

All code within the PIrANHA repository, including MAGNET v0.1.4 pipeline code, is available "AS IS" under a generous GNU license. See the [LICENSE](LICENSE) file for more information.

## CITATION

If you use scripts from this repository as part of your published research, then I require you to cite the PIrANHA repository and/or MAGNET package as follows (also see DOI information below):

  Bagley, J.C. 2017. PIrANHA. GitHub repository, Available at: http://github.com/justincbagley/PIrANHA/.
  
  Bagley, J.C. 2017. MAGNET. GitHub package, Available at: http://github.com/justincbagley/MAGNET. 

Alternatively, please provide the following link to this software program in your manuscript:

  http://github.com/justincbagley/MAGNET
  
**Example citations using the above URL:** 
	"We estimated a gene tree for each SNP locus in RAxML v8 (Stamatakis 2014) using 
	the MAGNET pipeline in the PIrANHA github repository (http://github.com/justincbagley/MAGNET).
	Each RAxML run specified the GTRGAMMA model and coestimated the maximum-likelihood phylogeny
	and bootstrap proportions from 500 bootstrap pseudoreplicates."

## DOI

The DOI for MAGNET, via Zenodo, is as follows:  [![DOI](https://zenodo.org/badge/66839898.svg)](https://zenodo.org/badge/latestdoi/66839898). Here is an example of citing MAGNET using the DOI: 
  
  Bagley, J.C. 2017. MAGNET. GitHub package, Available at: http://doi.org/10.5281/zenodo.166024.

## INTRODUCTION

The estimation of species-level phylogenies, or "species trees" is a fundamental goal in evolutionary biology. However, while "gene trees" estimated from different loci provide insight into the varying evolutionary histories of different parts of the genome, gene trees are random realizations of a stochastic evolutionary process. Thus, gene trees often exhibit conflicting topologies, being incongruent with each other and incongruent with the underlying species tree due to a variety of genetic and biological processes (e.g. gene flow, incomplete lineage sorting, introgression, selection). 

With the advent of recent advances in DNA sequencing technologies, biologists now commonly sequence data from multiple loci, and even hundreds to thousands of loci can quickly be sequenced using massively parallel sequencing on NGS sequencing platforms. Full-likelihood or Bayesian algorithms for inferring species trees and population-level parameters from multiple loci, such as \*BEAST and SNAPP, are computationally burdensome and may be difficult to apply to large amounts of data or distantly related taxa (or other cases that complicate obtaining MCMC convergence). By contrast, a number of recently developed and widely used "summary-statistics" approaches rely on sets of gene trees to infer a species tree for a set of taxa (reviewed by Chifman and Kubatko, 2014; Mirarab and Warnow, 2015). These methods are specifically designed to estimate gene trees or use gene trees input by the user, which are treated as observed data points analyzed in a distance-based or coalescent algorithm. Moreover, summary-statistics approaches to species tree inference tend to be accurate and typically much faster than full-data approaches (e.g. Mirarab et al., 2014;  Chifman and Kubatko, 2014). Examples of species tree software in this category include programs such as BUCKy (Larget et al., 2010), STEM (Liu et al., 2010), spedeSTEM, NJst (Liu and Yu, 2011), ASTRAL and ASTRAL-II (Mirarab and Warnow, 2015), and ASTRID (Vachaspati and Warnow, 2015). Phylogenetic network models implemented in recent software like SplitsTree and SNaQ also improve network and inference by analyzing sets of gene trees. 

Despite the importance of gene trees in species tree and network inference, few resources have been specifically designed to aid rapid estimation of gene trees for different loci. MAGNET (MAny GeNE Trees) is a shell script pipeline within the PIrANHA (PhylogenetIcs ANd PHylogeogrAphy) github repository (https://github.com/justincbagley/PIrANHA) that helps fill this gap by automating extracting individual loci from a multilocus sequence alignment file and inferring a maximum-likelihood (ML) gene tree for each locus. The MAGNET package was originally coded up to aid analyses of SNP loci generated by massively parallel sequencing of ddRAD-seq genomic libraries (Peterson et al. 2012) of freshwater fishes. However, MAGNET can be used to estimate gene trees for loci in other multilocus data types with the appropriate format using conversion scripts provided within the package (see below). 

## HARDWARE AND SETUP

:computer: MAGNET focuses on allowing users to automate the workflow necessary for quickly estimating many gene trees for many loci on their local machine. 

:thumbsup: No special hardware or setup is necessary, unless the user is interested in estimating gene trees on a remote supercomputing cluster. In that case, the user is referred to the c-MAGNET or "cluster MAGNET" software repository (https://github.com/justincbagley/c-MAGNET/; under development).


## SOFTWARE DEPENDENCIES

MAGNET v0.1.4 is composed of shell, R, and Perl scripts and also calls several software programs; thus, it relies on several software dependencies. These dependencies are described in some detail in README files for different scripts in the package. However, here I provide a list of them, with asterisks preceding those already included in the MAGNET distribution:

- Perl (available at: https://www.perl.org/get.html).
- \*Nayoki Takebayashi's file conversion Perl scripts (available at: http://raven.iab.alaska.edu/~ntakebay/teaching/programming/perl-scripts/perl-scripts.html).
- Python (available at: https://www.python.org/downloads/).
- bioscripts.convert v0.4 Python package (available at: https://pypi.python.org/pypi/bioscripts.convert/0.4; also see README for "NEXUS2gphocs.sh").
- RAxML, installed and running on local machine (available at: http://sco.h-its.org/exelixis/web/software/raxml/index.html).

Users must install all software not included in MAGNET, and ensure that it is available via the command line on their local machine. On the user's local machine, Perl should be available by simply typing "Perl" at the command line; Python should be available by typing "python" at the command line; and bioscripts.convert package should be available by typing "convbioseq" at the command line. Also, RAxML should be compiled using SSE3 install commands, so that RAxML can be called by simply typing "raxmlHPC-SSE3" on the command line. For detailed instructions for setting up RAxML this way, refer to the newest RAxML user manual (available at: http://sco.h-its.org/exelixis/resource/download/NewManual.pdf).

## INPUT FILE FORMAT

MAGNET assumes that you are starting from multilocus DNA sequence data in a single datafile in G-Phocs (Gronau et al. 2011) format, with the extension ".gphocs", or in NEXUS format with the extension ".nex". For genomic data such as RAD tags or other SNP data derived from genotyping-by-sequencing (GBS)-type methods, it is recommended that the user assemble the data, call SNPs, and output SNP data files in various formats including .gphocs format in pyRAD or ipyrad (Eaton 2014) prior to using MAGNET. However, this may not always be possible, and .gphocs format is not yet among the most popular file formats in phylogenomics/population genomics. Thus, I have added a "NEXUS2gphocs.sh" shell script utility within MAGNET (in the "shell" folder) that will convert a sequential NEXUS file into .gphocs format for you. An example NEXUS file "example.nex" is included in the distribution.

Feel free to use the NEXUS2gphocs.sh utility script independently of MAGNET to convert from .gphocs to NEXUS format. However, when doing this, _make sure to follow the usage guidelines below_.

## PIPELINE

Apart from input file conversion steps, the MAGNET pipeline works by calling three different scripts, in series, each designed to conduct a task that yields output that is processed in the next step of the pipeline. In STEP #1, the "gphocs2multiPhylip.sh" shell script is used to extract loci from the input file and place each locus in a Phylip-formatted file with extension ".phy". In STEP #2, a shell script named "MultiRAxMLPrepper.sh" is used to place the .phy files into separate folders ("run folders"), and prepare them to be run in RAxML. In STEP #3, a script named "RAxMLRunner.sh" is called to run RAxML on the contents of each run folder. In a "clean-up" step, MAGNET moves all .phy files files remaining in the working directory after STEP #3 to a new folder, "phylip_files", that is created in the working directory.

After running the MAGNET pipeline, the shell script "getGeneTrees.sh" automates post-processing of the output, including organizing all inferred gene trees into a single "gene_trees" folder in the working directory, and combining the individual 'best' gene trees resulting from each run into a single file named "besttrees.tre".

## USAGE

Additional input file and usage info is available in the usage or help texts. To get regular usage info for MAGNET, type ```$ ./MAGNET.sh```, ```$ ./MAGNET.sh -h .```, or ```./MAGNET.sh -help``` while in the MAGNET directory. However, it is more useful (particularly when running for the first time) to get _verbose usage info_ for MAGNET, including detailed descriptions of each option; do this by typing ```$ ./MAGNET.sh -H .``` or ```./MAGNET.sh -Help``` (capital "H" flag) at the command line while in the MAGNET directory. The verbose usage text is as follows:
```
$ ./MAGNET.sh
Usage: MAGNET.sh [Help: -h help H Help] [Options: -b r g m] inputNexus 
 ## Help:
  -h   help text (also: -help)
  -H   verbose help text (also: -Help)

 ## Options:
  -b   numBootstraps (def: 100) RAxML bootstrap pseudoreplicates
  -r   raxmlModel (def: GTRGAMMA; other: GTRGAMMAI, GTRCAT, GTRCATI)
  -g   gapThreshold (def: 0.001=essentially zero gaps allowed unless >1000 
       individuals; takes float proportion value)
  -m   indivMissingData (def: 1=allowed; 0=removed)

 OVERVIEW
 Reads in a single G-PhoCS ('*.gphocs') or NEXUS ('*.nex') datafile, splits each locus into 
 a separate phylip-formatted alignment file, and sets up and runs RAxML (Stamatakis 2014) to 
 infer gene trees for each locus. If a NEXUS datafile is supplied, it is converted into 
 G-PhoCS format (Gronau et al. 2011). Sequence names may not include hyphen characters, or 
 there will be issues. For info on various dependencies, see 'README.md' file in the 
 distribution folder; however, it is key that the dependencies are available from the command 
 line interface. 

 DETAILS
 The -b flag sets the number of boostrap pseudoreplicates for RAxML to perform while estimating 
 the gene tree for each locus. The default is 100; remove bootstrapping by setting to 0.

 The -r flag sets the RAxML model for each locus. This uses the full default GTRGAMMA model,
 and at present it is not possible to vary the model across loci. If you want to use HKY
 or K80, you will need to manually change the 'RAxMLRunner.sh' section of this script.

 The following options are available **ONLY** if you are starting from a NEXUS input file:

	The -g flag supplies a 'gap threshold' to an R script, which deletes all column sites in 
	the DNA alignment with a proportion of gap characters '-' at or above the threshold value. 
	If no gap threshold is specified, all sites with gaps are removed by default. If end goal
	is to produce a file for G-PhoCS, you  will want to leave gapThreshold at the default. 
	However, if the next step in your pipeline involves converting from .gphocs to other data 
	formats, you will likely want to set gapThreshold=1 (e.g. before converting to phylip 
	format for RAxML). 

	The -m flag allows users to choose their level of tolerance for individuals with missing
	data. The default is indivMissingData=1, allowing individuals with runs of 10 or more 
	missing nucleotide characters ('N') to be kept in the alignment. Alternatively, setting
	indivMissingData=0 removes all such individuals from each locus; thus, while the input
	file would have had the same number of individuals across loci, the resulting file could
	have varying numbers of individuals for different loci.

 CITATION
 Bagley, J.C. 2017. MAGNET. GitHub package, Available at: 
	<http://github.com/justincbagley/MAGNET>.
 or
 Bagley, J.C. 2017. MAGNET. GitHub package, Available at: 
	<http://doi.org/10.5281/zenodo.166024>.

 REFERENCES
 Gronau I, Hubisz MJ, Gulko B, Danko CG, Siepel A (2011) Bayesian inference of ancient human 
	demography from individual genome sequences. Nature Genetics, 43, 1031-1034.
 Stamatakis A (2014) RAxML version 8: a tool for phylogenetic analysis and post-analysis of 
	large phylogenies. Bioinformatics, 30, 1312-1313.
```

**_IMPORTANT NOTE on NEXUS2gphocs usage:_ In its current form, you must move NEXUS2gphocs.sh (out of the shell folder) _and_ rmGapSites.r (out of the R folder) into the MAGNET directory in order to run NEXUS2gphocs as a standalone script.** (This also assumes the target inputNexus is also located in the MAGNET dir.)

To get usage info for NEXUS2gphocs.sh, type the first line blow while in the MAGNET directory, and you will get the output that follows:

```
./NEXUS2gphocs.sh
Usage="Usage: NEXUS2gphocs.sh [Help: -h help H Help] [Options: -b r g m] inputNexus 
 ## Help:
  -h   help text (also: -help)
  -H   verbose help text (also: -Help)

 ## Options:
  -g   gapThreshold (def: $MY_GAP_THRESHOLD=essentially zero gaps allowed unless >1000 
       individuals; takes float proportion value)
  -m   indivMissingData (def: $MY_INDIV_MISSING_DATA=allowed; 0=removed)

 OVERVIEW
 Reads in a single NEXUS datafile and converts it to '.gphocs' format for G-PhoCS software
 (Gronau et al. 2011). Sequence names may not include hyphen characters, or there will be 
 issues. For best results, update to R v3.3.1 or higher.

 The -g flag supplies a 'gap threshold' to an R script, which deletes all column sites in 
 the DNA alignment with a proportion of gap characters '-' at or above the threshold value. 
 If no gap threshold is specified, all sites with gaps are removed by default. If end goal
 is to produce a file for G-PhoCS, you  will want to leave gapThreshold at the default. 
 However, if the next step in your pipeline involves converting from .gphocs to other data 
 formats, you will likely want to set gapThreshold=1 (e.g. before converting to phylip 
 format for RAxML). 

 The -m flag allows users to choose their level of tolerance for individuals with missing
 data. The default is indivMissingData=1, allowing individuals with runs of 10 or more 
 missing nucleotide characters ('N') to be kept in the alignment. Alternatively, setting
 indivMissingData=0 removes all such individuals from each locus; thus, while the input
 file would have had the same number of individuals across loci, the resulting file could
 have varying numbers of individuals for different loci.

 Dependencies: Perl; R; and Naoki Takebayashi Perl scripts 'fasta2phylip.pl' and 
 'selectSites.pl' in working directory or available from command line (in your path)."

 CITATION
 Bagley, J.C. 2017. MAGNET. GitHub package, Available at: 
	<http://github.com/justincbagley/MAGNET>.
 or
 Bagley, J.C. 2017. MAGNET. GitHub package, Available at: 
	<http://doi.org/10.5281/zenodo.166024>.
```

**Below I give some examples of how to use the software under the two most common scenarios:**

**SCENARIO 1.** If your data contain very little missing data and, in particular, they contain no individuals with all missing data for a locus, then it should be fine to run MAGNET using the default options (giving no flags) as follows:
```
##--Scenario 1, generic usage:
./MAGNET.sh inputNexus

##--Example:
cd ~/Downloads/MAGNET-master/
./MAGNET.sh example.nex
```
**SCENARIO 2.** If your data are relatively lower quality data (e.g. from NGS runs) and you have lots of missing data, including individuals with all missing data for a locus (as is common for RAD tag/SNP data), then RAxML will not run properly under the default MAGNET options. You will likely get up to ~10 messages like "ERROR: Sequence XXXXX consists entirely of undetermined values which will be treated as missing data", follwed by a summary like this: "ERROR: Found 10 sequences that consist entirely of undetermined values, exiting...", and RAxML will quit. The rest of the pipeline will be affected, for example the final summary gene tree file will make no sense because it will simply include a concatenation of all files in the working directory. 

To avoid the above issues caused by large amounts of missing data, you should run MAGNET while **setting the -m flag to 0** (indivMissingData=0) to specify that individuals with missing data are NOT allowed:
```
##--Scenario 2, all params except indivMissingData set to default options:
./MAGNET.sh -m0 inputNexus

##--Example:
cd ~/Downloads/MAGNET-master/
./MAGNET.sh -m0 example.nex
```

In addition to the above, here are illustrations of the **RAxML options**:
```
##--Scenario 1, GTRCAT model, instead of the default GTRGAMMA model:
./MAGNET.sh -rGTRCAT inputNexus

##--Scenario 2, 500 bootstrap reps per locus, instead of the default 100:
./MAGNET.sh -b500 -m0 inputNexus

##--Scenario 2, 0 bootstrap reps per locus:
./MAGNET.sh -b0 -m0 inputNexus
```

## ACKNOWLEDGEMENTS

I thank the Brigham Young University Fulton Supercomputing Lab (FSL) for providing computational resources used during the development of this software.

## REFERENCES

- Chifman J, Kubatko L (2014) Quartet inference from SNP data under the coalescent model. Bioinformatics, 30, pages 3317–3324.
- Eaton DAR (2014) PyRAD: assembly of de novo RADseq loci for phyloge-netic analyses. Bioinformatics, 30, 1844–1849.
- Gronau I, Hubisz MJ, Gulko B, Danko CG, Siepel A (2011) Bayesian inference of ancient human demography from individual genome sequences. Nature Genetics, 43, 1031-1034.
- Larget BR, Kotha SK, Dewey CN, Ané C (2010) BUCKy: gene tree/species tree reconciliation with Bayesian concordance analysis. Bioinformatics, 26(22):2910-2911.
- Liu L, Yu L, Edwards SV (2010) A maximum pseudo-likelihood approach for estimating species trees under the coalescent model. BMC Evol Biol, 10(1):302.
- Liu L, Yu L (2011) Estimating species trees from unrooted gene trees. Syst Biol, 60(5):661-667.
- Mirarab S, Warnow T (2015) ASTRAL-II: coalescent-based species tree estimation with many hundreds of taxa and thousands of genes. Bioinformatics, 30:44-52.
- Peterson BK, Weber JN, Kay EH, Fisher HS, Hoekstra HE (2012) Double digest RADseq: an inexpensive method for de novo SNP discovery and genotyping in model and non-model species. PLoS One, 7, e37135.
- Stamatakis A (2014) RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics, 30.9, 1312-1313.
- Vachaspati P, Warnow T (2015) ASTRID: Accurate Species TRees from Internode Distances. BMC Genomics, 16(Suppl 10):S3.


March 5, 2017
Justin C. Bagley, Tuscaloosa, AL, USA

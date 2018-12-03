# MAGNET Change Log

## MAGNET v0.1.6 (modified minor version release - several changes after official 0.1.6 release)
- **December 2018:** Added new MAGNET script updated to include --resume option, and to set raxml executable name one of two ways after detecting machine type (raxml on Mac, raxmlHPC-SS3 on Linux/supercomputer).
- **November 25 2018:** Added new 'RAxMLRunChecker.sh' script v1.0, which counts the number of completed RAxML runs during the course of, or after, a MAGNET pipeline run, and also collates information on the dataset (e.g. number of patterns) and run (e.g. run time, optimum likelihood) for each locus/partition.
- **November 20 2018:** Updated MAGNET with edited 'MAGNET.sh' (now v0.1.7) and 'NEXUS2gphocs.sh' (now v1.3) scripts containing an important bug fix and some new code checking for whether the NEXUS to fasta file conversion succeeded.

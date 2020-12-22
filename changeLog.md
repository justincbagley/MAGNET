# MAGNET Change Log

## MAGNET v1.2.0 (current official major version release)

-   **December 19-21, 2020:** Updated the standalone `MAGNET.sh` script and pipeline files (`shell` and `R` scripts) in this repository to version last released with the current head of PIrANHA (MAGNET v1.1.1), and then updated MAGNET to v1.2.0 here and in PIrANHA repo through various changes. Changes include a full rewrite of code for manually parsing the options, as an overhaul of the usage texts, and various minor edits to essentially all files (including some stylistic changes and some code improvements).

...

## MAGNET v0.1.9 (official minor version release - several changes after 0.1.6 release)

-   **February 2019:** Updated MAGNET script by adding a getBipartTrees function to the MAGNET pipeline, which organizes RAxML bipartitions trees for each locus (= best ML trees with bootstrap proportions along nodes the corresponding bootstrap searches search; resulting from ```-f a -x```options, which are included in all MAGNET calls to RAxML). Edited header and script banner to be unofficially prepped for future official release versioning 0.1.9.
-   **December 2018:** Added new MAGNET script updated to include --resume option, fix step info output to screen, and set raxml executable name one of two ways after detecting machine type (```raxml``` on Mac, ```raxmlHPC-SS3``` on Linux/supercomputer).
-   **November 25 2018:** Added new 'RAxMLRunChecker.sh' script v1.0, which counts the number of completed RAxML runs during the course of, or after, a MAGNET pipeline run, and also collates information on the dataset (e.g. number of patterns) and run (e.g. run time, optimum likelihood) for each locus/partition.
-   **November 20 2018:** Updated MAGNET with edited 'MAGNET.sh' (now v0.1.7) and 'NEXUS2gphocs.sh' (now v1.3) scripts containing an important bug fix and some new code checking for whether the NEXUS to fasta file conversion succeeded.

#!/bin/sh

##########################################################################################
#                                                                                        #
#                           phyNcharSumm v1.0.2, December 2020                           #
# Author: Justin C. Bagley                                                               #
# Date: Created by Justin Bagley on November 9, 2016.                                    #
# Last update: December 23, 2020                                                         #
# Copyright (c) 2016-2020 Justin C. Bagley. All rights reserved.                         #
# Please report bugs to <jbagley@jsu.edu>.                                               #
#                                                                                        #
# Description:                                                                           #
# SHELL SCRIPT THAT SUMMARIZES THE NUMBER OF CHARACTERS IN EACH OF N PHYLIP DNA SEQUENCE #
# ALIGNMENTS IN CURRENT WORKING DIRECTORY AND SAVES THIS INFORMATION TO FILE             #
#                                                                                        #
##########################################################################################

######################################## START ###########################################

echo "INFO      | $(date) |----------------------------------------------------------------"
echo "INFO      | $(date) | phyNcharSumm, v1.0.2 December 2020  (part of PIrANHA v0.4a4)   "
echo "INFO      | $(date) | Copyright (c) 2016-2020 Justin C. Bagley. All rights reserved. "
echo "INFO      | $(date) |----------------------------------------------------------------"

###### Starting from a folder containing multiple PHYLIP alignment files (e.g. as generated
## by MAGNET.sh for indiv. SNP loci), this script uses a for loop to echo all alignment
## names (containing locus or taxon information) to file, and then echo the number of 
## characters (Nchar) recursively to a file named "nchar.txt" in the working directory.

echo "INFO      | $(date) |          Saving number of characters for each alignment in file named './nchar.txt'. "
(
	for i in ./*.phy; do 
		echo "$i" >> phyalign_names.txt; 
		echo "$(head -n1 $i | awk -F"[0-9]*\ " '{print $NF}')" >> nchar.txt ; 
	done;
)

## I would like to build on this by calculating statistics from the nchar list from within
## the shell.

echo "-------------------------------------------------------------------------------------"

#
#
#
######################################### END ############################################

exit 0

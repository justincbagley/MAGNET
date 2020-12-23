#!/bin/sh

##########################################################################################
#  __  o  __   __   __  |__   __                                                         #
# |__) | |  ' (__( |  ) |  ) (__(                                                        # 
# |                                                                                      #
#                                                                                        #
# File: MAGNET.sh ~ MAny GeNE Trees, v1.2.0                                              #
  VERSION="v1.2.0"                                                                       #
# Author: Justin C. Bagley                                                               #
# Date: Created by Justin Bagley on Mon, Aug 29 13:12:45 2016 -0700.                     #
# Last update: December 23, 2020                                                         #
# Copyright (c) 2016-2020 Justin C. Bagley. All rights reserved.                         #
# Please report bugs to <jbagley@jsu.edu>.                                               #
#                                                                                        #
# Description:                                                                           #
# SHELL PIPELINE FOR AUTOMATING ESTIMATION OF A MAXIMUM-LIKELIHOOD (ML) GENE TREE IN     #
# RAxML FOR EACH OF MANY LOCI IN A RAD-seq, UCE, OR OTHER MULTILOCUS SEQUENCE DATASET    #
#                                                                                        #
##########################################################################################

# Short present working directory (PWD) functions
# ------------------------------------------------------
# New common functions (most scripts written/edited from
# July 2020 onwards) for echoing $PWD, but truncating 
# absolute path to precisely fit the Terminal window.
# Functions for two cases: 1) Regular PWD and 2) user-
# specified PWD.
# ------------------------------------------------------
function echoShortPWD () {
		MY_ABS_PATH_LENGTH="$(echo "$PWD" | wc -c | sed 's/\ //g')";
		MY_ABS_PATH_ECHO_LENGTH="$(calc "$MY_ABS_PATH_LENGTH"+43)";
		MY_BASH_WINDOW_COLS="$(tput cols | sed 's/\ //g')";
#
		if [[ "$MY_ABS_PATH_ECHO_LENGTH" -gt "85" ]] && [[ "$MY_ABS_PATH_ECHO_LENGTH" -gt "$MY_BASH_WINDOW_COLS" ]]; then
			MY_CORRECTION_LENGTH="$(calc "$MY_ABS_PATH_ECHO_LENGTH"-"$MY_BASH_WINDOW_COLS"+3)";
			MY_NUM_FINAL_PWD_CHARS="$(calc "$MY_ABS_PATH_LENGTH"-"$MY_CORRECTION_LENGTH")";
			SHORT_PWD="$(echo ${PWD:${#PWD}<$MY_NUM_FINAL_PWD_CHARS?0:-$MY_NUM_FINAL_PWD_CHARS})";   ## Get last $MY_NUM_FINAL_PWD_CHARS characters of $PWD
			echo "INFO      | $(date) | Starting input directory (using current dir): "
			echo "INFO      | $(date) | ...$SHORT_PWD"
		else
			echo "INFO      | $(date) | Starting input directory (using current dir): "
			echo "INFO      | $(date) | $PWD"	
		fi
}

# Working directory function
# ------------------------------------------------------
# Common function (most scripts) for echoing and cd'ing
# into user-specified working dir. Three cases: 1) cwd,
# 2) up one dir, and 3) some other dir.
# ------------------------------------------------------
function echoCDWorkingDir () {
if [[ -s "$USER_SPEC_PATH" ]]; then
	if [[ "$USER_SPEC_PATH" = "$(printf '%q\n' "$(pwd -P)")" ]] || [[ "$USER_SPEC_PATH" = "." ]]; then
		MY_CWD="$(printf '%q\n' "$(pwd -P)" | sed 's/\\//g')";
		echo "INFO      | $(date) | Setting working directory to:  "
		echo "INFO      | $(date) | $MY_CWD "
	elif [[ "$USER_SPEC_PATH" != "$(printf '%q\n' "$(pwd)")" ]]; then
		if [[ "$USER_SPEC_PATH" = ".." ]] || [[ "$USER_SPEC_PATH" = "../" ]] || [[ "$USER_SPEC_PATH" = "..;" ]] || [[ "$USER_SPEC_PATH" = "../;" ]]; then
			cd ..;
			MY_CWD="$(printf '%q\n' "$(pwd)" | sed 's/\\//g')";
			echo "INFO      | $(date) | Setting working directory to user-specified dir:  "	
			echo "INFO      | $(date) | $MY_CWD "
		else
			MY_CWD=$USER_SPEC_PATH
			cd "$MY_CWD";
			echo "INFO      | $(date) | Setting working directory to user-specified dir:  "	
			echo "INFO      | $(date) | $MY_CWD "
		fi
	else
		echo "WARNING   | $(date) | Null working directory path. Quitting... "
		exit 1
	fi
fi

}

# Calculator
# ------------------------------------------------------
# Function to send commands to Bash's arbitrary precision
# calculator language.
#
# A calculator function, which echos input to Bash bc. 
# Usage: calc <math_commands>, e.g. calc 1+2, or with 
# variables, e.g. calc $VAR1+$VAR2.
# ------------------------------------------------------
function calc () {
	bc -l <<< "$@" ;
}

# Check machine type
# ------------------------------------------------------
# Common function for checking machine type, setting related
# variables.
# ------------------------------------------------------
function checkMachineType () {
unameOut="$(uname -s)";
case "${unameOut}" in
	Linux*)     machine=Linux;;
	Darwin*)    machine=Mac;;
	CYGWIN*)    machine=Cygwin;;
	MINGW*)     machine=MinGw;;
	*)          export machine="UNKNOWN:${unameOut}"
esac;
}

# Import variables
# ------------------------------------------------------
export EP=$(printf '!')
export TAB=$(printf '\t')
export CR=$(printf '\r')




MAGNET () {

######################################## START ###########################################
##########################################################################################

echo "INFO      | $(date) |----------------------------------------------------------------"
echo "INFO      | $(date) | MAGNET, v1.2.0 December 2020                                   "
echo "INFO      | $(date) | Copyright (c) 2016-2020 Justin C. Bagley. All rights reserved. "
echo "INFO      | $(date) |----------------------------------------------------------------"


# Single NEXUS (or .gphocs) file run
# ------------------------------------------------------
# Run on single input NEXUS file ('.nex'; also accepts single G-PhoCS 
# formatted file, with extension '.gphocs') in current working directory, 
# specified with -i flag.
# ------------------------------------------------------

#######
if [[ "$STARTING_FILE_TYPE" = "1" ]] && [[ "$MY_NEXUS" != "NULL" ]]; then

	echo "INFO      | $(date) | Starting MAGNET pipeline... "
	echo "INFO      | $(date) | Running with the following options: "
	echo "INFO      | $(date) | - NEXUS file, <inputNEXUS> = ${MY_NEXUS} "
	echo "INFO      | $(date) | - Starting <fileType> = ${STARTING_FILE_TYPE} "
	echo "INFO      | $(date) | - RAxML <executable> = ${MY_RAXML_EXECUTABLE} "
	echo "INFO      | $(date) | - Bootstrap reps, <numBootstraps> = ${MY_NUM_BOOTREPS} "
	echo "INFO      | $(date) | - RAxML model, <raxmlModel> = ${MY_RAXML_MODEL} "
	echo "INFO      | $(date) | - <simpleModel> option = ${MY_SIMPLE_MODEL} "
	echo "INFO      | $(date) | - <gapThreshold> option = ${MY_GAP_THRESHOLD} "
	echo "INFO      | $(date) | - <indivMissingData> option = ${MY_INDIV_MISSING_DATA} "
	echo "INFO      | $(date) | - Outgroup taxon, <outgroup> = ${MY_OUTGROUP} "
	echo "INFO      | $(date) | - RAxML output name = ${MY_OUTPUT_NAME} "
	echo "INFO      | $(date) | - Resume switch (--resume) = ${MY_RESUME_SWITCH} "

# --------------------------------------------------
# -- STEP #1: SETUP.
# --------------------------------------------------
	echo "INFO      | $(date) | ----------------------------------- "
	echo "INFO      | $(date) | # Step #1: Set up workspace (e.g. functions, working directory) and check machine type. " # | tee -a "$MY_OUTPUT_FILE_SWITCH"
	echo "INFO      | $(date) | ----------------------------------- "

	# SET WORKING DIRECTORY AND CHECK MACHINE TYPE
	# --------------------------------------------------
	echoShortPWD
	export MY_WORKING_DIR="$(pwd)"
	checkMachineType

	# START DEBUG MODE, IF CALLED
	# --------------------------------------------------
	if [[ "$MY_DEBUG_MODE_SWITCH" != "0" ]]; then set -xv; fi

# --------------------------------------------------
# -- STEP #2: HANDLE INPUT FILE(S).
# --------------------------------------------------
	echo "INFO      | $(date) | ----------------------------------- "
	echo "INFO      | $(date) | # Step #2: Input single NEXUS (or G-PhoCS-formatted) file, or multiple PHYLIP files. "
	echo "INFO      | $(date) | ----------------------------------- "
	echo "INFO      | $(date) | For -f 1 or -f 2, if '.gphocs' input file present, continue; else convert NEXUS file "
	echo "INFO      | $(date) | to G-PhoCS format using NEXUS2gphocs code. If -f 3, then run multiple PHYLIP files in  "
	echo "INFO      | $(date) | RAxML."
	shopt -s nullglob
	if [[ -n $(find . -name "*.gphocs" -type f) ]]; then
		echo "INFO      | $(date) | Found '.gphocs' input file... "
	    MY_GPHOCS_DATA_FILE=./*.gphocs ;	 ## Assign G-PhoCS-formatted genomic/SNP data file (originally produced/output by pyRAD) in run directory to variable.
	else
	    echo "WARNING   | $(date) | No '.gphocs' input file in current working directory... "
	    echo "INFO      | $(date) | Attempting to convert NEXUS file, if present, to GPho-CS format... "
	fi

	# ---------------------------------------------------------- #
	# ---------------- NEXUS2gphocs.sh FUNCTION ---------------- #
	# ---------------------------------------------------------- #
	NEXUS2gphocs_function () {

	# GET NEXUS FILE & DATA CHARACTERISTICS, CONVERT NEXUS TO FASTA FORMAT
	# --------------------------------------------------
	# Extract charset info from sets block at end of NEXUS file: 
	# --------------------------------------------------
	MY_NEXUS_CHARSETS="$(egrep "charset|CHARSET" "$MY_NEXUS" | \
	awk -F"=" '{print $NF}' | sed 's/\;/\,/g' | \
	awk '{a[NR]=$0} END {for (i=1;i<NR;i++) print a[i];sub(/.$/,"",a[NR]);print a[NR]}' | \
	sed 's/\,/\,'$CR'/g' | sed 's/^\ //g')" ;
	
	# Count number of loci present in the NEXUS file, based on number of charsets defined.
	# Also get corrected count starting from 0 for numbering loci below...
	MY_NLOCI="$(echo "$MY_NEXUS_CHARSETS" | wc -l)";
	MY_CORR_NLOCI="$(calc "$MY_NLOCI" - 1)";
	
	# This is the base name of the original nexus file, so you have it. This will not work if NEXUS file name is written in all caps, ".NEX", in the file name.
	MY_NEXUS_BASENAME="$(echo "$MY_NEXUS" | sed 's/\.\///g; s/\.nex//g')";
	
	# Convert data file from NEXUS to FASTA format using bioscripts.convert v0.4 Python package:
	# However, if alignment is too long (>100,000 bp), then need to convert to FASTA using my 
	# script and then wrap to 60 characters with fold function (as suggested at stackexchange
	# post URL: https://unix.stackexchange.com/questions/25173/how-can-i-wrap-text-at-a-certain-column-size).
	# If this conversion failes because the alignment is too long, then the code to follow 
	# will have nothing to work with. So, I am here adding a conditional quit if the FASTA
	# file is not generated.

	#---------TODO: ADD IF/THEN CONDITIONAL AND MY OWN NEXUS2FASTA SCRIPT HERE!!!!----------#

	convbioseq fasta "$MY_NEXUS" > "$MY_NEXUS_BASENAME".fasta ;
	MY_FASTA="$(echo "$MY_NEXUS_BASENAME".fasta | sed 's/\.\///g; s/\.nex//g')";
	
	# The line above creates a file with the name basename.fasta, where basename is the base name of the original .nex file. For example, "hypostomus_str.nex" would be converted to "hypostomus_str.fasta".
	# Check to make sure the FASTA was created; if so, echo info, if not, echo warning and quit:
	if [[ -s "$MY_NEXUS_BASENAME".fasta ]]; then
		echo "INFO      | $(date) | Input NEXUS was successfully converted to FASTA format. Moving forward... "
	else
		echo "WARNING   | $(date) | NEXUS to FASTA file conversion FAILED. Quitting... "
		exit 1
	fi
	
	# PUT COMPONENTS OF ORIGINAL NEXUS FILE AND THE FASTA FILE TOGETHER TO MAKE A
	# A G-PhoCS-FORMATTED DATA FILE
	# --------------------------------------------------
	# Make top (first line) of the G-Phocs format file, which should
	# have the number of loci on the first line:
	# --------------------------------------------------
	echo "$MY_NLOCI" | sed 's/[\ ]*//g' > gphocs_top.txt ;
	
	echo "$MY_GAP_THRESHOLD" > ./gap_threshold.txt ;
	count=0
	(
		for j in $MY_NEXUS_CHARSETS; do
			echo "$j"
			charRange="$(echo "$j" | sed 's/\,//g')";
	        echo "$charRange"
	        export setLower="$(echo "$j" | sed 's/\-.*$//g')";
			export setUpper="$(echo "$j" | sed 's/[0-9]*\-//g' | sed 's/\,//g; s/\ //g')";
	
			**/selectSites.pl -s "$charRange" "$MY_FASTA" > ./sites.fasta ;
				
			**/fasta2phylip.pl ./sites.fasta > ./sites.phy ;

			# Need to make sure there is a space between the tip taxon name (10 characters as output
			# by the fasta2phylip.pl Perl script) and the corresponding sequence, for all tips. Use
			# a perl search and replace for this:

			perl -p -i -e 's/^([A-Za-z0-9\-\_\ ]{10})/$1\ /g' ./sites.phy ;

				# If .phy file from NEXUS charset $j has gaps in alignment, then call 
				# rmGapSites.R R script to remove all column positions with gaps from
				# alignment and output new, gapless PHYLIP file named "./sites_nogaps.phy". 
				# If charset $j does not have gaps, go to next line of loop. We do the 
				# above by first creating a temporary file containing all lines in
				# sites.phy with the gap character:
				grep -n "-" ./sites.phy > ./gaptest.tmp ;
				
				# Next, we test for nonzero testfile, indicating presence of gaps in $j, 
				# using UNIX test operator "-s" (returns true if file size is not zero). 
				# If fails, cat sites.phy into file with same name as nogaps file that
				# is output by rmGapSites.R and move forward:
				if [[ -s ./gaptest.tmp ]]; then
					echo "Removing column sites in locus${count} with gaps. "
					R CMD BATCH **/rmGapSites.R
				else
			   		echo ""
			   		cat ./sites.phy > ./sites_nogaps.phy ;
				fi
				
				export phylip_header="$(head -n1 ./sites_nogaps.phy)";
	        	locus_ntax="$(head -n1 ./sites_nogaps.phy | sed 's/[\ ]*[.0-9]*$//g')";
				locus_nchar="$(head -n1 ./sites_nogaps.phy | sed 's/[0-9]*\ //g')";
			
			
        		if [[ "$MY_INDIV_MISSING_DATA" = "0" ]]; then
					sed '1d' ./sites_nogaps.phy | egrep -v 'NNNNNNNNNN|nnnnnnnnnn' > ./cleanLocus.tmp ;
					cleanLocus_ntax="$(cat ./cleanLocus.tmp | wc -l)";
					echo locus"$((count++))" "$cleanLocus_ntax" "$locus_nchar" > ./locus_top.tmp ;
					cat ./locus_top.tmp ./cleanLocus.tmp >> ./gphocs_body.txt ;
				else
					echo locus"$((count++))" "$locus_ntax" "$locus_nchar" > ./locus_top.tmp ;
					cat ./locus_top.tmp ./sites_nogaps.phy >> ./gphocs_body.txt ;
				fi

			if [[ -s ./sites.fasta ]] && [[ -s ./sites.phy ]] && [[ ! -z ./*.tmp ]] && [[ -s ./sites_nogaps.phy ]]; then
				rm ./sites.fasta ./sites.phy ;
				if [[ "$(ls -1 ./*.tmp 2>/dev/null | wc -l | sed 's/\ //g')" != "0"  ]]; then 
					rm ./*.tmp ; 
				fi
				rm ./sites_nogaps.phy ;
			fi
		done
	)

	grep -v "^[0-9]*\ [0-9]*.*$" ./gphocs_body.txt > ./gphocs_body_fix.txt ;
	sed 's/locus/'$CR'locus/g' ./gphocs_body_fix.txt > ./gphocs_body_fix2.txt ;
	cat ./gphocs_top.txt ./gphocs_body_fix2.txt > "$MY_NEXUS_BASENAME".gphocs ;

	# CLEANUP: REMOVE UNNECESSARY FILES
	# --------------------------------------------------
	# Remove temporary header, threshold and body files.
	# --------------------------------------------------
	if [[ -s ./gphocs_top.txt ]]; then rm ./gphocs_top.txt ; fi ;
	if [[ -s ./gap_threshold.txt ]]; then rm ./gap_threshold.txt ; fi ;
	if [[ ! -z ./gphocs_body* ]]; then rm ./gphocs_body* ; fi ;

}
	shopt -s nullglob
	if [[ -n $(find . -name "*.nex" -type f) ]]; then
		# Run the function!
		NEXUS2gphocs_function
	else
		echo "INFO      | $(date) | No NEXUS files in current working directory. Continuing... "
	fi
	
	shopt -s nullglob
	if [[ -n $(find . -name "*.gphocs" -type f) ]]; then
		echo "INFO      | $(date) | MAGNET successfully created a '.gphocs' input file from the existing NEXUS file... "
	    MY_GPHOCS_DATA_FILE=./*.gphocs ;	 ## Assign G-PhoCS-formatted genomic/SNP data file (originally produced/output by pyRAD) in run directory to variable.
	else
	    echo "WARNING   | $(date) | Failed to convert NEXUS file into G-PhoCS format... "
	    echo "INFO      | $(date) | Quitting."
	    exit
	fi


	# ---------------------------------------------------------- #
	# ------------- gphocs2multiPhylip.sh FUNCTION ------------- #
	# ---------------------------------------------------------- #
	gphocs2multiPhylip_function () {

	MY_NLOCI="$(head -n1 "$MY_GPHOCS_DATA_FILE")";
	MY_CORR_NLOCI="$(calc "$MY_NLOCI" - 1)";

	# --------------------------------------------------
	# -- STEP #3: MAKE ALIGNMENTS FOR EACH LOCUS FROM GPHOCS FILE, IF STARTING FROM GPHOCS FILE.
	# --------------------------------------------------
	echo "INFO      | $(date) | ----------------------------------- "
	echo "INFO      | $(date) | # Step #3: Make alignments for each locus. "
	echo "INFO      | $(date) | ----------------------------------- "
	echo "INFO      | $(date) | In a single loop, using info from '.gphocs' file to split each locus block into a separate PHYLIP-formatted "
	echo "INFO      | $(date) | alignment file using gphocs2multiPhylip code... "
	(
		for (( i=0; i<=MY_CORR_NLOCI; i++ )); do
			echo "$i"
			MY_NTAX="$(grep -n "locus$i\ " "$MY_GPHOCS_DATA_FILE" | \
			awk -F"locus$i " '{print $NF}' | sed 's/\ [0-9]*//g')";			

			MY_NCHAR="$(grep -n "locus$i\ " "$MY_GPHOCS_DATA_FILE" | \
			awk -F"locus$i [0-9]*\ " '{print $NF}')";
		
			awk "/locus"$i"\ / {for(j=1; j<=MY_NTAX; j++) {getline; print}}" "$MY_GPHOCS_DATA_FILE" > ./locus"$i".tmp ;

			echo "$MY_NTAX $MY_NCHAR" > ./locus"$i"_header.tmp ;
				
			cat ./locus"$i"_header.tmp ./locus"$i".tmp > ./locus"$i".phy ;
		done
	)

	# CLEANUP: REMOVE UNNECESSARY OR TEMPORARY FILES
	# --------------------------------------------------
	# Remove temporary files generated during run.
	# --------------------------------------------------
	if [[ ! -z ./*.tmp ]]; then rm ./*.tmp ; fi
}
	# Run the function!
	gphocs2multiPhylip_function


	# ASSIGN PHYLIP FILES TO ENVIRONMENTAL VARIABLE
	# --------------------------------------------------
	# Assign PHYLIP-formatted genomic/SNP data files (e.g. output by gphocs2multiPhylip.sh 
	# shell script) in run directory to variable.
	# --------------------------------------------------	
	if [[ -n $(find . -name "*.phy" -type f) ]]; then
	    MY_PHYLIP_ALIGNMENTS=./*.phy ;		## Assign PHYLIP-formatted genomic/SNP data files (e.g. output by gphocs2multiPhylip.sh shell script) in run directory to variable.
	else
	    echo "..."
	fi


	# ---------------------------------------------------------- #
	# ------------- MultiRAxMLPrepper.sh FUNCTION -------------- #
	# ---------------------------------------------------------- #
	MultiRAxMLPrepper_function () {

	# --------------------------------------------------
	# -- STEP #4: MAKE RAxML RUN FOLDERS.
	# --------------------------------------------------
	echo "INFO      | $(date) | ----------------------------------- "
	echo "INFO      | $(date) | # Step #4: Make and check RAxML run folders. "
	echo "INFO      | $(date) | ----------------------------------- "

	if [[ "$MY_RESUME_SWITCH" = "0" ]]; then
	
		MY_N_PHYLIP_FILES="$(ls -1 ./*.phy 2>/dev/null | wc -l | sed 's/\ //g')";
	
		# MAKE RAXML RUN FOLDERS, ONE FOR EACH GENE/LOCUS
		# --------------------------------------------------
		# Loop through the input .phy files and do the following for each file: (A) generate one 
		# folder per .phy file with the same name as the file, only minus the extension, then 
		# (B) move input .phy file into corresponding folder.
		# --------------------------------------------------
		count=1
		(
			for i in $MY_PHYLIP_ALIGNMENTS; do
				MY_BASENAME="$(basename "$i" '.phy')";
				mkdir "$MY_BASENAME"/ ;
				cp "$i" "$MY_BASENAME"/ ;
				echo "$((count++))" >/dev/null 2>&1 ;
			done
		)
	
		# CHECK RUN FOLDERS
		# --------------------------------------------------
		# Setup and run check on the number of run folders created by the program:
		# --------------------------------------------------
		export MY_TOTAL_FILECOUNT="$(find . -type f | wc -l)";
		export MY_TOTAL_DIRCOUNT="$(find . -type d | wc -l)";
		#MY_NUM_RUN_FOLDERS="$(ls ./*/*.phy | wc -l | perl -pe 's/\t//g; s/\ //g')";
	
		echo "INFO      | $(date) | Number of run folders created: $count "
	
		if [[ "$count" = "$MY_N_PHYLIP_FILES" ]]; then
			echo "INFO      | $(date) | Folder check PASSED: number of run folders matches number of PHYLIP alignments. "
		else
			echo "WARNING   | $(date) | Folder check FAILED: number of run folders (${count}) does NOT match the number of PHYLIP alignments (${MY_N_PHYLIP_FILES}). This may cause errors. "
		fi
	
	elif [[ "$MY_RESUME_SWITCH" = "1" ]]; then
		if [[ "$count" = "$MY_N_PHYLIP_FILES" ]]; then
			echo "IMPORTANT${EP}| $(date) | Resuming a previous/existing run in current working dir. Skipping MultiRAxMLPrepper, using available run folders... "
			echo "INFO      | $(date) | Folder check PASSED: number of run folders matches number of PHYLIP alignments. "
		else
			echo "WARNING   | $(date) | Folder check FAILED: number of run folders does NOT match the number of PHYLIP alignments. There may be errors. "
		fi
	
	fi
}
	# Run the function!
	MultiRAxMLPrepper_function


	# ---------------------------------------------------------- #
	# ---------------- RAxMLRunner.sh FUNCTION ----------------- #
	# ---------------------------------------------------------- #
	RAxMLRunner_function () {
	
	if [[ "$MY_RESUME_SWITCH" = "0" ]]; then
	
	# --------------------------------------------------
	# -- STEP #5: ESTIMATE BEST ML GENE TREES.
	# --------------------------------------------------
	echo "INFO      | $(date) | ----------------------------------- "
	echo "INFO      | $(date) | # Step #5: Estimate best maximum-likelihood (ML) gene trees. "
	echo "INFO      | $(date) | ----------------------------------- "
	echo "INFO      | $(date) | Looping through and analyzing contents of each run folder in RAxML... "
		# Each folder is set with the locus name corresponding to the locus' position in the
		# original .gphocs alignment (which, if output by pyRAD, is simply in the order in which
		# the loci were logged to file by pyRAD, no special order). Also, each folder contains
		# one .phy file carrying the same basename as the folder name, e.g. "locus0.phy". So,
		# all we need to do here is loop through each folder and call RAxML to run using its
		# contents as the input file, as follows:
		(
			for i in ./*/; do
				if [[ "$i" != "./archive/" ]] && [[ "$i" != "./bad_genes/" ]] && [[ "$i" != "./R/" ]] && [[ "$i" != "./shell/" ]] && [[ "$i" != "./perl/" ]] && [[ "$i" != "./orig_phylip/" ]] && [[ "$i" != "./phylip/" ]] && [[ "$i" != "./orig_fasta/" ]] && [[ "$i" != "./fasta/" ]] && [[ "$i" != "./phylip_files/" ]]; then
					echo "$i"
					cd "$i";
					LOCUS_NAME="$(echo "$i" | sed 's/\.\///g; s/\/$//g')"; # NOTE: not currently using "$LOCUS_NAME" here, but leave for now, bc may need to use it later...
				#
					if [[ "$MY_OUTGROUP" = "NULL" ]] && [[ "$MY_SIMPLE_MODEL" = "NULL" ]]; then
						"$MY_RAXML_EXECUTABLE" -f a -x "$RANDOM""$RANDOM" -p "$RANDOM""$RANDOM" -# ${MY_NUM_BOOTREPS} -m ${MY_RAXML_MODEL} -s ./*.phy -n ${MY_OUTPUT_NAME}
					fi
				#
					if [[ "$MY_OUTGROUP" != "NULL" ]] && [[ "$MY_SIMPLE_MODEL" = "NULL" ]]; then
						"$MY_RAXML_EXECUTABLE" -f a -x "$RANDOM""$RANDOM" -p "$RANDOM""$RANDOM" -# ${MY_NUM_BOOTREPS} -m ${MY_RAXML_MODEL} -s ./*.phy -o ${MY_OUTGROUP} -n ${MY_OUTPUT_NAME}
					fi
				#
					if [[ "$MY_OUTGROUP" = "NULL" ]] && [[ "$MY_SIMPLE_MODEL" != "NULL" ]]; then
						"$MY_RAXML_EXECUTABLE" -f a -x "$RANDOM""$RANDOM" -p "$RANDOM""$RANDOM" -# ${MY_NUM_BOOTREPS} -m ${MY_RAXML_MODEL} -s ./*.phy --${MY_SIMPLE_MODEL} -n ${MY_OUTPUT_NAME}
					fi
				#
					if [[ "$MY_OUTGROUP" != "NULL" ]] && [[ "$MY_SIMPLE_MODEL" != "NULL" ]]; then
						"$MY_RAXML_EXECUTABLE" -f a -x "$RANDOM""$RANDOM" -p "$RANDOM""$RANDOM" -# ${MY_NUM_BOOTREPS} -m ${MY_RAXML_MODEL} -s ./*.phy --${MY_SIMPLE_MODEL} -o ${MY_OUTGROUP} -n ${MY_OUTPUT_NAME}
					fi
					cd ..;
				fi
			done
		)
	
		# Here: adding loop code to move all .phy files remaining in the current working 
		# directory, after Step #3 of the pipeline, to a new folder called "phylip_files". This
		# is done here because if the phylip_files folder is present at the end of Step #3,
		# then RAxML will also try to estimate a gene tree for .phy file(s) in this folder during
		# Step #5 of the pipeline above.
		mkdir ./phylip_files/ ;
		(
			for i in $MY_PHYLIP_ALIGNMENTS; do
				echo "$i"
				mv "$i" ./phylip_files/ ;
			done
		)
	
	elif [[ "$MY_RESUME_SWITCH" = "1" ]]; then
	
	# --------------------------------------------------
	# -- STEP #3 ALT: RESUME GENE TREE ESTIMATION FROM (EXPECTED) PREVIOUS RUN.
	# --------------------------------------------------
	echo "INFO      | $(date) | ----------------------------------- "
	echo "INFO      | $(date) | # Step #3: Resuming gene tree estimation.  "
	echo "INFO      | $(date) | ----------------------------------- "
	echo "INFO      | $(date) | Run on remaining/incomplete run folders, skip those with completed RAxML runs. "
		(
			for i in ./*/; do
				if [[ "$i" != "./archive/" ]] && [[ "$i" != "./bad_genes/" ]] && [[ "$i" != "./R/" ]] && [[ "$i" != "./shell/" ]] && [[ "$i" != "./perl/" ]] && [[ "$i" != "./orig_phylip/" ]] && [[ "$i" != "./phylip/" ]] && [[ "$i" != "./orig_fasta/" ]] && [[ "$i" != "./fasta/" ]] && [[ "$i" != "./phylip_files/" ]]; then
					cd "$i";
					LOCUS_NAME="$(echo "$i" | sed 's/\.\///g; s/\/$//g')";
				#
					if [[ "$MY_OUTPUT_NAME" = "raxml_out" ]] && [[ ! -s ./RAxML_info.raxml_out ]]; then
					echo "$i"
				#
						if [[ "$MY_OUTGROUP" = "NULL" ]] && [[ "$MY_SIMPLE_MODEL" = "NULL" ]]; then
							"$MY_RAXML_EXECUTABLE" -f a -x "$RANDOM""$RANDOM" -p "$RANDOM""$RANDOM" -# ${MY_NUM_BOOTREPS} -m ${MY_RAXML_MODEL} -s ./*.phy -n ${MY_OUTPUT_NAME}
						fi
				#
						if [[ "$MY_OUTGROUP" != "NULL" ]] && [[ "$MY_SIMPLE_MODEL" = "NULL" ]]; then
							"$MY_RAXML_EXECUTABLE" -f a -x "$RANDOM""$RANDOM" -p "$RANDOM""$RANDOM" -# ${MY_NUM_BOOTREPS} -m ${MY_RAXML_MODEL} -s ./*.phy -o ${MY_OUTGROUP} -n ${MY_OUTPUT_NAME}
						fi
				#
						if [[ "$MY_OUTGROUP" = "NULL" ]] && [[ "$MY_SIMPLE_MODEL" != "NULL" ]]; then
							"$MY_RAXML_EXECUTABLE" -f a -x "$RANDOM""$RANDOM" -p "$RANDOM""$RANDOM" -# ${MY_NUM_BOOTREPS} -m ${MY_RAXML_MODEL} -s ./*.phy --${MY_SIMPLE_MODEL} -n ${MY_OUTPUT_NAME}
						fi
				#
						if [[ "$MY_OUTGROUP" != "NULL" ]] && [[ "$MY_SIMPLE_MODEL" != "NULL" ]]; then
							"$MY_RAXML_EXECUTABLE" -f a -x "$RANDOM""$RANDOM" -p "$RANDOM""$RANDOM" -# ${MY_NUM_BOOTREPS} -m ${MY_RAXML_MODEL} -s ./*.phy --${MY_SIMPLE_MODEL} -o ${MY_OUTGROUP} -n ${MY_OUTPUT_NAME}
						fi
				#
					elif [[ "$MY_OUTPUT_NAME" != "raxml_out" ]] && [[ ! -s ./RAxML_info."$MY_OUTPUT_NAME" ]]; then
					echo "$i"
				#
						if [[ "$MY_OUTGROUP" = "NULL" ]] && [[ "$MY_SIMPLE_MODEL" = "NULL" ]]; then
							"$MY_RAXML_EXECUTABLE" -f a -x "$RANDOM""$RANDOM" -p "$RANDOM""$RANDOM" -# ${MY_NUM_BOOTREPS} -m ${MY_RAXML_MODEL} -s ./*.phy -n ${MY_OUTPUT_NAME}
						fi
				#
						if [[ "$MY_OUTGROUP" != "NULL" ]] && [[ "$MY_SIMPLE_MODEL" = "NULL" ]]; then
							"$MY_RAXML_EXECUTABLE" -f a -x "$RANDOM""$RANDOM" -p "$RANDOM""$RANDOM" -# ${MY_NUM_BOOTREPS} -m ${MY_RAXML_MODEL} -s ./*.phy -o ${MY_OUTGROUP} -n ${MY_OUTPUT_NAME}
						fi
				#
						if [[ "$MY_OUTGROUP" = "NULL" ]] && [[ "$MY_SIMPLE_MODEL" != "NULL" ]]; then
							"$MY_RAXML_EXECUTABLE" -f a -x "$RANDOM""$RANDOM" -p "$RANDOM""$RANDOM" -# ${MY_NUM_BOOTREPS} -m ${MY_RAXML_MODEL} -s ./*.phy --${MY_SIMPLE_MODEL} -n ${MY_OUTPUT_NAME}
						fi
				#
						if [[ "$MY_OUTGROUP" != "NULL" ]] && [[ "$MY_SIMPLE_MODEL" != "NULL" ]]; then
							"$MY_RAXML_EXECUTABLE" -f a -x "$RANDOM""$RANDOM" -p "$RANDOM""$RANDOM" -# ${MY_NUM_BOOTREPS} -m ${MY_RAXML_MODEL} -s ./*.phy --${MY_SIMPLE_MODEL} -o ${MY_OUTGROUP} -n ${MY_OUTPUT_NAME}
						fi
					fi
					cd ..;
				fi
			done
		)
	
		if [[ ! -s ./phylip_files/ ]]; then
			mkdir ./phylip_files/ ;
		fi
		(
			for i in $MY_PHYLIP_ALIGNMENTS; do
				echo "INFO      | $(date) | ${i}"
				mv "$i" ./phylip_files/ ;
			done
		)
	fi
}
	# Run the function!
	RAxMLRunner_function


	if [[ "$MY_RESUME_SWITCH" = "0" ]]; then
		# --------------------------------------------------
		# -- STEP #6: CONDUCT POST-PROCESSING OF RAxML RESULTS.
		# --------------------------------------------------
		echo "INFO      | $(date) | ----------------------------------- "
		echo "INFO      | $(date) | # Step #6: RAxML post-processing analyses. "
		echo "INFO      | $(date) | ----------------------------------- "
	elif [[ "$MY_RESUME_SWITCH" = "1" ]]; then
		# --------------------------------------------------
		# -- STEP #4 ALT: CONDUCT POST-PROCESSING OF RAxML RESULTS.
		# --------------------------------------------------
		echo "INFO      | $(date) | ----------------------------------- "
		echo "INFO      | $(date) | # Step #4: RAxML post-processing analyses. "
		echo "INFO      | $(date) | ----------------------------------- "
	fi

	# ---------------------------------------------------------- #
	# --------------- getGeneTrees.sh FUNCTION ----------------- #
	# ---------------------------------------------------------- #
	getGeneTrees_function () {

	echo "INFO      | $(date) | Organizing gene trees and making final output file containing all trees... "
	echo "INFO      | $(date) | Making list of ML gene trees generated by RAxML... "

	ls **/RAxML_bestTree.raxml_out > geneTrees.list ;

	# Assign gene tree list to variable
	MY_GENE_TREE_LIST="$(cat ./geneTrees.list)";

	# ORGANIZE GENE TREES INTO ONE LOCATION
	# --------------------------------------------------
	# Place all inferred gene trees into a single "gene_trees" folder in the current
	# working directory. However, all the gene tree files have the same name. So, in order
	# to do this, we have to give each gene tree a name that matches the corresponding run
	# folder, i.e. locus. We can rename each file right after downloading it.
	# --------------------------------------------------
	mkdir ./gene_trees/ ;

	echo "INFO      | $(date) | Copying *ALL* ML gene trees to 'gene_trees' folder in current directory for post-processing..."
	(
		for j in $MY_GENE_TREE_LIST; do
			echo "$j"
			cp "$j" ./gene_trees/ ;
			MY_LOCUS_NAME="$(echo "$j" | sed 's/\/[A-Za-z.\_\-]*//g')";
			cp ./gene_trees/RAxML_bestTree.raxml_out ./gene_trees/"$MY_LOCUS_NAME"_RAxML_best.tre ;
			if [[ -s ./gene_trees/RAxML_bestTree.raxml_out ]]; then rm ./gene_trees/RAxML_bestTree.raxml_out ; fi
		done
	)

	echo "INFO      | $(date) | Making final output file 'besttrees.tre' containing best ML trees from all runs/loci..."
	(
		for k in ./gene_trees/*; do
			echo "$k"
			cat "$k" >> ./besttrees.tre ;
		done
	)
}
	# Run the function!
	getGeneTrees_function


	# ---------------------------------------------------------- #
	# --------------- getBootTrees.sh FUNCTION ----------------- #
	# ---------------------------------------------------------- #
	getBootTrees_function () {

	echo "INFO      | $(date) | Organizing bootstrap trees and making final output file containing all trees... "
	echo "INFO      | $(date) | Making list of ML bootstrap trees generated by RAxML... "

	ls **/RAxML_bootstrap.raxml_out > bootTrees.list ;

	# Assign bootstrap tree list to variable
	MY_BOOT_TREE_LIST="$(cat ./bootTrees.list)";

	# ORGANIZE BOOTSTRAP TREES INTO ONE LOCATION
	# --------------------------------------------------
	# Place all inferred bootstrap tree files into a single "bootstrap_trees" folder in 
	# working directory. However, all the boot tree files have the same name. So, in order
	# to do this, we have to give each boot tree file a name that matches the corresponding
	# run folder, i.e. locus. We can rename each file right after downloading it.
	# --------------------------------------------------
	mkdir ./bootstrap_trees ;

	echo "INFO      | $(date) | Copying *ALL* ML bootstrap trees to 'bootstrap_trees' folder in current directory for post-processing..."
	(
		for l in $MY_BOOT_TREE_LIST; do
			echo "$l"
			cp "$l" ./bootstrap_trees/ ;
			MY_LOCUS_NAME="$(echo "$l" | sed 's/\/[A-Za-z.\_\-]*//g')";
			cp ./bootstrap_trees/RAxML_bootstrap.raxml_out ./bootstrap_trees/"$MY_LOCUS_NAME"_RAxML_boot.tre ;
			if [[ -s ./bootstrap_trees/RAxML_bootstrap.raxml_out ]]; then rm ./bootstrap_trees/RAxML_bootstrap.raxml_out ; fi
		done
	)

	echo "INFO      | $(date) | Making final output file 'boottrees.tre' containing best ML trees from all runs/loci..."
	(
		for m in ./bootstrap_trees/*; do
			echo "$m"
			cat "$m" >> ./boottrees.tre ;
		done
	)

	echo "INFO      | $(date) | Making final list of ML bootstrap trees ('final_bootTrees.list') in bootstrap_trees directory..."
	ls ./bootstrap_trees/*.tre > final_bootTrees.list ;
}
	# Run the function!
	getBootTrees_function


	# ---------------------------------------------------------- #
	# -------------- getBipartTrees.sh FUNCTION ---------------- #
	# ---------------------------------------------------------- #
	getBipartTrees_function () {

	echo "INFO      | $(date) | Organizing bipartitions trees (with bootstrap proportion labels) and making final output file containing all bipartitions trees... "
	ls **/RAxML_bipartitions.raxml_out > bipartTrees.list ;

	# Assign bootstrap tree list to variable
	MY_BIPART_TREE_LIST="$(cat ./bipartTrees.list)";

	# ORGANIZE BIPARTITIONS TREES INTO ONE LOCATION
	# --------------------------------------------------
	mkdir ./bipartitions_trees

	echo "INFO      | $(date) | Copying *ALL* RAxML bootstrap bipartitions trees to 'bipartitions_trees' folder in current directory for post-processing..."
	(
		for l in ${MY_BIPART_TREE_LIST}; do
			echo "$l"
			cp "$l" ./bipartitions_trees/ ;
			MY_LOCUS_NAME="$(echo $l | sed 's/\/[A-Za-z.\_\-]*//g')";
			cp ./bipartitions_trees/RAxML_bipartitions.raxml_out ./bipartitions_trees/"$MY_LOCUS_NAME"_RAxML_bipartitions.tre ;
			if [[ -s ./bipartitions_trees/RAxML_bipartitions.raxml_out ]]; then rm ./bipartitions_trees/RAxML_bipartitions.raxml_out ; fi
		done
	)

	echo "INFO      | $(date) | Making final output file 'biparttrees.tre' containing RAxML bipartitions trees from all runs/loci..."
	(
		for m in ./bipartitions_trees/*; do
			echo "$m"
			cat "$m" >> ./biparttrees.tre ;
		done
	)

	echo "INFO      | $(date) | Making final list of RAxML bipartitions trees ('final_bipartTrees.list') in bipartitions_trees directory..."
	ls ./bipartitions_trees/*.tre > final_bipartTrees.list ;
}
	# Run the function!
	getBipartTrees_function


fi
#######


# Multi-PHYLIP run
# ------------------------------------------------------
# Run on multiple PHYLIP files (extension '.phy') in the 
# current working directory.
# ------------------------------------------------------

#######
if [[ "$STARTING_FILE_TYPE" = "2" ]]; then

	echo "INFO      | $(date) | Starting MAGNET pipeline... "
	echo "INFO      | $(date) | Running with the following options: "
	echo "INFO      | $(date) | - NEXUS file, <inputNEXUS> = ${MY_NEXUS} "
	echo "INFO      | $(date) | - Starting <fileType> = ${STARTING_FILE_TYPE} "
	echo "INFO      | $(date) | - RAxML <executable> = ${MY_RAXML_EXECUTABLE} "
	echo "INFO      | $(date) | - Bootstrap reps, <numBootstraps> = ${MY_NUM_BOOTREPS} "
	echo "INFO      | $(date) | - RAxML model, <raxmlModel> = ${MY_RAXML_MODEL} "
	echo "INFO      | $(date) | - <simpleModel> option = ${MY_SIMPLE_MODEL} "
	echo "INFO      | $(date) | - <gapThreshold> option = ${MY_GAP_THRESHOLD} "
	echo "INFO      | $(date) | - <indivMissingData> option = ${MY_INDIV_MISSING_DATA} "
	echo "INFO      | $(date) | - Outgroup taxon, <outgroup> = ${MY_OUTGROUP} "
	echo "INFO      | $(date) | - RAxML output name = ${MY_OUTPUT_NAME} "
	echo "INFO      | $(date) | - Resume switch (--resume) = ${MY_RESUME_SWITCH} "

# --------------------------------------------------
# -- STEP #1: SETUP.
# --------------------------------------------------
	echo "INFO      | $(date) | ----------------------------------- "
	echo "INFO      | $(date) | # Step #1: Set up workspace (e.g. functions, working directory) and check machine type. " # | tee -a "$MY_OUTPUT_FILE_SWITCH"
	echo "INFO      | $(date) | ----------------------------------- "

	############ SET WORKING DIRECTORY AND CHECK MACHINE TYPE
	USER_SPEC_PATH="$(printf '%q\n' "$(pwd)")";
	echoCDWorkingDir
	MY_WORKING_DIR="$(pwd)"
	checkMachineType

# --------------------------------------------------
# -- STEP #2: HANDLE INPUT FILE(S).
# --------------------------------------------------
	echo "INFO      | $(date) | ----------------------------------- "
	echo "INFO      | $(date) | # Step #2: Input single NEXUS (or G-PhoCS-formatted) file, or multiple PHYLIP files. "
	echo "INFO      | $(date) | ----------------------------------- "
	echo "INFO      | $(date) | For -f 1 or -f 2, if '.gphocs' input file present, continue; else convert NEXUS file "
	echo "INFO      | $(date) | to G-PhoCS format using NEXUS2gphocs code. If -f 3, then run multiple PHYLIP files in  "
	echo "INFO      | $(date) | RAxML."

	MY_PHYLIP_ALIGNMENTS=./*.phy ;		## Assign PHYLIP-formatted multilocus gene / genomic/SNP / RAD locus sequence alignment files (e.g. output by gphocs2multiPhylip.sh shell script) in run directory to variable.

	# ---------------------------------------------------------- #
	# ------------- MultiRAxMLPrepper.sh FUNCTION -------------- #
	# ---------------------------------------------------------- #
	MultiRAxMLPrepper_function () {

	if [[ "$MY_RESUME_SWITCH" = "0" ]]; then
	
		# --------------------------------------------------
		# -- STEP #3: MAKE RAxML RUN FOLDERS.
		# --------------------------------------------------
		echo "INFO      | $(date) | ----------------------------------- "
		echo "INFO      | $(date) | # Step #3: Make and check RAxML run folders. "
		echo "INFO      | $(date) | ----------------------------------- "
	
		MY_N_PHYLIP_FILES="$(ls -1 ./*.phy 2>/dev/null | wc -l | sed 's/\ //g')";
	
		# MAKE RAXML RUN FOLDERS, ONE FOR EACH GENE/LOCUS
		# --------------------------------------------------
		# Loop through the input .phy files and do the following for each file: (A) generate one 
		# folder per .phy file with the same name as the file, only minus the extension, then 
		# (B) move input .phy file into corresponding folder.
		# --------------------------------------------------
		count=1
		(
			for i in $MY_PHYLIP_ALIGNMENTS; do
				MY_BASENAME="$(basename "$i" '.phy')";
				mkdir "$MY_BASENAME"/ ;
				cp "$i" "$MY_BASENAME"/ ;
				echo "$((count++))" >/dev/null 2>&1 ;
			done
		)
	
		# CHECK RUN FOLDERS
		# --------------------------------------------------
		# Setup and run check on the number of run folders created by the program:
		# --------------------------------------------------
		export MY_TOTAL_FILECOUNT="$(find . -type f | wc -l)";
		export MY_TOTAL_DIRCOUNT="$(find . -type d | wc -l)";
		#MY_NUM_RUN_FOLDERS="$(ls ./*/*.phy | wc -l | perl -pe 's/\t//g; s/\ //g')";
	
		echo "INFO      | $(date) | Number of run folders created: ${count} "
	
		if [[ "$count" = "$MY_N_PHYLIP_FILES" ]]; then
			echo "INFO      | $(date) | Folder check PASSED: number of run folders matches number of PHYLIP alignments. "
		else
			echo "WARNING   | $(date) | Folder check FAILED: number of run folders (${count}) does NOT match the number of PHYLIP alignments (${MY_N_PHYLIP_FILES}). This may cause errors. "
		fi
	
	elif [[ "$MY_RESUME_SWITCH" = "1" ]]; then
		if [[ "$count" = "$MY_N_PHYLIP_FILES" ]]; then
			echo "IMPORTANT${EP}| $(date) | Resuming a previous/existing run in current working dir. Skipping MultiRAxMLPrepper, using available run folders... "
			echo "INFO      | $(date) | Folder check PASSED: number of run folders matches number of PHYLIP alignments. "
		else
			echo "WARNING   | $(date) | Folder check FAILED: number of run folders does NOT match the number of PHYLIP alignments. There may be errors. "
		fi
	
	fi
}
	# Run the function!
	MultiRAxMLPrepper_function


	# ---------------------------------------------------------- #
	# ---------------- RAxMLRunner.sh FUNCTION ----------------- #
	# ---------------------------------------------------------- #
	RAxMLRunner_function () {
	
	if [[ "$MY_RESUME_SWITCH" = "0" ]]; then
	
		# --------------------------------------------------
		# -- STEP #4: ESTIMATE BEST ML GENE TREES.
		# --------------------------------------------------
		echo "INFO      | $(date) | ----------------------------------- "
		echo "INFO      | $(date) | # Step #4: Estimate best maximum-likelihood (ML) gene trees. "
		echo "INFO      | $(date) | ----------------------------------- "
		echo "INFO      | $(date) | Looping through and analyzing contents of each run folder in RAxML... "
		echo "INFO      | $(date) | ----------------------------------- "
		# Each folder is set with the locus name corresponding to the locus' position in the
		# original .gphocs alignment (which, if output by pyRAD, is simply in the order in which
		# the loci were logged to file by pyRAD, no special order). Also, each folder contains
		# one .phy file carrying the same basename as the folder name, e.g. "locus0.phy". So,
		# all we need to do here is loop through each folder and call RAxML to run using its
		# contents as the input file, as follows:
		(
			for i in ./*/; do
				if [[ "$i" != "./archive/" ]] && [[ "$i" != "./bad_genes/" ]] && [[ "$i" != "./R/" ]] && [[ "$i" != "./shell/" ]] && [[ "$i" != "./perl/" ]] && [[ "$i" != "./orig_phylip/" ]] && [[ "$i" != "./phylip/" ]] && [[ "$i" != "./orig_fasta/" ]] && [[ "$i" != "./fasta/" ]] && [[ "$i" != "./phylip_files/" ]]; then
					echo "$i"
					cd "$i";
					LOCUS_NAME="$(echo "$i" | sed 's/\.\///g; s/\/$//g')";  # NOTE: not currently using "$LOCUS_NAME" here, but leave for now, bc may need to use it later...
				#
					if [[ "$MY_OUTGROUP" = "NULL" ]] && [[ "$MY_SIMPLE_MODEL" = "NULL" ]]; then
						"$MY_RAXML_EXECUTABLE" -f a -x "$RANDOM""$RANDOM" -p "$RANDOM""$RANDOM" -# ${MY_NUM_BOOTREPS} -m ${MY_RAXML_MODEL} -s ./*.phy -n ${MY_OUTPUT_NAME}
					fi
				#
					if [[ "$MY_OUTGROUP" != "NULL" ]] && [[ "$MY_SIMPLE_MODEL" = "NULL" ]]; then
						"$MY_RAXML_EXECUTABLE" -f a -x "$RANDOM""$RANDOM" -p "$RANDOM""$RANDOM" -# ${MY_NUM_BOOTREPS} -m ${MY_RAXML_MODEL} -s ./*.phy -o ${MY_OUTGROUP} -n ${MY_OUTPUT_NAME}
					fi
				#
					if [[ "$MY_OUTGROUP" = "NULL" ]] && [[ "$MY_SIMPLE_MODEL" != "NULL" ]]; then
						"$MY_RAXML_EXECUTABLE" -f a -x "$RANDOM""$RANDOM" -p "$RANDOM""$RANDOM" -# ${MY_NUM_BOOTREPS} -m ${MY_RAXML_MODEL} -s ./*.phy --${MY_SIMPLE_MODEL} -n ${MY_OUTPUT_NAME}
					fi
				#
					if [[ "$MY_OUTGROUP" != "NULL" ]] && [[ "$MY_SIMPLE_MODEL" != "NULL" ]]; then
						"$MY_RAXML_EXECUTABLE" -f a -x "$RANDOM""$RANDOM" -p "$RANDOM""$RANDOM" -# ${MY_NUM_BOOTREPS} -m ${MY_RAXML_MODEL} -s ./*.phy --${MY_SIMPLE_MODEL} -o ${MY_OUTGROUP} -n ${MY_OUTPUT_NAME}
					fi
					cd ..;
				fi
			done
		)
	
		# Here: adding loop code to move all .phy files remaining in the current working 
		# directory, after Step #3 of the pipeline, to a new folder called "phylip_files". This
		# is done here because if the phylip_files folder is present at the end of Step #3,
		# then RAxML will also try to estimate a gene tree for .phy file(s) in this folder during
		# Step #5 of the pipeline above.
		mkdir ./phylip_files
		(
			for i in $MY_PHYLIP_ALIGNMENTS; do
				echo "$i"
				mv "$i" ./phylip_files/ ;
			done
		)
	
	elif [[ "$MY_RESUME_SWITCH" = "1" ]]; then
	
	echo "INFO      | $(date) | Step #3: Resuming gene tree estimation. Run on remaining/incomplete run folders, skip those with completed RAxML runs. "
	
		(
			for i in ./*/; do
				if [[ "$i" != "./archive/" ]] && [[ "$i" != "./bad_genes/" ]] && [[ "$i" != "./R/" ]] && [[ "$i" != "./shell/" ]] && [[ "$i" != "./perl/" ]] && [[ "$i" != "./orig_phylip/" ]] && [[ "$i" != "./phylip/" ]] && [[ "$i" != "./orig_fasta/" ]] && [[ "$i" != "./fasta/" ]] && [[ "$i" != "./phylip_files/" ]]; then
					cd "$i";
					LOCUS_NAME="$(echo "$i" | sed 's/\.\///g; s/\/$//g')"; # NOTE: not currently using $LOCUS_NAME here, but leave for now, bc may need to use it later...
				#
					if [[ ! -s ./RAxML_info.raxml_out ]]; then
						echo "$i"
						if [[ "$MY_OUTGROUP" = "NULL" ]] && [[ "$MY_SIMPLE_MODEL" = "NULL" ]]; then
							"$MY_RAXML_EXECUTABLE" -f a -x "$RANDOM""$RANDOM" -p "$RANDOM""$RANDOM" -# ${MY_NUM_BOOTREPS} -m ${MY_RAXML_MODEL} -s ./*.phy -n ${MY_OUTPUT_NAME}
						fi
					#
						if [[ "$MY_OUTGROUP" != "NULL" ]] && [[ "$MY_SIMPLE_MODEL" = "NULL" ]]; then
							"$MY_RAXML_EXECUTABLE" -f a -x "$RANDOM""$RANDOM" -p "$RANDOM""$RANDOM" -# ${MY_NUM_BOOTREPS} -m ${MY_RAXML_MODEL} -s ./*.phy -o ${MY_OUTGROUP} -n ${MY_OUTPUT_NAME}
						fi
					#
						if [[ "$MY_OUTGROUP" = "NULL" ]] && [[ "$MY_SIMPLE_MODEL" != "NULL" ]]; then
							"$MY_RAXML_EXECUTABLE" -f a -x "$RANDOM""$RANDOM" -p "$RANDOM""$RANDOM" -# ${MY_NUM_BOOTREPS} -m ${MY_RAXML_MODEL} -s ./*.phy --${MY_SIMPLE_MODEL} -n ${MY_OUTPUT_NAME}
						fi
					#
						if [[ "$MY_OUTGROUP" != "NULL" ]] && [[ "$MY_SIMPLE_MODEL" != "NULL" ]]; then
							"$MY_RAXML_EXECUTABLE" -f a -x "$RANDOM""$RANDOM" -p "$RANDOM""$RANDOM" -# ${MY_NUM_BOOTREPS} -m ${MY_RAXML_MODEL} -s ./*.phy --${MY_SIMPLE_MODEL} -o ${MY_OUTGROUP} -n ${MY_OUTPUT_NAME}
						fi
					fi
					cd ..;
				fi
			done
		)
	
		if [[ ! -s ./phylip_files/ ]]; then
			mkdir ./phylip_files/ ;
		fi
		(
			for i in $MY_PHYLIP_ALIGNMENTS; do
				echo "$i"
				mv "$i" ./phylip_files/ ;
			done
		)
	fi
}
	# Run the function!
	RAxMLRunner_function


	if [[ "$MY_RESUME_SWITCH" = "0" ]]; then
		# --------------------------------------------------
		# -- STEP #5: CONDUCT POST-PROCESSING OF RAxML RESULTS.
		# --------------------------------------------------
		echo "INFO      | $(date) | ----------------------------------- "
		echo "INFO      | $(date) | # Step #5: RAxML post-processing analyses. "
		echo "INFO      | $(date) | ----------------------------------- "
	elif [[ "$MY_RESUME_SWITCH" = "1" ]]; then
		# --------------------------------------------------
		# -- STEP #4 ALT: CONDUCT POST-PROCESSING OF RAxML RESULTS.
		# --------------------------------------------------
		echo "INFO      | $(date) | ----------------------------------- "
		echo "INFO      | $(date) | # Step #4: RAxML post-processing analyses. "
		echo "INFO      | $(date) | ----------------------------------- "
	fi


	# ---------------------------------------------------------- #
	# --------------- getGeneTrees.sh FUNCTION ----------------- #
	# ---------------------------------------------------------- #
	getGeneTrees_function () {

	echo "INFO      | $(date) | Organizing gene trees and making final output file containing all trees... "
	echo "INFO      | $(date) | Making list of ML gene trees generated by RAxML... "

	ls **/RAxML_bestTree.raxml_out > geneTrees.list ;

	# Assign gene tree list to variable
	MY_GENE_TREE_LIST="$(cat ./geneTrees.list)";

	# ORGANIZE GENE TREES INTO ONE LOCATION
	# --------------------------------------------------
	# Place all inferred gene trees into a single "gene_trees" folder in the current
	# working directory. However, all the gene tree files have the same name. So, in order
	# to do this, we have to give each gene tree a name that matches the corresponding run
	# folder, i.e. locus. We can rename each file right after downloading it.
	# --------------------------------------------------
	mkdir ./gene_trees/ ;

	echo "INFO      | $(date) | Copying *ALL* ML gene trees to 'gene_trees' folder in current directory for post-processing..."
	(
		for j in $MY_GENE_TREE_LIST; do
			echo "$j"
			cp "$j" ./gene_trees/ ;
			MY_LOCUS_NAME="$(echo "$j" | sed 's/\/[A-Za-z.\_\-]*//g')";
			cp ./gene_trees/RAxML_bestTree.raxml_out ./gene_trees/"$MY_LOCUS_NAME"_RAxML_best.tre ;
			if [[ -s ./gene_trees/RAxML_bestTree.raxml_out ]]; then rm ./gene_trees/RAxML_bestTree.raxml_out ; fi
		done
	)

	echo "INFO      | $(date) | Making final output file 'besttrees.tre' containing best ML trees from all runs/loci..."
	(
		for k in ./gene_trees/*; do
			echo "$k"
			cat "$k" >> ./besttrees.tre ;
		done
	)
}
	# Run the function!
	getGeneTrees_function


	# ---------------------------------------------------------- #
	# --------------- getBootTrees.sh FUNCTION ----------------- #
	# ---------------------------------------------------------- #
	getBootTrees_function () {

	echo "INFO      | $(date) | Organizing bootstrap trees and making final output file containing all trees... "
	echo "INFO      | $(date) | Making list of ML bootstrap trees generated by RAxML... "

	ls **/RAxML_bootstrap.raxml_out > bootTrees.list ;

	# Assign bootstrap tree list to variable
	MY_BOOT_TREE_LIST="$(cat ./bootTrees.list)";

	# ORGANIZE BOOTSTRAP TREES INTO ONE LOCATION
	# --------------------------------------------------
	# Place all inferred bootstrap tree files into a single "bootstrap_trees" folder in 
	# working directory. However, all the boot tree files have the same name. So, in order
	# to do this, we have to give each boot tree file a name that matches the corresponding
	# run folder, i.e. locus. We can rename each file right after downloading it.
	# --------------------------------------------------
	mkdir ./bootstrap_trees/ ;

	echo "INFO      | $(date) | Copying *ALL* ML bootstrap trees to 'bootstrap_trees' folder in current directory for post-processing..."
	(
		for l in ${MY_BOOT_TREE_LIST}; do
			echo "$l"
			cp "$l" ./bootstrap_trees/ ;
			MY_LOCUS_NAME="$(echo "$l" | sed 's/\/[A-Za-z.\_\-]*//g')";
			cp ./bootstrap_trees/RAxML_bootstrap.raxml_out ./bootstrap_trees/"$MY_LOCUS_NAME"_RAxML_boot.tre ;
			if [[ -s ./bootstrap_trees/RAxML_bootstrap.raxml_out ]]; then rm ./bootstrap_trees/RAxML_bootstrap.raxml_out ; fi
		done
	)

	echo "INFO      | $(date) | Making final output file 'boottrees.tre' containing best ML trees from all runs/loci..."
	(
		for m in ./bootstrap_trees/*; do
			echo "$m"
			cat "$m" >> ./boottrees.tre ;
		done
	)

	echo "INFO      | $(date) | Making final list of ML bootstrap trees ('final_bootTrees.list') in bootstrap_trees directory..."
	ls ./bootstrap_trees/*.tre > final_bootTrees.list ;
}
	# Run the function!
	getBootTrees_function


	# ---------------------------------------------------------- #
	# -------------- getBipartTrees.sh FUNCTION ---------------- #
	# ---------------------------------------------------------- #
	getBipartTrees_function () {

	echo "INFO      | $(date) | Organizing bipartitions trees (with bootstrap proportion labels) and making final output file containing all bipartitions trees... "
	ls **/RAxML_bipartitions.raxml_out > bipartTrees.list ;

	# Assign bootstrap tree list to variable
	MY_BIPART_TREE_LIST="$(cat ./bipartTrees.list)";

	############ ORGANIZE BIPARTITIONS TREES INTO ONE LOCATION
	mkdir ./bipartitions_trees

	echo "INFO      | $(date) | Copying *ALL* RAxML bootstrap bipartitions trees to 'bipartitions_trees' folder in current directory for post-processing..."
	(
		for l in ${MY_BIPART_TREE_LIST}; do
			echo "$l"
			cp "$l" ./bipartitions_trees/ ;
			MY_LOCUS_NAME="$(echo "$l" | sed 's/\/[A-Za-z.\_\-]*//g')";
			cp ./bipartitions_trees/RAxML_bipartitions.raxml_out ./bipartitions_trees/"$MY_LOCUS_NAME"_RAxML_bipartitions.tre ;
			if [[ -s ./bipartitions_trees/RAxML_bipartitions.raxml_out ]]; then rm ./bipartitions_trees/RAxML_bipartitions.raxml_out ; fi
		done
	)

	echo "INFO      | $(date) | Making final output file 'biparttrees.tre' containing RAxML bipartitions trees from all runs/loci..."
	(
		for m in ./bipartitions_trees/*; do
			echo "$m"
			cat "$m" >> ./biparttrees.tre ;
		done
	)

	echo "INFO      | $(date) | Making final list of RAxML bipartitions trees ('final_bipartTrees.list') in bipartitions_trees directory..."
	ls ./bipartitions_trees/*.tre > final_bipartTrees.list ;
}
	# Run the function!
	getBipartTrees_function

fi
#######


	if [[ "$STARTING_FILE_TYPE" = "1" ]]; then
		if [[ "$MY_RESUME_SWITCH" = "0" ]]; then
			# --------------------------------------------------
			# -- STEP #7: CLEAN UP WORKSPACE
			# --------------------------------------------------
			# Clean up workspace by removing temporary files generated during run. 
			# --------------------------------------------------
			echo "INFO      | $(date) | ----------------------------------- "
			echo "INFO      | $(date) | # Step #7: Clean up workspace by removing temporary files generated during run. "
			echo "INFO      | $(date) | ----------------------------------- "
			## Remove arguments file generated when parsing the options:
			echo "INFO      | $(date) | Removing arguments file generated when parsing the options..."
			if [[ -s ./args.txt ]]; then rm ./args.txt ; fi
		elif [[ "$MY_RESUME_SWITCH" = "1" ]]; then
			# --------------------------------------------------
			# -- STEP #5: CLEAN UP WORKSPACE
			# --------------------------------------------------
			# Clean up workspace by removing temporary files generated during run. 
			# --------------------------------------------------
			echo "INFO      | $(date) | ----------------------------------- "
			echo "INFO      | $(date) | # Step #5: Clean up workspace by removing temporary files generated during run. "
			echo "INFO      | $(date) | ----------------------------------- "
			## Remove arguments file generated when parsing the options:
			echo "INFO      | $(date) | Removing arguments file generated when parsing the options..."
			if [[ -s ./args.txt ]]; then rm ./args.txt ; fi
		fi
	fi
	if [[ "$STARTING_FILE_TYPE" = "2" ]]; then
		if [[ "$MY_RESUME_SWITCH" = "0" ]]; then
			# --------------------------------------------------
			# -- STEP #6: CLEAN UP WORKSPACE
			# --------------------------------------------------
			# Clean up workspace by removing temporary files generated during run. 
			# --------------------------------------------------
			echo "INFO      | $(date) | ----------------------------------- "
			echo "INFO      | $(date) | # Step #6: Clean up workspace by removing temporary files generated during run. "
			echo "INFO      | $(date) | ----------------------------------- "
			## Remove arguments file generated when parsing the options:
			echo "INFO      | $(date) | Removing arguments file generated when parsing the options..."
			if [[ -s ./args.txt ]]; then rm ./args.txt ; fi
		elif [[ "$MY_RESUME_SWITCH" = "1" ]]; then
			# --------------------------------------------------
			# -- STEP #5: CLEAN UP WORKSPACE
			# --------------------------------------------------
			# Clean up workspace by removing temporary files generated during run. 
			# --------------------------------------------------
			echo "INFO      | $(date) | ----------------------------------- "
			echo "INFO      | $(date) | # Step #5: Clean up workspace by removing temporary files generated during run. "
			echo "INFO      | $(date) | ----------------------------------- "
			## Remove arguments file generated when parsing the options:
			echo "INFO      | $(date) | Removing arguments file generated when parsing the options..."
			if [[ -s ./args.txt ]]; then rm ./args.txt ; fi
		fi
	fi
	

	echo "INFO      | $(date) | Done."
	echo "----------------------------------------------------------------------------------------------------------"
	echo ""

	# END DEBUG MODE
	# --------------------------------------------------
	if [[ "$MY_DEBUG_MODE_SWITCH" != "0" ]]; then set +xv; fi

##########################################################################################
######################################### END ############################################

}



############ CREATE USAGE & HELP TEXTS
USAGE="
Usage: $(basename "$0") [OPTION]...

 ${bold}Options:${reset}
  -f, --filetype     fileType (def: 1; also: 2) starting file type; if 1, script expects as
                     stdin a single input NEXUS file in the current directory; if 2, then
                     script expects multiple input PHYLIP files in current directory
  -i, --input        inputNEXUS (def: NULL) input NEXUS file (mandatory for -f 1)
  -e, --exec         executable (def: $MY_RAXML_EXECUTABLE) name of RAxML executable available
                     from user's command line interface
  -b, --boot         numBootstraps (def: $MY_NUM_BOOTREPS) RAxML bootstrap pseudoreplicates
  -r, --raxmlmodel   raxmlModel (def: $MY_RAXML_MODEL; other: GTRGAMMAI, GTRCAT, GTRCATI)
  -s, --simplemodel  simpleModel (def: $MY_SIMPLE_MODEL; other: JC69, K80, HKY85) specifies 
                     simple DNA substitution model that will override any other model (even 
                     across partitions)
  -g, --gapthresh    gapThreshold (def: $MY_GAP_THRESHOLD=essentially zero gaps allowed unless 
                     >1000 individuals; takes float proportion value) gap threshold value
  -m, --missing      indivMissingData (def: $MY_INDIV_MISSING_DATA=allowed; 0=removed) missing  
                     data setting
  -o, --outgroup     outgroup (def: NULL) outgroup given as single taxon name (tip label) or
                     comma-separted list
  -h, --help         echo this help text and exit
  -H, --Help         echo verbose help text and exit
  -V, --version      echo version and exit
  -R, --resume       resume (def: 0, off; 1, on) option allowing user to resume a previous 
                     MAGNET run in the current working directory
  -d, --debug        debug (def: 0, off; 1, on) run function in Bash debug mode

 ${bold}OVERVIEW${reset}
 The goal of MAGNET is to infer a maximum-likelihood (ML) gene tree in RAxML for each of 
 multiple loci, starting from one or multiple DNA sequence alignment input files. If supplied
 with a single G-PhoCS ('*.gphocs') or NEXUS ('*.nex') data file (using -f1 or -i <inputNEXUS> 
 -f1 options), then each locus is split into a separate PHYLIP alignment file, and RAxML 
 (Stamatakis 2014) is run to infer gene trees for each locus. If a NEXUS datafile is supplied, 
 it is converted into G-PhoCS format (Gronau et al. 2011) while splitting loci into separate 
 interleaved sequence blocks based on information provided in a sets block at the end of the 
 NEXUS file (e.g. defined using 'charset' commands), which is mandatory. However, if -f2, then 
 the program will run in current directory, assuming it contains multiple PHYLIP-formatted 
 alignment files. Under this scenario, MAGNET will skip directly to running the PHYLIP files 
 in RAxML using user-specified options. 
	Sequence names may not include hyphen characters, or there could be issues. For detailed 
 information on MAGNET and its various dependencies, see 'README.md' file in the distribution 
 folder; however, it is key that dependencies are available from the command line interface. 
 Among the most important options is <resume> (-r|--resume, off by default), which tells the
 program to resume a previous MAGNET run in current directory, including detecting incomplete 
 RAxML run folders, and running RAxML without overwriting results from the previous run(s).

 ${bold}Usage examples:${reset}

    ./MAGNET -f 2 -b 100 -g 1 -m 1                        Run MAGNET with 100 bootstrap pseudo-
                                                          replicates, gaps allowed, missing 
                                                          data allowed, and the GTRGAMMA model
    ./MAGNET -f 2 -b 100 -s HKY85 -g 1 -m 1               Same as above, but using the simpler
                                                          HKY85 substitution model for all loci    
    ./MAGNET -f 2 -e raxmlHPC -b 100 -s HKY85 -g 1 -m 1   Same as above, but using raxmlHPC 
                                                          executable
    ./MAGNET -H                                           Show this help text and exit

 ${bold}CITATION${reset}
 Bagley, J.C. 2020. PIrANHA v0.4a4. GitHub repository, Available at:
	<https://github.com/justincbagley/piranha>.

 ${bold}REFERENCES${reset}
 Gronau, I., Hubisz, M.J., Gulko, B., Danko, C.G., Siepel, A. 2011. Bayesian inference of 
	ancient human demography from individual genome sequences. Nature Genetics, 43, 1031-1034.
 Stamatakis, A. 2014. RAxML version 8: a tool for phylogenetic analysis and post-analysis of 
	large phylogenies. Bioinformatics, 30, 1312-1313.

 Created by Justin Bagley on/before Aug 29 13:12:45 2016 -0700.
 Copyright (c) 2016-2020 Justin C. Bagley. All rights reserved.
"

VERBOSE_USAGE="
Usage: $(basename "$0") [OPTION]...

 ${bold}Options:${reset}
  -f, --filetype     fileType (def: 1; also: 2) starting file type; if 1, script expects as
                     stdin a single input NEXUS file in the current directory; if 2, then
                     script expects multiple input PHYLIP files in current directory
  -i, --input        inputNEXUS (def: NULL) input NEXUS file (mandatory for -f 1)
  -e, --exec         executable (def: $MY_RAXML_EXECUTABLE) name of RAxML executable available
                     from user's command line interface
  -b, --boot         numBootstraps (def: $MY_NUM_BOOTREPS) RAxML bootstrap pseudoreplicates
  -r, --raxmlmodel   raxmlModel (def: $MY_RAXML_MODEL; other: GTRGAMMAI, GTRCAT, GTRCATI)
  -s, --simplemodel  simpleModel (def: $MY_SIMPLE_MODEL; other: JC69, K80, HKY85) specifies 
                     simple DNA substitution model that will override any other model (even 
                     across partitions)
  -g, --gapthresh    gapThreshold (def: $MY_GAP_THRESHOLD=essentially zero gaps allowed unless 
                     >1000 individuals; takes float proportion value) gap threshold value
  -m, --missing      indivMissingData (def: $MY_INDIV_MISSING_DATA=allowed; 0=removed) missing  
                     data setting
  -o, --outgroup     outgroup (def: NULL) outgroup given as single taxon name (tip label) or
                     comma-separted list
  -h, --help         echo this help text and exit
  -H, --Help         echo verbose help text and exit
  -V, --version      echo version and exit
  -R, --resume       resume (def: 0, off; 1, on) option allowing user to resume a previous 
                     MAGNET run in the current working directory
  -d, --debug        debug (def: 0, off; 1, on) run function in Bash debug mode

 ${bold}OVERVIEW${reset}
 The goal of MAGNET is to infer a maximum-likelihood (ML) gene tree in RAxML for each of 
 multiple loci, starting from one or multiple DNA sequence alignment input files. If supplied
 with a single G-PhoCS ('*.gphocs') or NEXUS ('*.nex') data file (using -f1 or -i <inputNEXUS> 
 -f1 options), then each locus is split into a separate PHYLIP alignment file, and RAxML 
 (Stamatakis 2014) is run to infer gene trees for each locus. If a NEXUS datafile is supplied, 
 it is converted into G-PhoCS format (Gronau et al. 2011) while splitting loci into separate 
 interleaved sequence blocks based on information provided in a sets block at the end of the 
 NEXUS file (e.g. defined using 'charset' commands), which is mandatory. However, if -f2, then 
 the program will run in current directory, assuming it contains multiple PHYLIP-formatted 
 alignment files. Under this scenario, MAGNET will skip directly to running the PHYLIP files 
 in RAxML using user-specified options. 
	Sequence names may not include hyphen characters, or there could be issues. For detailed 
 information on MAGNET and its various dependencies, see 'README.md' file in the distribution 
 folder; however, it is key that dependencies are available from the command line interface. 
 Among the most important options is <resume> (-r|--resume, off by default), which tells the
 program to resume a previous MAGNET run in current directory, including detecting incomplete 
 RAxML run folders, and running RAxML without overwriting results from the previous run(s).

 ${bold}DETAILS${reset}
 The -f flag (also --filetype) specifies the starting fileType. If -f 1, then the mandatory 
	input is the name or path to the corresponding <inputNEXUS> starting file, which is 
	passed using the -i|--input flag. If -f 2, then mandatory input is the name or path to 
	the working directory (type '.' for current directory, or supply a relative or absolute 
	path).
 
 The -i flag (also --input) passess the name of the input NEXUS file, <inputNEXUS> parameter, 
    to the program.

 The -e flag (also --exec) sets the name of the RAxML executable that will be called. The 
    default executable name is 'raxml', but the user may wish to change this to something 
    specific to their install or parallelization needs (e.g. 'raxmlHPC-PTHREADS-SSE3'). The
    default setting should work on local machine or supercomputing cluster installs. However, 
    this should be tested beforehand by entering 'raxml' at the command prompt. On some 
    version fo Linux this yields the following error message:

	   'raxml: error while loading shared libraries: libmpi.so.12: cannot open shared object 
	    file: No such file or directory'.

	If this occurs, then Open MPI related libraries are installed in a non-standard location 
	and you will need to add this location to your LD_LIBRARY_PATH, e.g.:

	'export LD_LIBRARY_PATH=/usr/local/openmpi-1.8.1/lib:$LD_LIBRARY_PATH'
 
	See the following URL: for more insight into this problem: https://stackoverflow.com/
    questions/14769599/mpi-error-loading-shared-libraries. However, simply using a different 
    raxml executable that does not rely on these libararies will also immediately solve the 
    problem. In my experience, just setting MAGNET to call the 'raxmlHPC' executable immed-
    iately solves this issue on Mac and Linux (so also try simply running MAGNET with '-e 
    raxmlHPC' or '--exec raxmlHPC').
  
 The -b flag sets the number of boostrap pseudoreplicates for RAxML to perform while 
    estimating the gene tree for each locus. The default is 100; remove bootstrapping by 
    setting to 0.

 The -r flag sets the RAxML model for each locus. This uses the full default GTRGAMMA model,
	and at present it is not possible to vary the model across loci. If you want to use HKY
	or K80, you will need to use the -s flag (below).

 The -s flag sets a simple RAxML model for each locus/partition, which will override any
	model set using the -r flag above and apply to all partitions. In the current version of 
	RAxML, it is possible to specify the JC69, K80, and HKY85 models as overrides. By default,
	this option is turned off and the model set under the -r flag is used instead.

 The following two options are available **ONLY** if you are starting from a NEXUS input file:

	The -g flag supplies a 'gap threshold' to an R script, which deletes all column sites in 
	the DNA alignment with a proportion of gap characters '-' at or above the threshold value. 
	If no gap threshold is specified, all sites with gaps are removed by default. If end goal
	is to produce a file for G-PhoCS, you  will want to leave <gapThreshold> at the default. 
	However, if the next step in your pipeline involves converting from .gphocs to other data 
	formats, you will likely want to set <gapThreshold> = 1 (e.g. before converting to PHYLIP 
	format for RAxML). 

	The -m flag allows users to choose their level of tolerance for individuals with missing
	data. The default is <indivMissingData> = 1, allowing individuals with runs of 10 or more 
	missing nucleotide characters ('N') to be kept in the alignment. Alternatively, setting
	<indivMissingData> = 0 removes all such individuals from each locus; thus, while the input
	file would have had the same number of individuals across loci, the resulting file could
	have varying numbers of individuals for different loci.

 The -o flag sets the outgroup exactly the same way as that described in the RAxML v8 user's
	manual, as a single name or as a comma-separated list with no spaces between taxon names. 
	The first name in the list is prioritized, e.g. when members of the list are not mono-
    phyletic.

 -R | --resume is among the most important options available in MAGNET because it tells the 
	program to resume a previous run in current directory, including to detect incomplete run 
	subfolders and run RAxML there without overwriting results from run folders with finished 
	runs. The default setting is to run without this option.

 The -d flag runs this function in Bash debug mode (set -xv), which is intended for debugging
	for development purposes. If you find a bug, please contact the author at jbagley@jsu.edu.

 ${bold}Usage examples:${reset}

    ./MAGNET -f 2 -b 100 -g 1 -m 1                        Run MAGNET with 100 bootstrap pseudo-
                                                          replicates, gaps allowed, missing 
                                                          data allowed, and the GTRGAMMA model
    ./MAGNET -f 2 -b 100 -s HKY85 -g 1 -m 1               Same as above, but using the simpler
                                                          HKY85 substitution model for all loci    
    ./MAGNET -f 2 -e raxmlHPC -b 100 -s HKY85 -g 1 -m 1   Same as above, but using raxmlHPC 
                                                          executable
    ./MAGNET -H                                           Show this help text and exit

 ${bold}CITATION${reset}
 Bagley, J.C. 2020. PIrANHA v0.4a4. GitHub repository, Available at:
	<https://github.com/justincbagley/piranha>.

 ${bold}REFERENCES${reset}
 Gronau, I., Hubisz, M.J., Gulko, B., Danko, C.G., Siepel, A. 2011. Bayesian inference of 
	ancient human demography from individual genome sequences. Nature Genetics, 43, 1031-1034.
 Stamatakis, A. 2014. RAxML version 8: a tool for phylogenetic analysis and post-analysis of 
	large phylogenies. Bioinformatics, 30, 1312-1313.

 Created by Justin Bagley on/before Aug 29 13:12:45 2016 -0700.
 Copyright (c) 2016-2020 Justin C. Bagley. All rights reserved.
"

if [[ -z "$*" ]]; then
	echo "$USAGE"
	exit
fi

if [[ "$1" == "-h" ]] || [[ "$1" == "-help" ]]; then
	echo "$USAGE"
	exit
fi

if [[ "$1" == "-H" ]] || [[ "$1" == "-Help" ]]; then
	echo "$VERBOSE_USAGE"
	exit
fi

if [[ "$1" == "-v" ]] || [[ "$1" == "-V" ]] || [[ "$1" == "--version" ]]; then
	echo "$(basename "$0") $VERSION";
	exit
fi

############ CHECK ARGUMENTS
	# echo "$@"; echo "$#"; echo "$1" 
	# for i in "$@"; do
	# 	echo "$i";
	# done
	# MY_ARGS="$(echo "$@" | perl -pe $'s/\ /\n/')"
	# echo "$MY_ARGS"


############ CLEAN WORKING DIR, CAPTURE ARGUMENTS, SEND TO FILE FOR PARSING
	if [[ -s ./args.tmp ]]; then rm ./args.tmp ; fi ;
	if [[ -s ./args.txt ]]; then rm ./args.txt ; fi ;
	ALL_MY_ARGUMENTS="$(echo "$@")"
	echo "$ALL_MY_ARGUMENTS" > ./args.txt ;
	perl -p -i -e $'s/\-/\n\-/g' ./args.txt ;
	perl -p -i -e $'s/\-filetype/\-\-filetype/g' ./args.txt ;
	perl -p -i -e $'s/\-input/\-\-input/g' ./args.txt ;
	perl -p -i -e $'s/\-exec/\-\-exec/g' ./args.txt ;
#	perl -p -i -e $'s/\-part/\-\-part/g' ./args.txt ;
	perl -p -i -e $'s/\-boot/\-\-boot/g' ./args.txt ;
	perl -p -i -e $'s/\-raxmlmodel/\-\-raxmlmodel/g' ./args.txt ;
	perl -p -i -e $'s/\-simplemodel/\-\-simplemodel/g' ./args.txt ;
	perl -p -i -e $'s/\-outgroup/\-\-outgroup/g' ./args.txt ;
	perl -p -i -e $'s/\-name/\-\-name/g' ./args.txt ;
	perl -p -i -e $'s/\-resume/\-\-resume/g' ./args.txt ;
	perl -p -i -e $'s/\-debug/\-\-debug/g' ./args.txt ;


############ MANUALLY PARSE THE OPTIONS FROM ARGS

### SET OPTIONS TO DEFAULT VALUES, EXCEPT WHERE VALUES WERE READ IN FROM USER ARGS
	if [[  "$(grep -h '\-f' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]] && [[  "$(grep -h '\-\-filetype' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]]; then
		STARTING_FILE_TYPE=1 ;
	elif [[  "$(grep -h '\-f' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]] && [[  "$(grep -h '\-\-filetype' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]]; then
		MY_ARG="$(grep -h '\-f' ./args.txt | perl -pe 's/\-f//g' | perl -pe 's/\ //g')";
		STARTING_FILE_TYPE="$MY_ARG" ;
	elif [[  "$(grep -h '\-f' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]] && [[  "$(grep -h '\-\-filetype' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]]; then
		MY_ARG="$(grep -h '\-\-filetype' ./args.txt | perl -pe 's/\-\-filetype//g' | perl -pe 's/\ //g')";
		STARTING_FILE_TYPE="$MY_ARG" ;
	fi
#
	if [[  "$(grep -h '\-i' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]] && [[  "$(grep -h '\-\-input' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]]; then
		MY_NEXUS=NULL ;
	elif [[  "$(grep -h '\-i' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]] && [[  "$(grep -h '\-\-input' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]]; then
		MY_ARG="$(grep -h '\-i' ./args.txt | perl -pe 's/\-i//g' | perl -pe 's/\ //g')";
		MY_NEXUS="$MY_ARG" ;
	elif [[  "$(grep -h '\-i' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]] && [[  "$(grep -h '\-\-input' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]]; then
		MY_ARG="$(grep -h '\-\-input' ./args.txt | perl -pe 's/\-\-input//g' | perl -pe 's/\ //g')";
		MY_NEXUS="$MY_ARG" ;
	fi
#
	if [[  "$(grep -h '\-e' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]] && [[  "$(grep -h '\-\-exec' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]]; then
		MY_RAXML_EXECUTABLE=raxml ;
	elif [[  "$(grep -h '\-e' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]] && [[  "$(grep -h '\-\-exec' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]]; then
		MY_ARG="$(grep -h '\-e' ./args.txt | perl -pe 's/\-e//g' | perl -pe 's/\ //g')";
		MY_RAXML_EXECUTABLE="$MY_ARG" ;
	elif [[  "$(grep -h '\-e' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]] && [[  "$(grep -h '\-\-exec' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]]; then
		MY_ARG="$(grep -h '\-\-exec' ./args.txt | perl -pe 's/\-\-exec//g' | perl -pe 's/\ //g')";
		MY_RAXML_EXECUTABLE="$MY_ARG" ;
	fi
#
# 	if [[  "$(grep -h '\-p' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]] && [[  "$(grep -h '\-\-part' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]]; then
# 		MY_PARTITIONS_FILE=partitions.txt ;
# 	elif [[  "$(grep -h '\-p' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]] && [[  "$(grep -h '\-\-part' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]]; then
# 		MY_ARG="$(grep -h '\-p' ./args.txt | perl -pe 's/\-p//g' | perl -pe 's/\ //g')";
# 		MY_PARTITIONS_FILE="$MY_ARG" ;
# 	elif [[  "$(grep -h '\-p' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]] && [[  "$(grep -h '\-\-part' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]]; then
# 		MY_ARG="$(grep -h '\-\-part' ./args.txt | perl -pe 's/\-\-part//g' | perl -pe 's/\ //g')";
# 		MY_PARTITIONS_FILE="$MY_ARG" ;
# 	fi
# #
	if [[  "$(grep -h '\-b' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]] && [[  "$(grep -h '\-\-boot' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]]; then
		MY_NUM_BOOTREPS=100 ;
	elif [[  "$(grep -h '\-b' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]] && [[  "$(grep -h '\-\-boot' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]]; then
		MY_ARG="$(grep -h '\-b' ./args.txt | perl -pe 's/\-b//g' | perl -pe 's/\ //g')";
		MY_NUM_BOOTREPS="$MY_ARG" ;
	elif [[  "$(grep -h '\-b' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]] && [[  "$(grep -h '\-\-boot' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]]; then
		MY_ARG="$(grep -h '\-\-boot' ./args.txt | perl -pe 's/\-\-boot//g' | perl -pe 's/\ //g')";
		MY_NUM_BOOTREPS="$MY_ARG" ;
	fi
#
	if [[  "$(grep -h '\-r' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]] && [[  "$(grep -h '\-\-raxmlmodel' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]]; then
		MY_RAXML_MODEL=GTRGAMMA ;
	elif [[  "$(grep -h '\-r' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]] && [[  "$(grep -h '\-\-raxmlmodel' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]]; then
		MY_ARG="$(grep -h '\-r' ./args.txt | perl -pe 's/\-r//g' | perl -pe 's/\ //g')";
		MY_RAXML_MODEL="$MY_ARG" ;
	elif [[  "$(grep -h '\-r' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]] && [[  "$(grep -h '\-\-raxmlmodel' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]]; then
		MY_ARG="$(grep -h '\-\-raxmlmodel' ./args.txt | perl -pe 's/\-\-raxmlmodel//g' | perl -pe 's/\ //g')";
		MY_RAXML_MODEL="$MY_ARG" ;
	fi
#
	if [[  "$(grep -h '\-s' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]] && [[  "$(grep -h '\-\-simplemodel' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]]; then
		MY_SIMPLE_MODEL=NULL ;
	elif [[  "$(grep -h '\-s' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]] && [[  "$(grep -h '\-\-simplemodel' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]]; then
		MY_ARG="$(grep -h '\-s' ./args.txt | perl -pe 's/\-s//g' | perl -pe 's/\ //g')";
		MY_SIMPLE_MODEL="$MY_ARG" ;
	elif [[  "$(grep -h '\-s' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]] && [[  "$(grep -h '\-\-simplemodel' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]]; then
		MY_ARG="$(grep -h '\-\-simplemodel' ./args.txt | perl -pe 's/\-\-simplemodel//g' | perl -pe 's/\ //g')";
		MY_SIMPLE_MODEL="$MY_ARG" ;
	fi
#
	if [[  "$(grep -h '\-g' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]] && [[  "$(grep -h '\-\-gapthresh' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]]; then
		MY_GAP_THRESHOLD=0.001 ;
	elif [[  "$(grep -h '\-g' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]] && [[  "$(grep -h '\-\-gapthresh' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]]; then
		MY_ARG="$(grep -h '\-g' ./args.txt | perl -pe 's/\-g//g' | perl -pe 's/\ //g')";
		MY_GAP_THRESHOLD="$MY_ARG" ;
	elif [[  "$(grep -h '\-g' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]] && [[  "$(grep -h '\-\-gapthresh' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]]; then
		MY_ARG="$(grep -h '\-\-gapthresh' ./args.txt | perl -pe 's/\-\-gapthresh//g' | perl -pe 's/\ //g')";
		MY_GAP_THRESHOLD="$MY_ARG" ;
	fi
#
	if [[  "$(grep -h '\-m' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]] && [[  "$(grep -h '\-\-missing' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]]; then
		MY_INDIV_MISSING_DATA=1 ;
	elif [[  "$(grep -h '\-m' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]] && [[  "$(grep -h '\-\-missing' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]]; then
		MY_ARG="$(grep -h '\-m' ./args.txt | perl -pe 's/\-m//g' | perl -pe 's/\ //g')";
		MY_INDIV_MISSING_DATA="$MY_ARG" ;
	elif [[  "$(grep -h '\-m' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]] && [[  "$(grep -h '\-\-missing' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]]; then
		MY_ARG="$(grep -h '\-\-missing' ./args.txt | perl -pe 's/\-\-missing//g' | perl -pe 's/\ //g')";
		MY_INDIV_MISSING_DATA="$MY_ARG" ;
		if [[ -z "$MY_INDIV_MISSING_DATA" ]] && [[ "$MY_INDIV_MISSING_DATA" != "1" ]] && [[ "$MY_INDIV_MISSING_DATA" != "0" ]]; then MY_INDIV_MISSING_DATA=0 ; fi
	fi
#
	if [[  "$(grep -h '\-o' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]] && [[  "$(grep -h '\-\-outgroup' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]]; then
		MY_OUTGROUP=NULL ;
	elif [[  "$(grep -h '\-o' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]] && [[  "$(grep -h '\-\-outgroup' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]]; then
		MY_ARG="$(grep -h '\-o' ./args.txt | perl -pe 's/\-o//g' | perl -pe 's/\ //g')";
		MY_OUTGROUP="$MY_ARG" ;
	elif [[  "$(grep -h '\-o' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]] && [[  "$(grep -h '\-\-outgroup' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]]; then
		MY_ARG="$(grep -h '\-\-outgroup' ./args.txt | perl -pe 's/\-\-outgroup//g' | perl -pe 's/\ //g')";
		MY_OUTGROUP="$MY_ARG" ;
	fi
#
	if [[  "$(grep -h '\-n' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]] && [[  "$(grep -h '\-\-name' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]]; then
		MY_OUTPUT_NAME=raxml_out ;
	elif [[  "$(grep -h '\-n' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]] && [[  "$(grep -h '\-\-name' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]]; then
		MY_ARG="$(grep -h '\-n' ./args.txt | perl -pe 's/\-n//g' | perl -pe 's/\ //g')";
		MY_OUTPUT_NAME="$MY_ARG" ;
	elif [[  "$(grep -h '\-n' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]] && [[  "$(grep -h '\-\-name' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]]; then
		MY_ARG="$(grep -h '\-\-name' ./args.txt | perl -pe 's/\-\-name//g' | perl -pe 's/\ //g')";
		MY_OUTPUT_NAME="$MY_ARG" ;
	fi
#
	if [[  "$(grep -h '\-r' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]] && [[  "$(grep -h '\-\-resume' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]]; then
		MY_RESUME_SWITCH=0 ;
	elif [[  "$(grep -h '\-r' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]] && [[  "$(grep -h '\-\-resume' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]]; then
		MY_ARG="$(grep -h '\-r' ./args.txt | perl -pe 's/\-r//g' | perl -pe 's/\ //g')";
		MY_RESUME_SWITCH="$MY_ARG" ;
	elif [[  "$(grep -h '\-r' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]] && [[  "$(grep -h '\-\-resume' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]]; then
		MY_ARG="$(grep -h '\-\-resume' ./args.txt | perl -pe 's/\-\-resume//g' | perl -pe 's/\ //g')";
		MY_RESUME_SWITCH="$MY_ARG" ;
		if [[ -z "$MY_RESUME_SWITCH" ]] && [[ "$MY_RESUME_SWITCH" != "0" ]] && [[ "$MY_RESUME_SWITCH" != "1" ]]; then MY_RESUME_SWITCH=1 ; fi
	fi
#
	if [[  "$(grep -h '\-d' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]] && [[  "$(grep -h '\-\-debug' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]]; then
		MY_DEBUG_MODE_SWITCH=0 ;
	elif [[  "$(grep -h '\-d' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]] && [[  "$(grep -h '\-\-debug' ./args.txt | wc -l | perl -pe 's/\ //g')" = "0" ]]; then
		MY_ARG="$(grep -h '\-d' ./args.txt | perl -pe 's/\-d//g' | perl -pe 's/\ //g')";
		MY_DEBUG_MODE_SWITCH="$MY_ARG" ;
	elif [[  "$(grep -h '\-d' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]] && [[  "$(grep -h '\-\-debug' ./args.txt | wc -l | perl -pe 's/\ //g')" != "0" ]]; then
		MY_ARG="$(grep -h '\-\-debug' ./args.txt | perl -pe 's/\-\-debug//g' | perl -pe 's/\ //g')";
		MY_DEBUG_MODE_SWITCH="$MY_ARG" ;
		if [[ -z "$MY_DEBUG_MODE_SWITCH" ]] && [[ "$MY_DEBUG_MODE_SWITCH" != "0" ]] && [[ "$MY_DEBUG_MODE_SWITCH" != "1" ]]; then MY_DEBUG_MODE_SWITCH=1 ; fi
	fi
#


# ############# ############# #############
# ##       TIME TO RUN THE SCRIPT        ##
# ##                                     ##
# ## You shouldn't need to edit anything ##
# ## beneath this line                   ##
# ##                                     ##
# ############# ############# #############

# Trap bad exits with your cleanup function
# trap trapCleanup EXIT INT TERM

# Set IFS to preferred implementation
IFS=$'\n\t'

# Exit on error. Append '||true' when you run the script if you expect an error.
set -o errexit

# Run in debug mode, if set
# if ${debug}; then set -x ; fi

# Exit on empty variable
# if ${strict}; then set -o nounset ; fi

# Bash will remember & return the highest exitcode in a chain of pipes.
# This way you can catch the error in case mysqldump fails in `mysqldump |gzip`, for example.
set -o pipefail

# Invoke the checkDependenices function to test for Bash packages.  Uncomment if needed.
# checkDependencies

# Run the script
MAGNET

exit 0

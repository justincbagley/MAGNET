#!/bin/sh

##########################################################################################
#                          RAxMLRunChecker v1.3.2, December 2020                         #
# Author: Justin C. Bagley                                                               #
# Date: Created by Justin Bagley on/before November 29, 2018.                            #
# Last update: December 23, 2020                                                         #
# Copyright (c) 2018-2020 Justin C. Bagley. All rights reserved.                         #
# Please report bugs to <jbagley@jsu.edu>.                                               #
#                                                                                        #
# Description:                                                                           #
# SHELL SCRIPT THAT COUNTS NUMBER OF LOCI/PARTITIONS WITH COMPLETED RAxML RUNS DURING    #
# OR AFTER A RUN OF THE MAGNET PIPELINE, AND COLLATES RUN INFORMATION                    #
#                                                                                        #
##########################################################################################
# USAGE                                                                                  #
# $ ./RAxMLRunChecker.sh workingDir                                                      #
# $ ./shell/RAxMLRunChecker.sh workingDir                                                #
#                                                                                        #
# Examples                                                                               #
# e.g. run in current working directory (cwd), where MAGNET pipeline has been run, or is #
# currently running, by entering the following from the command line from within cwd:    #
# $ ./RAxMLRunChecker.sh .                                                               #
##########################################################################################

######################################## START ###########################################

echo "INFO      | $(date) |----------------------------------------------------------------"
echo "INFO      | $(date) | RAxMLRunChecker, v1.3.2 December 2020  (part of PIrANHA v0.4a4)"
echo "INFO      | $(date) | Copyright (c) 2018-2020 Justin C. Bagley. All rights reserved. "
echo "INFO      | $(date) |----------------------------------------------------------------"
echo "INFO      | $(date) | Starting RAxMLRunChecker pipeline... "

############ Check for mandatory positional parameters
if [ $# -lt 1 ]; then
  echo "WARNING${EP}  | $(date) |          Missing argument for working directory path. Quitting... "
  exit 1
fi
USER_SPEC_PATH="$1"


echo "INFO      | $(date) | Step #1: Set up workspace. "
############ Set workingDir
if [[ "$USER_SPEC_PATH" = "$(printf '%q\n' "$(pwd)")" ]] || [[ "$USER_SPEC_PATH" = "." ]]; then
	#MY_CWD=`pwd -P`
	MY_CWD="$(printf '%q\n' "$(pwd)" | sed 's/\\//g')"
	echo "INFO      | $(date) |          Setting working directory to:  "
	echo "$MY_CWD "
elif [[ "$USER_SPEC_PATH" != "$(printf '%q\n' "$(pwd)")" ]]; then
	if [[ "$USER_SPEC_PATH" = ".." ]] || [[ "$USER_SPEC_PATH" = "../" ]] || [[ "$USER_SPEC_PATH" = "..;" ]] || [[ "$USER_SPEC_PATH" = "../;" ]]; then
		cd ..;
		MY_CWD="$(printf '%q\n' "$(pwd)" | sed 's/\\//g')"
	else
		MY_CWD=$USER_SPEC_PATH
		echo "INFO      | $(date) |          Setting working directory to user-specified dir:  "	
		echo "$MY_CWD "
		cd "$MY_CWD"
	fi
else
	echo "WARNING${EP}  | $(date) |          Null working directory path. Quitting... "
	exit 1
fi

	# Import variables and functions
	# ------------------------------------------------------
	EP=$(printf '!')
	CR=$(printf '\r')
	TAB=$(printf '\t'); 
	calc () {
	bc -l <<< "$@"
}

echo "INFO      | $(date) | Step #2: Check RAxML runs in subfolders in current directory (assumed to be a MAGNET run folder). "

	MY_N_LOCI_FOLD="$(ls -d ./locus*/ | wc -l | sed 's/^[\ ]*//g')";
	MY_N_COMPLETED="$(ls ./locus*/RAxML_info.raxml_out | wc -l | sed 's/^[\ ]*//g')";
	MY_N_REMAINING="$(calc "$MY_N_LOCI_FOLD" - "$MY_N_COMPLETED")";

	echo "INFO      | $(date) |          Total no. RAxML runs: $TAB$TAB$MY_N_LOCI_FOLD "
	echo "INFO      | $(date) |          No. completed RAxML runs: $TAB$MY_N_COMPLETED "
	echo "INFO      | $(date) |          No. remaining RAxML runs:    $TAB$MY_N_REMAINING "

	if [[ -s ./completed_run_info.tmp ]]; then
		rm ./completed_run_info.tmp ;
	fi
	if [[ -s ./completed_run_info.txt ]]; then
		rm ./completed_run_info.txt ;
	fi
	if [[ -s ./remaining_run_info.tmp ]]; then
		rm ./remaining_run_info.tmp ;
	fi
	if [[ -s ./remaining_run_info.txt ]]; then
		rm ./remaining_run_info.txt ;
	fi

	echo "INFO      | $(date) |          Saving RAxML run info to file... "

	count=1
	echo "INFO      | $(date) |          ...  $count / $MY_N_LOCI_FOLD ..."
(
	for i in ./locus*/; do 
		MY_LOCUS="$(echo "$i" | sed 's/\.\///g; s/\///g; s/\ //g')"; 
		MY_COUNT_HUND_CHECK="$(calc $count / 100 | sed 's/^[0-9]*\.//g; s/^[0]\{1\}//g')"
		if [[ "$MY_COUNT_HUND_CHECK" -eq "0" ]]; then
			echo "INFO      | $(date) |          ...  $count / $MY_N_LOCI_FOLD ..."
		fi
		if [[ "$count" -eq "$MY_N_LOCI_FOLD" ]]; then
			echo "INFO      | $(date) |          ...  $MY_N_LOCI_FOLD / $MY_N_LOCI_FOLD ..."
		fi
		cd "$i"; 
			if [[ -s RAxML_bipartitions.raxml_out ]]; then 

				MY_ALIGN_PATT="$(grep -h '^Alignment\ Patterns\:\ ' ./RAxML_info.raxml_out | sed 's/^.*\:\ //g')"
				MY_SUBST_MODEL="$(grep -h '^Substitution\ Matrix\:\ ' ./RAxML_info.raxml_out | sed 's/^.*\:\ //g')"
				MY_OPTIM_LIKE="$(grep -h 'Final\ ML\ Optimization\ Likelihood\:\ ' ./RAxML_info.raxml_out | sed 's/^.*\:\ //g')"
				MY_ML_RUN_TIME="$(grep -h 'Overall\ execution\ time\ ' ./RAxML_info.raxml_out | sed 's/^Overall\ execution\ time\ [A-Za-z\ ]*\:\ //g; s/or\ .*//g')"
				
				echo "$count$TAB$MY_LOCUS$TAB$MY_ALIGN_PATT$TAB$MY_SUBST_MODEL$TAB$MY_OPTIM_LIKE$TAB$MY_ML_RUN_TIME$TAB complete" >> ../completed_run_info.tmp; 
			fi
			if [[ ! -s RAxML_bipartitions.raxml_out ]]; then 
				
				echo "$count$TAB$MY_LOCUS$TAB$MY_ALIGN_PATT$TAB$MY_SUBST_MODEL$TAB incomplete" >> ../remaining_run_info.tmp ; 
				
			fi
		cd ..; 
		echo "$((count++))" > count.tmp ;
	done
)


	echo "No$TAB Locus$TAB No. Patterns$TAB Subst. Model$TAB Likelihood$TAB ML Run Time$TAB Status" > ./header.tmp ;
	echo "No$TAB Locus$TAB No. Patterns$TAB Subst. Model$TAB Status" > ./rem_header.tmp ;
	cat ./header.tmp ./completed_run_info.tmp > ./completed_run_info.txt ;
	cat ./rem_header.tmp ./remaining_run_info.tmp > ./remaining_run_info.txt ;

		# Check machine type and delete spaces from run info file using slightly different 
		# sed according to machine type:
		unameOut="$(uname -s)"
		case "${unameOut}" in
		    Linux*)     machine=Linux;;
		    Darwin*)    machine=Mac;;
		    CYGWIN*)    machine=Cygwin;;
		    MINGW*)     machine=MinGw;;
		    *)          export machine="UNKNOWN:${unameOut}"
		esac
		# echo "INFO      | $(date) |          System: ${machine}"

	echo "INFO      | $(date) |          Editing final RAxML run information files... "
		if [[ "${machine}" = "Mac" ]]; then
			sed -i '' 's/\ //g' ./completed_run_info.txt ;
			sed -i '' 's/\ //g' ./remaining_run_info.txt ;
		fi
		if [[ "${machine}" = "Linux" ]]; then
			sed -i 's/\ //g' ./completed_run_info.txt ;
			sed -i 's/\ //g' ./remaining_run_info.txt ;
		fi


############ III. CLEAN UP WORKSPACE BY REMOVING TEMPORARY FILES.
echo "INFO      | $(date) | Step #3: Clean up workspace. "
echo "INFO      | $(date) |          Cleaning up workspace by removing temporary files generated during run... "

	if [[ "$(ls -1 ./*.tmp 2>/dev/null | wc -l | sed 's/\ //g')" != "0"  ]]; then 
		rm ./*.tmp ; 
	fi

echo "-------------------------------------------------------------------------------------"
echo "output file(s): ./completed_run_info.txt "
echo "                ./remaining_run_info.txt "
echo ""

#
#
#
######################################### END ############################################

exit 0


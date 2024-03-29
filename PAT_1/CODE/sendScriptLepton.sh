#!/bin/bash

# Submission script for GridEngine (GE). Each job will 
# be executed via the jobScript.sh
# This jobScript supports up to 7 parameters. Edit 
# the user specific part of the script according to 
# your program.
#
# Input to the script is a filelist with 1 file per line.
# For each file a job is started. With the parameter 
# nFilesPerJob a comma separated filelist will be 
# generated and handed to the job script. This feature
# is usefull when running many small jobs. Each
# job has its own logfile. All needed directories for the
# logfiles will be created if non existing.
#
# IMPORTANT: the hera/prometheus cluster jobs will only
# see the /hera file system. All needed scripts, programs
# and parameters have to be located on /hera or the job
# will crash. This script syncs your working dir to the submission
# dir on /hera . Make sure your scripts use the submission dir!
# Software should be taken from /cvmfs/hades.gsi.de/install/
#
# job log files will be named like inputfiles. If nFilesPerJob > 1
# the log files will contain the partnumber.
#
######################################################################
#   CONFIGURATION

user=$(whoami)
currentDir=$(pwd)
submmissionbase=/hera/hades/user/${user}
submissiondir=${submmissionbase}/PAT64
nFilesPerJob=10                                 # number of files to be analyzed by 1 job (default==1)
jobscript=${submissiondir}/jobScript.sh      # exec script (full path, call without dot, set it executable!)
outputdir=/hera/hades/user/${user}/apr06/pat_lepton/NEW_delta_L075_PENA # outputdir for files AND logFiles
pathoutputlog=${outputdir}/out                   # protocol from batch farm for each file
par1=/cvmfs/hades.gsi.de/install/5.34.01/old/hydra1-17012013/defall.sh # optional par1 : environment script
par2=${submissiondir}/anapp             # optional par2 : executable
par3=""                                          # optional par3 : input file list
par4=""                                # optional par4 : outputfile (part number will be appended (_num.root))
par5=""                                                         # optional par5 : number of events
par6=""                                        # optional par6
par7="no"                                        # optional par7
resources="-l h_rt=10:0:0,h_vmem=2G"             # runtime < 10h, mem < 2GB

filelist=${currentDir}/delta_Dalitz_PENA.list  # file list in local dir! not in submissiondir!!!
######################################################################




nFiles=$( cat $filelist | wc -l)

#---------------------------------------------------------------------
# create needed dirs
if [ ! -d $submmissionbase ]
then
    echo "===> CREATE SUBMISSIONBASEDIR : $submmissionbase"
    mkdir -p $submmissionbase
else
    echo "===> USE SUBMISSIONBASEDIR : $submmissionbase"
fi

#---------------------------------------------------------------------
# output dirs

if [ ! -d $pathoutputlog ]
then
   echo "===> CREATE LOGDIR : $pathoutputlog"
   mkdir -p $pathoutputlog
else
   echo "===> USE LOGDIR : $pathoutputlog"
fi

if [ ! -d $outputdir ]
then
   echo "===> CREATE OUTPUTDIR : $outputdir"
   mkdir -p $outputdir
   mkdir -p $outputdir/crash
else
   echo "===> USE OUTPUTDIR : $outputdir"
fi
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# sync the local modified stuff 
# to the submission dir
echo "===> SYNC CURENTDIR TO SUBMISSIONDIR : rsync  -vHaz $currentDir ${submmissionbase}"
rsync  -vHaz $currentDir ${submmissionbase}/

syncStat=$?

if [ ! $syncStat -eq 0 ]
then
     echo "===> ERROR : SYNCHRONIZATION ENCOUNTERED PROBLEMS"
fi

echo "-------------------------------------------------"
#---------------------------------------------------------------------



ctF=0          # counter for file number
ctJ=0          # counter for job number
partNumber=0   # counter for part number


#---------------------------------------------------------------------
# read the files list into an job array
declare -a jobarray
ct1=0
for file in $(cat $filelist)
do
   jobarray[$ct1]=$file
   ((ct1+=1))
done
#---------------------------------------------------------------------


#---------------------------------------------------------------------
# loop over the job array and submit parts with
# nFilesPerJob to GE

while ((ctF<$nFiles))
do
     #---------------------------------------------------------------------
     # build comma separated file list
     # per job
     if [ $nFilesPerJob -gt 1 ]
     then
        infileList=${jobarray[${ctF}]}
        ((ctF+=1))
        for (( ctList=1;ctList<$nFilesPerJob; ctList++ ))
        do   	
            if [ $ctF -lt ${nFiles} ]
            then
               infileList="${infileList},${jobarray[${ctF}]}"
               ((ctF+=1))
            fi
        done
     else 
        infileList=${jobarray[${ctF}]}
        ((ctF+=1))
     fi
     #---------------------------------------------------------------------
     
     ((partNumber+=1))

     logfile="${pathoutputlog}/${partNumber}.log"

     if [ $nFilesPerJob -eq 1 ]
     then
        file=$(basename ${infileList})
        logfile="${pathoutputlog}/${file}.log"
     fi
     
     
     if [ -f ${logfile} ]
     then
        rm -f ${logfile}
     fi
     
     
     echo "-----------------------------------------------------------------------------"
     echo "add part ${partNumber}  last file ${ctF} of $nFiles ====> add new job ${infileList}"

     ######################################################################
     #  SEND NEW JOB (USER SPECIFIC)
     
     par3=${infileList}


     command="-j y -wd ${submissiondir} ${resources} -o ${logfile} \
     ${jobscript}  ${par1}   ${par2} ${par3}  ${par4}  ${par5}  ${par6} ${par7}"
     
     echo qsub ${command}

     if [ ! $syncStat -eq 0 ]
     then
        echo "===> ERROR : SYNCHRONIZATION ENCOUNTERED PROBLEMS"
     else
         # echo ""
         qsub ${command}
         sleep 0
     fi
     ######################################################################
     
done
#---------------------------------------------------------------------

echo "${nFiles} jobs submitted"


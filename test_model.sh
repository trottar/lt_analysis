#! /bin/bash

# Runs script in the ltsep python package that grabs current path enviroment
if [[ ${HOSTNAME} = *"cdaq"* ]]; then
    PATHFILE_INFO=`python3 /home/cdaq/pionLT-2021/hallc_replay_lt/UTIL_PION/bin/python/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
elif [[ ${HOSTNAME} = *"farm"* ]]; then
    PATHFILE_INFO=`python3 $replay_lt_env/lib/python3.9/site-packages/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
fi

# Split the string we get to individual variables, easier for printing and use later
VOLATILEPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f1` # Cut the string on , delimitter, select field (f) 1, set variable to output of command
ANALYSISPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f2`
HCANAPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f3`
REPLAYPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f4`
UTILPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f5`
PACKAGEPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f6`
OUTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f7`
ROOTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f8`
REPORTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f9`
CUTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f10`
PARAMPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f11`
SCRIPTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f12`
ANATYPE=`echo ${PATHFILE_INFO} | cut -d ','  -f13`
USER=`echo ${PATHFILE_INFO} | cut -d ','  -f14`
HOST=`echo ${PATHFILE_INFO} | cut -d ','  -f15`
SIMCPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f16`
LTANAPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f17`

# Define global variables for lt_analysis scripts
#Q2="4p4"
#W="2p74"
#TMIN=0.400
#TMAX=0.750
#ParticleType="kaon"
#
#Q2="1p6"
#W="2p22"
#TMIN=0.001
#TMAX=0.300
#
Q2="2p4"
W="2p22"
TMIN=0.100
TMAX=0.600
ParticleType="pion"

POL="+1" # Positive polarity

DEBUG="False" # Flag for no plot splash
#DEBUG="True" # Flag for plot splash

KIN="Q${Q2}W${W}"
OutFilename="Testing_${KIN}"
# - `date` is a command that prints or sets the system date and time.
# - `+%H` extracts the hour in 24-hour format.
# - `+%M` extracts the minute.
# - `+%S` extracts the second.
# - `+%Y` extracts the year.
# - `+%B` extracts the full month name.
# - `%d` extracts the day of the month.
formatted_date=$(date +%Y%B%d_H%HM%MS%S) 

echo
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "Running fitting algorithm for ${ParticleType} at Kinematics: ${KIN} with POL: ${POL}"
echo "Output will be saved to: ${OutFilename}"
echo "Current time and date: ${formatted_date}"
echo "DEBUG mode is: ${DEBUG}"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo

cd "${LTANAPATH}/testing"
python3 main.py ${KIN} ${W} ${Q2} ${TMIN} ${TMAX} ${ParticleType} ${POL} ${OutFilename} ${formatted_date} ${DEBUG}
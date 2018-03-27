#!/bin/bash
#PBS -e ERRFILE
#PBS -o LOGFILE

START_TIME=$(date) 
START_TIME_S=$(date +"%s") 
echo " "
echo "========================="
echo "       Beginning job     "
echo "                         "
echo "Start time : $START_TIME "
echo "Hostname   : $HOSTNAME   "
lsb_release -a
echo "========================="
echo " "

######################
# Set up environment #
######################
echo " "
echo "> Setting up environment"
echo "(" `date` ")"
echo " "

#export LD_PRELOAD=/usr/lib64/libglobus_gssapi_gsi.so.4
export V0_CMS_SW_DIR=/cvmfs/cms.cern.ch/
source $V0_CMS_SW_DIR/cmsset_default.sh

# Grid/DPM proxy
# --------------

export X509_USER_PROXY=/home-pbs/$USER/proxy/x509up_u9029

# LD_LIBRARY_PATH
# ---------------

cd /home-pbs/mjansova/CMSSW_10_0_2/src/CalibTracker/SiStripCommon/test
eval `scramv1 runtime -sh`
#cd -

$LD_LIBRARY_PATH



############################################
# Move to working area and launch analysis #
############################################

echo " "
echo "> Moving to working area"
echo "(" `date` ")"
echo " "

which root
echo $ROOTSYS
root -l -b --version -q


cd /home-pbs/mjansova/CMSSW_10_0_2/src/CalibTracker/SiStripCommon/test

echo " "
echo "> Launch the analysis"
echo "(" `date` ")"
echo " "

cmsRun run_simulation_step2_CrossTalk_ePerAdc_2017Coll.py /opt/sbg/scratch1/cms/mjansova/store/tmp/MCtuning/FNAME EPERADC X1 X2 X3
cmsRun run_simulation_step3_CrossTalk_ePerAdc_2017Coll.py /opt/sbg/scratch1/cms/mjansova/store/tmp/MCtuning/FNAME EPERADC X1 X2 X3
cmsRun trackCluster_testCollMCValidation.py /opt/sbg/scratch1/cms/mjansova/store/tmp/MCtuning/FNAMEstep3.root none 
rm /opt/sbg/scratch1/cms/mjansova/store/tmp/MCtuning/FNAME.root
rm /home-pbs/mjansova/CMSSW_10_0_2/src/CalibTracker/SiStripCommon/test/JOBN

######################
# End of job message #
######################

END_TIME=$(date) 
END_TIME_S=$(date +"%s") 
DURATION=$(($END_TIME_S - $START_TIME_S))
echo " "
echo "======================="
echo "       Job ending      "
echo "                       "
echo "End time : $END_TIME   "
echo "Duration : $(($DURATION / 60)) min $(($DURATION % 60)) s "
echo "======================="
echo " "

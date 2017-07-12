#!/bin/env zsh

for i in `seq 1 20`
do
    par2=$((0.07+($i*0.0005)))
    par3=0.02
    par1=$((1-(2*${par2})-(2*${par3})))
    echo "cmsRun run_simulation_step2_CrossTalk_2016.py CT3fix_"${par1}" "${par1}" "${par2}" "${par3} 
    echo "cmsRun run_simulation_step3_2016.py CT3fix_"${par1}
    echo "cmsRun trackCluster_test.py CT3fix_"${par1}"step3.root"
    echo "./distributionsChi2Test noPU2016Data_v2_fullInfo.root test_shallowTrackAndClusterFullInfo.rootCT3fix_"${par1}"step3.root "${par2}" "${par3}" TOB 0 5 yes CTtuningFix3" 
done

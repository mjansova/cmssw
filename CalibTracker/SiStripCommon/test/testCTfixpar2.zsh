#!/bin/env zsh

for i in `seq 1 20`
do
    par2=0.0745
    par3=$((0.01+($i*0.0015)))
    par1=$((1-(2*${par2})-(2*${par3})))
    echo "cmsRun run_simulation_step2_CrossTalk_2016.py CT2fix_"${par1}" "${par1}" "${par2}" "${par3} 
    echo "cmsRun run_simulation_step3_2016.py CT2fix_"${par1}
    echo "cmsRun trackCluster_test.py CT2fix_"${par1}"step3.root"
    echo "./distributionsChi2Test noPU2016Data_v2_fullInfo.root test_shallowTrackAndClusterFullInfo.rootCT2fix_"${par1}"step3.root "${par2}" "${par3}" TOB 0 5 yes CTtuningFix2" 
done

cmsRun ../python/UndergroundCosmicMu_cfi_py_GEN_SIM_2017Cos.py
cmsRun run_simulation_step2_CrossTalk_2017CRAFT.py CTcraft_0.805 0.805 0.0775 0.02
cmsRun run_simulation_step3_2017Cos.py CTcraft_0.805
cmsRun trackCluster_testCosmicsCRAFT.py  CTcraft_0.805step3.root
cmsRun run_simulation_Default_2017CRAFT.py CTcraft_def
cmsRun run_simulation_step3_2017Cos.py CTcraft_def
cmsRun trackCluster_testCosmicsCRAFT.py  CTcraft_defstep3.root

#cmsRun trackCluster_testCosmicsCRUZET.py  CTcruzet_0.805step3.root


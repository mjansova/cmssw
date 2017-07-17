cmsRun ../python/UndergroundCosmicMu_cfi_py_GEN_SIM_2017CosCRUZET.py
cmsRun run_simulation_step2_CrossTalk_2017CRUZET.py CTcruzet_0.805 0.805 0.0775 0.02
cmsRun run_simulation_step3_2017CosCRUZET.py CTcruzet_0.805
cmsRun run_simulation_Default_2017CRUZET.py CTcruzet_def
cmsRun run_simulation_step3_2017CosCRUZET.py CTcruzet_def

#cmsRun trackCluster_testCosmicsCRUZET.py  CTcruzet_0.805step3.root


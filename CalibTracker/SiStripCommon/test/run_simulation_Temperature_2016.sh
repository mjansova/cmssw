cmsRun ../python/MinBias_13TeV_pythia8_TuneCUETP8M1_cfi_GEN_SIM_2016.py #runs MC generation and geant4 simulation

cmsRun run_simulation_step2_Temperature_2016.py Temperature_min20_218.4 218.4 #runs digitization with temperature parameter changed to 218.4 
cmsRun run_simulation_step3_2016.py Temperature_min20_218.4 #runs reconstruction of the new simulation
cmsRun trackCluster_test.py  Temperature_min20_218.4step3.root #gets the on track clusters properties

cmsRun run_simulation_step2_Temperature_2016.py Temperature_min10_245.7 245.7
cmsRun run_simulation_step3_2016.py Temperature_min10_245.7
cmsRun trackCluster_test.py  Temperature_min10_245.7step3.root

cmsRun run_simulation_step2_Temperature_2016.py Temperature_def_273 273.0
cmsRun run_simulation_step3_2016.py Temperature_def_273
cmsRun trackCluster_test.py  Temperature_def_273step3.root

cmsRun run_simulation_step2_Temperature_2016.py Temperature_pl10_300.3 300.3
cmsRun run_simulation_step3_2016.py Temperature_pl10_300.3
cmsRun trackCluster_test.py  Temperature_pl10_300.3step3.root

cmsRun run_simulation_step2_Temperature_2016.py Temperature_pl20_327.6 327.6
cmsRun run_simulation_step3_2016.py Temperature_pl20_327.6
cmsRun trackCluster_test.py  Temperature_pl20_327.6step3.root

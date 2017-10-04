#cmsRun ../python/MinBias_13TeV_pythia8_TuneCUETP8M1_cfi_GEN_SIM_2016.py #runs MC generation and geant4 simulation

#cmsRun run_simulation_step2_electronPerAdc_2016.py electronPerAdc_min20_198.0 198.0 #runs digitization with temperature parameter changed to 218.4 
#cmsRun run_simulation_step3_2016.py electronPerAdc_min20_198.0 #runs reconstruction of the new simulation
cmsRun trackCluster_test.py  electronPerAdc_min20_198.0step3.root no #gets the on track clusters propertiesip

#cmsRun run_simulation_step2_electronPerAdc_2016.py electronPerAdc_min10_222.0 222.0 #runs digitization with temperature parameter changed to 218.4 
#cmsRun run_simulation_step3_2016.py electronPerAdc_min10_222.0 #runs reconstruction of the new simulation
cmsRun trackCluster_test.py  electronPerAdc_min10_222.0step3.root no #gets the on track clusters properties

#cmsRun run_simulation_step2_electronPerAdc_2016.py electronPerAdc_pl10_272.0 272.0 #runs digitization with temperature parameter changed to 218.4 
#cmsRun run_simulation_step3_2016.py electronPerAdc_pl10_272.0 #runs reconstruction of the new simulation
cmsRun trackCluster_test.py  electronPerAdc_pl10_272.0step3.root no #gets the on track clusters properties

#cmsRun run_simulation_step2_electronPerAdc_2016.py electronPerAdc_pl20_296.0 296.0 #runs digitization with temperature parameter changed to 218.4 
#cmsRun run_simulation_step3_2016.py electronPerAdc_pl20_296.0 #runs reconstruction of the new simulation
cmsRun trackCluster_test.py  electronPerAdc_pl20_296.0step3.root no #gets the on track clusters properties

#cmsRun run_simulation_step2_electronPerAdc_2016.py electronPerAdc_def_247.0 247.0 #runs digitization with temperature parameter changed to 218.4 
#cmsRun run_simulation_step3_2016.py electronPerAdc_def_247.0 #runs reconstruction of the new simulation
cmsRun trackCluster_test.py  electronPerAdc_def_247.0step3.root no #gets the on track clusters properties

#!/usr/bin/env python

import subprocess
import os
import regtools
from regtools import RegArgs
    
def main():
    run_step1 = False
    run_step2 = True
    run_step3 = True

    #modify the parameters as you wish and then re-run
    input_ideal_ic  = "/mercury/data1/harper/EgRegsNtups/DoubleElectron_FlatPt-1To300_2017ConditionsFlatPU0to70ECALGT_105X_mc2017_realistic_IdealEcalIC_v5-v2_AODSIM_EgRegTreeV1_extraVars.root"
    input_real_ic = "/mercury/data1/harper/EgRegsNtups/DoubleElectron_FlatPt-1To300_2017ConditionsFlatPU0to70_105X_mc2017_realistic_v5-v2_AODSIM_EgRegTreeV1_4.root"
    regArgs = RegArgs()
 #   regArgs.input_training = "/mercury/data1/harper/EgRegsNtups/DoubleElectron_FlatPt-1To300_2017ConditionsFlatPU0to70ECALGT_105X_mc2017_realistic_IdealEcalIC_v5-v2_AODSIM_EgRegTreeV1_extraVars.root"
#    regArgs.input_testing = "/mercury/data1/harper/EgRegsNtups/DoubleElectron_FlatPt-1To300_2017ConditionsFlatPU0to70ECALGT_105X_mc2017_realistic_IdealEcalIC_v5-v2_AODSIM_EgRegTreeV1_extraVars.root"
    regArgs.input_training =  str(input_ideal_ic)
    regArgs.input_testing = str(input_ideal_ic)

    regArgs.cfg_dir = "configs"
    regArgs.out_dir = "results"
    regArgs.cuts_name = "stdCuts" 
    regArgs.base_name = "scRegUL_1050"
    regArgs.ntrees = 1500
 
    if run_step1: regArgs.run_eb_and_ee()
    
    regArgs.do_eb = True
    forest_eb_file = regArgs.output_name()
    regArgs.do_eb = False
    forest_ee_file = regArgs.output_name()
    
    

    regArgs.base_name = "scRegUL_1050_realIC_IdealTraining"
#    regArgs.base_name = "scRegUL_1050_idealIC_IdealTraining"
    input_for_res_training = str(regArgs.applied_name()) #save the output name before we change it
    input_for_input_for_res_training = str(input_real_ic)
  #  input_for_input_for_res_training = str(input_ideal_ic)
    
    if run_step2: subprocess.Popen(["bin/slc6_amd64_gcc700/RegressionApplierExe",input_for_input_for_res_training,input_for_res_training,"--gbrForestFileEE",forest_ee_file,"--gbrForestFileEB",forest_eb_file,"--nrThreads","4","--treeName",regArgs.tree_name,"--writeFullTree","1","--regOutTag","Ideal"]).communicate()

 #   regArgs.base_name = "scRegUL_1050_idealIC_TestResTraining"
    regArgs.base_name = "scRegUL_1050_realIC_RealTraining"
    regArgs.input_training = input_for_res_training
    regArgs.input_testing = input_for_res_training
    regArgs.target = "mc.energy/(sc.rawEnergy*regIdealMean)"
    regArgs.mean_min = 0.2
    regArgs.mean_max = 2.0
    regArgs.fix_mean = True
    if run_step3: regArgs.run_eb_and_ee()

if __name__ =='__main__':
    main()



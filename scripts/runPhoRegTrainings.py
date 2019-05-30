#!/usr/bin/env python

import subprocess
import os
import regtools
from regtools import RegArgs
    

def main():

    #step 1, run calo only regression on the ideal IC to get the mean
    #step 2, apply the mean to the real IC sample and save the result in a tree
    #step 3, retrain the resolution for the real IC on the corrected energy
    #step 4, run trk-calo regression using the real IC corrections as inputs 

    #event split: ECAL Ideal IC train = eventnr%10=0
    #             ECAL Real IC train = eventnr%10=1
    run_step1 = True
    run_step2 = True
    run_step3 = True
    
    base_pho_cuts = "(mc.energy>0 && ssFrac.sigmaIEtaIEta>0 && ssFrac.sigmaIPhiIPhi>0 && pho.et>0 && {extra_cuts})"

    input_ideal_ic  = "/mercury/data1/harper/EgRegsNtups/DoublePhoton_FlatPt-5To300_2017ConditionsFlatPU0to70ECALGT_105X_mc2017_realistic_IdealEcalIC_v5-v2_AODSIM_EgRegTreeV5Refined.root"
    input_real_ic = "/mercury/data1/harper/EgRegsNtups/DoublePhoton_FlatPt-5To300_2017ConditionsFlatPU0to70_105X_mc2017_realistic_v5-v2_AODSIM_EgRegTreeV5Refined.root"
    #step1 train the calo only regression using IDEAL intercalibration constants
    print "starting step1"
    regArgs = RegArgs()
    regArgs.input_training = str(input_ideal_ic)
    regArgs.input_testing = str(input_ideal_ic)  
    regArgs.set_phoecal_default()
    regArgs.cuts_base = base_pho_cuts.format(extra_cuts = "evt.eventnr%10==0")
    regArgs.cfg_dir = "configs"
    regArgs.out_dir = "resultsPhoV5" 
    regArgs.ntrees = 1500  
    regArgs.base_name = "regPhoEcal2017UL_IdealIC_IdealTraining"
    if run_step1: regArgs.run_eb_and_ee()
    
    #step2 now we run over the REAL intercalibration constant data and make a rew tree with this regression included
    print "starting step2"
    regArgs.do_eb = True
    forest_eb_file = regArgs.output_name()
    regArgs.do_eb = False
    forest_ee_file = regArgs.output_name()

    regArgs.base_name = "regPhoEcal2017UL_RealIC_IdealTraining"
    input_for_res_training = str(regArgs.applied_name()) #save the output name before we change it
    if run_step2: subprocess.Popen(["bin/slc6_amd64_gcc700/RegressionApplierExe",input_real_ic,input_for_res_training,"--gbrForestFileEE",forest_ee_file,"--gbrForestFileEB",forest_eb_file,"--nrThreads","4","--treeName",regArgs.tree_name,"--writeFullTree","1","--regOutTag","Ideal"]).communicate()
    
    #step3 we now run over re-train with the REAL sample for the sigma, changing the target to have the correction applied 
    print "starting step3"
    regArgs.base_name = "regPhoEcal2017UL_RealIC_RealTraining"
    regArgs.input_training = input_for_res_training
    regArgs.input_testing = input_for_res_training
    regArgs.target = "mc.energy/((sc.rawEnergy+sc.rawESEnergy)*regIdealMean)"
    regArgs.fix_mean = True
    regArgs.write_full_tree = "1"
    regArgs.reg_out_tag = "Real"
    regArgs.cuts_base = base_pho_cuts.format(extra_cuts = "evt.eventnr%10==1")
    if run_step3: regArgs.run_eb_and_ee()

    
          
    
if __name__ =='__main__':
    main()



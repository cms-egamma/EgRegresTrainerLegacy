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
    #             ECAL ECAL-Trk IC train = eventnr%10=2
    run_step1 = False
    run_step2 = False
    run_step3 = True
    run_step4 = True
    run_step4_extra = True
    
    base_ele_cuts = "(mc.energy>0 && ssFrac.sigmaIEtaIEta>0 && ssFrac.sigmaIPhiIPhi>0 && ele.et>0 && {extra_cuts})"

    #input_ideal_ic  = "/eos/cms/store/group/phys_egamma/EgRegression/DoubleElectron_FlatPt-1To300/2018-Prompt-IDEALIC/EgRegTreeV1/DoubleElectron_FlatPt-1To300_ntuples_FlatPU0to70IdealECAL_102X_upgrade2018_realistic_forECAL_v15-v1_EgRegTreeV1_1.root"
    #input_real_ic = "/eos/cms/store/group/phys_egamma/EgRegression/DoubleElectron_FlatPt-1To300/2018-Prompt/EgRegTreeV1/DoubleElectron_FlatPt-1To300_ntuples_FlatPU0to70RAW_102X_upgrade2018_realistic_v15-v1_EgRegTreeV1_1.root"
    input_ideal_ic  = "/mercury/data1/harper/EgRegsNtups/DoubleElectron_FlatPt-1To300_2017ConditionsFlatPU0to70ECALGT_105X_mc2017_realistic_IdealEcalIC_v5-v2_AODSIM_EgRegTreeV2Refined.root"
    input_real_ic = "/mercury/data1/harper/EgRegsNtups/DoubleElectron_FlatPt-1To300_2017ConditionsFlatPU0to70_105X_mc2017_realistic_v5-v2_AODSIM_EgRegTreeV2Refined.root"
    #step1 train the calo only regression using IDEAL intercalibration constants
    print "starting step1"
    regArgs = RegArgs()
    regArgs.input_training = str(input_ideal_ic)
    regArgs.input_testing = str(input_ideal_ic)  
    regArgs.set_ecal_default()
    regArgs.cuts_base = base_ele_cuts.format(extra_cuts = "evt.eventnr%10==0")
    regArgs.cfg_dir = "configs"
    regArgs.out_dir = "resultsEle" 
    regArgs.ntrees = 1500  
    regArgs.base_name = "regEle2017UL_IdealIC_IdealTraining"
    if run_step1: regArgs.run_eb_and_ee()
    
    #step2 now we run over the REAL intercalibration constant data and make a rew tree with this regression included
    print "starting step2"
    regArgs.do_eb = True
    forest_eb_file = regArgs.output_name()
    regArgs.do_eb = False
    forest_ee_file = regArgs.output_name()

    regArgs.base_name = "regEleEcal2017UL_RealIC_IdealTraining"
    input_for_res_training = str(regArgs.applied_name()) #save the output name before we change it
    if run_step2: subprocess.Popen(["bin/slc6_amd64_gcc700/RegressionApplierExe",input_real_ic,input_for_res_training,"--gbrForestFileEE",forest_ee_file,"--gbrForestFileEB",forest_eb_file,"--nrThreads","4","--treeName",regArgs.tree_name,"--writeFullTree","1","--regOutTag","Ideal"]).communicate()
    
    #step3 we now run over re-train with the REAL sample for the sigma, changing the target to have the correction applied 
    print "starting step3"
    regArgs.base_name = "regEleEcal2017UL_RealIC_RealTraining"
    regArgs.input_training = input_for_res_training
    regArgs.input_testing = input_for_res_training
    regArgs.target = "mc.energy/((sc.rawEnergy+sc.rawESEnergy)*regIdealMean)"
    regArgs.fix_mean = True
    regArgs.write_full_tree = "1"
    regArgs.reg_out_tag = "Real"
    regArgs.cuts_base = base_ele_cuts.format(extra_cuts = "evt.eventnr%10==1")
    if run_step3: regArgs.run_eb_and_ee()

    
    #step4 do the E/p low combination
    #remember we use the Ideal Mean but Real Sigma (real mean is 1 by construction)
    print "starting step4"
    input_for_comb = str(regArgs.applied_name())

    regArgs.base_name = "regEleEcalTrk2017UL_RealIC"
    regArgs.var_eb =":".join(["(sc.rawEnergy+sc.rawESEnergy)*regIdealMean","regRealSigma/regIdealMean","ele.trkPModeErr/ele.trkPMode","(sc.rawEnergy+sc.rawESEnergy)*regIdealMean/ele.trkPMode","ele.ecalDrivenSeed","ssFull.e3x3/sc.rawEnergy","ele.fbrem","ele.trkEtaMode","ele.trkPhiMode"])
    regArgs.var_ee =":".join(["(sc.rawEnergy+sc.rawESEnergy)*regIdealMean","regRealSigma/regIdealMean","ele.trkPModeErr/ele.trkPMode","(sc.rawEnergy+sc.rawESEnergy)*regIdealMean/ele.trkPMode","ele.ecalDrivenSeed","ssFull.e3x3/sc.rawEnergy","ele.fbrem","ele.trkEtaMode","ele.trkPhiMode"])
    regArgs.target = "(mc.energy * (ele.trkPModeErr*ele.trkPModeErr + (sc.rawEnergy+sc.rawESEnergy)*(sc.rawEnergy+sc.rawESEnergy)*regRealSigma*regRealSigma) / ( (sc.rawEnergy+sc.rawESEnergy)*regIdealMean*ele.trkPModeErr*ele.trkPModeErr + ele.trkPMode*(sc.rawEnergy+sc.rawESEnergy)*(sc.rawEnergy+sc.rawESEnergy)*regRealSigma*regRealSigma ))"
    regArgs.input_training = input_for_comb
    regArgs.input_testing = input_for_comb
    regArgs.write_full_tree = "0"  
    regArgs.fix_mean = False
    regArgs.reg_out_tag = "EcalTrk"
    regArgs.cuts_base = base_ele_cuts.format(extra_cuts = "evt.eventnr%10==2")
    if run_step4: 
        regArgs.run_eb_and_ee()
    if run_step4_extra:
        regArgs.base_name = "regEleEcalTrkLowPt2017UL_RealIC"
        regArgs.cuts_base = base_ele_cuts.format(extra_cuts = "evt.eventnr%10==2 && mc.pt<50")
        forest_eb,forest_ee = regArgs.forest_filenames()
        regArgs.run_eb_and_ee()
        regArgs.base_name = "regEleEcalTrkHighPt2017UL_RealIC"
        regArgs.cuts_base = base_ele_cuts.format(extra_cuts = "evt.eventnr%10==2 && mc.pt>50 && mc.pt<200")
        forest_eb_highpt,forest_ee_highpt = regArgs.forest_filenames()
        regArgs.run_eb_and_ee()
        regArgs.base_name = "regEleEcalTrkLowHighPt2017UL_RealIC"
        subprocess.Popen(["bin/slc6_amd64_gcc700/RegressionApplierExe",regArgs.input_testing,regArgs.applied_name(),"--gbrForestFileEB",forest_eb,"--gbrForestFileEE",forest_ee,"--gbrForestFileEBHighEt",forest_eb_highpt,"--gbrForestFileEEHighEt",forest_ee_highpt,"--highEtThres","50.","--nrThreads","4","--treeName",regArgs.tree_name,"--writeFullTree","0"]).communicate()
    
        
    
if __name__ =='__main__':
    main()



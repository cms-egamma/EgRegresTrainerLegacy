#!/usr/bin/env python

import subprocess
import os
try:
    from regtools import RegArgs
except ImportError:
    raise ImportError('module regtools not found, please add "python" to PYTHON27PATH or PYTHONPATH as appropriate eg\n  export PYTHON27PATH=$PYTHON27PATH:python\nand try again')
import time
import argparse
def main():  

    parser = argparse.ArgumentParser(description='runs the SC regression trainings')
    parser.add_argument('--era',required=True,help='year to produce for, 2016, 2017, 2018 are the options')
    parser.add_argument('--input_dir','-i',default='/home/hep/wrtabb/Egamma/input_trees',help='input directory with the ntuples')
    parser.add_argument('--output_dir','-o',default="results",help='output dir')
    args = parser.parse_args()

    #step 1, run calo only regression on the ideal IC to get the mean
    #step 2, apply the mean to the real IC sample and save the result in a tree
    #step 3, retrain the resolution for the real IC on the corrected energy
    run_step1 = True
    run_step2 = True
    run_step3 = True

    #setup the selection (event number cuts come later)
    cuts_name = "stdCuts" 
    base_ele_cuts = "(mc.energy>0 && ssFrac.sigmaIEtaIEta>0 && ssFrac.sigmaIPhiIPhi>0 && {extra_cuts})"
    
    #prefixes all the regressions produced
    if args.era=='2016':
        base_reg_name = "scReg2016UL"
        input_ideal_ic  = "{}/DoubleElectron_FlatPt-1To300_2016ConditionsFlatPU0to70ECALGT_105X_mcRun2_asymptotic_IdealEcalIC_newSR_v2-v2_AODSIM_EgRegTreeV5.root".format(args.input_dir)
        input_real_ic = "{}/DoubleElectron_FlatPt-1To300_2016ConditionsFlatPU0to70RAW_105X_mcRun2_asymptotic_newECALSR_v2-v2_AODSIM_EgRegTreeV5.root".format(args.input_dir)    
        ideal_eventnr_cut = "evt.eventnr%5==0"  #4million electrons
        real_eventnr_cut = "evt.eventnr%5==1" #4million electrons
    elif args.era=='2017':
        base_reg_name = "scReg2017UL"    
        input_ideal_ic  = "{}/DoubleElectron_FlatPt-1To300_2017ConditionsFlatPU0to70ECALGT_105X_mc2017_realistic_IdealEcalIC_v5-v2_AODSIM_EgRegTreeV1_extraVars.root".format(args.input_dir)
        input_real_ic = "{}/DoubleElectron_FlatPt-1To300_2017ConditionsFlatPU0to70_105X_mc2017_realistic_v5-v2_AODSIM_EgRegTreeV1_4.root".format(args.input_dir)   
        ideal_eventnr_cut = "evt.eventnr%2==0"
        real_eventnr_cut = "evt.eventnr%2==0" #events in the ntuple are different so can get away with this
    elif args.era=='2018':
        base_reg_name = "scReg2018UL"    
        input_ideal_ic  = "{}/DoubleElectron_FlatPt-1To300_2018ConditionsFlatPU0to70ECALGT_105X_upgrade2018_realistic_IdealEcalIC_v4-v1_AODSIM_EgRegTreeV5_partStatsV2.root".format(args.input_dir)
        input_real_ic = "{}/DoubleElectron_FlatPt-1To300_2018ConditionsFlatPU0to70RAW_105X_upgrade2018_realistic_v4-v1_AODSIM_EgRegTreeV5_partStatsV2.root".format(args.input_dir)    
        ideal_eventnr_cut = "evt.eventnr%5==0"  #4million electrons
        real_eventnr_cut = "evt.eventnr%5==1" #4million electrons
    else:
        raise ValueError("era {} is invalid, options are 2016/2017/2018".format(era))

    
    regArgs = RegArgs()
    regArgs.input_training =  str(input_ideal_ic)
    regArgs.input_testing = str(input_ideal_ic)
    regArgs.set_sc_default()
    regArgs.cfg_dir = "configs"
    regArgs.out_dir = args.output_dir
    regArgs.cuts_name = cuts_name
    regArgs.base_name = "{}_IdealIC_IdealTraining".format(base_reg_name)
    regArgs.cuts_base = base_ele_cuts.format(extra_cuts = ideal_eventnr_cut)
    regArgs.ntrees = 1500
 
    print """about to run the supercluster regression with: 
    name: {name}
    ideal ic input: {ideal_ic}
    real ic input: {real_ic}
    output dir: {out_dir}
steps to be run:
    step 1: ideal training for mean       = {step1}
    step 2: apply ideal training to real  = {step2}
    step 3: real training for sigma       = {step3}""".format(name=base_reg_name,ideal_ic=input_ideal_ic,real_ic=input_real_ic,out_dir=args.output_dir,step1=run_step1,step2=run_step2,step3=run_step3)
    time.sleep(20)

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    if run_step1: regArgs.run_eb_and_ee()
    
    regArgs.do_eb = True
    forest_eb_file = regArgs.output_name()
    regArgs.do_eb = False
    forest_ee_file = regArgs.output_name()

    regArgs.base_name = "{}_RealIC_IdealTraining".format(base_reg_name)
    input_for_res_training = str(regArgs.applied_name()) #save the output name before we change it
    input_for_input_for_res_training = str(input_real_ic)
    
    if run_step2: subprocess.Popen(["bin/slc6_amd64_gcc700/RegressionApplierExe",input_for_input_for_res_training,input_for_res_training,"--gbrForestFileEE",forest_ee_file,"--gbrForestFileEB",forest_eb_file,"--nrThreads","4","--treeName",regArgs.tree_name,"--writeFullTree","1","--regOutTag","Ideal"]).communicate()

    regArgs.base_name = "{}_RealIC_RealTraining".format(base_reg_name)
    regArgs.input_training = input_for_res_training
    regArgs.input_testing = input_for_res_training
    regArgs.target = "mc.energy/(sc.rawEnergy*regIdealMean)"
    regArgs.fix_mean = True
    regArgs.cuts_base = base_ele_cuts.format(extra_cuts = real_eventnr_cut)
    if run_step3: regArgs.run_eb_and_ee()

if __name__ =='__main__':
    main()



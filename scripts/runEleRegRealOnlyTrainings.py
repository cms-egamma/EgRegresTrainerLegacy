#!/usr/bin/env python

import subprocess
import os
import regtools
from regtools import RegArgs
import argparse

def main():
    parser = argparse.ArgumentParser(description='runs the SC regression trainings')
    parser.add_argument('--era',required=True,help='year to produce for, 2016, 2017, 2018 are the options')
    parser.add_argument('--input_dir','-i',default='/mercury/data1/harper/EgRegsNtups',help='input directory with the ntuples')
    parser.add_argument('--output_dir','-o',default="results",help='output dir')
    args = parser.parse_args()

    #this skips ideal/real ICs and only trains on real
    #for specialised trainings when there is no ideal
    #step 1, run calo only regression on the ideal IC to get the mean
    #step 2, run trk-calo regression using the real IC corrections as inputs 

    #event sample breakdown, 33% for ECAL, 33% for track training, 33% for testing

    run_step1 = True
    run_step2 = True
    
    base_ele_cuts = "(mc.energy>0 && ssFrac.sigmaIEtaIEta>0 && ssFrac.sigmaIPhiIPhi>0 && ele.et>0 && {extra_cuts})"

    if args.era=='2018':
        era_name = "2018Run3Proj"
        input_real_ic = "{}/DoubleElectron_FlatPt-1To100_2018ConditionsFlatPU0to70RAW_105X_upgrade2018_realistic_v4-v1_AODSIM_EgRegTreeV5Refined2018Reg_1MEles.root".format(args.input_dir)    
        ecal_eventnr_cut = "evt.eventnr%3==1"
        ep_eventnr_cut = "evt.eventnr%3==2" 
    elif args.era=='2023':
        era_name = "2023Run3Proj"
        input_real_ic = "{}/DoubleElectron_FlatPt-1To100_2023ScenarioFlatPU0to80RAW_106X_mcRun3_2023_realistic_v3_ext1-v2_AODSIM_EgRegTreeV5Refined2023Reg2018.root".format(args.input_dir)    
        ecal_eventnr_cut = "evt.eventnr%3==1"
        ep_eventnr_cut = "evt.eventnr%3==2" 
    else:
        raise ValueError("era {} is invalid, options are 2016/2017/2018".format(era))

    
    #step1 train the calo only regression using IDEAL intercalibration constants
    print "starting step1"
    regArgs = RegArgs()
    regArgs.input_training = str(input_real_ic)
    regArgs.input_testing = str(input_real_ic)  
    regArgs.set_ecal_default()
    regArgs.cuts_base = base_ele_cuts.format(extra_cuts = ecal_eventnr_cut)
    regArgs.cuts_name = "stdCuts"
    regArgs.cfg_dir = "configs"
    regArgs.out_dir = args.output_dir
    regArgs.ntrees = 1500  
    regArgs.write_full_tree = "1"
    regArgs.reg_out_tag = "Real"
    regArgs.base_name = "regEleEcal{era_name}_RealIC".format(era_name=era_name)
    if run_step1: regArgs.run_eb_and_ee()
    
    
    #step2 do the E/p combination
    print "starting step2"
    input_for_comb = str(regArgs.applied_name())

    regArgs.base_name = "regEleEcalTrk{era_name}_RealIC".format(era_name=era_name)
    regArgs.var_eb =":".join(["(sc.rawEnergy+sc.rawESEnergy)*regRealMean","regRealSigma/regRealMean","ele.trkPModeErr/ele.trkPMode","(sc.rawEnergy+sc.rawESEnergy)*regRealMean/ele.trkPMode","ele.ecalDrivenSeed","ssFull.e3x3/sc.rawEnergy","ele.fbrem","ele.trkEtaMode","ele.trkPhiMode"])
    regArgs.var_ee =":".join(["(sc.rawEnergy+sc.rawESEnergy)*regRealMean","regRealSigma/regRealMean","ele.trkPModeErr/ele.trkPMode","(sc.rawEnergy+sc.rawESEnergy)*regRealMean/ele.trkPMode","ele.ecalDrivenSeed","ssFull.e3x3/sc.rawEnergy","ele.fbrem","ele.trkEtaMode","ele.trkPhiMode"])
    regArgs.target = "(mc.energy * (ele.trkPModeErr*ele.trkPModeErr + (sc.rawEnergy+sc.rawESEnergy)*(sc.rawEnergy+sc.rawESEnergy)*regRealSigma*regRealSigma) / ( (sc.rawEnergy+sc.rawESEnergy)*regRealMean*ele.trkPModeErr*ele.trkPModeErr + ele.trkPMode*(sc.rawEnergy+sc.rawESEnergy)*(sc.rawEnergy+sc.rawESEnergy)*regRealSigma*regRealSigma ))"
    regArgs.input_training = input_for_comb
    regArgs.input_testing = input_for_comb
    regArgs.write_full_tree = "0"  
    regArgs.fix_mean = False
    regArgs.reg_out_tag = "EcalTrk"
    regArgs.cuts_base = base_ele_cuts.format(extra_cuts = ep_eventnr_cut)

    if run_step2: regArgs.run_eb_and_ee()        
    
if __name__ =='__main__':
    main()



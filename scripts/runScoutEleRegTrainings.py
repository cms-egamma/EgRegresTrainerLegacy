#!/usr/bin/env python

import subprocess
import os
import regtools_scouting as regtools
from regtools_scouting import RegArgs
import argparse

def main():
    parser = argparse.ArgumentParser(description='runs the SC regression trainings')
    parser.add_argument('--era',required=True,help='year to produce for, 2016, 2017, 2018 are the options')
    parser.add_argument('--input_dir','-i',default='/mercury/data1/harper/EgRegsNtups',help='input directory with the ntuples')
    parser.add_argument('--output_dir','-o',default="results",help='output dir')
    args = parser.parse_args()

    #step 1, ecal only
    #step 2, trk
    #step 3, ecal-trk 


    run_step1 = False
    run_step2 = False
    run_step3 = True
    
    
    base_ele_cuts = "(mc.energy>0 && ele.pt>0 && {extra_cuts})"


    
    era_name = "Scout2023"
    input_real_ic  = "{}/DoubleElectron_Pt-1To300_gun_130X_mcRun3_2023_realistic_postBPix_v2-v2_AODSIM_EGRegSC_500k.root".format(args.input_dir)

    input_real_ic  = os.path.join(args.input_dir,"DoubleElectron_Pt-1To300_gun_130X_mcRun3_2023_realistic_postBPix_v2-v2_AODSIM_EGRegSC_remade.root")

    input_real_ic = os.path.join(args.input_dir,"DoubleElectron_Pt-1To300_gun_130X_mcRun3_2023_realistic_postBPix_v2-v2_AODSIM_EGRegSC_WithOffline.root")

    ecal_eventnr_cut = "evt.eventnr%10==0" 
    trk_eventnr_cut = "evt.eventnr%10==1" 
    ep_eventnr_cut = "evt.eventnr%10==2" 
    
    #step1 train the calo only regression using IDEAL intercalibration constants
    print("starting step1")
    regArgs = RegArgs()
    regArgs.input_training = str(input_real_ic)
    regArgs.input_testing = str(input_real_ic)  

    regArgs.cuts_base = base_ele_cuts.format(extra_cuts = ecal_eventnr_cut)
    regArgs.cuts_name = "stdCuts"
    regArgs.cfg_dir = "configs"
    regArgs.out_dir = "resultsEleTrkTrainOff" 
    regArgs.ntrees = 1500  
    regArgs.base_name = "regEleEcal{era_name}".format(era_name=era_name)
    regArgs.reg_out_tag = "Ecal"

    if run_step1:
        regArgs.run_eb_and_ee()

    #step2 this is the trk p training
    print("starting step2")
    input_for_trk_training = str(regArgs.applied_name())
    regArgs.base_name = "regEleTrk{era_name}".format(era_name=era_name)
    regArgs.set_trk_vars()
    regArgs.input_training = input_for_trk_training
    regArgs.input_testing = input_for_trk_training
    regArgs.write_full_tree = "1"
    regArgs.reg_out_tag = "Trk"
    regArgs.cuts_base = base_ele_cuts.format(extra_cuts = trk_eventnr_cut)
    if run_step2:
        regArgs.run_eb_and_ee()
    
    #step3 do the E/p combination
    #remember we use the Ideal Mean but Real Sigma (real mean is 1 by construction)
    print("starting step3")
    input_for_comb = str(regArgs.applied_name())

    regArgs.base_name = "regEleEcalTrkTrain{era_name}".format(era_name=era_name)
    regArgs.set_elecomb_trktrain()
    regArgs.input_training = input_for_comb
    regArgs.input_testing = input_for_comb
    regArgs.write_full_tree = "0"  
    regArgs.fix_mean = False
    regArgs.reg_out_tag = "EcalTrkTrain"
    regArgs.cuts_base = base_ele_cuts.format(extra_cuts = ep_eventnr_cut)
    if run_step3: 
        regArgs.run_eb_and_ee()
    
        
    
if __name__ =='__main__':
    main()



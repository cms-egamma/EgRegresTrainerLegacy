#!/usr/bin/env python

import subprocess
import os
import regtools
from regtools import RegArgs
import argparse

def main():
    parser = argparse.ArgumentParser(description='runs the SC regression trainings')
    
    parser.add_argument('--input_dir','-i',default='/mercury/data3/harper/workexpfiles/',help='input directory with the ntuples')
    parser.add_argument('--output_dir','-o',default="results",help='output dir')
    args = parser.parse_args()

    #step 1, ecal only
    #step 2, trk
    #step 3, ecal-trk 


    run_step1 = True
    run_step2 = False
    run_step3 = False
    
    
    era_name = "Retrain2026_FullPosEtaZQ_EtaOrPhiGap1Clus"
    input_real_ic  = "{}/DoubleElectron_FlatPT-1t500_skimmedntup.root".format(args.input_dir)

    base_ele_cuts = "(mc.energy>0 && ssFrac.sigmaIEtaIEta>0 && ssFrac.sigmaIPhiIPhi>0 && ele.et>0 && {extra_cuts})"

    base_ele_cuts = "(mc.energy>0 && ssFrac.sigmaIEtaIEta>0 && ssFrac.sigmaIPhiIPhi>0 && ele.et>0 && sc.scEta>0 && eleExtra.q>0 && eleExtra.z0>0 && ((sc.etaGapCode!=0 || sc.phiGapCode!=0) || !sc.isEB) && sc.numberOfClusters==1 && {extra_cuts})"

    ecal_eventnr_cut = "evt.eventnr%2==0" 
    trk_eventnr_cut = "evt.eventnr%10==1" 
    ep_eventnr_cut = "evt.eventnr%10==2" 
    
    #step1 train the calo only regression using IDEAL intercalibration constants
    print("starting step1")
    regArgs = RegArgs()
    regArgs.set_ecal_default()
    regArgs.input_training = str(input_real_ic)
    regArgs.input_testing = str(input_real_ic)  


    #regArgs.var_eb = ':'.join(["sc.rawEnergy","sc.etaWidth","sc.phiWidth","sc.seedClusEnergy/sc.rawEnergy","ssFull.e5x5/sc.rawEnergy",
    #                          "sc.scEta",
    #                          "rho","ssFull.sigmaIEtaIEta","ssFull.sigmaIEtaIPhi","ssFull.sigmaIPhiIPhi",#12
    #                         ])
    regArgs.var_eb+=":eleExtra.dEtaInSeed*abs(sc.scEta)/sc.scEta"
    regArgs.var_eb+=":eleExtra.dPhiIn*eleExtra.q"
    regArgs.var_eb+=":eleExtra.z0"
    regArgs.var_eb+=":eleExtra.q"


    regArgs.cuts_base = base_ele_cuts.format(extra_cuts = ecal_eventnr_cut)
    regArgs.cuts_name = "stdCuts"
    regArgs.cfg_dir = "configs"
    regArgs.out_dir = "resultsEle2" 
    regArgs.ntrees = 1500 
    regArgs.min_events = 300
    regArgs.base_name = "regEleEcal{era_name}".format(era_name=era_name)
    regArgs.reg_out_tag = "Ecal"

    if run_step1:
        regArgs.run_eb_and_ee()

 
    
        
    
if __name__ =='__main__':
    main()



#!/usr/bin/env python

import subprocess
import os

class RegArgs:
    def set_defaults(self):
        self.base_name = "reg_eleBParking"
        self.cuts_name = "stdCuts"
        self.vars_name = "stdVar"  
        self.cfg_dir = "configs"
        self.out_dir = "results"
        self.tree_name = "egRegTree"
        self.write_full_tree = "0"
        self.min_events = 300
        self.shrinkage = 0.15
        self.min_significance = 5.0
        self.event_weight = 1.
        self.input_testing = "test.root"
        self.input_training = "train.root"
        self.target = "mc.energy/(sc.rawEnergy + sc.rawESEnergy)"
        self.var_eb = ':'.join(["sc.rawEnergy","sc.etaWidth","sc.phiWidth","sc.seedClusEnergy/sc.rawEnergy","ssFull.e5x5/sc.rawEnergy","ele.hademTow",
                                "rho","sc.dEtaSeedSC","sc.dPhiSeedSC","ssFull.e3x3/sc.rawEnergy","ssFull.sigmaIEtaIEta","ssFull.sigmaIEtaIPhi","ssFull.sigmaIPhiIPhi",#12
                                "ssFull.eMax/ssFull.e5x5","ssFull.e2nd/ssFull.e5x5","ssFull.eTop/ssFull.e5x5","ssFull.eBottom/ssFull.e5x5","ssFull.eLeft/ssFull.e5x5","ssFull.eRight/ssFull.e5x5",
                                "ssFull.e2x5Max/ssFull.e5x5","ssFull.e2x5Left/ssFull.e5x5","ssFull.e2x5Right/ssFull.e5x5","ssFull.e2x5Top/ssFull.e5x5","ssFull.e2x5Bottom/ssFull.e5x5",
                                "ele.nrSatCrys","sc.numberOfClusters","sc.iEtaOrX","sc.iPhiOrY","sc.iEtaMod5","sc.iPhiMod2","sc.iEtaMod20","sc.iPhiMod20"])

        self.var_ee = ':'.join(["sc.rawEnergy","sc.etaWidth","sc.phiWidth","sc.seedClusEnergy/sc.rawEnergy","ssFull.e5x5/sc.rawEnergy","ele.hademTow", #5
                                "rho","sc.dEtaSeedSC","sc.dPhiSeedSC","ssFull.e3x3/sc.rawEnergy","ssFull.sigmaIEtaIEta","ssFull.sigmaIEtaIPhi","ssFull.sigmaIPhiIPhi", #12
                                "ssFull.eMax/ssFull.e5x5","ssFull.e2nd/ssFull.e5x5","ssFull.eTop/ssFull.e5x5","ssFull.eBottom/ssFull.e5x5","ssFull.eLeft/ssFull.e5x5","ssFull.eRight/ssFull.e5x5",
                                "ssFull.e2x5Max/ssFull.e5x5","ssFull.e2x5Left/ssFull.e5x5","ssFull.e2x5Right/ssFull.e5x5","ssFull.e2x5Top/ssFull.e5x5","ssFull.e2x5Bottom/ssFull.e5x5",
                                "ele.nrSatCrys","sc.numberOfClusters","sc.iEtaOrX","sc.iPhiOrY","sc.rawESEnergy/sc.rawEnergy"])

        self.cuts_base = "(mc.energy>0 && ssFull.sigmaIEtaIEta>0 && ssFull.sigmaIPhiIPhi>0 && ele.et>0 && eventnr%2==0)"
        self.ntrees = 1500
        self.do_eb = True

    def __init__(self):
        self.set_defaults()

    def name(self):
        if self.do_eb: region = "EB"
        else: region = "EE"
        return "{args.base_name}_{args.vars_name}_{args.cuts_name}_{region}_ntrees{args.ntrees}".format(args=self,region=region)

    def applied_name(self):
        return "{args.out_dir}/{args.base_name}_{args.vars_name}_{args.cuts_name}_ntrees{args.ntrees}_applied.root".format(args=self)
    
    def cfg_name(self):
        return "{}/{}.config".format(self.cfg_dir,self.name())

    def output_name(self):
        return "{}/{}_results.root".format(self.out_dir,self.name())
        
    
def make_cfg(args):
    base_cfg = """
Trainer: GBRLikelihoodTrain
NumberOfRegressions: 1
TMVAFactoryOptions: !V:!Silent:!Color:!DrawProgressBar
OutputDirectory: {args.out_dir}
Regression.1.Name: {name}
Regression.1.InputFiles: {args.input_training}
Regression.1.Tree: {args.tree_name}
Regression.1.trainingOptions: SplitMode=random:!V
Regression.1.Options: MinEvents={args.min_events}:Shrinkage={args.shrinkage}:NTrees={args.ntrees}:MinSignificance={args.min_significance}:EventWeight={args.event_weight}
Regression.1.DoCombine: False
Regression.1.DoEB: {args.do_eb}
Regression.1.VariablesEB: {args.var_eb}
Regression.1.VariablesEE: {args.var_ee}
Regression.1.Target: {args.target}
Regression.1.CutBase: {args.cuts_base} 
Regression.1.CutEB: sc.isEB
Regression.1.CutEE: !sc.isEB

""".format(args=args,name=args.name())
    if not os.path.isdir(args.cfg_dir):
        os.mkdir(args.cfg_dir)
    with open(args.cfg_name(),"w") as f:
        f.write(base_cfg)
    

def run_eb_and_ee(regArgs):  

    if not os.path.isdir(regArgs.out_dir):
        os.mkdir(regArgs.out_dir)

    regArgs.do_eb = True
    make_cfg(args=regArgs)
    print "starting: {}".format(regArgs.name())
    subprocess.Popen(["bin/slc6_amd64_gcc700/RegressionTrainerExe",regArgs.cfg_name()]).communicate()
    forest_eb_file = regArgs.output_name()
    
    regArgs.do_eb = False
    make_cfg(args=regArgs)
    print "starting: {}".format(regArgs.name())
    subprocess.Popen(["bin/slc6_amd64_gcc700/RegressionTrainerExe",regArgs.cfg_name()]).communicate()
    forest_ee_file = regArgs.output_name()

    subprocess.Popen(["bin/slc6_amd64_gcc700/RegressionApplierExe",regArgs.input_testing,regArgs.applied_name(),"--gbrForestFileEE",forest_ee_file,"--gbrForestFileEB",forest_eb_file,"--nrThreads","4","--treeName",regArgs.tree_name,"--writeFullTree",regArgs.write_full_tree]).communicate()
    
    print "made ",regArgs.applied_name()

def main():

    #step 1, run calo only regression and stick it into a tree so it can be used for ecal-trk combination
    #step 2, put apply the regresion to the real IC and save the result in a tree
    #step 3, run trk-calo regression
    run_step1 = False
    run_step2 = True
    run_step3 = True

    #step1 train the calo only regression using IDEAL intercalibration constants
    regArgs = RegArgs()
    regArgs.input_training = "/mercury/data1/harper/BParking_pt1To30Ntup/doubleElectron_FlatPt-1To300_FlatPU0to70RAWIDEALIC_102X_upgrade2018_realistic_v15-v1_BParkRECO_ntupV3.root"
    regArgs.input_testing = "/mercury/data1/harper/BParking_pt1To30Ntup/doubleElectron_FlatPt-1To300_FlatPU0to70RAWIDEALIC_102X_upgrade2018_realistic_v15-v1_BParkRECO_ntupV3.root"
    regArgs.target = "mc.energy/(sc.rawEnergy + sc.rawESEnergy)"
    regArgs.cfg_dir = "configs"
    regArgs.out_dir = "results" 
    regArgs.ntrees = 1500  
    regArgs.base_name = "regEleBParkingIDEALECAL"
    if run_step1: run_eb_and_ee(regArgs=regArgs)
    
    #step2 now we run over the REAL intercalibration constant data and make a rew tree with this regression included
    regArgs.do_eb = True
    forest_eb_file = regArgs.output_name()
    regArgs.do_eb = False
    forest_ee_file = regArgs.output_name()

    input_real_ic = "/mercury/data1/harper/BParking_pt1To30Ntup/doubleElectron_FlatPt-1To300_FlatPU0to70RAW_102X_upgrade2018_realistic_v15-v1_BParkRECO_ntupV3.root"
    regArgs.base_name = "regEleBParkingECAL_IDEALTraining"
    input_for_comb = str(regArgs.applied_name()) #save the output name before we change it
    if run_step2: subprocess.Popen(["bin/slc6_amd64_gcc700/RegressionApplierExe",input_real_ic,input_for_comb,"--gbrForestFileEE",forest_ee_file,"--gbrForestFileEB",forest_eb_file,"--nrThreads","4","--treeName",regArgs.tree_name,"--writeFullTree","1","--regOutTag","Ecal"]).communicate()
    

    #step3 do the E/p combination
    regArgs.base_name = "regEleBParkingIDEALECALTrk"
    regArgs.var_eb =":".join(["(sc.rawEnergy+sc.rawESEnergy)*regEcalMean","regEcalMean/regEcalSigma","ele.trkPModeErr/ele.trkPMode","(sc.rawEnergy+sc.rawESEnergy)*regEcalMean/ele.trkPMode","ele.ecalDrivenSeed","ssFull.e3x3/sc.rawEnergy","ele.fbrem","ele.trkEtaMode","ele.trkPhiMode"])
    regArgs.var_ee =":".join(["(sc.rawEnergy+sc.rawESEnergy)*regEcalMean","regEcalMean/regEcalSigma","ele.trkPModeErr/ele.trkPMode","(sc.rawEnergy+sc.rawESEnergy)*regEcalMean/ele.trkPMode","ele.ecalDrivenSeed","ssFull.e3x3/sc.rawEnergy","ele.fbrem","ele.trkEtaMode","ele.trkPhiMode"])
    regArgs.target = "(mc.energy * (ele.trkPModeErr*ele.trkPModeErr + (sc.rawEnergy+sc.rawESEnergy)*(sc.rawEnergy+sc.rawESEnergy)*regEcalSigma*regEcalSigma) / ( (sc.rawEnergy+sc.rawESEnergy)*regEcalMean*ele.trkPModeErr*ele.trkPModeErr + ele.trkPMode*(sc.rawEnergy+sc.rawESEnergy)*(sc.rawEnergy+sc.rawESEnergy)*regEcalSigma*regEcalSigma ))"
    regArgs.input_training = input_for_comb
    regArgs.input_testing = input_for_comb
    regArgs.write_full_tree = "1"
    if run_step3: run_eb_and_ee(regArgs=regArgs)
    
    
if __name__ =='__main__':
    main()



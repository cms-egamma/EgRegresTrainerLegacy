#!/usr/bin/env python

import subprocess
import os

class RegArgs:
    def set_defaults(self):
        self.base_name = "reg_sc"
        self.cuts_name = "stdCuts"
        self.vars_name = "stdVar"  
        self.cfg_dir = "configs"
        self.out_dir = "results"
        self.input_testing = "test.root"
        self.input_training = "train.root"
        self.var_eb = "nrVert:sc.rawEnergy:sc.etaWidth:sc.phiWidth:ssFrac.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:ssFrac.eMax/sc.rawEnergy:ssFrac.e2nd/sc.rawEnergy:ssFrac.eLeftRightDiffSumRatio:ssFrac.eTopBottomDiffSumRatio:ssFrac.sigmaIEtaIEta:ssFrac.sigmaIEtaIPhi:ssFrac.sigmaIPhiIPhi:sc.numberOfClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.iEtaOrX:sc.iPhiOrY"  
        self.var_ee = "nrVert:sc.rawEnergy:sc.etaWidth:sc.phiWidth:ssFrac.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:ssFrac.eMax/sc.rawEnergy:ssFrac.e2nd/sc.rawEnergy:ssFrac.eLeftRightDiffSumRatio:ssFrac.eTopBottomDiffSumRatio:ssFrac.sigmaIEtaIEta:ssFrac.sigmaIEtaIPhi:ssFrac.sigmaIPhiIPhi:sc.numberOfClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.iEtaOrX:sc.iPhiOrY:sc.seedEta"
        self.cuts_base = "(mc.energy>0 && ssFrac.sigmaIEtaIEta>0 && ssFrac.sigmaIPhiIPhi>0 && evt.eventnr%2==0)"
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
Regression.1.Tree: egRegTree
Regression.1.Method: BDT
Regression.1.trainingOptions: SplitMode=random:!V
Regression.1.Options: MinEvents=300:Shrinkage=0.15:NTrees={args.ntrees}:MinSignificance=5.0:EventWeight=1
Regression.1.DoErrors: True
Regression.1.DoCombine: False
Regression.1.DoEB: {args.do_eb}
Regression.1.VariablesEB: {args.var_eb}
Regression.1.VariablesEE: {args.var_ee}
Regression.1.Target: mc.energy / ( sc.rawEnergy  )
Regression.1.TargetError: 1.253*abs( BDTresponse - mc.energy / ( sc.rawEnergy ) )
Regression.1.HistoConfig: jobs/dummy_Histo.config
Regression.1.CutBase: {args.cuts_base} 
Regression.1.CutEB: sc.isEB
Regression.1.CutEE: !sc.isEB
Regression.1.CutError: {args.cuts_base}

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

    subprocess.Popen(["bin/slc6_amd64_gcc700/RegressionApplierExe",regArgs.input_testing,regArgs.applied_name(),"--gbrForestFileEE",forest_ee_file,"--gbrForestFileEB",forest_eb_file,"--nrThreads","4"]).communicate()
    
    print "made ",regArgs.applied_name()
    
def main():

    #modify the parameters as you wish and then re-run

    regArgs = RegArgs()
    regArgs.input_training = "/mercury/data1/harper/EgRegsNtups/DoubleElectron_FlatPt-1To300_2017ConditionsFlatPU0to70ECALGT_105X_mc2017_realistic_IdealEcalIC_v5-v2_AODSIM_EgRegTreeV1_extraVars.root"
    regArgs.input_testing = "/mercury/data1/harper/EgRegsNtups/DoubleElectron_FlatPt-1To300_2017ConditionsFlatPU0to70ECALGT_105X_mc2017_realistic_IdealEcalIC_v5-v2_AODSIM_EgRegTreeV1_extraVars.root"
    regArgs.cfg_dir = "configs"
    regArgs.out_dir = "results"
    regArgs.cuts_name = "stdCuts" 
    regArgs.base_name = "scReg_1050_std"
    regArgs.cuts_base = "(mc.energy>0 && ssFrac.sigmaIEtaIEta>0 && ssFrac.sigmaIPhiIPhi>0 && evt.eventnr%2==0)"
    regArgs.ntrees = 1500

    regArgs.vars_name = "stdVarsNoWidth"
    regArgs.var_eb = "nrVert:sc.rawEnergy:ssFrac.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:ssFrac.eMax/sc.rawEnergy:ssFrac.e2nd/sc.rawEnergy:ssFrac.eLeftRightDiffSumRatio:ssFrac.eTopBottomDiffSumRatio:ssFrac.sigmaIEtaIEta:ssFrac.sigmaIEtaIPhi:ssFrac.sigmaIPhiIPhi:sc.numberOfClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.iEtaOrX:sc.iPhiOrY"  
    regArgs.var_ee = "nrVert:sc.rawEnergy:ssFrac.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:ssFrac.eMax/sc.rawEnergy:ssFrac.e2nd/sc.rawEnergy:ssFrac.eLeftRightDiffSumRatio:ssFrac.eTopBottomDiffSumRatio:ssFrac.sigmaIEtaIEta:ssFrac.sigmaIEtaIPhi:ssFrac.sigmaIPhiIPhi:sc.numberOfClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.iEtaOrX:sc.iPhiOrY:sc.seedEta"

    regArgs.vars_name = "stdVarsNoWidthFull"
    regArgs.var_eb = "nrVert:sc.rawEnergy:ssFull.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:ssFull.eMax/sc.rawEnergy:ssFull.e2nd/sc.rawEnergy:ssFull.eLeftRightDiffSumRatio:ssFull.eTopBottomDiffSumRatio:ssFull.sigmaIEtaIEta:ssFull.sigmaIEtaIPhi:ssFull.sigmaIPhiIPhi:sc.numberOfClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.iEtaOrX:sc.iPhiOrY"  
    regArgs.var_ee = "nrVert:sc.rawEnergy:ssFull.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:ssFull.eMax/sc.rawEnergy:ssFull.e2nd/sc.rawEnergy:ssFull.eLeftRightDiffSumRatio:ssFull.eTopBottomDiffSumRatio:ssFull.sigmaIEtaIEta:ssFull.sigmaIEtaIPhi:ssFull.sigmaIPhiIPhi:sc.numberOfClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.iEtaOrX:sc.iPhiOrY:sc.seedEta"

    good_ele = False
    if good_ele:
        regArgs.vars_name = "goodEleVars"
        regArgs.cuts_name = "goodEleCuts"
        regArgs.cuts_base = "(mc.energy>0 && ssFrac.sigmaIEtaIEta>0 && ssFrac.sigmaIPhiIPhi>0 && evt.eventnr%2==0 && (sc.etaGapCode==0 && sc.phiGapCode==0 && sc.nearbyChanStatus==0))"
        regArgs.var_eb = "nrVert:sc.rawEnergy:ssFrac.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:ssFrac.eMax/sc.rawEnergy:ssFrac.e2nd/sc.rawEnergy:ssFrac.eLeftRightDiffSumRatio:ssFrac.eTopBottomDiffSumRatio:ssFrac.sigmaIEtaIEta:ssFrac.sigmaIEtaIPhi:ssFrac.sigmaIPhiIPhi:sc.numberOfClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.seedEta"
        regArgs.var_ee = "nrVert:sc.rawEnergy:ssFrac.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:ssFrac.eMax/sc.rawEnergy:ssFrac.e2nd/sc.rawEnergy:ssFrac.eLeftRightDiffSumRatio:ssFrac.eTopBottomDiffSumRatio:ssFrac.sigmaIEtaIEta:ssFrac.sigmaIEtaIPhi:ssFrac.sigmaIPhiIPhi:sc.numberOfClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.seedEta"
    else:
        regArgs.vars_name = "badEleVars"
        #regArgs.cuts_name = "badEleCuts"
        #regArgs.cuts_base = "(mc.energy>0 && ssFrac.sigmaIEtaIEta>0 && ssFrac.sigmaIPhiIPhi>0 && evt.eventnr%2==0 && !(sc.etaGapCode==0 && sc.phiGapCode==0 && sc.nearbyChanStatus==0))"
        regArgs.var_eb = "nrVert:sc.rawEnergy:ssFrac.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:ssFrac.eMax/sc.rawEnergy:ssFrac.e2nd/sc.rawEnergy:ssFrac.eLeftRightDiffSumRatio:ssFrac.eTopBottomDiffSumRatio:ssFrac.sigmaIEtaIEta:ssFrac.sigmaIEtaIPhi:ssFrac.sigmaIPhiIPhi:sc.numberOfClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.seedEta:sc.etaGapCode:sc.phiGapCode:sc.nearbyChanStatus"
        regArgs.var_ee = "nrVert:sc.rawEnergy:ssFrac.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:ssFrac.eMax/sc.rawEnergy:ssFrac.e2nd/sc.rawEnergy:ssFrac.eLeftRightDiffSumRatio:ssFrac.eTopBottomDiffSumRatio:ssFrac.sigmaIEtaIEta:ssFrac.sigmaIEtaIPhi:ssFrac.sigmaIPhiIPhi:sc.numberOfClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.seedEta:sc.etaGapCode:sc.phiGapCode:sc.nearbyChanStatus" 

    
    regArgs.vars_name = "stdVarsNoWidthIValues"
    regArgs.var_eb = "nrVert:sc.rawEnergy:ssFrac.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:ssFrac.eMax/sc.rawEnergy:ssFrac.e2nd/sc.rawEnergy:ssFrac.eLeftRightDiffSumRatio:ssFrac.eTopBottomDiffSumRatio:ssFrac.sigmaIEtaIEta:ssFrac.sigmaIEtaIPhi:ssFrac.sigmaIPhiIPhi:sc.numberOfClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.seedEta:sc.etaGapCode:sc.phiGapCode:sc.nearbyChanStatus" 
    regArgs.var_ee = "nrVert:sc.rawEnergy:ssFrac.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:ssFrac.eMax/sc.rawEnergy:ssFrac.e2nd/sc.rawEnergy:ssFrac.eLeftRightDiffSumRatio:ssFrac.eTopBottomDiffSumRatio:ssFrac.sigmaIEtaIEta:ssFrac.sigmaIEtaIPhi:ssFrac.sigmaIPhiIPhi:sc.numberOfClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.seedEta:sc.etaGapCode:sc.phiGapCode:sc.nearbyChanStatus" 
    regArgs.vars_name = "stdVarsNoWidthNoGapDead"
    regArgs.var_eb = "nrVert:sc.rawEnergy:ssFrac.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:ssFrac.eMax/sc.rawEnergy:ssFrac.e2nd/sc.rawEnergy:ssFrac.eLeftRightDiffSumRatio:ssFrac.eTopBottomDiffSumRatio:ssFrac.sigmaIEtaIEta:ssFrac.sigmaIEtaIPhi:ssFrac.sigmaIPhiIPhi:sc.numberOfClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.seedEta" 
    regArgs.var_ee = "nrVert:sc.rawEnergy:ssFrac.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:ssFrac.eMax/sc.rawEnergy:ssFrac.e2nd/sc.rawEnergy:ssFrac.eLeftRightDiffSumRatio:ssFrac.eTopBottomDiffSumRatio:ssFrac.sigmaIEtaIEta:ssFrac.sigmaIEtaIPhi:ssFrac.sigmaIPhiIPhi:sc.numberOfClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.seedEta" 

    regArgs.vars_name = "stdVarsOld"
    regArgs.var_eb = "nrVert:sc.rawEnergy:sc.etaWidth:sc.phiWidth:ssFrac.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:ssFrac.eMax/sc.rawEnergy:ssFrac.e2nd/sc.rawEnergy:ssFrac.eLeftRightDiffSumRatio:ssFrac.eTopBottomDiffSumRatio:ssFrac.sigmaIEtaIEta:ssFrac.sigmaIEtaIPhi:ssFrac.sigmaIPhiIPhi:sc.numberOfClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.iEtaOrX:sc.iPhiOrY"  
    regArgs.var_ee = "nrVert:sc.rawEnergy:sc.etaWidth:sc.phiWidth:ssFrac.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:ssFrac.eMax/sc.rawEnergy:ssFrac.e2nd/sc.rawEnergy:ssFrac.eLeftRightDiffSumRatio:ssFrac.eTopBottomDiffSumRatio:ssFrac.sigmaIEtaIEta:ssFrac.sigmaIEtaIPhi:ssFrac.sigmaIPhiIPhi:sc.numberOfClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.iEtaOrX:sc.iPhiOrY"
    run_eb_and_ee(regArgs=regArgs)
    
if __name__ =='__main__':
    main()



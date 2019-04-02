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
        self.var_eb = "nrVert:sc.rawEnergy:sc.etaWidth:sc.phiWidth:ssFrac.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:ssFrac.eMax/sc.rawEnergy:ssFrac.e2nd/sc.rawEnergy:ssFrac.eLeftRightDiffSumRatio:ssFrac.eTopBottomDiffSumRatio:ssFrac.sigmaIEtaIEta:ssFrac.sigmaIEtaIPhi:ssFrac.sigmaIPhiIPhi:sc.numberOfSubClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.iEtaOrX:sc.iPhiOrY"
        self.var_ee = "nrVert:sc.rawEnergy:sc.etaWidth:sc.phiWidth:ssFrac.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:ssFrac.eMax/sc.rawEnergy:ssFrac.e2nd/sc.rawEnergy:ssFrac.eLeftRightDiffSumRatio:ssFrac.eTopBottomDiffSumRatio:ssFrac.sigmaIEtaIEta:ssFrac.sigmaIEtaIPhi:ssFrac.sigmaIPhiIPhi:sc.numberOfSubClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.iEtaOrX:sc.iPhiOrY:sc.seedEta"

        self.cuts_base = "(mc.energy>0 && ssFrac.sigmaIEtaIEta>0 && ssFrac.sigmaIPhiIPhi>0)"
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

    subprocess.Popen(["bin/slc6_amd64_gcc700/RegressionApplierExe",regArgs.input_testing,regArgs.applied_name(),"--gbrForestFileEE",forest_ee_file,"--gbrForestFileEB",forest_eb_file,"--nrThreads","4","--treeName","egRegTree"]).communicate()
    
    print "made ",regArgs.applied_name()
    
def run_job(scenario="AC1Sigma"):
    #modify the parameters as you wish and then re-run
    regArgs = RegArgs()
    regArgs.input_training = "/mercury/data1/harper/mcFiles/pfRecHitValid/DoubleElePt1To100_EGReg/doubleElePt1To100_EGRegV2_{scenario}.root".format(scenario=scenario)
    regArgs.input_testing = "/mercury/data1/harper/mcFiles/pfRecHitValid/DoubleElePt1To100_EGReg/doubleElePt1To100_EGRegV2_{scenario}.root".format(scenario=scenario)
    regArgs.cfg_dir = "configs"
    regArgs.out_dir = "results"
    regArgs.base_name = "scReg_mustFixed_{scenario}".format(scenario=scenario)
   # regArgs.cuts_name = "stdCutsAllEvts" 
   # regArgs.cuts_base = "(mc.energy>0 && sc.sigmaIEtaIEta>0 && sc.sigmaIPhiIPhi>0)"
 #   regArgs.ntrees = 1500
    #regArgs.vars_name = "stdVarsNoWidth"
    #regArgs.var_eb = "nrVert:sc.rawEnergy:sc.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:sc.eMax/sc.rawEnergy:sc.e2nd/sc.rawEnergy:sc.eLeftRightDiffSumRatio:sc.eTopBottomDiffSumRatio:sc.sigmaIEtaIEta:sc.sigmaIEtaIPhi:sc.sigmaIPhiIPhi:sc.numberOfClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.iEtaOrX:sc.iPhiOrY"
    #regArgs.var_ee = "nrVert:sc.rawEnergy:sc.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:sc.eMax/sc.rawEnergy:sc.e2nd/sc.rawEnergy:sc.eLeftRightDiffSumRatio:sc.eTopBottomDiffSumRatio:sc.sigmaIEtaIEta:sc.sigmaIEtaIPhi:sc.sigmaIPhiIPhi:sc.numberOfClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.iEtaOrX:sc.iPhiOrY:sc.seedEta"
   # regArgs.vars_name = "stdVarsNoWidthNoSigma"
   # regArgs.var_eb = "nrVert:sc.rawEnergy:sc.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:sc.eMax/sc.rawEnergy:sc.e2nd/sc.rawEnergy:sc.eLeftRightDiffSumRatio:sc.eTopBottomDiffSumRatio:sc.numberOfClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.iEtaOrX:sc.iPhiOrY"
    #regArgs.var_ee = "nrVert:sc.rawEnergy:sc.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:sc.eMax/sc.rawEnergy:sc.e2nd/sc.rawEnergy:sc.eLeftRightDiffSumRatio:sc.eTopBottomDiffSumRatio:sc.numberOfClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.iEtaOrX:sc.iPhiOrY:seedEta"

  #  regArgs.vars_name = "stdVarsNoWidthNoSigma"
   # regArgs.var_eb = "nrVert:sc.rawEnergy:sc.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:sc.eMax/sc.rawEnergy:sc.e2nd/sc.rawEnergy:sc.eLeftRightDiffSumRatio:sc.eTopBottomDiffSumRatio:sc.numberOfClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.iEtaOrX:sc.iPhiOrY"
  #  regArgs.var_ee = "nrVert:sc.rawEnergy:sc.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:sc.eMax/sc.rawEnergy:sc.e2nd/sc.rawEnergy:sc.eLeftRightDiffSumRatio:sc.eTopBottomDiffSumRatio:sc.numberOfClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.iEtaOrX:sc.iPhiOrY:seedEta"
    run_eb_and_ee(regArgs=regArgs)
    
def main():
#    scenarios = ["AB1Sigma","AC1Sigma","AC2Sigma","AC3Sigma","NoThres"]
    scenarios = ["ACMixedSigma","AC2Sigma","AC3Sigma"]
    scenarios = ["ACMixedSigmaIdealIC"]
    for scenario in scenarios:
        run_job(scenario=scenario)

if __name__ =='__main__':
    main()



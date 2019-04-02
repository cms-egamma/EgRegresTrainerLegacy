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
        self.input_testing = "test.root"
        self.input_training = "train.root"
        self.var_eb = ':'.join(["sc.rawEnergy","sc.etaWidth","sc.phiWidth","sc.seedClusEnergy/sc.rawEnergy","ssFull.e5x5/sc.rawEnergy","ele.hademTow",
                                "rho","sc.dEtaSeedSC","sc.dPhiSeedSC","ssFull.e3x3/sc.rawEnergy","ssFull.sigmaIEtaIEta","ssFull.sigmaIEtaIPhi","ssFull.sigmaIPhiIPhi",
                                "ssFull.eMax/ssFull.e5x5","ssFull.e2nd/ssFull.e5x5","ssFull.eTop/ssFull.e5x5","ssFull.eBottom/ssFull.e5x5","ssFull.eLeft*ssFull.e5x5","ssFull.eRight*ssFull.e5x5",
                                "ssFull.e2x5Max/ssFull.e5x5","ssFull.e2x5Left/ssFull.e5x5","ssFull.e2x5Right/ssFull.e5x5","ssFull.e2x5Top/ssFull.e5x5","ssFull.e2x5Bottom/ssFull.e5x5",
                                "sc.numberOfClusters","sc.iEtaOrX","sc.iPhiOrY","sc.iEtaMod5","sc.iPhiMod2","sc.iEtaMod20","sc.iPhiMod20"])

        self.var_ee = ':'.join(["sc.rawEnergy","sc.etaWidth","sc.phiWidth","sc.seedClusEnergy/sc.rawEnergy","ssFull.e5x5/sc.rawEnergy","ele.hademTow",
                                "rho","sc.dEtaSeedSC","sc.dPhiSeedSC","ssFull.e3x3/sc.rawEnergy","ssFull.sigmaIEtaIEta","ssFull.sigmaIEtaIPhi","ssFull.sigmaIPhiIPhi",
                                "ssFull.eMax/ssFull.e5x5","ssFull.e2nd/ssFull.e5x5","ssFull.eTop/ssFull.e5x5","ssFull.eBottom/ssFull.e5x5","ssFull.eLeft*ssFull.e5x5","ssFull.eRight*ssFull.e5x5",
                                "ssFull.e2x5Max/ssFull.e5x5","ssFull.e2x5Left/ssFull.e5x5","ssFull.e2x5Right/ssFull.e5x5","ssFull.e2x5Top/ssFull.e5x5","ssFull.e2x5Bottom/ssFull.e5x5",
                                "sc.numberOfClusters","sc.iEtaOrX","sc.iPhiOrY","sc.rawESEnergy/sc.rawEnergy"])

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
Regression.1.Tree: egRegTree
Regression.1.Method: BDT
Regression.1.trainingOptions: SplitMode=random:!V
Regression.1.Options: MinEvents=300:Shrinkage=0.15:NTrees={args.ntrees}:MinSignificance=5.0:EventWeight=1
Regression.1.DoErrors: True
Regression.1.DoCombine: False
Regression.1.DoEB: {args.do_eb}
Regression.1.VariablesEB: {args.var_eb}
Regression.1.VariablesEE: {args.var_ee}
Regression.1.TargetError: 1.253*abs( BDTresponse - (mc.energy * (ele.trkPModeErr*ele.trkPModeErr + (sc.rawEnergy+sc.rawESEnergy)*(sc.rawEnergy+sc.rawESEnergy)*regSigma*regSigma) / ( (sc.rawEnergy+sc.rawESEnergy)*regMean*ele.trkPModeErr*ele.trkPModeErr + ele.trkPMode*(sc.rawEnergy+sc.rawESEnergy)*(sc.rawEnergy+sc.rawESEnergy)*regSigma*regSigma ) ) )
Regression.1.Target: (mc.energy * (ele.trkPModeErr*ele.trkPModeErr + (sc.rawEnergy+sc.rawESEnergy)*(sc.rawEnergy+sc.rawESEnergy)*regSigma*regSigma) / ( (sc.rawEnergy+sc.rawESEnergy)*regMean*ele.trkPModeErr*ele.trkPModeErr + ele.trkPMode*(sc.rawEnergy+sc.rawESEnergy)*(sc.rawEnergy+sc.rawESEnergy)*regSigma*regSigma ))
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
    
def main():

    #modify the parameters as you wish and then re-run

    regArgs = RegArgs()
    regArgs.input_training = "test.root"
    regArgs.input_testing = "test.root"
    regArgs.cfg_dir = "configs"
    regArgs.out_dir = "results"
    regArgs.ntrees = 1500
    regArgs.vars_name = "trkComb"
    regArgs.var_eb ="(sc.rawEnergy+sc.rawESEnergy)*regMean:regMean/regSigma:ele.trkPModeErr/ele.trkPMode:(sc.rawEnergy+sc.rawESEnergy)*regMean/ele.trkPMode:ele.fbrem:clus1.clusterRawEnergy/sc.rawEnergy:ele.trkEtaMode:trkPhiMode"
    regArgs.var_ee ="(sc.rawEnergy+sc.rawESEnergy)*regMean:regMean/regSigma:ele.trkPModeErr/ele.trkPMode:(sc.rawEnergy+sc.rawESEnergy)*regMean/ele.trkPMode:ele.fbrem:clus1.clusterRawEnergy/sc.rawEnergy:ele.trkEtaMode:trkPhiMode"

    
 #   regArgs.var_eb = "nrVert:sc.rawEnergy:sc.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:sc.eMax/sc.rawEnergy:sc.e2nd/sc.rawEnergy:sc.eLeftRightDiffSumRatio:sc.eTopBottomDiffSumRatio:sc.sigmaIEtaIEta:sc.sigmaIEtaIPhi:sc.sigmaIPhiIPhi:sc.numberOfClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.iEtaOrX:sc.iPhiOrY"
#    regArgs.var_ee = "nrVert:sc.rawEnergy:sc.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:sc.eMax/sc.rawEnergy:sc.e2nd/sc.rawEnergy:sc.eLeftRightDiffSumRatio:sc.eTopBottomDiffSumRatio:sc.sigmaIEtaIEta:sc.sigmaIEtaIPhi:sc.sigmaIPhiIPhi:sc.numberOfClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.iEtaOrX:sc.iPhiOrY:seedEta"
   # regArgs.vars_name = "stdVarsNoWidthNoSigma"
   # regArgs.var_eb = "nrVert:sc.rawEnergy:sc.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:sc.eMax/sc.rawEnergy:sc.e2nd/sc.rawEnergy:sc.eLeftRightDiffSumRatio:sc.eTopBottomDiffSumRatio:sc.numberOfClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.iEtaOrX:sc.iPhiOrY"
    #regArgs.var_ee = "nrVert:sc.rawEnergy:sc.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:sc.eMax/sc.rawEnergy:sc.e2nd/sc.rawEnergy:sc.eLeftRightDiffSumRatio:sc.eTopBottomDiffSumRatio:sc.numberOfClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.iEtaOrX:sc.iPhiOrY:seedEta"

  #  regArgs.vars_name = "stdVarsNoWidthNoSigma"
   # regArgs.var_eb = "nrVert:sc.rawEnergy:sc.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:sc.eMax/sc.rawEnergy:sc.e2nd/sc.rawEnergy:sc.eLeftRightDiffSumRatio:sc.eTopBottomDiffSumRatio:sc.numberOfClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.iEtaOrX:sc.iPhiOrY"
  #  regArgs.var_ee = "nrVert:sc.rawEnergy:sc.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:sc.eMax/sc.rawEnergy:sc.e2nd/sc.rawEnergy:sc.eLeftRightDiffSumRatio:sc.eTopBottomDiffSumRatio:sc.numberOfClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.iEtaOrX:sc.iPhiOrY:seedEta"
    run_eb_and_ee(regArgs=regArgs)
    
if __name__ =='__main__':
    main()



"""
regtools is a python module to collect various functions for running the regression
"""
import os
import subprocess

class RegArgs:
    def set_defaults(self):
        self.base_name = "reg_sc"
        self.cuts_name = "stdCuts"
        self.vars_name = "stdVar"  
        self.cfg_dir = "configs"
        self.out_dir = "results" 
        self.tree_name = "egRegTree"
        self.write_full_tree = "0"
        self.reg_out_tag = ""
        self.min_events = 300
        self.shrinkage = 0.15
        self.min_significance = 5.0
        self.event_weight = 1.
        self.mean_min = 0.2
        self.mean_max = 2.0
        self.fix_mean = False
        self.input_testing = "test.root"
        self.input_training = "train.root"
        self.target = "mc.energy/(sc.rawEnergy)"
        self.var_eb = "nrVert:sc.rawEnergy:sc.etaWidth:sc.phiWidth:ssFrac.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:ssFrac.eMax/sc.rawEnergy:ssFrac.e2nd/sc.rawEnergy:ssFrac.eLeftRightDiffSumRatio:ssFrac.eTopBottomDiffSumRatio:ssFrac.sigmaIEtaIEta:ssFrac.sigmaIEtaIPhi:ssFrac.sigmaIPhiIPhi:sc.numberOfSubClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.iEtaOrX:sc.iPhiOrY"  
        self.var_ee = "nrVert:sc.rawEnergy:sc.etaWidth:sc.phiWidth:ssFrac.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:ssFrac.eMax/sc.rawEnergy:ssFrac.e2nd/sc.rawEnergy:ssFrac.eLeftRightDiffSumRatio:ssFrac.eTopBottomDiffSumRatio:ssFrac.sigmaIEtaIEta:ssFrac.sigmaIEtaIPhi:ssFrac.sigmaIPhiIPhi:sc.numberOfSubClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.iEtaOrX:sc.iPhiOrY:sc.seedEta"
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


    def set_sc_default(self):
        self.target = "mc.energy/(sc.rawEnergy)"
        self.var_eb = "nrVert:sc.rawEnergy:sc.etaWidth:sc.phiWidth:ssFrac.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:ssFrac.eMax/sc.rawEnergy:ssFrac.e2nd/sc.rawEnergy:ssFrac.eLeftRightDiffSumRatio:ssFrac.eTopBottomDiffSumRatio:ssFrac.sigmaIEtaIEta:ssFrac.sigmaIEtaIPhi:ssFrac.sigmaIPhiIPhi:sc.numberOfSubClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.iEtaOrX:sc.iPhiOrY"  
        self.var_ee = "nrVert:sc.rawEnergy:sc.etaWidth:sc.phiWidth:ssFrac.e3x3/sc.rawEnergy:sc.seedClusEnergy/sc.rawEnergy:ssFrac.eMax/sc.rawEnergy:ssFrac.e2nd/sc.rawEnergy:ssFrac.eLeftRightDiffSumRatio:ssFrac.eTopBottomDiffSumRatio:ssFrac.sigmaIEtaIEta:ssFrac.sigmaIEtaIPhi:ssFrac.sigmaIPhiIPhi:sc.numberOfSubClusters:sc.clusterMaxDR:sc.clusterMaxDRDPhi:sc.clusterMaxDRDEta:sc.clusterMaxDRRawEnergy/sc.rawEnergy:clus1.clusterRawEnergy/sc.rawEnergy:clus2.clusterRawEnergy/sc.rawEnergy:clus3.clusterRawEnergy/sc.rawEnergy:clus1.clusterDPhiToSeed:clus2.clusterDPhiToSeed:clus3.clusterDPhiToSeed:clus1.clusterDEtaToSeed:clus2.clusterDEtaToSeed:clus3.clusterDEtaToSeed:sc.iEtaOrX:sc.iPhiOrY:sc.seedEta"
    
    def set_ecal_default(self):
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

    def set_phoecal_default(self):
        #note photon uses cone based H/E rather than tower based H/E as electrons do
        #also note, photon sigmaIEtaIPhi is different to the rest
        self.target = "mc.energy/(sc.rawEnergy + sc.rawESEnergy)"
        self.var_eb = ':'.join(["sc.rawEnergy","sc.etaWidth","sc.phiWidth","sc.seedClusEnergy/sc.rawEnergy","ssFull.e5x5/sc.rawEnergy","pho.hademCone",
                                "rho","sc.dEtaSeedSC","sc.dPhiSeedSC","ssFull.e3x3/sc.rawEnergy","ssFull.sigmaIEtaIEta","phoSSFull.sigmaIEtaIPhi","ssFull.sigmaIPhiIPhi",#12
                                "ssFull.eMax/ssFull.e5x5","ssFull.e2nd/ssFull.e5x5","ssFull.eTop/ssFull.e5x5","ssFull.eBottom/ssFull.e5x5","ssFull.eLeft/ssFull.e5x5","ssFull.eRight/ssFull.e5x5",
                                "ssFull.e2x5Max/ssFull.e5x5","ssFull.e2x5Left/ssFull.e5x5","ssFull.e2x5Right/ssFull.e5x5","ssFull.e2x5Top/ssFull.e5x5","ssFull.e2x5Bottom/ssFull.e5x5",
                                "pho.nrSatCrys","sc.numberOfClusters","sc.iEtaOrX","sc.iPhiOrY","sc.iEtaMod5","sc.iPhiMod2","sc.iEtaMod20","sc.iPhiMod20"])

        self.var_ee = ':'.join(["sc.rawEnergy","sc.etaWidth","sc.phiWidth","sc.seedClusEnergy/sc.rawEnergy","ssFull.e5x5/sc.rawEnergy","pho.hademCone", #5
                                "rho","sc.dEtaSeedSC","sc.dPhiSeedSC","ssFull.e3x3/sc.rawEnergy","ssFull.sigmaIEtaIEta","phoSSFull.sigmaIEtaIPhi","ssFull.sigmaIPhiIPhi", #12
                                "ssFull.eMax/ssFull.e5x5","ssFull.e2nd/ssFull.e5x5","ssFull.eTop/ssFull.e5x5","ssFull.eBottom/ssFull.e5x5","ssFull.eLeft/ssFull.e5x5","ssFull.eRight/ssFull.e5x5",
                                "ssFull.e2x5Max/ssFull.e5x5","ssFull.e2x5Left/ssFull.e5x5","ssFull.e2x5Right/ssFull.e5x5","ssFull.e2x5Top/ssFull.e5x5","ssFull.e2x5Bottom/ssFull.e5x5",
                                "pho.nrSatCrys","sc.numberOfClusters","sc.iEtaOrX","sc.iPhiOrY","sc.rawESEnergy/sc.rawEnergy"])

    def set_elecomb_default(self):
        self.var_eb =":".join(["(sc.rawEnergy+sc.rawESEnergy)*regEcalMean","regEcalSigma/regEcalMean","ele.trkPModeErr/ele.trkPMode","(sc.rawEnergy+sc.rawESEnergy)*regEcalMean/ele.trkPMode","ele.ecalDrivenSeed","ssFull.e3x3/sc.rawEnergy","ele.fbrem","ele.trkEtaMode","ele.trkPhiMode"])
        self.var_ee =":".join(["(sc.rawEnergy+sc.rawESEnergy)*regEcalMean","regEcalSigma/regEcalMean","ele.trkPModeErr/ele.trkPMode","(sc.rawEnergy+sc.rawESEnergy)*regEcalMean/ele.trkPMode","ele.ecalDrivenSeed","ssFull.e3x3/sc.rawEnergy","ele.fbrem","ele.trkEtaMode","ele.trkPhiMode"])
        self.target = "(mc.energy * (ele.trkPModeErr*ele.trkPModeErr + (sc.rawEnergy+sc.rawESEnergy)*(sc.rawEnergy+sc.rawESEnergy)*regEcalSigma*regEcalSigma) / ( (sc.rawEnergy+sc.rawESEnergy)*regEcalMean*ele.trkPModeErr*ele.trkPModeErr + ele.trkPMode*(sc.rawEnergy+sc.rawESEnergy)*(sc.rawEnergy+sc.rawESEnergy)*regEcalSigma*regEcalSigma ))"

    def make_cfg(self):
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
Regression.1.MeanMin: {args.mean_min}
Regression.1.MeanMax: {args.mean_max}
Regression.1.FixMean: {args.fix_mean}

""".format(args=self,name=self.name())
        if not os.path.isdir(self.cfg_dir):
            os.mkdir(self.cfg_dir)
        with open(self.cfg_name(),"w") as f:
            f.write(base_cfg)

    def run_eb_and_ee(self):  

        if not os.path.isdir(self.out_dir):
            os.mkdir(self.out_dir)

        self.do_eb = True
        self.make_cfg()

	# Scram arch has to be set manually
	# This also must be done in the training script in scripts/
	arch = "slc7_amd64_gcc700"

        print "starting: {}".format(self.name())
        subprocess.Popen(["bin/"+arch+"/RegressionTrainerExe",self.cfg_name()]).communicate()
        forest_eb_file = self.output_name()
    
        self.do_eb = False
        self.make_cfg()
        print "starting: {}".format(self.name())
        subprocess.Popen(["bin/"+arch+"/RegressionTrainerExe",self.cfg_name()]).communicate()
        forest_ee_file = self.output_name()

        
        subprocess.Popen(["bin/"+arch+"/RegressionApplierExe",self.input_testing,self.applied_name(),"--gbrForestFileEE",forest_ee_file,"--gbrForestFileEB",forest_eb_file,"--nrThreads","4","--treeName",self.tree_name,"--writeFullTree",self.write_full_tree,"--regOutTag",self.reg_out_tag]).communicate()

        print "made ",self.applied_name()


    def forest_filenames(self):
        do_eb_org = self.do_eb
        self.do_eb = True
        forest_eb_file = self.output_name()
        self.do_eb = False
        forest_ee_file = self.output_name()
        self.do_eb = do_eb_org
        return forest_eb_file,forest_ee_file

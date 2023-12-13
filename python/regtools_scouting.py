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
        self.tree_name = "egScoutRegTree"
        self.write_full_tree = "1"
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
        self.target = "mc.energy/ele.energy"
        self.var_eb = "rho:ele.energy:ele.eta:ele.phi:ele.sigmaIetaIeta:ele.r9:ele.sMin:ele.sMaj:ele.rechitZeroSuppression:ele.iEtaOrIX:ele.iPhiOrIY";        
        self.var_ee = "rho:ele.energy:ele.eta:ele.phi:ele.sigmaIetaIeta:ele.r9:ele.sMin:ele.sMaj:ele.rechitZeroSuppression:ele.iEtaOrIX:ele.iPhiOrIY";
        self.cuts_base = "(mc.energy>0 && ele.pt>0 && evt.eventnr%2==0)"
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

    def set_vars_withtrk(self):
        self.var_eb = "rho:ele.pt:ele.eta:ele.phi:ele.trkdz[ele.bestTrkIndx]:ele.trkpt[ele.bestTrkIndx]:ele.trketa[ele.bestTrkIndx]:ele.trkchi2overndf[ele.bestTrkIndx]:ele.sigmaIetaIeta:ele.r9:ele.sMin:ele.sMaj:ele.rechitZeroSuppression:ele.iEtaOrIX:ele.iPhiOrIY";        
        self.var_ee = "rho:ele.pt:ele.eta:ele.phi:ele.trkdz[ele.bestTrkIndx]:ele.trkpt[ele.bestTrkIndx]:ele.trketa[ele.bestTrkIndx]:ele.trkchi2overndf[ele.bestTrkIndx]:ele.sigmaIetaIeta:ele.r9:ele.sMin:ele.sMaj:ele.rechitZeroSuppression:ele.iEtaOrIX:ele.iPhiOrIY";
        
    def set_trk_vars(self):
        self.var_eb = "ele.trkdz[ele.bestTrkIndx]:ele.trkpt[ele.bestTrkIndx]:ele.trketa[ele.bestTrkIndx]:ele.trkphi[ele.bestTrkIndx]:ele.trkchi2overndf[ele.bestTrkIndx]";        
        self.var_ee = "ele.trkdz[ele.bestTrkIndx]:ele.trkpt[ele.bestTrkIndx]:ele.trketa[ele.bestTrkIndx]:ele.trkphi[ele.bestTrkIndx]:ele.trkchi2overndf[ele.bestTrkIndx]";  
        self.target = "mc.energy/ele.trkp[ele.bestTrkIndx]"

        
    def set_elecomb_default(self):
        self.var_eb =":".join(["ele.pt*regEcalMean","regEcalSigma/regEcalMean","ele.trkd0[ele.bestTrkIndx]","ele.pt*regEcalMean/ele.trkpt[ele.bestTrkIndx]","ele.trkchi2overndf[ele.bestTrkIndx]","ele.nrTrks"])
        self.var_ee =":".join(["ele.pt*regEcalMean","regEcalSigma/regEcalMean","ele.trkd0[ele.bestTrkIndx]","ele.pt*regEcalMean/ele.trkpt[ele.bestTrkIndx]","ele.trkchi2overndf[ele.bestTrkIndx]","ele.nrTrks"])
        self.target = "mc.energy * (0.0036*ele.trkp[ele.bestTrkIndx]*ele.trkp[ele.bestTrkIndx] + ele.energy*ele.energy*regEcalSigma*regEcalSigma) / (ele.energy*regEcalMean*ele.trkp[ele.bestTrkIndx]*ele.trkp[ele.bestTrkIndx]*0.0036 + ele.trkp[ele.bestTrkIndx]*ele.energy*ele.energy*regEcalSigma*regEcalSigma) "
        self.mean_min = 0.2
        self.mean_max = 3.0

    def set_elecomb_trktrain(self):
        self.var_eb =":".join(["ele.pt*regEcalMean","regEcalSigma/regEcalMean","regTrkSigma/regTrkMean","ele.trkd0[ele.bestTrkIndx]","ele.pt*regEcalMean/ele.trkpt[ele.bestTrkIndx]/regTrkMean","ele.trkchi2overndf[ele.bestTrkIndx]","ele.nrTrks"])
        self.var_ee =":".join(["ele.pt*regEcalMean","regEcalSigma/regEcalMean","regTrkSigma/regTrkMean","ele.trkd0[ele.bestTrkIndx]","ele.pt*regEcalMean/ele.trkpt[ele.bestTrkIndx]/regTrkMean","ele.trkchi2overndf[ele.bestTrkIndx]","ele.nrTrks"])
        self.target = "mc.energy * (regTrkSigma*regTrkSigma*ele.trkp[ele.bestTrkIndx]*ele.trkp[ele.bestTrkIndx] + ele.energy*ele.energy*regEcalSigma*regEcalSigma) / (ele.energy*regEcalMean*ele.trkp[ele.bestTrkIndx]*ele.trkp[ele.bestTrkIndx]*regTrkSigma*regTrkSigma + ele.trkp[ele.bestTrkIndx]*ele.energy*ele.energy*regEcalSigma*regEcalSigma) "
        self.mean_min = 0.2
        self.mean_max = 3.0

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
Regression.1.CutEB: ele.isEB
Regression.1.CutEE: !ele.isEB
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
        print("starting: {}".format(self.name()))
        print("cfg name",self.cfg_name())
        subprocess.Popen(["bin/el8_amd64_gcc11/RegressionTrainerExe",self.cfg_name()]).communicate()
        forest_eb_file = self.output_name()
    
        self.do_eb = False
        self.make_cfg()
        print("starting: {}".format(self.name()))
        subprocess.Popen(["bin/el8_amd64_gcc11/RegressionTrainerExe",self.cfg_name()]).communicate()
        forest_ee_file = self.output_name()

        
        subprocess.Popen(["bin/el8_amd64_gcc11/RegressionApplierExe",self.input_testing,self.applied_name(),"--gbrForestFileEE",forest_ee_file,"--gbrForestFileEB",forest_eb_file,"--nrThreads","4","--treeName",self.tree_name,"--writeFullTree",self.write_full_tree,"--regOutTag",self.reg_out_tag,"--meanLow",str(self.mean_min),"--meanHigh",str(self.mean_max)]).communicate()

        print("made ",self.applied_name())


    def forest_filenames(self):
        do_eb_org = self.do_eb
        self.do_eb = True
        forest_eb_file = self.output_name()
        self.do_eb = False
        forest_ee_file = self.output_name()
        self.do_eb = do_eb_org
        return forest_eb_file,forest_ee_file

import ROOT
import argparse
import numpy as np

class BDTTransformer:
    def __init__(self,range_min,range_max):
        self.range_min = range_min
        self.range_max = range_max
        self.offset = self.range_min + 0.5 * (self.range_max - self.range_min)
        self.scale = 0.5 * (self.range_max - self.range_min)
    
    def transform(self,raw_value):
        #features_ptr = features.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        return self.offset + self.scale * np.sin(raw_value)

class RegressionContainer:
    def __init__(self,base_file,range_mean_min,range_mean_max,range_sigma_min,range_sigma_max):
        
        
        self.root_file_eb = ROOT.TFile(base_file.format(region="EB"),"READ")
        self.root_file_ee = ROOT.TFile(base_file.format(region="EE"),"READ")
        print(self.root_file_eb.GetName())
        print(self.root_file_ee.GetName())

        self.forest_mean_eb = getattr(self.root_file_eb,f"EBCorrection")
        self.forest_mean_ee = getattr(self.root_file_ee,f"EECorrection")
        
        self.forest_sigma_eb = getattr(self.root_file_eb,f"EBUncertainty")
        self.forest_sigma_ee = getattr(self.root_file_ee,f"EEUncertainty")

        self.transformer_mean = BDTTransformer(range_mean_min,range_mean_max)
        self.transformer_sigma = BDTTransformer(range_sigma_min,range_sigma_max)
       
        
    def get_meansigma(self,features,isEB):
        forest_mean = self.forest_mean_eb if isEB else self.forest_mean_ee
        forest_sigma = self.forest_sigma_eb if isEB else self.forest_sigma_ee
        
        raw_mean = forest_mean.GetResponse(features.data())
        raw_sigma = forest_sigma.GetResponse(features.data()) 
        return self.transformer_mean.transform(raw_mean),self.transformer_sigma.transform(raw_sigma)

def safe_divide(numer,denom,default_val =0):
    if denom == 0:
        return default_val
    return numer/denom

def ptToP(pt,eta):
    return safe_divide(pt,np.sin(ROOT.MathFuncs.etaToTheta(eta)),0)
    

def get_best_trk_indx(calo_pt,calo_eta,trk_pts,trk_etas):
    """
    Returns the index of the track that is closest to the calo pt
    """
    best_trk_indx = 0
    best_abs_epm1 = 9999
    calo_energy = ptToP(calo_pt,calo_eta)
    for trk_indx in range(0,trk_pts.size()):
        trk_p = ptToP(trk_pts[trk_indx],trk_etas[trk_indx])
        abs_epm1 = abs(trk_p/calo_energy -1)
        if abs_epm1 < best_abs_epm1:
            best_abs_epm1 = abs_epm1
            best_trk_indx = trk_indx

    return best_trk_indx

def get_best_trk_indx_tree(tree,ele_nr):
    return get_best_trk_indx(tree.Electron_pt[ele_nr],tree.Electron_eta[ele_nr],tree.Electron_trkpt[ele_nr],tree.Electron_trketa[ele_nr])


def get_ietaiphi(seedid):
    detid = ROOT.DetId(seedid)
    if detid.subdetId() == 1:
        ebdetid = ROOT.EBDetId(detid)
        return ebdetid.ieta(),ebdetid.iphi()
    else:
        eedetid = ROOT.EEDetId(detid)
        return eedetid.ix(),eedetid.iy()

def isEB(seedid):
    detid = ROOT.DetId(seedid)
    return detid.subdetId() == 1

def get_features_calo(tree):
    features = []
    for ele_nr in range(0,tree.n_ele):
        ele_features = ROOT.std.vector('float')(11)
        
        ele_features[0] = tree.rho[0]
        ele_features[1] = ptToP(tree.Electron_pt[ele_nr],tree.Electron_eta[ele_nr])
        ele_features[2] = tree.Electron_eta[ele_nr]
        ele_features[3] = tree.Electron_phi[ele_nr]
        ele_features[4] = tree.Electron_sigmaietaieta[ele_nr]
        ele_features[5] = tree.Electron_r9[ele_nr]
        ele_features[6] = tree.Electron_smin[ele_nr]
        ele_features[7] = tree.Electron_smaj[ele_nr]
        ele_features[8] = tree.Electron_rechitzerosuppression[ele_nr]
        ele_features[9],ele_features[10] = get_ietaiphi(tree.Electron_seedid[ele_nr])
        features.append({"features" : ele_features,"isEB" : isEB(tree.Electron_seedid[ele_nr])})
    return features
        
def get_features_trk(tree):
    features = []
    for ele_nr in range(0,tree.n_ele):
        ele_features = ROOT.std.vector('float')(5)
        trk_indx = get_best_trk_indx(tree.Electron_pt[ele_nr],tree.Electron_eta[ele_nr],tree.Electron_trkpt[ele_nr],tree.Electron_trketa[ele_nr])
        ele_features[0] = tree.Electron_trkdz[ele_nr][trk_indx]
        ele_features[1] = tree.Electron_trkpt[ele_nr][trk_indx]
        ele_features[2] = tree.Electron_trketa[ele_nr][trk_indx]
        ele_features[3] = tree.Electron_trkphi[ele_nr][trk_indx]        
        ele_features[4] = tree.Electron_trkrchi2[ele_nr][trk_indx]
        features.append({"features" : ele_features,"isEB" : isEB(tree.Electron_seedid[ele_nr])})
    return features

def get_get_raw_comb(tree,ecal_meansigmas,trk_meansigmas):
    raw_comb = []
    for ele_nr in range(0,tree.n_ele):            
        trk_indx = get_best_trk_indx(tree.Electron_pt[ele_nr],tree.Electron_eta[ele_nr],tree.Electron_trkpt[ele_nr],tree.Electron_trketa[ele_nr])
        trk_p = ptToP(tree.Electron_trkpt[ele_nr][trk_indx],tree.Electron_trketa[ele_nr][trk_indx])
        calo_e = ptToP(tree.Electron_pt[ele_nr],tree.Electron_eta[ele_nr])
        calo_e_corr = calo_e * ecal_meansigmas[ele_nr][0]
        trk_p_corr = trk_p * trk_meansigmas[ele_nr][0]
        trk_p_err = trk_p * trk_meansigmas[ele_nr][1]
        calo_e_err = calo_e * ecal_meansigmas[ele_nr][1]
        #small bug it should be trp_p_corr**2
        comb = (calo_e_corr * trk_p_err**2 + trk_p *calo_e_err**2) / (calo_e_err**2 + trk_p_err**2)
        raw_comb.append(comb)
    return raw_comb

def get_features_comb(tree,ecal_meansigmas,trk_meansigmas):
    features = []
    for ele_nr in range(0,tree.n_ele):
        ele_features = ROOT.std.vector('float')(7)
        trk_indx = get_best_trk_indx(tree.Electron_pt[ele_nr],tree.Electron_eta[ele_nr],tree.Electron_trkpt[ele_nr],tree.Electron_trketa[ele_nr])

        ele_features[0] = tree.Electron_pt[ele_nr]*ecal_meansigmas[ele_nr][0]
        ele_features[1] = ecal_meansigmas[ele_nr][1]/ecal_meansigmas[ele_nr][0]
        ele_features[2] = trk_meansigmas[ele_nr][1]/trk_meansigmas[ele_nr][0]
        ele_features[3] = tree.Electron_trkd0[ele_nr][trk_indx]
        ele_features[4] = tree.Electron_pt[ele_nr]*ecal_meansigmas[ele_nr][0]/(tree.Electron_trkpt[ele_nr][trk_indx]*trk_meansigmas[ele_nr][0])
        ele_features[5] = tree.Electron_trkrchi2[ele_nr][trk_indx]
        ele_features[6] = tree.Electron_trkpt.size()
        features.append({"features" : ele_features,"isEB" : isEB(tree.Electron_seedid[ele_nr])})
    return features

if __name__ == "__main__":
    ROOT.gSystem.Load("./libs/$SCRAM_ARCH/libUtility.so")
    parser = argparse.ArgumentParser(description='Add in the regressed energy to the tree')
    parser.add_argument('-i', '--inputfile', help='Input file', required=True)
    parser.add_argument('-o', '--outputfile', help='Output file', required=True)
    args = parser.parse_args()

    input_file = ROOT.TFile(args.inputfile, "READ")

    input_tree = input_file.Get("mmtree/tree")

    output_file = ROOT.TFile(args.outputfile, "RECREATE")
    output_file.mkdir("mmtree")
    output_file.cd("mmtree")
    output_tree = input_tree.CloneTree(0)
    
    energy_corr = ROOT.std.vector('float')()
    trkp_corr = ROOT.std.vector('float')()
    energy_comb = ROOT.std.vector('float')()

    output_tree.Branch("Electron_ecorr", energy_corr)
    output_tree.Branch("Electron_trkpcorr", trkp_corr)
    output_tree.Branch("Electron_energycomb", energy_comb)
    
    ecal_reg = RegressionContainer("resultsEleTrkTrainOff/regEleEcalScout2023_stdVar_stdCuts_{region}_ntrees1500_results.root",0.2,2,0.0002,0.5)
    trk_reg = RegressionContainer("resultsEleTrkTrainOff/regEleTrkScout2023_stdVar_stdCuts_{region}_ntrees1500_results.root",0.2,2,0.0002,0.5)
    comb_reg = RegressionContainer("resultsEleTrkTrainOff/regEleEcalTrkScout2023_stdVar_stdCuts_{region}_ntrees1500_results.root",0.2,3,0.0002,0.5)

    n_entries = input_tree.GetEntries()
    for entrynr in range(0,n_entries):
        input_tree.GetEntry(entrynr)
        if entrynr % 1000 == 0:
            print(f"Processing entry {entrynr}/{n_entries}")
#        if entrynr > 100:
#            break
        energy_corr.clear()
        trkp_corr.clear()
        energy_comb.clear()
        energy_corr.resize(input_tree.n_ele)
        trkp_corr.resize(input_tree.n_ele)
        energy_comb.resize(input_tree.n_ele)        
        
        features_ecal = get_features_calo(input_tree)
        features_trk = get_features_trk(input_tree)
        

        ecal_meansigmas = [ecal_reg.get_meansigma(features["features"],features["isEB"]) for features in features_ecal]
        trk_meansigmas = [trk_reg.get_meansigma(features["features"],features["isEB"]) for features in features_trk]

        features_comb = get_features_comb(input_tree,ecal_meansigmas,trk_meansigmas)
        comb_meansigmas = [comb_reg.get_meansigma(features["features"],features["isEB"]) for features in features_comb]

        raw_comb = get_get_raw_comb(input_tree,ecal_meansigmas,trk_meansigmas)

        for elenr in range(0,input_tree.n_ele):
            trk_indx = get_best_trk_indx_tree(input_tree,elenr)
            trkp = ptToP(input_tree.Electron_trkpt[elenr][trk_indx],input_tree.Electron_trketa[elenr][trk_indx])
            energy = ptToP(input_tree.Electron_pt[elenr],input_tree.Electron_eta[elenr])

            energy_corr[elenr] = energy * ecal_meansigmas[elenr][0]
            trkp_corr[elenr] = trkp * trk_meansigmas[elenr][0]
            energy_comb[elenr] = raw_comb[elenr] * comb_meansigmas[elenr][0]
            
            
        
        output_tree.Fill()

    output_file.Write()


    
    

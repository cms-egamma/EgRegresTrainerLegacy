import ROOT
import argparse
import numpy as np
import itertools

class Cut:
    def __init__(self, func,cut_eb,cut_ee):
        self.func = func
        self.cut_eb = cut_eb
        self.cut_ee = cut_ee

    def __call__(self, tree, ele_nr):
        cut = self.cut_eb if abs(tree.Electron_eta[ele_nr]) < 1.5 else self.cut_ee
        return self.func(tree,ele_nr) < cut


selection = [
    Cut(lambda tree,ele_nr : tree.Electron_sigmaietaieta[ele_nr],0.0105,0.035),
    #Cut(lambda tree,ele_nr : tree.Electron_hoe[ele_nr],0.1,0.1),
    Cut(lambda tree,ele_nr : tree.Electron_hcaliso[ele_nr]/tree.Electron_pt[ele_nr],0.15,0.1),
    Cut(lambda tree,ele_nr : tree.Electron_ooemoop[ele_nr],0.015,0.01),
    Cut(lambda tree,ele_nr : tree.Electron_detain[ele_nr],0.005,0.01),
    Cut(lambda tree,ele_nr : tree.Electron_dphiin[ele_nr],0.03,0.05),
    Cut(lambda tree,ele_nr : tree.Electron_tkiso[ele_nr],5,5),
]

def get_eles_passing(tree):
    result = []
    for ele_nr in range(0,tree.n_ele):
        if all([cut(tree,ele_nr) for cut in selection]):
            result.append(ele_nr)
    return result

if __name__ == "__main__":
    ROOT.gSystem.Load("./libs/$SCRAM_ARCH/libUtility.so")
    parser = argparse.ArgumentParser(description='Add in the regressed energy to the tree')
    parser.add_argument('-i', '--inputfile', help='Input file', required=True)
    parser.add_argument('-o', '--outputfile', help='Output file', required=True)
    args = parser.parse_args()

    input_file = ROOT.TFile(args.inputfile, "READ")

    input_tree = input_file.Get("mmtree/tree")

    output_file = ROOT.TFile(args.outputfile, "RECREATE")
    mass_hist = ROOT.TH1D("massHist","massHist",150,0,150)
    comb_mass_hist = ROOT.TH1D("combMassHist","massHist",150,0,150)
    
    n_entries = input_tree.GetEntries()
    for entrynr in range(0,n_entries):
        input_tree.GetEntry(entrynr)
        if entrynr % 1000 == 0:
            print(f"Processing entry {entrynr}/{n_entries}")
#        if entrynr > 100:
#            break


        
        eles_passing = get_eles_passing(input_tree)
        if len(eles_passing) < 1:
            continue

        for ele1_nr,ele2_nr in itertools.combinations(eles_passing,2):
            ele1_p4 = ROOT.Math.LorentzVector("ROOT::Math::PtEtaPhiMVector")(input_tree.Electron_pt[ele1_nr],input_tree.Electron_eta[ele1_nr],input_tree.Electron_phi[ele1_nr],0.0)
            ele2_p4 = ROOT.Math.LorentzVector("ROOT::Math::PtEtaPhiMVector")(input_tree.Electron_pt[ele2_nr],input_tree.Electron_eta[ele2_nr],input_tree.Electron_phi[ele2_nr],0.0)
            
            mass = (ele1_p4+ele2_p4).M()

            ele1_p4_comb = ROOT.Math.LorentzVector("ROOT::Math::PtEtaPhiMVector")(input_tree.Electron_energycomb[ele1_nr]*input_tree.Electron_sintheta_trk[ele1_nr],input_tree.Electron_besttrketa[ele1_nr],input_tree.Electron_besttrkphi[ele1_nr],0.0)
            ele2_p4_comb = ROOT.Math.LorentzVector("ROOT::Math::PtEtaPhiMVector")(input_tree.Electron_energycomb[ele2_nr]*input_tree.Electron_sintheta_trk[ele2_nr],input_tree.Electron_besttrketa[ele2_nr],input_tree.Electron_besttrkphi[ele2_nr],0.0)

            mass_comb = (ele1_p4_comb+ele2_p4_comb).M()

            mass_hist.Fill(mass)
            comb_mass_hist.Fill(mass_comb)

    output_file.Write()


    
    

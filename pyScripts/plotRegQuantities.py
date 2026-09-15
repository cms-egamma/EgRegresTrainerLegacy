import ROOT
import argparse
import os
import re
import numpy as np
import regtools_scouting as regtools


def sanitise(expr):
    """turns a draw expression into something safe to use in histogram / file names"""
    return re.sub(r'\W+','_',expr).strip('_')


def get_auto_binning(tree,expr,cut="",frac=0.95,nr_bins=100,pad_frac=0.05,max_int_bins=1000,nr_entries=None):
    """
    returns (nr_bins,xmin,xmax) where the x-axis range contains the central frac of the
    values of expr passing cut, padded by pad_frac on each side
    integer valued quantities instead get unit width bins centred on the integers
    (unless that needs more than max_int_bins bins)
    returns None if no values pass the cut
    """
    if nr_entries is None: nr_entries = ROOT.TTree.kMaxEntries
    nr_sel = tree.Draw(expr,cut,"goff",nr_entries)
    if nr_sel<=0:
        return None
    if nr_sel>=tree.GetEstimate():
        tree.SetEstimate(nr_sel+1)
        nr_sel = tree.Draw(expr,cut,"goff",nr_entries)
    vals = np.frombuffer(tree.GetV1(),dtype=np.float64,count=nr_sel)
    xmin,xmax = np.quantile(vals,[(1.-frac)/2.,(1.+frac)/2.])
    if np.all(vals==np.floor(vals)) and xmax-xmin<max_int_bins:
        xmin,xmax = np.floor(xmin),np.ceil(xmax)
        return int(xmax-xmin)+1,xmin-0.5,xmax+0.5
    if xmax==xmin:
        return nr_bins,xmin-1.,xmax+1.
    pad = (xmax-xmin)*pad_frac
    return nr_bins,xmin-pad,xmax+pad


def plot_features(tree,features,cut,plot_dir,plot_name,frac=0.95,nr_entries=None):
    """
    plots each feature expression with an auto-ranged x-axis containing the central
    frac of the values passing cut, writing a png per feature and a multi-page pdf
    named plot_name.pdf to plot_dir
    """
    if nr_entries is None: nr_entries = ROOT.TTree.kMaxEntries
    os.makedirs(plot_dir,exist_ok=True)
    canvas = ROOT.TCanvas("plot_canvas","",900,700)
    pdf_name = os.path.join(plot_dir,plot_name+".pdf")
    canvas.Print(pdf_name+"[")
    hists = []
    for feature in features:
        binning = get_auto_binning(tree,feature,cut,frac=frac,nr_entries=nr_entries)
        if binning is None:
            print("no entries selected for {}, skipping".format(feature))
            continue
        nr_bins,xmin,xmax = binning
        hist_name = "{}_{}".format(plot_name,sanitise(feature))
        hist = ROOT.TH1D(hist_name,"{feature} {tag};{feature};# entries".format(feature=feature,tag=plot_name),nr_bins,xmin,xmax)
        tree.Draw("{}>>{}".format(feature,hist_name),cut,"goff",nr_entries)
        hist.SetLineWidth(2)
        hist.Draw("HIST")
        canvas.Print(os.path.join(plot_dir,hist_name+".png"))
        canvas.Print(pdf_name)
        hists.append(hist)
    canvas.Print(pdf_name+"]")


def main():
    parser = argparse.ArgumentParser(description='plots the features and target of the regression training')
    parser.add_argument('--input_file','-i',required=True,help='input ntuple with the regression tree')
    parser.add_argument('--output_dir','-o',default="plots",help='output dir for the plots')
    parser.add_argument('--frac',type=float,default=0.95,help='fraction of events the auto ranged x-axis contains')
    parser.add_argument('--nr_entries',type=int,default=None,help='number of tree entries to use (default all)')
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(111111) #with under/overflows as the auto range clips the tails

    regargs = regtools.RegArgs()

    tree = ROOT.TChain(regargs.tree_name)
    tree.Add(args.input_file)
    regargs.set_elecomb_default(use_existing_energy=False) #to get the correct feature names for the combined regression

    #region cuts mirror CutEB/CutEE in the trainer cfg
    #(regtools_scouting uses ele.isEB, regtools uses sc.isEB)
    regions = [("EB",regargs.var_eb,"ele.isEB"),
               ("EE",regargs.var_ee,"!ele.isEB")]
    for region_name,region_vars,region_cut in regions:
        features = list(dict.fromkeys(region_vars.split(":")+[regargs.target]))
        cut = "({}) && {}".format(regargs.cuts_base,region_cut)
        plot_features(tree,features,cut,args.output_dir,"regFeats_{}".format(region_name),
                      frac=args.frac,nr_entries=args.nr_entries)


if __name__ == "__main__":
    main()

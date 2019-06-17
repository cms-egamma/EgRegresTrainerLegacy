#!/usr/bin/env python

import argparse
import subprocess
    
if __name__ =='__main__':
    
    parser = argparse.ArgumentParser(description='makes a list of files to run over')
    parser.add_argument('input_file',help='input filename')
    parser.add_argument('output_file',help='output filename')
    parser.add_argument('--ideal',help='ideal training files with {region} instead of EB,EE')
    parser.add_argument('--real',help='real training files with {region} instead of EB,EE')
    parser.add_argument('--ecaltrk',help='ecaltrk training files with {region} instead of EB,EE')
    args = parser.parse_args()

    
    base_cmd = "./bin/slc6_amd64_gcc700/RegressionApplierExe {input_file} {output_file} --gbrForestFileEB {gbrEB} --gbrForestFileEE {gbrEE} --nrThreads 4 --writeFullTree 1 --regOutTag {reg_out_tag}"
        
    ecal_ideal_file = "ecalIdealTmp.root"
    ecal_real_file = "ecalRealTmp.root"

    #apply the ideal training
    
    ideal_args = {}
    ideal_args['input_file'] = args.input_file
    ideal_args['output_file'] = ecal_ideal_file
    ideal_args['gbrEB'] = args.ideal.format(region="EB")
    ideal_args['gbrEE'] = args.ideal.format(region="EE")
    ideal_args['reg_out_tag'] = "Ideal"
    cmd = base_cmd.format(**ideal_args)
    print cmd
    subprocess.Popen(cmd.split()).communicate()
    
    real_args = {}
    real_args['input_file'] = ecal_ideal_file
    real_args['output_file'] = ecal_real_file
    real_args['gbrEB'] = args.real.format(region="EB")
    real_args['gbrEE'] = args.real.format(region="EE")
    real_args['reg_out_tag'] = "Real"    
    cmd = base_cmd.format(**real_args)
    subprocess.Popen(cmd.split()).communicate()
    

    ecaltrk_args = {}
    ecaltrk_args['input_file'] = ecal_real_file
    ecaltrk_args['output_file'] = args.output_file
    ecaltrk_args['gbrEB'] = args.ecaltrk.format(region="EB")
    ecaltrk_args['gbrEE'] = args.ecaltrk.format(region="EE")
    ecaltrk_args['reg_out_tag'] = "EcalTrk"
    cmd = base_cmd.format(**ecaltrk_args)
    subprocess.Popen(cmd.split()).communicate()
    
    

import sys
from mhfp.encoder import MHFPEncoder
from rdkit import Chem
import numpy as np
import pickle
import argparse

def get_MACCS(mol,output,fpSize=2048):
    fp =  np.frombuffer(Chem.MACCSkeys.GenMACCSKeys(mol).ToBitString().encode(), 'u1') - ord('0')
    with open(output+"MACCS.bin","ab") as file:
        pickle.dump(fp,file)
    return 0

def get_Morgan(mol,output,radius=2, fpSize=1024):
    fp =  np.frombuffer(Chem.AllChem.GetMorganFingerprintAsBitVect(mol,radius,fpSize).ToBitString().encode(), 'u1') - ord('0')
    with open(output+f"ECFP{radius*2}.bin","ab") as file:
        pickle.dump(fp,file)
    return 0

def get_MHFP(mol, mhfp,output, fpSize = 2048,radius=3):
    fp =  mhfp.secfp_from_mol(mol,fpSize,radius)
    with open(output+f"MHFP{radius*2}.bin","ab") as file:
        pickle.dump(fp,file)
    return 0

parser = argparse.ArgumentParser(prog='FPython',
     description='''FingerPrint Calculator''',
     epilog=''' ''',
        formatter_class=argparse.MetavarTypeHelpFormatter)
requiredNamed = parser.add_argument_group('Required named arguments')
group1 = parser.add_argument_group('FP', 'FingerPrints to be calculated.')

requiredNamed.add_argument('-i', '--input',type=str, action='store', required=True,help='The path to the input file. Must be files with just smiles.')
parser.add_argument('-o', '--output',type=str, action='store',default="./FP_",help='The path to the output files. One .bin file will be created for each fingerprint.')
group1.add_argument('-e2',  default=False,action='store_true',help='ECFP2 FingerPrint Calculation')
group1.add_argument('-e4',  default=False,action='store_true',help='ECFP4 FingerPrint Calculation')
group1.add_argument('-e6',  default=False,action='store_true',help='ECFP6 FingerPrint Calculation')
group1.add_argument('-M',  default=False,action='store_true',help='MACCS FingerPrint Calculation')
group1.add_argument('-F',  default=False,action='store_true',help='MHFP FingerPrint Calculation')

if __name__=="__main__":
    if len(sys.argv) < 2:
        parser.print_usage()
        sys.exit(1)
    parsed_args = parser.parse_args()
    try:
        with open(parsed_args.input) as f:
            if parsed_args.F:
                mhfp_obj = MHFPEncoder()
            index = 0
            for line in f:
                try:
                    molecule = Chem.MolFromSmiles(line)
                except Exception as e:
                    sys.stderr.write(f"Something went really wrong. Exception: {e}")
                if molecule is None:
                    sys.stderr.write(f"Wrong molecule at line {index}.\n")
                if parsed_args.e2:
                    get_Morgan(molecule,parsed_args.output,radius=1)
                if parsed_args.e4:
                    get_Morgan(molecule,parsed_args.output,radius=2)
                if parsed_args.e6:
                    get_Morgan(molecule,parsed_args.output,radius=3)
                if parsed_args.M:
                    get_MACCS(molecule,parsed_args.output)
                if parsed_args.F:
                    get_MHFP(molecule,mhfp_obj,parsed_args.output)
                if index % 10000 == 0:
                    sys.stdout.write(f"Molecule {index}.\n")
                index += 1
    except FileNotFoundError:
        parser.print_help()
        sys.stderr.write("\n\nError Opening File.\n\n")
        sys.exit(1)

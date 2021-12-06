import numpy as np
import pandas as pd

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MACCSkeys, Descriptors, Descriptors3D, Draw, rdMolDescriptors, Draw, PandasTools
from rdkit.DataManip.Metric.rdMetricMatrixCalc import GetTanimotoSimMat, GetTanimotoDistMat
from rdkit import RDConfig

#Utilities 
import os 
from tqdm import tqdm



def load_data(file_path:str): 
    """Loads Data from SDF file

    Args:
        file_path (string): Path to file 
    """
    
    sdfFile = os.path.join(file_path)
    mols=Chem.SDMolSupplier(sdfFile)
    frame = PandasTools.LoadSDF(sdfFile,smilesName='SMILES',molColName='Molecule',
           includeFingerprints=True, removeHs=False, strictParsing=True)
    frame = calculate_descriptors(data_set=frame, mols=mols)
    return frame

def calculate_descriptors(data_set:pd.DataFrame, mols:pd.Series):
    """Calculate Descriptors using RDKit, needs a molecular descriptor class from RDKit 

    Args:
        data_set (pd.DataFrame): [description]

    Returns:
        [type]: [description]
    """
    
     
    for i,mol in tqdm(enumerate(mols)):
        Chem.SanitizeMol(mol)
        data_set.loc[i,'Molecule']=mol
        data_set.loc[i,'MolWt']=Descriptors.MolWt(mol)
        data_set.loc[i,'LogP']=Descriptors.MolLogP(mol)
        data_set.loc[i,'NumHAcceptors']=Descriptors.NumHAcceptors(mol)
        data_set.loc[i,'NumHDonors']=Descriptors.NumHDonors(mol)
        data_set.loc[i,'NumHeteroatoms']=Descriptors.NumHeteroatoms(mol)
        data_set.loc[i,'NumRotatableBonds']=Descriptors.NumRotatableBonds(mol)
        data_set.loc[i,'NumHeavyAtoms']=Descriptors.HeavyAtomCount (mol)
        data_set.loc[i,'NumAliphaticCarbocycles']=Descriptors.NumAliphaticCarbocycles(mol)
        data_set.loc[i,'NumAliphaticHeterocycles']=Descriptors.NumAliphaticHeterocycles(mol)
        data_set.loc[i,'NumAliphaticRings']=Descriptors.NumAliphaticRings(mol)
        data_set.loc[i,'NumAromaticCarbocycles']=Descriptors.NumAromaticCarbocycles(mol)
        data_set.loc[i,'NumAromaticHeterocycles']=Descriptors.NumAromaticHeterocycles(mol)
        data_set.loc[i,'NumAromaticRings']=Descriptors.NumAromaticRings(mol)
        data_set.loc[i,'RingCount']=Descriptors.RingCount(mol)
        data_set.loc[i,'FractionCSP3']=Descriptors.FractionCSP3(mol)

        data_set.loc[i,'TPSA']=Descriptors.TPSA(mol)
        data_set.loc[i,'NPR1']=rdMolDescriptors.CalcNPR1(mol)
        data_set.loc[i,'NPR2']=rdMolDescriptors.CalcNPR2(mol)
        data_set.loc[i,'InertialShapeFactor']=Descriptors3D.InertialShapeFactor(mol)
        data_set.loc[i,'RadiusOfGyration']=Descriptors3D.RadiusOfGyration(mol)
    return data_set

if __name__ == "__main__": 
    from utils import save_csv
    data = load_data("tests/data/Biodeg.sdf")
    save_csv(data, "converted_data.csv")
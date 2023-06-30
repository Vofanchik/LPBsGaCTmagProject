from rdkit import Chem, DataStructs
from rdkit.Chem import rdDepictor, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdMolDescriptors
from rdkit import DataStructs

import pubchempy


def getMolSvg(smile='NC(C(=O)O)CS'):
    mol = Chem.MolFromSmiles(smile)
    mc = Chem.Mol(mol.ToBinary())

    try:
        Chem.Kekulize(mc)
    except:
        mc = Chem.Mol(mol.ToBinary())

    drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('opacity:1.0', 'opacity:0.0')
    return svg


def mol_mass_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mass = rdMolDescriptors.CalcExactMolWt(mol)
    return mass


def iupac_from_smiles(smiles):
    return pubchempy.get_compounds(smiles, namespace='smiles')[0].iupac_name


def return_morganfp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    fp = AllChem.GetMorganFingerprint(mol, 2)
    return fp


def similiaryty_list_return(smile, list_of_id_name_smiles):
    targ = Chem.MolFromSmiles(smile)
    ms = [[i[0], i[1], Chem.MolFromSmiles(i[2])] for i in filter(lambda x: x[2] != smile, list_of_id_name_smiles)]
    fpgen = AllChem.GetRDKitFPGenerator()
    fps = [[x[0], x[1], fpgen.GetFingerprint(x[2])] for x in ms]
    tfp = fpgen.GetFingerprint(targ)
    r = [[i[0], i[1], DataStructs.TanimotoSimilarity(tfp, [i][0][2])] for i in fps]
    r.sort(key=lambda x: x[2], reverse=True)
    return r

def similiaryty_map_rerurn()


if __name__ == '__main__':
    m1 = return_morganfp('ClC1=COCNC1')
    m2 = return_morganfp('IC1=COCNC1')
    print(DataStructs.DiceSimilarity(m1, m2))

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdMolDescriptors
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
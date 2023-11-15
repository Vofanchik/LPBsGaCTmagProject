from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D, SimilarityMaps

import pubchempy
from rdkit.Chem.rdFingerprintGenerator import GetRDKitFPGenerator

from Classes.Testin import timeit

@timeit
def getMolSvg(smile='NC(C(=O)O)CS'):
    mol = Chem.MolFromSmiles(smile)
    mc = Chem.Mol(mol.ToBinary())
    print(mc)
    try:
        Chem.Kekulize(mc)
    except:
        pass

    drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('opacity:1.0', 'opacity:0.0')
    return svg

@timeit
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

@timeit
def similiaryty_list_return(smile, list_of_id_name_smiles):
    targ = Chem.MolFromSmiles(smile)
    ms = [[i[0], i[1], Chem.MolFromSmiles(i[2])] for i in filter(lambda x: x[2] != smile, list_of_id_name_smiles)]
    fpgen = GetRDKitFPGenerator()
    fps = [[x[0], x[1], fpgen.GetFingerprint(x[2])] for x in ms]
    tfp = fpgen.GetFingerprint(targ)
    r = [[i[0], i[1], DataStructs.TanimotoSimilarity(tfp, [i][0][2])] for i in fps]
    r.sort(key=lambda x: x[2], reverse=True)
    return r

@timeit
def similiaryty_map_rerurn(target, referent):
    targetmol = Chem.MolFromSmiles(target)
    refmol = Chem.MolFromSmiles(referent)

    d = Draw.MolDraw2DSVG(400, 400)
    d.ClearDrawing()
    target_mol_simi_fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(
        refmol,
        targetmol,
        lambda m, i: SimilarityMaps.GetMorganFingerprint(
            m, i, radius=2, fpType='bv'),
        draw2d=d,
    )
    d.FinishDrawing()
    return d.GetDrawingText().replace('opacity:1.0', 'opacity:0.0')

@timeit
def similiaryty_map_rerurn_png(target, referent):
    targetmol = Chem.MolFromSmiles(target)
    refmol = Chem.MolFromSmiles(referent)

    d = Draw.MolDraw2DCairo(400, 400)
    d.ClearDrawing()
    SimilarityMaps.GetSimilarityMapForFingerprint(
        refmol,
        targetmol,
        lambda m, i: SimilarityMaps.GetMorganFingerprint(
            m, i, radius=2, fpType='bv'),
        draw2d=d,
    )
    d.FinishDrawing()

    return d.GetDrawingText()


if __name__ == '__main__':
    m1 = return_morganfp('ClC1=COCNC1')
    m2 = return_morganfp('IC1=COCNC1')
    print(DataStructs.DiceSimilarity(m1, m2))

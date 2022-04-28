from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.Chem.rdchem import BondType

def bond_insensitive(mol):
    molrw = Chem.RWMol(mol)
    Chem.Kekulize(molrw, clearAromaticFlags=True)
    for b in molrw.GetBonds():
        b.SetBondType(BondType.SINGLE)
    return molrw

def getConformers(file_rec_lig_sdf, file_warhead_sdf, protac_smi, file_docked, file_out):
    rmsList = []
    #rdkit might meet some errors for some ligands
    e3_ligand = Chem.SDMolSupplier(file_rec_lig_sdf)[0]
    warhead = Chem.SDMolSupplier(file_warhead_sdf)[0]
    docked_head = Chem.SDMolSupplier(file_docked)[0]
    with open(protac_smi,'rb') as protac_smi_input:
        protac_smi = protac_smi_input.read().splitlines()[0]
    protac = Chem.MolFromSmiles(protac_smi)
    Chem.AddHs(protac)
    Chem.AddHs(docked_head)
    docked_head_bi = bond_insensitive(docked_head)
    protac_bi = bond_insensitive(protac)
    e3_ligand_bi = bond_insensitive(e3_ligand)
    warhead_bi = bond_insensitive(warhead)
    docked_e3 = docked_head_bi.GetSubstructMatch(e3_ligand_bi)
    docked_warhead = docked_head_bi.GetSubstructMatch(warhead_bi)
    protac_e3 = protac_bi.GetSubstructMatch(e3_ligand_bi)
    protac_warhead = protac_bi.GetSubstructMatch(warhead_bi)
    if len(docked_warhead) == 0 or len(docked_e3) == 0 or len(protac_e3) == 0 or len(protac_warhead) == 0:
        print("The smiles of PROTAC doesn't match the structures of ligands of target or receptor proteins.")
    protac_align_id = list(protac_e3)+list(protac_warhead)
    docked_align_id = list(docked_e3)+list(docked_warhead)
    if not (len(docked_e3) == 0 or len(docked_warhead) == 0):
        docked_pos = {protac_e3[j]: docked_head.GetConformer().GetAtomPosition(docked_e3[j]) for j in range(len(docked_e3))}
        warhead_pos = {protac_warhead[j]: docked_head.GetConformer().GetAtomPosition(docked_warhead[j]) for j in range(len(docked_warhead))}
        cmap = {}
        cmap.update(docked_pos)
        cmap.update(warhead_pos)
        cids = AllChem.EmbedMultipleConfs(protac, numConfs=100, coordMap=cmap, maxAttempts=1000, numThreads=1, ignoreSmoothingFailures=True)
        if len(cids) > 0:
            writer = Chem.SDWriter(file_out)
            for i in range(len(cids)):
                rms = rdMolAlign.AlignMol(protac, docked_head, prbCid=i,atomMap=zip(protac_align_id,docked_align_id))
                if rms < 0.5:
                    rmsList.append(rms)
                    writer.write(protac, confId=i)
    return len(rmsList)
    
if __name__ == "__main__":
    print(getConformers('rec_lig.sdf', 'target_lig.sdf', 'protac.smi', 'lig_1.sdf', 'protacs.sdf'))

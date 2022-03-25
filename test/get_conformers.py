import os,re
from pathlib2 import Path
from string import digits
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign

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
    docked_e3 = docked_head.GetSubstructMatch(e3_ligand)
    docked_warhead = docked_head.GetSubstructMatch(warhead)
    protac_e3 = protac.GetSubstructMatch(e3_ligand)
    protac_warhead = protac.GetSubstructMatch(warhead)
    if len(docked_warhead) == 0 or len(docked_e3) == 0 or len(protac_e3) == 0 or len(protac_warhead) == 0:
        print "The smiles of PROTAC doesn't match the structures of ligands of target or receptor proteins."
    print docked_e3
    print docked_warhead
    print protac_e3
    print protac_warhead
    protac_align_id = list(protac_e3)+list(protac_warhead)
    docked_align_id = list(docked_e3)+list(docked_warhead)
    if not (len(docked_e3) == 0 or len(docked_warhead) == 0):
        cmap = {protac_e3[j]: docked_head.GetConformer().GetAtomPosition(docked_e3[j]) for j in range(len(docked_e3))}
        cmap.update({protac_warhead[j]: docked_head.GetConformer().GetAtomPosition(docked_warhead[j]) for j in range(len(docked_warhead))})
        cids = AllChem.EmbedMultipleConfs(protac, numConfs=100, coordMap=cmap, maxAttempts=1000, numThreads=1)
        if len(cids) > 0:
            writer = Chem.SDWriter(file_out)
            for i in range(len(cids)):
                rms = rdMolAlign.AlignMol(protac, docked_head, prbCid=i,atomMap=zip(protac_align_id,docked_align_id))
                if rms < 0.5:
                    rmsList.append(rms)
                    writer.write(protac, confId=i)
    return len(rmsList)
    
if __name__ == "__main__":
    print(getConformers('rec_lig.sdf', 'target_lig.sdf', 'protac.smi', 'lig_1.sdf', 'protac.sdf'))
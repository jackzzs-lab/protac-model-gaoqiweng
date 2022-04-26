from pathlib2 import Path
import shutil

if __name__ == "__main__":
    
    with open("protacs.smi", "r") as f:
        for line in f.readlines():
            ligand_dir = Path(line.split()[1])
            shutil.rmtree(str(ligand_dir), ignore_errors=True)
            ligand_dir.mkdir()
            with open(str(ligand_dir / 'protac.smi'), "w+") as f:
                f.write(line.split()[0])
            for pdb in Path.cwd().glob('*.pdb'):
                if ligand_dir.name.startswith(pdb.stem):
                    (ligand_dir / 'target.pdb').symlink_to(pdb)
                elif ligand_dir.name.endswith(pdb.stem):
                    (ligand_dir / 'receptor.pdb').symlink_to(pdb)
            for sdf in Path.cwd().glob('*.sdf'):
                if ligand_dir.name.endswith(sdf.stem.split("_")[0]):
                    (ligand_dir / "_".join(sdf.name.split("_")[1:])).symlink_to(sdf)
            (ligand_dir / 'site_info.txt').symlink_to(Path.cwd() / 'site_info.txt')
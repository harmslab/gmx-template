import setup_gromacs_md as sgm

print("*** Protein alone ***")
sgm.setup_run(output_dir="example-outputs/prot-alone",
              ff_source_dir="ff/charmm36-feb2021.ff/",
              mdp_files_dir="ff/mdp-files/",
              protein_pdb_file="example-inputs/tlr4-md2.pdb")

print("*** Custom ligand alone, not posed ***")
sgm.setup_run(output_dir="example-outputs/lig-alone",
              ff_source_dir="ff/charmm36-feb2021.ff/",
              mdp_files_dir="ff/mdp-files/",
              charmm_gui_gromacs_dir="ff/e-coli_l-i-o-o1/")

print("*** Custom ligand alone, posed ***")
sgm.setup_run(output_dir="example-outputs/lig-alone_posed",
              ff_source_dir="ff/charmm36-feb2021.ff/",
              mdp_files_dir="ff/mdp-files/",
              charmm_gui_gromacs_dir="ff/e-coli_l-i-o-o1/",
              ligand_pose_files="example-inputs/e-coli_l-i-o-o1_pose1.pdb")

print("*** Protein with single pose of custom ligand ***")
sgm.setup_run(output_dir="example-outputs/prot-one-lig_posed",
              ff_source_dir="ff/charmm36-feb2021.ff/",
              mdp_files_dir="ff/mdp-files/",
              charmm_gui_gromacs_dir="ff/e-coli_l-i-o-o1/",
              protein_pdb_file="example-inputs/tlr4-md2.pdb",
              ligand_pose_files="example-inputs/e-coli_l-i-o-o1_pose1.pdb")

print("*** Protein with two poses of custom ligand ***")
sgm.setup_run(output_dir="example-outputs/prot-two-lig_posed",
              ff_source_dir="ff/charmm36-feb2021.ff/",
              mdp_files_dir="ff/mdp-files/",
              charmm_gui_gromacs_dir="ff/e-coli_l-i-o-o1/",
              protein_pdb_file="example-inputs/tlr4-md2.pdb",
              ligand_pose_files=["example-inputs/e-coli_l-i-o-o1_pose1.pdb",
                                 "example-inputs/e-coli_l-i-o-o1_pose2.pdb"])

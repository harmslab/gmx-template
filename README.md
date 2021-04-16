# gmx-template

Template for setting up MD simulations in the Harms lab. 

## Package contents

+ `example.py`: example python file showing different ways to set up md runs. 
+ `ff/`: directory with forcefield parameters. 
  + Models of various LPS and LPS-like molecules. Generated using the
    charmm-gui LPS modeler tool. These have names like `e-coli_l-i-o-o1`, which mean the LPS is taken from *E. coli* and has a lipid A (first `l`), inner core (`i`), outer core (`o`), and one O-antigen (`o1`). 
  + `charmm36-feb2021.ff`. The CHARMM36m forcefield, ported to GROMACS. 
    http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs
  + `mdp-files`. Files used for setting up a GROMACS run. 
+ `example-inputs/`: directories with inputs used to show different ways to run `example.py`. 
+ `example-outputs/`: directory where outputs are written from example runs of `example.py`. 

## Build initial model

### Prepping your protein pdb file

+ Delete all HETATM except metal ions at key sites (e.g. CA in an S100). 
  This includes deleting solvent. 
+ Renumber each chain so residues are different numbers in each chain.  
  Not all MD analysis code knows about chains, meaning if you have chain A
  and B, each with residue 1-100, it can get confusing.  I would make 
  chain A go from 1-100 and chain B go from 1001-1100.  You can do this in
  pymol (https://pymolwiki.org/index.php?title=Alter&redirect=no)
+ Make sure TER entries break each chain. 
+ If you are docking a ligand, you might have to move sidechains using the pymol "Mutagenesis" wizard. 
+ If you are combining several structures into one, pose them relative to one another.  I like the "Actions->Drag coordinates" feature in pymol. 

### Prepping LPS (or other custom ligands)

1. Generate MD parameters for LPS of interest by going to charmm-gui.org and using the "LPS Modeler Tool."  
2. Download the resulting forcefield parameters.
3. Convert the forcefield parameters to gromacs format by going to charmm-gui.org and using the "Force field converter" tool. Upload the `.psf`, `.crd`, and all `.str` files.  
4. Download the zipped file and uncompress it.  Find the  `gromacs` directory (which will
   have a bunch of files like `.gro`, `.itp`, and a `toppar` folder).  Note where 
   this ends up on your computer. 
5. Obtain the matched forcefield in gromacs format (such as `charmm36-feb2021.ff`, which is included in this repo)
6. In pymol, open the `.gro` file from the gromacs directory pose it where you want in the protein structure.  
7. Save out the posed LPS molecule as a .pdb file.  **MAKE SURE TO CHECK** "Original atom order" and "Retain atom ids" when writing out.   
8. If you are docking two or more LPS molecules, save out each pose as its
   own `.pdb` file.  

## Configure MD simulation 

### Setting up a run

1. Use `scp` to copy the following to talapas.  Make sure these end up in `/projects/harmslab/USERNAME`. 
   + The cleaned protein pdb file.
   + Any posed LPS pdb files.
   + The appropriate LPS gromacs directory 
   + The `gmx-template` directory
2. On talapas, modify the copy and modify the `gmx-template/example.py` script appropriately for your system. This should be run within the `gmx-template/` directory. 
3. Start an interactive session on talapas (`qsub -I -A harmslab`). 
4. Load gromacs:

```
module load intel/19
module load gromacs/2019.4
```

5. Run the `copied_example.py` script:

```
cd gmx-template
python copied_example.py
```

6. This script will generate an output directory in `gmx-template`.  Let's say you said `output_dir="my_output_dir"` in the `copied_example.py ` script.  This will make a new directory called `my_output_dir` with several sub directories.  Copy `my_output_dir/final` to some convenient location. I usually make a directory for each system I plan to simulate (`MD-2-with-LPS`, `MD-2-apo`, etc.) and then have replicate directories within each one (`MD-2-with-LPS/run000`, `MD2-with-LPS/run001`). The `final` directory you are copying would be `run000`, as this will be where the MD simulation actually runs.  If you plan to do several replicate runs, you can copy `final` to `run000`, `run001`, etc.  An independent simulation can then be run in each directory. 

### Starting a run

1. Go into your run directory (the copy of the `final` directory generated by `copied_example.py`, e.g. `run000`) and open the `run_md.srun` file using a text editor .  Make sure the `slurm` bits at the top are correct--especially these bits.

```
#SBATCH --partition=longfat
#SBATCH --time=07-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
```

 If you're not sure how to decide what these should be, talk to Mike. 

2. In the run directory, open the `md.mdp` file. Most of this file should not be changed, however, check the following line:

```
nsteps                  = 50000000  ; 2 * 50000000 = 1000 ps (100 ns)
```



This sets how long the simulation will run.  Assuming a `2 fs` step (`dt` line below in mdp file):

$$run\ time = 2 \times 10^{-15} s \times nsteps$$

This means the default (`50000000`) will run for 100 ns:

$$100\ ns = 1 \times 10^{-7}\ s = 2 \times 10^{-15} s \times 50,000,000\ steps$$

If you don't know how long you need to do your run for, talk to Mike. 

3. Start the run by typing:

```
qsub run_md.srun
```

4. Once the run starts, it will sequentially make the following log files: `em.log`, `nvt.log`, `npt.log`, and `md_0_1.log`.  You can follow the progress of the run by typing `tail -f BLAH.log`.  (You'll have to hit `CTRL+C` to exit the `tail` call).  The vast majority of the run will be occupied by the `md` run part of the script.  
5. The run will go until it either completes all steps or you run out of run-time, at which point `slurm` will kill your job. 

## After the run

### Extracting your run into useful compressed file

On talapas, navigate to the directory you ran the simulation in and type `qsub extract_md.srun md_0_1 npt.gro output`.  This will make a directory called `output` that has three files: `traj.gro`, `traj.tpr`, and `traj.xtc`.  These contain the main outputs from your run for further analysis. Copy this whole directory down to `spock` or your personal computer for analyses.  You can run `extract_md.srun` this while the run is in progress.  It will simply extract up to however far the run has progressed.  

### Extending a run

If you want to restart a stopped job--whether because it ran out of cluster time or because it finished all specified steps--edit the `run_md.srun` file.  There are a few lines commented at the bottom of the file that give instructions for how to do this. 

### Visualizing results

Unless your runs are short and/or your computer is powerful, you'll likely need look at the runs sitting at spock.  You'll want to use the program `vmd`. Navigate to wherever you download the simulation and type `vmd traj.gro traj.xtc`.  A vmd tutorial is here: http://www.ks.uiuc.edu/Training/Tutorials/vmd/tutorial-html/index.html
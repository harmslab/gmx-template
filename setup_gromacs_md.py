__description__ = \
"""
Small library for generating GROMACS simulation inputs for proteins, small
molecules with custom forcefields, or both together. 
"""
__author__ = "Michael J. Harms"
__date__ = "2021-04-16"

import re, glob, os, subprocess, shutil, random, string, sys

def _run_cmd(cmd_list,base=None,run_in_dir=None,input=None):
    """
    Run a bash command, doing error checking. 
    
    cmd_list: bash command, with each argument an individual entry in the list.
    base: base to write out logs.  If none, do not write. 
    run_in_dir: directory in which to run the command. 
    """

    # Go into specified working directory
    cwd = os.getcwd()
    if run_in_dir is not None:
        os.chdir(run_in_dir)
    
    # Run subprocess
    s = subprocess.run(cmd_list,capture_output=True,input=input)
    
    # Check for success
    if s.returncode != 0:
        err = ["Command failed.  Command was:\n\n"]
        err.append(" ".join(cmd_list))
        err.append("\n\n")
        err.append("standard error was:\n\n")

        stderr = s.stderr.decode("utf-8")
        
        err.append(stderr)
        
        os.chdir(cwd)
        
        raise RuntimeError("".join(err))
    
    # Write out logs
    if base:
        f = open(f"{base}.stdout","w")
        f.write(s.stdout.decode("utf-8"))
        f.close()
        
        f = open(f"{base}.stderr","w")
        f.write(s.stderr.decode("utf-8"))
        f.close()
    
    os.chdir(cwd)
    
def _generic_handler(line,num_for_key=None):
    """
    Handle a gromacs itp-style line.   
    
    num_for_key: number of fields to use for the key, starting from the first 
                 field.  This means an atom field (1,) or bond field (1,2), or
                 angle field (1,2,3), etc. can use the same parsing function
    """
    
    # Get rid of everything after comment
    line = line.split(";")[0].strip()
    
    # If line only has comments or is blank, return nothing
    if line == "":
        return None
    else:
        
        # If no num_for_key is specified, use the whole line (less the comment)
        # as a key
        if num_for_key is None:
            key = line
            
        # If num_for_key is specified, grab fields 0-num_for_key and use as key
        else:
            key = tuple(line.split()[0:num_for_key])
        value = line
        
        return key, value

def _moleculetype_handler(line):
    """moleculetype handler: field 1 as key"""
    return _generic_handler(line,1)
    
def _atoms_handler(line):
    """atoms handler: field 1 as key"""
    return _generic_handler(line,1)
    
def _bonds_handler(line):
    """atoms handler: fields 1-2 as key"""
    return _generic_handler(line,2)

def _pairs_handler(line):
    """pairs handler: fields 1-2 as key"""
    return _generic_handler(line,2)

def _angles_handler(line):
    """angles handler: fields 1-3 as key"""
    return _generic_handler(line,3)
    
def _dihedrals_handler(line):
    """dihedrals handler: fields 1-4 as key"""
    return _generic_handler(line,4)
    
def _load_itp(itp_file,filter_set=None):
    """
    Load an itp file into a dictionary or load an itp file and filter based on 
    values in filter_set. 
    
    itp_file: a gromacs itp file
    filter_set: dictionary of key (e.g. "atoms") and sub-keys (e.g. (1,3,4))
    
    If only itp_file is stored, return a dictionary with key/sub-keys 
    representing the itp file.
    
    If itp_file and filter_set are specified, return the text of the itp file
    filtered such that it only has the key/sub-key/value information from 
    filter_set. 
    """
    
    # Dictionary to hold itp data
    itp_dict = {}
    
    # Dictionary mapping keys from itp data to handler functions.
    handler_dict = {}
    handler_dict["moleculetype"] = _moleculetype_handler
    handler_dict["atoms"] = _atoms_handler
    handler_dict["atomtypes"] = _atoms_handler
    handler_dict["bonds"] = _bonds_handler
    handler_dict["bondtypes"] = _bonds_handler
    handler_dict["pairs"] = _pairs_handler
    handler_dict["pairtypes"] = _pairs_handler
    handler_dict["angles"] = _angles_handler
    handler_dict["angletypes"] = _angles_handler
    handler_dict["dihedrals"] = _dihedrals_handler
    handler_dict["dihedraltypes"] = _dihedrals_handler
    
    # handler key is current key (e.g. "angles"). current handler is the
    # handler function in use. 
    handler_key = None
    current_handler = None
    
    # Go through the line
    out_lines = []
    with open(itp_file) as f:
        for line in f:
            
            # If a new key
            if re.match("\[ .*? \]",line):
                
                # Get the handler_key and current_handler.  (The 
                # _generic_handler uses the whole line as a key and value). 
                handler_key = line.strip()[1:-1].strip()
                try:
                    current_handler = handler_dict[handler_key]
                except KeyError:
                    current_handler = _generic_handler
                
                # make sure itp_dict has a place for the key. 
                try:
                    itp_dict[handler_key]
                except KeyError:
                    itp_dict[handler_key] = {}

                if filter_set is not None:
                    out_lines.append(line)
                    
                continue
            
            # If we have a real handler...
            if current_handler is not None:
                
                # Get the value from the line.
                value = current_handler(line)
                
                # If we got a value, store in itp_dict
                if value is not None:
                    itp_dict[handler_key][value[0]] = value[1]
                    
                    # Store reverse key too (the angle (a,b,c) should be the same
                    # as (c,b,a)). 
                    if len(value[0]) > 1:
                        itp_dict[handler_key][tuple(value[0][::-1])] = value[1]
                    
                # If we have a filter set, decide what to do with the this line
                if filter_set is not None:
                    
                    # If we didn't get anything from the line, pass it along
                    # without changing it
                    if value is None:
                        out_lines.append(line)
                    else:
                        
                        # If we got a key/value pair and it is in the fitler_set,
                        # we should grab it. 
                        try:
                            filter_set[handler_key][value[0]]
                            out_lines.append(line)
                        except KeyError:
                            pass
            
            # If we don't have a real handler...
            else:
                
                # If a filter_set is specified, write the line without changing
                # it. 
                if filter_set is not None:
                    out_lines.append(line)

    # If we have a a filter_set, return the filtered lines rather than the itp_dict
    if filter_set is not None:
        return "".join(out_lines)
                    
    return itp_dict

def _load_itp_into_ff(new_itp,ff_dir,never_load=["defaults"]):
    """
    Take some new_itp file (say, for a ligand) and load this into a ff_dir, 
    stripping all duplicate definitions from the new_itp file.  This is useful 
    because something like charmm-gui will spit out an itp file for a ligand 
    that includes all atoms.  If we want to combine that ligand with a protein;
    however, we need to use the generic forcefield. The generic forcefield will
    include duplicate definitions of the atoms.  This script whacks those 
    duplicates from the new itp out. 
        
    WARNING: this function will modify the contents of the forcefield 
    directory.  Make sure this is a local copy you do not mind changing. 
    
    new_itp: itp file for the ligand.
    ff_dir: gromacs-style forcefield dir (generally downloaded from 
            http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs)
    never_load: fields that will never be loaded from the new_itp file.
    """
    
    # Load new itp in to an ff
    new_ff = _load_itp(new_itp)

    # Load itp files into dictionaries
    itp_list = glob.glob(f"{ff_dir}/*.itp")
    
    ff_itp = [_load_itp(itp) for itp in itp_list]

    # Make a set of all keys seen in itp files, except new one (stuff like
    # 'atomtypes', 'dihedraltypes', etc.)
    all_existing_key_types = []
    for ff in ff_itp:
        all_existing_key_types.extend(ff.keys())
    all_existing_key_types = set(all_existing_key_types)

    # Go over all existing key types
    for key_type in all_existing_key_types:

        # Get entries with key_type from new ff. If no entries, move along
        try:
            new_ff_entries = new_ff[key_type]
        except KeyError:
            continue
            
        # Make so certain key types will never load
        if key_type in never_load:
            new_ff[key_type] = {}
            continue

        # This will hold sub-keys that are in the new ff but not the old ones
        new_dict = {}

        # Go through all of the entries under the key in the new ff
        for k in new_ff_entries.keys():

            # See if this entry is found 
            entry_found = False
            for ff in ff_itp:

                try:
                    ff[key_type][k]
                    entry_found = True
                    break
                except KeyError:
                    pass

            # if the entry wasn't in any of the itp files, include it in the new
            # itp 
            if not entry_found:
                new_dict[k] = new_ff_entries[k]
                
        # Update the new_ff with a (likely) truncated set of entires
        new_ff[key_type] = new_dict

    # Now, reload the itp file, filtering by entries in new_ff
    new_file_contents = _load_itp(new_itp,new_ff)
    
    f = open(os.path.join(ff_dir,"custom-ligand.itp"),'w')
    f.write(new_file_contents)
    f.close()
    
    append_next = False
    out_lines = []
    with open(os.path.join(ff_dir,"forcefield.itp")) as f:
        for line in f:
            if append_next:
                out_lines.append("#include \"custom-ligand.itp\"\n")
                append_next = False
            
            if line.startswith("#include \"ffbonded.itp\""):
                append_next = True
            out_lines.append(line)
            
    f = open(os.path.join(ff_dir,"forcefield.itp"),"w")
    f.write("".join(out_lines))
    f.close()

    
def _load_pose(pose_file,tmp_gro,output_file):
    """
    load the coordinates from a pose pdb file into a gro file.  (This is useful
    because the pdb format truncates the atom names, but pymol only saves out
    pdb files...  The work flow is to load the .gro file into pymol, tweak the
    coordinates of the file, and then save them out as a .pdb file.)
    
    pose_file: pdb file with new atom coordinates.  these atoms must be in the 
               same order and have the same (first four character) names as 
               the gro file.  This is what is written out by pymol if you 
               load in from the .gro file and then export with the 'Original
               atom order' and 'Retain atom ids' options checked. 
    tmp_gro: gro file that will have data loaded into it.
    output_file: gro file for writing out.
    """
    
    # Read pdb file, grabbing atom names and coordinates in nm
    atoms = []
    with open(pose_file) as f:
        for line in f:
            if line[:6] in ["ATOM  ","HETATM"]:
                atom_name = line[12:16].strip()
                coord = [float(line[30+8*i:(38+8*i)])/10.0 for i in range(3)]
                atoms.append((atom_name,coord))
        
    # Go through gro file and replace atom coordinates with those from the pdb
    # file
    out_lines = []
    with open(tmp_gro) as f:
        for i, line in enumerate(f):
            
            # Skip two header lines
            if i < 2:
                out_lines.append(line)
                continue
            
            # Last line is not atom in gro file
            if i - 2 == len(atoms):
                out_lines.append(line)
                continue
            
            # Make sure atoms match
            atom = line[10:15].strip()
            if atom[:4] != atoms[i-2][0]:
                err = "apparent mismatch in atoms in gro and pdb file\n"
                err += f"   {atom[:4]} vs {atoms[i-2][0]}\n"
                raise ValueError(err)
            
            # Load in new coordinates
            front = line[:20]
            coord = "".join(["{:>8.3f}".format(c) for c in atoms[i-2][1]])
            back = line[44:]
            
            out_lines.append(f"{front}{coord}{back}")
            
    # Write to output file
    f = open(output_file,'w')
    f.write("".join(out_lines))
    f.close()

def _load_ligands_into_system(new_gro_files,
                              new_itp_files,
                              system_gro=None,
                              system_top=None):
    """
    new_gro_files: list of new gro files to load into the system gro file
    new_itp_files: list of itp files matched to gro files to load into the 
                   system top file
    system_gro: system gro file
    system_top: system topology file
    """
    
    # sanity check
    if len(new_gro_files) != len(new_itp_files):
        err = "gro files and itp files must be matched\n"
        raise ValueError(err)

    # Read current system gro file
    if system_gro is not None:
        
        f = open(system_gro,'r')
        gro_lines = f.readlines()
        f.close()
        out_lines = gro_lines[2:-1]
    else:
        gro_lines = ["gmx\n"]
        out_lines = []
        
    # Go through new gro files
    itp_to_import = {}
    molecule_counts = {}
    for i in range(len(new_gro_files)):
        
        # Read lines from new gro files
        f = open(new_gro_files[i])
        gro_file = f.readlines()[2:-1]
        f.close()
        
        # Figure out residue offset we need to add to new gro file
        this_start = int(gro_file[0][:5])
        if len(out_lines) == 0:
            last_resid = 0
        else:
            last_resid = int(out_lines[-1][:5])
        offset = last_resid - this_start + 1
        
        # Go through lines from new gro file and stick into system gro file
        for g in gro_file:
            v = int(g[:5]) + offset
            out_lines.append("{:5d}{}".format(v,g[5:]))
            
        # If we've already seen this itp file, increment counter
        if new_itp_files[i] in itp_to_import:
            molecule_counts[itp_to_import[new_itp_files[i]]] += 1
            
        # If we haven't seen this itp file, load it in to our counter
        else:
            itp_dict = _load_itp(new_itp_files[i])
            mtype = list(itp_dict["moleculetype"].keys())[0][0]
            try:
                molecule_counts[mtype] 
                err = "molecule with same name specified in two different itp files\n"
                raise ValueError(err)
            except KeyError:
                itp_to_import[new_itp_files[i]] = mtype
                molecule_counts[mtype] = 1
    
    # Construct itp inclusion lines
    itp_inclusion_string = []
    for itp in itp_to_import:
        local_name = os.path.split(itp)[-1]
        itp_inclusion_string.append(f"#include \"{local_name}\"")
    itp_inclusion_string = "\n".join(itp_inclusion_string)
                                    
    # Construct molecule counts data
    m_counts_string = []
    for m in molecule_counts:
        m_counts_string.append(f"{m}            {molecule_counts[m]}")
    m_counts_string.append("\n")
    m_counts_string = "\n".join(m_counts_string)
    
    # Construct finalized gro file
    final_gro_contents = [gro_lines[0]]
    final_gro_contents.append("{:d}\n".format(len(out_lines)))
    final_gro_contents.extend(out_lines)
    final_gro_contents.append(gro_lines[-1])
    
    if system_gro is None:
        system_gro = os.path.join(os.path.dirname(system_top),"system.gro")
        
    f = open(system_gro,'w')
    f.write("".join(final_gro_contents))
    f.close()
    
    # Construct finalized top file
    looking_for_new_molecule = False
    final_top = []
    with open(system_top) as f:
        for line in f:

            if line.startswith("[ system ]"):
                final_top.append(itp_inclusion_string)
                final_top.append("\n\n")
                
            if line.startswith("[ molecules ]"):
                looking_for_new_molecule = True
            
            if looking_for_new_molecule and line.strip() == "":
                final_top.append(m_counts_string)
                looking_for_new_molecule = False
                
            final_top.append(line)
                
    if looking_for_new_molecule:
        final_top.append(m_counts_string)
        
    f = open(system_top,'w')
    f.write("".join(final_top))
    f.close()
    
    
def setup_run(output_dir,
              ff_source_dir,
              mdp_files_dir,
              charmm_gui_gromacs_dir=None,
              protein_pdb_file=None,
              ligand_pose_files=None,
              overwrite=False,
              box_size=1.0,
              ion_conc=0.1):
    """
    Set up an MD run, possibly with a ligand that has a custom forcefield
    generated by charmm-gui. This will add a solvent box and ions to neutralize
    the system.

    output_dir: directory to write files to
    ff_source_dir: directory of forcefield in gromacs format (e.g. CHARMM36m.ff)
    mdp_files_dir: directory with mdp files to set up run
    charmm_gui_gromacs_dir: gromacs dir spit out by charmm-gui
    protein_pdb_file: protein pdb file, posed as desired. if not specified,
                      will build simulation without protein.
    ligand_pose_files: ligand pdb file (or list of files) with pose(s). 
    overwrite: overwrite output_dir (default is false)
    box_size: length of edge of box away from molecule in all directions (nm)
    ion_conc: ion concentration in M.  will add NA and CL ions to this conc
              in numbers to neutralize system. 
              
    WARNING: this does *not* check for forcefield consistency.  You could 
    generate the ligand ff on charmm-gui using CHARMM36m and then stick it into
    the CHARMM22 forcefield.  This script may work, without error, but could 
    lead to wacky outcomes in the simulations. 

    HACKS:
        * This excludes the position_restraints field from the ligand .itp
          file.  This was a single, screwed up, line for LPS as generated
          by charmm-gui, so we simply dropped.
        * Always uses TIP3P waters. 
    """

    # Raise error if insufficient input specified
    if charmm_gui_gromacs_dir is None and \
       protein_pdb_file is None and \
       ligand_pose_files is None:
        err = "You must specify at least one of charmm_gui_gromacs dir,\n"
        err += "protein_pdb_file, or ligand_pose_files\n"
        raise ValueError(err)
        
    if ligand_pose_files is not None and charmm_gui_gromacs_dir is None:
        err = "ligand_pose_files require charmm_gui_gromacs_dir is given\n"
        raise ValueError(err)
    
    # Random 10-letter string to append to temporary files
    tmp_base = "".join([random.choice(string.ascii_letters) for _ in range(10)])
    
    # Create output directory
    if os.path.isdir(output_dir):
        if overwrite:
            shutil.rmtree(output_dir)
        else:        
            err = f"'{output_dir}' already exists!\n"
            raise FileExistsError(err)
    os.mkdir(output_dir)
    
    # Create a source directory to store the raw ff and charmm-gui output
    src_dir = os.path.join(output_dir,"src")
    os.mkdir(src_dir)
    
    # Create working directory
    working_dir = os.path.join(output_dir,"working")
    os.mkdir(working_dir)
    
    # Create final direcotry
    final_dir = os.path.join(output_dir,"final")
    os.mkdir(final_dir)
    
    # Copy in the forcefield directory
    ff_source_dir_base = os.path.split(os.path.dirname(ff_source_dir))[-1]
    shutil.copytree(ff_source_dir,os.path.join(src_dir,ff_source_dir_base))
    shutil.copytree(ff_source_dir,os.path.join(working_dir,ff_source_dir_base))
    ff_gromacs_dir = os.path.join(working_dir,ff_source_dir_base)
    
    # Copy in the mdp files
    mdp_files_dir_base = os.path.split(os.path.dirname(mdp_files_dir))[-1]
    shutil.copytree(mdp_files_dir,os.path.join(src_dir,mdp_files_dir_base))
    mdp_files = glob.glob(os.path.join(mdp_files_dir,"*.*"))
    for m in mdp_files:
        shutil.copy(m,working_dir)

    # If custom gromacs dir specified
    if charmm_gui_gromacs_dir is not None:
        
        # Copy in the charmm-gui gromacs directory
        charmm_gui_gromacs_dir_base = os.path.split(os.path.dirname(charmm_gui_gromacs_dir))[-1]    
        shutil.copytree(charmm_gui_gromacs_dir,os.path.join(src_dir,charmm_gui_gromacs_dir_base))
        charmm_gui_gromacs_dir = os.path.join(src_dir,charmm_gui_gromacs_dir_base)    

        # Look for gro file for ligand
        gro_files = glob.glob(os.path.join(charmm_gui_gromacs_dir,"*.gro"))
        if len(gro_files) != 1:
            err = f"There should be a single .gro file in {charmm_gui_gromacs_dir}\n"
            raise ValueError(err)
        gro_file = gro_files[0]

        # Look for itp file for ligand
        itp_files = glob.glob(os.path.join(charmm_gui_gromacs_dir,"toppar","*.itp"))
        itp_files = [itp for itp in itp_files if os.path.split(itp)[-1] != "forcefield.itp"]
        if len(itp_files) != 1:
            err = f"There should be a single .itp file besides forcefield.itp "
            err += f"in {charmm_gui_gromacs_dir}/toppar/\n"
            raise ValueError(err)
        itp_file = itp_files[0]

        # Look for the topol file for the ligand
        topol_file = os.path.join(charmm_gui_gromacs_dir,"topol.top")

        # Look for the forcefield itp file
        ff_itp_file = os.path.join(charmm_gui_gromacs_dir,"toppar","forcefield.itp")

        # Look for unique names in gro_file
        f = open(gro_file)
        lines = f.readlines()[2:-1]
        f.close()
        unique_resid_in_gro = set([l[5:10] for l in lines])

        # Look for unique residue names in itp file
        itp_atoms = _load_itp(itp_file)["atoms"]
        unique_resid_in_itp = []
        for k in itp_atoms:
            unique_resid_in_itp.append(itp_atoms[k].split()[3])
        unique_resid_in_itp = set(unique_resid_in_itp)

        # Find stuff in gro file that does not match itp file.  Clean up.
        gro_mismatches = unique_resid_in_gro.difference(unique_resid_in_itp)
        itp_mismatches = unique_resid_in_itp.difference(unique_resid_in_gro)
        to_gro_name = {}
        for g in gro_mismatches:
            for p in itp_mismatches:
                if p[:len(g)] == g:

                    try:
                        to_gro_name[p]
                        print(f"itp residue {p} ambiguous. maps to {to_gro_name[p]} and {g}\n")
                        raise ValueError(err)
                    except KeyError:
                        to_gro_name[p] = g + " "*(len(p) - len(g))

        out_lines = []
        in_atoms = False
        in_pos_restraints = False
        with open(itp_file,'r') as f:
            for line in f:
                if line.startswith("[ atoms ]"):
                    in_atoms = True
                    out_lines.append(line)
                    continue

                if in_atoms:
                    if line.startswith("["):
                        in_atoms = False
                    else:
                        for p in to_gro_name:
                            line = re.sub(p,to_gro_name[p],line)

                # Remove posres
                if line.startswith("#ifdef POSRES"):
                    in_pos_restraints = True
                    continue

                if in_pos_restraints:
                    if line.startswith("#endif"):
                        in_pos_restraints = False
                    continue 

                out_lines.append(line)

        # Create temporary itp file
        tmp_itp = os.path.join(working_dir,f"{tmp_base}_ligand.itp")
        f = open(tmp_itp,'w')
        f.write("".join(out_lines))
        f.close()

        # Create temporary gro file
        tmp_gro = os.path.join(working_dir,f"{tmp_base}_ligand.gro")
        shutil.copy(gro_file,tmp_gro)

        # Create temporary ligand topology file
        tmp_top = os.path.join(working_dir,f"{tmp_base}_ligand-topol.top")
        shutil.copy(topol_file,tmp_top)

        # Create temporary ligand_ff.itp
        tmp_ff_itp = os.path.join(working_dir,f"{tmp_base}_ligand_ff.itp")
        shutil.copy(ff_itp_file,tmp_ff_itp)

        # Load the tmp forcefield itp file into the gromacs forcefield directory
        _load_itp_into_ff(tmp_ff_itp,ff_gromacs_dir,never_load=["defaults","position_restraints"])
        
    # If ligand poses are specified, load them into .gro files in the 
    # working directory
    if charmm_gui_gromacs_dir is not None:
        ligand_gro_files = []
        if ligand_pose_files is not None:

            if type(ligand_pose_files) is str:
                ligand_pose_files = [ligand_pose_files]

            for i, p in enumerate(ligand_pose_files):
                p = os.path.abspath(p)
                p1 = os.path.split(p)[-1]
                pose_file = os.path.join(working_dir,f"{tmp_base}_{p1}")
                shutil.copy(p,pose_file)

                output_file = os.path.join(working_dir,f"{tmp_base}_pose_{i}.gro")
                _load_pose(pose_file,tmp_gro,output_file)
                ligand_gro_files.append(output_file)

        # If no poses are specified, use the simple ligand.gro file from the initial topology
        else:
            output_file = os.path.join(working_dir,f"{tmp_base}_pose_0.gro")
            shutil.copy(tmp_gro,output_file)
            ligand_gro_files.append(output_file)

    print("Creating initial topology")
    sys.stdout.flush()
        
    # If there is a protein pdb file, run pdb2gmx
    if protein_pdb_file is not None:
        
        tmp_protein_pdb = os.path.join(working_dir,f"{tmp_base}_protein.pdb")
        shutil.copy(protein_pdb_file,tmp_protein_pdb)
        
        local_ff = ff_source_dir_base.split(".ff")[0].strip()

        _run_cmd(["gmx","pdb2gmx",
                  "-f",os.path.split(tmp_protein_pdb)[-1],
                  "-o","system.gro",
                  "-p","topol.top",
                  "-ff",local_ff,
                  "-water","tip3p",
                  "-cmap"],
                  base=f"{tmp_base}_pdb2gmx",
                  run_in_dir=working_dir)
        
        system_gro = os.path.join(working_dir,"system.gro")
        system_topol = os.path.join(working_dir,"topol.top")
    
    # Otherwise, create a starting topology file
    else:
        
        system_gro = None
        system_topol = os.path.join(working_dir,"topol.top")
        
        f = open(system_topol,"w")
        f.write(f"#include \"./{ff_source_dir_base}/forcefield.itp\"\n")
        f.write(f"#include \"./{ff_source_dir_base}/tip3p.itp\"\n")
        f.write("\n\n")
        f.write("#ifdef POSRES_WATER\n")
        f.write("[ position_restraints ]\n")
        f.write("1    1       1000       1000       1000\n")
        f.write("#endif\n\n")
        f.write(f"#include \"./{ff_source_dir_base}/ions.itp\"\n")
        f.write("\n\n")
        f.write("[ system ]\n; Name\ngmx\n\n[ molecules ]\n\n")
        f.close()
        
    if charmm_gui_gromacs_dir is not None:
        
        # Load ligands into the system
        output_itp = os.path.join(working_dir,"custom-lig_local.itp")
        shutil.copy(tmp_itp,output_itp)
        new_itp_files = [output_itp for _ in ligand_gro_files]
        _load_ligands_into_system(ligand_gro_files,
                                  new_itp_files,
                                  system_gro,
                                  system_topol)
    
    
    # When we reach this point, we should have system.gro, topol.top, the 
    # local forcefield, and a bunch of autogenerated .itp files that fully 
    # specify the system.  This is not solvated or charge neutralized. 

    print("Create box")
    sys.stdout.flush()
    _run_cmd(["gmx","editconf",
             "-f","system.gro",
             "-o","system.gro",
             "-c",
             "-d","{:.1f}".format(box_size),
             "-bt","cubic"],
            base=f"{tmp_base}_editconf",
                  run_in_dir=working_dir)
    
    print("Solvating")
    sys.stdout.flush()
    _run_cmd(["gmx","solvate",
             "-cp","system.gro",
             "-cs","spc216.gro",
             "-o","system.gro",
             "-p","topol.top"],
            base=f"{tmp_base}_solvate",
                  run_in_dir=working_dir)
    
    print("Ionizing")
    sys.stdout.flush()
    
    _run_cmd(["gmx","grompp",
             "-f","ions.mdp",
             "-c","system.gro",
             "-p","topol.top",
             "-o","ions.tpr"],
            base=f"{tmp_base}_grompp_ion",
                  run_in_dir=working_dir)
            
    _run_cmd(["gmx","genion",
             "-s","ions.tpr",
             "-o","system.gro",
             "-p","topol.top",
             "-pname","NA",
             "-nname","CL",
             "-neutral",
             "-conc","{:.3}".format(ion_conc)],
            base=f"{tmp_base}_genion",
            run_in_dir=working_dir,
            input=b"SOL")
            
        
    # Copy out to the final directory
    non_tmp = [c for c in os.listdir(working_dir) if not c.startswith(tmp_base)]
    non_tmp = [c for c in non_tmp if not c.startswith("#")]
    for c in non_tmp:        
        source = os.path.join(working_dir,c)
        target = os.path.join(final_dir,c)
        if os.path.isdir(source):
            shutil.copytree(source,target)
        else:
            shutil.copy(source,target)
    

'''
Functions in this file will serve the purpose to perform alignments so that conserved networks can be generated

'''

#from modeller import *

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, PDBIO, Select
from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq
import os



##First function could be revised to use salign instead of multiple instances of malign

def perform_structure_alignment(pdb_dir, same_chain='A',out_dir='aligned_pdbs'):
    #Optional structural alignment -- alignment would be required, but this method uses modeller in particular


    env = Environ()
    aln = Alignment(env)
    pdbs = os.listdir(pdb_dir) 
    pdbs.sort()  

    if type(same_chain) != list:  #provide ability for a list of chains to be given, but not required, by default will take chain A from everything
        same_chain = [same_chain]*len(pdbs)

 
    for (pdb, chain) in zip(pdbs, same_chain):
        print(pdb, chain)
        pdb_file = os.path.join(pdb_dir, pdb)
        m = Model(env, file=pdb_file, model_segment=('FIRST:'+chain, 'LAST:'+chain))
        aln.append_model(m, atom_files=os.path.join(pdb_dir, pdb), align_codes=pdb)

    aln.malign()
    aln.salign()
    aln.compare_structures()
    
    mdl = model(env)
    mdl.read(file=aln[pdbs[0]].atom_file, model_segment=aln[pdbs[0]].range)
    mdl.color(aln=aln)
    mdl.write(file=f"{out_dir}/{pdbs[0].split('.pdb')[0]}.aln.pdb")

    
    for pdb in pdbs[1:]:
        mdl2 = model(env)
        mdl2.read(file=aln[pdb].atom_file, model_segment=aln[pdb].range)
        sel = selection(mdl).only_atom_types('CA')
        sel.superpose(mdl2, aln)
        mdl2.write(file=f"{out_dir}/{pdb.split('.pdb')[0]}.aln.pdb")

def msa_with_modeller(alignment_file, combined_fasta):
    """
    Perform MSA using modeller

    Note: If your proteins are not closely related, it may be better to use a more sophisticated alignment method

    Paramenters:
    - alignment_file: Name of file to write to
    - combined_fasta: Collection of sequences with names of all proteins in fasta format

    Returns:
    None
    """

    log.verbose()

    #Initialize environment
    env = environ()
    env.io.atom_files_directory='./'
    
    #Perform alignment, using default modeller parameters
    aln = alignment(env, file=combined_fasta, alignment_format='FASTA')
    aln.salign(rr_file='$(LIB)/as1.sim.mat',  # Substitution matrix used
            output='',
            max_gap_length=20,
            gap_function=False,              # If False then align2d not done
                feature_weights=(1., 0., 0., 0., 0., 0.),
            gap_penalties_1d=(-100, 0),
            output_weights_file='saligni1d.mtx',
            similarity_flag=True)   # Ensuring that the dynamic programming
                                    # matrix is not scaled to a
                                    # difference matrix

    #Write final alignment
    aln.write(file=alignment_file, alignment_format='PIR')


def perform_structure_alignment(pdb_dir, same_chain='A',out_dir='aligned_pdbs'):
    """
    Option for a build-in structural alignment. Will align all pdbs to the first pdb when pdbs_dir is sorted and output corresponding trans/rot matrices

    Parameters:
    - pdb_dir: directory containing all pdbs
    - same_chain: str or list containing all chains of interest, in order of pdbs when sorted

    Returns:
    Dictionary containing rotation and translation matrices for each pdb
    """

    #Initialize environment
    env = Environ()

    #Sort pdbs by name
    pdbs = os.listdir(pdb_dir) 
    pdbs.sort()  

    #Initialize dictionary to contain rotation and translation matrices
    rotation_information = {'Rot': [], 'Trans': []}


    #provide ability for a list of chains to be given, but not required, by default will take chain A from everything
    if type(same_chain) != list:  
        same_chain = [same_chain]*len(pdbs)

    #Exclude pdbs[0], since that will be aligned to
    for i, ref_pdb in enumerate(pdbs[1:]):

        aln = Alignment(env)

        m = Model(env, file=os.path.join(pdb_dir,pdbs[0]), model_segment=('FIRST:'+same_chain[0], 'LAST:'+same_chain[0]))
        aln.append_model(m, atom_files=os.path.join(pdb_dir, pdbs[0]), align_codes=pdbs[0])

        #Select chain of interest
        chain = same_chain[i]
        pdb_file = os.path.join(pdb_dir, ref_pdb)

        m = Model(env, file=pdb_file, model_segment=('FIRST:'+chain, 'LAST:'+chain))
        aln.append_model(m, atom_files=os.path.join(pdb_dir, ref_pdb), align_codes=ref_pdb)
    
        aln.malign()
        aln.malign3d()
        aln.compare_structures()
        
        mdl = model(env)
        mdl.read(file=aln[pdbs[0]].atom_file, model_segment=aln[pdbs[0]].range)
        mdl.color(aln=aln)
        mdl.write(file=f"{out_dir}/{pdbs[0].split('.pdb')[0]}.aln.pdb")
    
        
        mdl2 = model(env)
        mdl2.read(file=aln[ref_pdb].atom_file, model_segment=aln[ref_pdb].range)
        sel = selection(mdl).only_atom_types('CA')
        s = sel.superpose(mdl2, aln)
        mdl2.write(file=f"{out_dir}/{ref_pdb.split('.pdb')[0]}.aln.pdb")
        
        #Append to dictionary
        rotation_information['Rot'].append(np.array(s.rotation))   
        rotation_information['Trans'].append(np.array(s.translation))   

    return(rotation_information)


class ChainAndNonProteinSelect(Select):
    def accept_chain(self, chain):
        return chain.id == 'A'
    
    def accept_residue(self, residue):
        # Accept residue if it belongs to chain A or if it's a non-protein residue
        if residue.parent.id == 'A' or not is_aa(residue):
            return True
        return False

def align_with_waters(pdb_dir, rotation_matrices, translation_vectors, out_dir='aligned_pdbs_with_water'):
    """
    Use trans/rot matrices from structural alignment to translate waters

    Parameters: 
    - pdb_dir: Directory of given pdbs of interest
    - rotation_matrices, translation_vectors: Outputted from perform_structural_alignment

    Returns:
    None
    """
    def apply_transformation(coordinates, rotation_matrix, translation_vector):
        """
        Internal function to transform coordinates

        Parameters:
        - coordinates: coordinates of one atom
        - rotation_matrix, translation_vector: Outputted from perform_structural_alignment

        Returns:
        Numpy array of transformed coordinates
        """
        transformed_coords = np.dot(coordinates, rotation_matrix.T) + translation_vector
        return transformed_coords

    parser = PDBParser()
    pdbs = os.listdir(pdb_dir)
    pdbs.sort()
    for i, pdb in enumerate(pdbs):
        structure = parser.get_structure(pdb.split('.')[0], os.path.join(pdb_dir, pdb))
        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        atom_coord = atom.get_coord()
                        if i > 0:
                            transformed_coord = apply_transformation(atom_coord, rotation_matrices[i-1], translation_vectors[i-1])
                            atom.set_coord(transformed_coord)
                        else:
                            atom.set_coord(atom_coord)

        io = PDBIO()
        io.set_structure(structure)
        io.save(f"{out_dir}/{pdb.split('.')[0]}_aligned.pdb",ChainAndNonProteinSelect())


def seq_similarity(seq1, seq2):
    """
    Compute sequence similarity between two sequences

    Method taken from https://www.kaggle.com/code/stpeteishii/cafa-5-calculate-sequence-similarity

    Parameters:
    -seq1, seq2: sequences to compare
    
    Returns:
    float
    """

    alignments = pairwise2.align.globalxx(seq1,seq2)
    best_alignment = alignments[0]
    aligned_seq1 = best_alignment[0]
    aligned_seq2 = best_alignment[1]
    num_matches = 0
    num_mismatches = 0
    for i in range(len(aligned_seq1)):
        if aligned_seq1[i] == aligned_seq2[i]:
            num_matches += 1
        else:
            num_mismatches += 1
    similarity = num_matches / (num_matches + num_mismatches)
    return similarity   


def pdb_to_fastas(pdb_file, fasta_out, name='STATE', custom_residues=None):
    """
    Convert pdb to fasta, includes nonstandard residue names

    Parameters:
    - pdb_file: input
    - fasta_out: output directory
    - name: name (.fa extension added automatically)

    Returns:
    None
    """
    os.makedirs(fasta_out, exist_ok=True)

    amino_acid_dict = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
    'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'HIP': 'H', 'HID': 'H', 'HIE': 'H', 'HISD': 'H',
    'HISE': 'H', 'HISP': 'H', 'AS4': 'D','ASH': 'D',
    'GL4': 'E', 'GLH': 'E', 'ARN': 'R', 'LYN': 'K',
    'CYX': 'C', 'CYM': 'C', 'CSP': 'C', 'SEP': 'S', 'ASX': 'D'
}
    if custom_residues is not None:
        for key, val in custom_residues.items():
            amino_acid_dict[key] = val
        
    sequence = []
    with open(pdb_file, 'r') as PDB:
        pdb_lines = PDB.readlines()

    for line in pdb_lines:
        if ('ATOM' in line) and ('CA' in line):
            aa = line.split()[3]
            sequence.append(amino_acid_dict[aa])


    with open(os.path.join(fasta_out,f"{name}.fa"), 'w') as FASTA:
        FASTA.write(f'>{name}\n')
        FASTA.write(''.join(sequence))


def parse_fasta(fasta):
    """
    Get sequence information from fasta

    Returns:
    Seq objects
    """
    records = SeqIO.parse(fasta, 'fasta')
    df = pd.DataFrame(columns=['sequence', 'EntryID'])
    for record in records:
        df_new_row = pd.DataFrame({'sequence': str(record.seq), 'EntryID': record.id}, index=[0])
        df = pd.concat([df, df_new_row])

    seq1 = Seq(df.iloc[0,0])
    seq2 = Seq(df.iloc[1,0])
    return seq1, seq2

def generate_msa_alignment(alignment_file, combined_fasta, fasta_individual):
    """
    Renumber residue indices based on MSA alignment

    Parameters:
    - alignment_file: name of MSA alignment file -- will generate if the file does not exist
    - combined_fasta: Fasta file with all sequences
    - fasta_cur: Particular fasta file with sequence of interest

    Returns:
    List of msa_indices
    """
    #Make alignment file with modeller if name does not exist
    if not os.path.exists(alignment_file):
        msa_with_modeller(alignment_file, combined_fasta)

    #Open alignment and fasta files
    with open(alignment_file, 'r') as ALIGN:
        align_data = ALIGN.readlines()
    with open(fasta_individual, 'r') as FASTA:
        fasta_cur = FASTA.readlines()

    #Name of sequence
    names = [f_ind for (f_ind, f) in enumerate(align_data) if fasta_cur[0].replace('>','').replace('\n','') in f]
    end_index = [f_ind for (f_ind, f) in enumerate(align_data) if (f_ind > names[0] and '*' in f)][0]

    #Isolate sequence
    sequence = ''.join(align_data[names[0]+2:end_index+1]).replace('*','').replace('\n','')

    #Initiate empty list for indices
    msa_indices = []

    #Renumber based on alignment
    for val_cur, letter in enumerate(sequence):
        if letter != '-':
            msa_indices.append(val_cur+1)
            
    return msa_indices


def convert_msa_to_individual(msa_indices, msa_indices_ref, resids, resid_sequence_ref, resid_individual_ref):
    """
    Get residue indices from MSA using known reference

    msa_indices: List of MSA indices for particular sequence
    msa_indices_ref: List of MSA indices for reference sequence
    resids: List of residues for particular sequence
    resid_sequence_ref: List of residues for reference
    resid_individual_ref: Reference residue for one structure

    Returns: (Int) Desired resid for paricular structure
    """

    # Find the index in the reference MSA
    reference_msa = msa_indices_ref[np.where(np.array([int(f) for f in resid_sequence_ref]) == int(resid_individual_ref))[0][0]]

    # Find the corresponding index in the individual MSA
    #TEMPORARY
    # Find the closest index in msa_indices to reference_msa
    closest_index = np.argmin(np.abs(msa_indices - np.array(reference_msa)))

    try:
        # Use the closest index to find the desired residue ID
        desired_resid = int(resids[closest_index])
    except:
        print('Closest index', closest_index)
        print('Reference MSA', reference_msa)
        print('Calculated MSA', msa_indices[closest_index])
        print(len(resids), len(msa_indices))
        desired_resid=None
        print('Failed')

    #individual_index = np.where(msa_indices == reference_msa)[0][0]
    #print(individual_index)

    # Extract the desired residue ID
    #desired_resid = int(resids[individual_index[0]]) 

    return desired_resid


if __name__ == '__main__':
    fasta = 'fasta/all_seqs.fa'
    seq1, seq2 = parse_fasta(fasta) 
    similarity = seq_similarity(seq1, seq2)
    rotation_information = perform_structure_alignment('PTP_crystals_example/clean_pdbs')

    align_with_waters('PTP_crystals_example/clean_pdbs', rotation_information['Rot'], rotation_information['Trans'], out_dir='PTP_crystals_example/aligned_with_waters')

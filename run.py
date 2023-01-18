import argparse
import sys
import os
import glob
import logging
import time

from config import *
sys.path.append(lib_dir)
sys.path.append(scripts_dir)
from my_log import *
from validators import *
from utils import *
from prepare_loops import *
from alignment_generator import *
from superimposition_generator import *
from partial_pdb_generator import *
from image_helper import *
from validators import *

from relational_graph_generator import *

from cif import *

import collections

def get_residue_count(d):
    return sum(map(lambda x: len(d[x]), d))

def main():

    process_start_time = time.time()
    parser = argparse.ArgumentParser(description='Prepare input for RNAMotifContrast')
    parser.add_argument('-i', required=True, help='Input file containing motifs')
    parser.add_argument('-t', nargs='?', default='ScanX', const='ScanX', help='Input file containing motifs')

    try:
        args = parser.parse_args()
    except Exception as e:
        parser.print_help()
        sys.exit()

    user_input_fname = args.i
    
    global alignment_tool
    alignment_tool = args.t

    if alignment_tool.lower() == 'scanx':
        alignment_tool = 'ScanX'
    elif alignment_tool.lower() == 'tmalign' or alignment_tool.lower() == 'rnaalign':
        alignment_tool = 'TMalign'
    else:
        logger.warning('Invalid alignment tool selected. Continuing with ScanX.')
        alignment_tool = 'ScanX'

    
    loop_cif_extension = 0
    partial_pdbx_dir = os.path.join(data_dir, 'pdbx_extracted_ext' + str(loop_cif_extension))
    graphs_and_pickles_dir = os.path.join(data_dir, 'graphs_and_pickles')

    global alignment_dir_auto

    alignment_dir = os.path.join(alignment_dir_auto, alignment_tool + '_alignments')
    
    input_fname = os.path.join(data_dir, user_input_fname)
    input_fname_base = os.path.basename(input_fname)

    output_dir = os.path.join(root_dir, 'output', 'RNAalign' if alignment_tool == 'TMalign' else alignment_tool)

    prepare_executables()
    validate_all(input_fname)
    # create_required_directories(alignment_dir, graphs_and_pickles_dir)
    create_directory(pdbx_dir)
    create_directory(fasta_dir)
    create_directory(pdb_fasta_mapping_dir)
    create_directory(annotation_dir)
    create_directory(loop_dir)
    create_directory(alignment_dir)
    create_directory(graphs_and_pickles_dir)
    create_directory(output_dir)

    print('')
    logger.info('Reading input from ' + input_fname[base_path_len:])
    print('')

    families = {}
    fp_input = open(input_fname)
    loop_list = csv_to_list(fp_input.readlines())
    fp_input.close()

    loop_count = 0
    for item in loop_list:
        if len(item) > 2:
            # families[item[0]] = map(lambda x: str(strToNode(x)), item[1:]) # item[1:]
            families[item[0]] = item[1:]

    prepare_data(families)
    if input_index_type == 'pdb':
        families = convert_a_cluster_from_PDB_to_FASTA(families)
        for family_id in families:
            families[family_id] = list(map(lambda x: str(strToNode(x)), families[family_id]))
    
    loop_count = 0
    loop_node_list_str = []
    for family_id in families:
        loops = families[family_id]
        loop_count += len(loops)
        for loop in loops:
            loop = str(strToNode(loop))
            loop_node_list_str.append(loop)
    
    duplicates = [item for item, count in collections.Counter(loop_node_list_str).items() if count > 1]
    if len(duplicates) > 0:
        print('duplicates:')
        print(duplicates)

    loop_node_list_str = sorted(list(set(loop_node_list_str)))

    logger.info(str(loop_count) + ' loops (' + str(len(loop_node_list_str)) + ' unique) found in ' + str(len(families)) + ' famil' + ('ies' if len(families) > 1 else 'y') + '.')
    print('')

    # get pdb and fasta files 
    # and
    # generate loop.cif files
    # and annotation files
    
    previous_graph_file_reused = True
    all_alignment_files_found = False
    
    # prepare_data(families)
    prepare_loop_files(loop_node_list_str)    #chkd
    prepare_partial_pdbs(partial_pdbx_dir, loop_node_list_str, loop_cif_extension)

    # filtering for HETATM to use TM-Align, might not work if loop is extended
    families, filtered_count = filter_loops_containing_HETATM(families, partial_pdbx_dir)
    logger.info(str(filtered_count) + ' loops are filtered due to having HETATM.')
    
    alignment_data = None
    rmsd_data_dict = None
    
    attempt = 1
    while attempt <= 5:
        attempt += 1
        # generate_best_alignment_data function generates .graph file
        previous_graph_file_reused, all_alignment_files_found = generate_best_alignment_data(alignment_tool, input_fname_base, graphs_and_pickles_dir, alignment_dir, families)
        if all_alignment_files_found == False:
            get_alignment_files(alignment_tool, alignment_dir, families, loop_node_list_str, False)
            continue
        # family-wise
        alignment_data, rmsd_data_dict = load_alignment_and_rmsd_data(families, loop_node_list_str, input_fname_base, alignment_tool, partial_pdbx_dir, alignment_dir, graphs_and_pickles_dir, previous_graph_file_reused)
        if alignment_data == None:
            delete_graph_file(input_fname_base, graphs_and_pickles_dir) # as alignments need to be generated, it seems existing graph file is invalid
            get_alignment_files(alignment_tool, alignment_dir, families, loop_node_list_str, False)
            continue
        break

    # print('Check:')
    # print(alignment_data['E-loop'][strToNode('4V88_A6:626-630_967-971')][strToNode('1MFQ_A:78-82_93-97')])
    # print(alignment_data['E-loop'][strToNode('1MFQ_A:78-82_93-97')][strToNode('4V88_A6:626-630_967-971')])
    # sys.exit()
    # print('Checking alignment_data')
    # print(len(alignment_data['E-loop'][strToNode('4V88_A6:626-630_967-971')]))
    # print(len(alignment_data['Tandem-shear'][strToNode('5J7L_AA:1415-1418_1480-1483')]))

    # print('Checking rmsd_data')
    # ky = list(rmsd_data_dict['E-loop'][1].keys())[0]
    # print(ky)
    # print(len(rmsd_data_dict['E-loop'][1][ky][1]))
    # ky = list(rmsd_data_dict['Tandem-shear'][1].keys())[0]
    # print(ky)
    # print(len(rmsd_data_dict['Tandem-shear'][1][ky][1]))

    # sys.exit()
    generate_relative_graph_among_motif_families(families, loop_node_list_str, alignment_data, rmsd_data_dict, user_input_fname, alignment_tool, alignment_dir, graphs_and_pickles_dir, partial_pdbx_dir, output_dir)

    logger.info('\nTotal time taken: ' + str(round((time.time() - process_start_time), 3)) + ' seconds.\n')
        
if __name__ == '__main__':
    main()

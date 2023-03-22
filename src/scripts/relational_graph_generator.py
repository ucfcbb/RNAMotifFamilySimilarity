import sys
import os
import glob
import logging
import operator
import time
import pickle
import multiprocessing as mp
import numpy as np
from datetime import datetime
from Bio.PDB import *
import copy
import gc
import collections

import networkx as nx
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont

# python 3 compatibility
from functools import reduce
from past.builtins import map

sys.path.append('../../')
from config import *
sys.path.append(scripts_dir)
from utils import *
from pymol_helper import *
from classes import *
from validators import *
# from image_helper import *

from superimposition_generator import *

def load_inter_cluster_alignment_data(clusters, corresponding_loop_data):
    cluster_alignment_data = {}

    for c_id in clusters:
        cluster_alignment_data[c_id] = {}

        for i in range(len(clusters[c_id])):

            node1 = strToNode(clusters[c_id][i])

            if node1 not in cluster_alignment_data[c_id]:
                cluster_alignment_data[c_id][node1] = {}

            # for j in range(len(clusters[c_id])):
            #     if i == j:
            #         continue
            for j in range(len(corresponding_loop_data[c_id])):

                # node2 = strToNode(clusters[c_id][j])
                node2 = strToNode(corresponding_loop_data[c_id][j])

                if node2 not in cluster_alignment_data[c_id]:
                    cluster_alignment_data[c_id][node2] = {}

                # if node1 not in alignment_data or node2 not in alignment_data[node1]:
                #     logger.error('ERROR: (' + str(node1) + ', ' + str(node2) + ') pair not found in graph file!')
                #     sys.exit()

                # cluster_alignment_data[c_id][node1][node2] = alignment_data[node1][node2]
                # cluster_alignment_data[c_id][node2][node1] = alignment_data[node2][node1]

                cluster_alignment_data[c_id][node1][node2] = ('', '', 0., '', '', '', '', 0.)
                cluster_alignment_data[c_id][node2][node1] = ('', '', 0., '', '', '', '', 0.)

    return cluster_alignment_data, get_loops_in_cluster(clusters)

# def load_inter_family_alignment_data(input_fname_base, alignment_tool, alignment_dir, graphs_and_pickles_dir, alignment_fname, clusters, corresponding_loop_data, previous_graph_file_reused):
def load_inter_family_alignment_data(input_fname_base, alignment_tool, alignment_dir, graphs_and_pickles_dir, clusters, corresponding_loop_data):
    alignment_data_fname = os.path.join(graphs_and_pickles_dir, 'inter_family_alignment_data_' + input_fname_base + '_' + alignment_tool + '.pickle2')
    if (sys.version_info >= (3, 0)):
        alignment_data_fname = os.path.join(graphs_and_pickles_dir, 'inter_family_alignment_data_' + input_fname_base + '_' + alignment_tool + '.pickle3')

    # if previous_graph_file_reused == True:
        # if is_valid_pickle(alignment_data_fname, clusters) == True:
        # if True:
    if os.path.exists(alignment_data_fname):
        logger.info('Loading saved alignment data from previous run (In ' + alignment_data_fname[base_path_len:] + ') ...\n')
        f = open(alignment_data_fname, 'rb')
        cluster_alignment_data = pickle.load(f)
        f.close()
        return cluster_alignment_data

    logger.info('Loading inter-family alignment data from alignment files.')
    start_time = time.time()

    # alignment_data = load_graph_data(alignment_fname)

    # inter_cluster_alignment_data, loops_node_list = load_inter_cluster_alignment_data(alignment_data, clusters, corresponding_loop_data)
    inter_cluster_alignment_data, loops_node_list = load_inter_cluster_alignment_data(clusters, corresponding_loop_data)

    # print(inter_family_alignment_data)
    # sys.exit()


    # cluster_alignment_data, loops_node_list = load_inter_cluster_alignment_data(alignment_data, clusters)

    # print('Number of loops in cluster file: ' + str(len(loops_in_cluster)))

    file_counter = 0
    for node in loops_node_list:
        for r1 in get_all_loop_combination(str(node)):
            node1 = strToNode(r1)
    # for fn in glob.glob(os.path.join(alignment_dir, '*.aln')):
            fn = os.path.join(alignment_dir, r1 + '.aln')

            if not os.path.isfile(fn):
                return None

            # stime = datetime.now()
            file_counter += 1
            print('Processsing ' + fn[base_path_len:] + ' ... (' + str(file_counter) + ')')
            # sys.stdout.flush()
            # r1 = os.path.basename(fn)[:-4]
            # node1 = strToNode(r1)

            # if node1 not in loops_in_cluster:
            #     continue

            cid_nodelist_pair = find_nodes_in_cluster(node1, inter_cluster_alignment_data)

            fp = open(fn)
            lines = fp.readlines()
            fp.close()
            # flag = 0
            # print('No. of lines: ' + str(len(lines)))
            # print('Connected nodes: ' + str(len(node_dict)))
            line_index = 0
            while line_index < len(lines):
                # print 'Reading line ' + str(line_index)
                # sys.exit()
                if lines[line_index].startswith('#  Aligning'):
                    test_r1 = lines[line_index].split('::')[1].split(' and ')[0].strip().strip(':')
                    if r1 != test_r1:
                        logger.error('ERROR: filename and loop mismatch. r1(filename): ' + r1 + ', r1(loop): ' + test_r1)
                        sys.exit()

                    r2 = lines[line_index].split('::')[1].split(' and ')[1].strip().strip(':')
                    node2 = strToNode(r2)

                    for c_id, node_dict in cid_nodelist_pair:
                        if node2 in node_dict:

                            # t1 = inter_cluster_alignment_data[c_id][node1][node2][0]
                            # t2 = inter_cluster_alignment_data[c_id][node1][node2][1]
                            # # make sure the best alignment ordering is same as the ordering of loop in current file
                            # if not ((r1 == t1 and r2 == t2) or (r1 == t2 and r2 == t1)):
                            #     # line_index += 12
                            #     continue

                            # #safety condition for reverse order (might be redundant)
                            # if (r1 == t2 and r2 == t1):
                            #     t1 = r2
                            #     t2 = r1

                            score_text = lines[line_index+1].split(':')[1].strip()
                            if score_text == '':
                                score = -50.
                                logger.error('ERROR: No alignment score found for: ' + r1 + ' and ' + r2)
                                sys.exit()
                            else:
                                score = float(score_text)

                            text = lines[line_index+3].split(':')[1].strip()
                            cr1 = get_local_alignment_index(r1, text)

                            text = lines[line_index+4].split(':')[1].strip()
                            cr2 = get_local_alignment_index(r2, text)

                            aln1 = lines[line_index+6].strip()
                            aln2 = lines[line_index+7].strip()

                            if len(aln1) == 0 or len(aln2) == 0:
                                # set dummy alignment
                                aln1 = 'A'
                                aln2 = 'A'
                                temp_cr1 = cr1.split(':')
                                dummy_index = temp_cr1[1].split('_')[0].split('-')[0]
                                cr1 = temp_cr1[0] + ':' + dummy_index + '-' + dummy_index
                                temp_cr2 = cr2.split(':')
                                dummy_index = temp_cr2[1].split('_')[0].split('-')[0]
                                cr2 = temp_cr2[0] + ':' + dummy_index + '-' + dummy_index

                            zscore = inter_cluster_alignment_data[c_id][node1][node2][2]

                            prev_score = inter_cluster_alignment_data[c_id][node1][node2][7]
                            if score > prev_score:
                                inter_cluster_alignment_data[c_id][node1][node2] = (r1, r2, zscore, cr1, cr2, aln1, aln2, score)
                                inter_cluster_alignment_data[c_id][node2][node1] = (r2, r1, zscore, cr2, cr1, aln2, aln1, score)
                            # break
                    # we can skip at least 11 lines from current alignment data
                    line_index += 11
                line_index += 1
            # while end
            
            # etime = datetime.now()
            # difftime = (etime - stime).total_seconds()
            # print('Time taken: ')
            # print(difftime)
            # break

    f = open(alignment_data_fname,"wb")
    pickle.dump(inter_cluster_alignment_data, f)
    f.close()

    logger.info('Done')
    logger.info('Time taken: ' + str(round((time.time() - start_time), 3)) + ' seconds.\n')

    return inter_cluster_alignment_data

# def load_inter_family_TM_alignment_data(input_fname_base, alignment_tool, alignment_dir, graphs_and_pickles_dir, alignment_fname, clusters, corresponding_loop_data):
def load_inter_family_TM_alignment_data(input_fname_base, alignment_tool, alignment_dir, graphs_and_pickles_dir, clusters, corresponding_loop_data):
    alignment_data_fname = os.path.join(graphs_and_pickles_dir, 'inter_family_alignment_data_' + input_fname_base + '_' + alignment_tool + '.pickle2')
    if (sys.version_info >= (3, 0)):
        alignment_data_fname = os.path.join(graphs_and_pickles_dir, 'inter_family_alignment_data_' + input_fname_base + '_' + alignment_tool + '.pickle3')

    # if previous_graph_file_reused == True:
        # if is_valid_pickle(alignment_data_fname, clusters) == True:
        # if True:
    if os.path.exists(alignment_data_fname):
        logger.info('Loading saved alignment data from previous run (In ' + alignment_data_fname[base_path_len:] + ') ...\n')
        f = open(alignment_data_fname, 'rb')
        cluster_alignment_data = pickle.load(f)
        f.close()
        return cluster_alignment_data

    logger.info('Loading inter-family alignment data from alignment files.')
    start_time = time.time()

    # alignment_data = load_graph_data(alignment_fname)

    # inter_cluster_alignment_data, loops_node_list = load_inter_cluster_alignment_data(alignment_data, clusters, corresponding_loop_data)
    inter_cluster_alignment_data, loops_node_list = load_inter_cluster_alignment_data(clusters, corresponding_loop_data)

    # print(inter_cluster_alignment_data)
    # sys.exit()


    # cluster_alignment_data, loops_node_list = load_inter_cluster_alignment_data(alignment_data, clusters)

    # print('Number of loops in cluster file: ' + str(len(loops_in_cluster)))

    file_counter = 0
    # print(loops_node_list)
    # sys.exit()
    for node in loops_node_list:
        print(node)
        # for r1 in get_all_loop_combination(str(node)):
        for r1 in [str(node)]:
            node1 = strToNode(r1)
    # for fn in glob.glob(os.path.join(alignment_dir, '*.aln')):
            fn = os.path.join(alignment_dir, r1 + '.aln')

            if not os.path.isfile(fn):
                return None

            # stime = datetime.now()
            file_counter += 1
            print('Processsing ' + fn[base_path_len:] + ' ... (' + str(file_counter) + ')')
            # sys.stdout.flush()
            # r1 = os.path.basename(fn)[:-4]
            # node1 = strToNode(r1)

            # if node1 not in loops_in_cluster:
            #     continue

            cid_nodelist_pair = find_nodes_in_cluster(node1, inter_cluster_alignment_data)

            fp = open(fn)
            lines = fp.readlines()
            fp.close()
            # flag = 0
            # print('No. of lines: ' + str(len(lines)))
            # print('Connected nodes: ' + str(len(node_dict)))
            line_index = 0
            while line_index < len(lines):
                # print 'Reading line ' + str(line_index)
                # sys.exit()
                line = lines[line_index]
                # if lines[line_index].startswith('#  Aligning'):
                if line.startswith('Name of Chain_1:'):
                    r1_line = line
                    r2_line = lines[line_index+1]
                    score_line = lines[line_index+8]

                    test_r1 = os.path.basename(r1_line.strip().split(' ')[3].strip())[:-4]
                    if r1 != test_r1:
                        logger.error('ERROR: filename and loop mismatch. r1(filename): ' + r1 + ', r1(loop): ' + test_r1)
                        sys.exit()

                    r2 = os.path.basename(r2_line.strip().split(' ')[3].strip())[:-4]
                    node2 = strToNode(r2)

                    for c_id, node_dict in cid_nodelist_pair:
                        if node2 in node_dict:

                            # t1 = inter_cluster_alignment_data[c_id][node1][node2][0]
                            # t2 = inter_cluster_alignment_data[c_id][node1][node2][1]
                            # # make sure the best alignment ordering is same as the ordering of loop in current file
                            # if not ((r1 == t1 and r2 == t2) or (r1 == t2 and r2 == t1)):
                            #     # line_index += 12
                            #     continue

                            # #safety condition for reverse order (might be redundant)
                            # if (r1 == t2 and r2 == t1):
                            #     t1 = r2
                            #     t2 = r1

                            score_line = lines[line_index+8]
                            score_text = score_line.strip().split('(')[0].strip().split('=')[1].strip()
                            if score_text == '':
                                score = -50.
                                logger.error('ERROR: No alignment score found for: ' + r1 + ' and ' + r2)
                                sys.exit()
                            else:
                                score = float(score_text)

                            # text = lines[line_index+3].split(':')[1].strip()
                            # cr1 = get_local_alignment_index(r1, text)

                            # text = lines[line_index+4].split(':')[1].strip()
                            # cr2 = get_local_alignment_index(r2, text)

                            cr1 = r1
                            cr2 = r2

                            # aln1 = lines[line_index+6].strip()
                            # aln2 = lines[line_index+7].strip()

                            aln1 = lines[line_index + 12].strip()
                            aln2 = lines[line_index + 14].strip()

                            aln1, aln2 = adjust_rotated_tmalign_alignment_result_v2(r1, r2, aln1, aln2)

                            if '*' in aln1:
                                aln1, aln2 = adjust_rotated_tmalign_alignment_result(r1, r2, aln1, aln2)
                                # print(aln1, aln2)
                                # sys.exit()
                                star_ind = aln1.index('*')
                                aln1 = aln1[:star_ind] + aln1[star_ind+1:]
                                aln2 = aln2[:star_ind] + aln2[star_ind+1:]
                                t1 = r1.strip().split(':')[0] + ':' + '_'.join(rotate_l(r1.strip().split(':')[1].strip().split('_'), 1))
                                t2 = r2
                                # print(aln1, aln2)
                                # print(t1, t2)
                                # sys.exit()

                            if len(aln1) == 0 or len(aln2) == 0:
                                # set dummy alignment
                                aln1 = 'A'
                                aln2 = 'A'
                                temp_cr1 = cr1.split(':')
                                dummy_index = temp_cr1[1].split('_')[0].split('-')[0]
                                cr1 = temp_cr1[0] + ':' + dummy_index + '-' + dummy_index
                                temp_cr2 = cr2.split(':')
                                dummy_index = temp_cr2[1].split('_')[0].split('-')[0]
                                cr2 = temp_cr2[0] + ':' + dummy_index + '-' + dummy_index

                            zscore = inter_cluster_alignment_data[c_id][node1][node2][2]

                            prev_score = inter_cluster_alignment_data[c_id][node1][node2][7]
                            if score > prev_score:
                                inter_cluster_alignment_data[c_id][node1][node2] = (r1, r2, zscore, cr1, cr2, aln1, aln2, score)
                                inter_cluster_alignment_data[c_id][node2][node1] = (r2, r1, zscore, cr2, cr1, aln2, aln1, score)
                            # break
                    # we can skip at least 11 lines from current alignment data
                    line_index += 11
                line_index += 1
            # while end
            # print('ok')
            # etime = datetime.now()
            # difftime = (etime - stime).total_seconds()
            # print('Time taken: ')
            # print(difftime)
            # break

    f = open(alignment_data_fname,"wb")
    pickle.dump(inter_cluster_alignment_data, f)
    f.close()

    logger.info('Done')
    logger.info('Time taken: ' + str(round((time.time() - start_time), 3)) + ' seconds.\n')

    return inter_cluster_alignment_data

def generate_inter_family_rmsd_data(input_fname_base, alignment_tool, partial_pdbx_dir, graphs_and_pickles_dir, alignment_data, clusters, loop_list):
# def generate_inter_family_rmsd_data(input_fname_base, alignment_tool, partial_pdbx_dir, graphs_and_pickles_dir, alignment_data, clusters, loop_list, previous_graph_file_reused):
# def generate_rmsd_data(input_fname_base, partial_pdbx_dir, graphs_and_pickles_dir, alignment_data, clusters, loop_list, previous_graph_file_reused):

    rmsd_data_fname = os.path.join(graphs_and_pickles_dir, 'inter_family_rmsd_data_' + input_fname_base + '_' + alignment_tool + '.pickle2')
    if (sys.version_info >= (3, 0)):
        rmsd_data_fname = os.path.join(graphs_and_pickles_dir, 'inter_family_rmsd_data_' + input_fname_base + '_' + alignment_tool + '.pickle3')
     
    # if previous_graph_file_reused == True:
        # if is_valid_pickle(rmsd_data_fname, clusters) == True:
        # if True:
    if os.path.exists(rmsd_data_fname):
        logger.info('Loading saved RMSD data from previous run (In ' + rmsd_data_fname[base_path_len:] + ') ...\n')
        f = open(rmsd_data_fname, 'rb')
        rmsd_data_dict = pickle.load(f)
        f.close()
        return rmsd_data_dict

    # print(alignment_data['GNGA'][strToNode('4V9F_0:2621-2626')][strToNode('4V9F_0:460-465')])
    # print(alignment_data['GNAA'][strToNode('4V9F_0:2621-2626')][strToNode('4V9F_0:460-465')])
    # sys.exit()
    logger.info('Generating inter-family RMSD data.')
    start_time = time.time()

    rmsd_data_dict = {}
    coord_dict = {}
    pdb_structure = None
    prev_pdb_chain = ''
    structure_counter = 0

    # for lp in loop_list:
    #     pdb_chain, regions = lp.split(':')
    #     pdb = pdb_chain.split('_')[0]
    #     pdb_pm = get_pdb_index_list(lp)
    #     if prev_pdb_chain != pdb_chain:
    #         pdb_structure = None
    #     else:
    #         structure_counter += 1
    #         if structure_counter % 500 == 0:
    #             gc.collect()
    #     coord_backbone, coord_sugar, pdb_structure = get_atom_coordinate(os.path.join(pdbx_dir, pdb+'.cif'), pdb_pm, pdb_structure)
    #     coord_dict[lp] = (coord_backbone, coord_sugar)
    #     prev_pdb_chain = pdb_chain

    for lp in loop_list:
        pdb_pm = get_pdb_index_list(lp)
        coord_backbone, coord_sugar, pdb_structure = get_atom_coordinate(os.path.join(partial_pdbx_dir, lp + '.cif'), pdb_pm)
        coord_dict[lp] = (coord_backbone, coord_sugar)

    # time_align_residue = 0
    # time_get_coordinate = 0
    # time_rmsd = 0
    # time_start = time.time()
    for cluster_id in alignment_data:
        # if cluster_id not in clusters:
        #     continue
        sum_of_avg_rmsd_for_c = 0.
        rmsd_data_list_dict = {}
        index_dict = {}
        i = 0
        # print(loop_list)
        for l1 in alignment_data[cluster_id]:
            index_dict[l1] = i
            i += 1

        # pdb_chain_dict = {}
        # pdb_res_mapping_dict = {}
        # fasta_seq_dict = {}
        # for l1 in alignment_data[cluster_id]:
        #     pdb_chain, _ = str(l1).strip().split(':')
            
        #     pdb_id, chain_id = pdb_chain.strip().split('_')
        #     if pdb_id not in pdb_chain_dict:
        #         pdb_chain_dict[pdb_id] = []
        #     pdb_chain_dict[pdb_id].append(chain_id)

        #     if pdb_chain not in pdb_res_mapping_dict:
        #         pdb_res_mapping_dict[pdb_chain] = load_pdb_res_map(pdb_chain)

        # for pdb_id in pdb_chain_dict:
        #     fasta_seq_dict.update(load_fasta_seq(pdb_id, pdb_chain_dict[pdb_id]))

        pdb_res_mapping_dict, fasta_seq_dict = load_pdb_fasta_mapping_and_fasta_seq_dict(cluster_id, alignment_data)

        # i = 0
        for l1 in alignment_data[cluster_id]:
            # if str(l1) not in loop_list:
            #     continue
            fit_ret = []
            sum_of_rmsd_for_l1 = 0.
            # rmsd_data_list_item = {}
            # j = 0
            for l2 in alignment_data[cluster_id][l1]:
                # if str(l2) not in loop_list:
                #     continue
                if l1 != l2:
                    # print(cluster_id)
                    r1, r2, zscore, cr1, cr2, aln_1, aln_2, score = alignment_data[cluster_id][l1][l2]
                    # print('r1, r2: ')
                    # print(r1, r2)
                    # print(r1, r2, zscore, cr1, cr2, aln_1, aln_2, score)
                    # if output_env == 'local':
                    # time_s = time.time()
                    pdb1_pm, pdb2_pm, i1_pm, i2_pm = aln_residue_temp(pdb_res_mapping_dict, fasta_seq_dict, r1, r2, cr1, cr2, aln_1, aln_2, 0, len(aln_1)-1, 0)
                    # time_align_residue += time.time() - time_s
                    # else:
                    #     pdb1_pm, pdb2_pm, i1_pm, i2_pm = aln_residue(r1, r2, cr1, cr2, aln_1, aln_2, 0, len(aln_1)-1, 0)

                    pdb_chain1, _ = r1.split(':')
                    pdb_chain2, _ = r2.split(':')
                    pdb1 = pdb_chain1.split('_')[0]
                    pdb2 = pdb_chain2.split('_')[0]

                    # structures = {}
                    #returns centroid of backbone atoms
                    # time_s = time.time()
                    coord1 = extract_atom_coordinate(coord_dict[str(l1)], pdb1_pm, pdb1)
                    coord2 = extract_atom_coordinate(coord_dict[str(l2)], pdb2_pm, pdb2)
                    # time_get_coordinate += time.time() - time_s
                    # print(coord1, coord2)

                    X, Y = convert_array(coord1, coord2)

                    if len(X) != len(Y):
                        logger.warning('WARNING: Corresponding co-ordinates for alignments not found! rmsd = 20 assigned.')
                        rmsd = 20.
                    elif len(X) == 0:
                        logger.warning('WARNING: Co-ordinates for alignments not found! rmsd = 20 assigned.')
                        rmsd = 20.
                    else:
                        XC = sum(X)/len(X)
                        YC = sum(Y)/len(Y)
                        # calculating relative co-ordinate using mean as reference
                        X -= XC
                        Y -= YC
                        # time_s = time.time()
                        rmsd = kabsch_rmsd(X, Y)
                        # time_rmsd += time.time() - time_s

                    sum_of_rmsd_for_l1 += rmsd
                    # X.shape[0] represents the number of aligned nucleotides
                    # fit_ret.append((index_dict[l2], str(l2), rmsd, X.shape[0]))
                    fit_ret.append((index_dict[l2], str(l2), rmsd, len(pdb1_pm)))
                    # end of if
                # j += 1
                # end of l2 for loop

            # time_diff = time.time() - time_start
            # if time_diff > 30:
            #     print(cluster_id, l1)
            #     print('time_align_residue', time_align_residue)
            #     print('time_get_coordinate', time_get_coordinate)
            #     print('time_rmsd', time_rmsd)
            #     time_start = time.time()
            #     print('')

            avg_of_rmsd_for_l1 = 0.0
            if (len(clusters[cluster_id]) - 1) > 0:
                avg_of_rmsd_for_l1 = sum_of_rmsd_for_l1 / (len(clusters[cluster_id]) - 1)
            sum_of_avg_rmsd_for_c += avg_of_rmsd_for_l1
            # print str(i)+','+l1+'\t'+'|'.join(map(lambda x: str(x[0])+','+x[1]+','+str(x[2])+','+str(x[3]), sorted(fit_ret, key=lambda x: x[2])))
            rmsd_data_list_dict[(index_dict[l1], str(l1))] = (avg_of_rmsd_for_l1, sorted(fit_ret, key=lambda x: x[2]))
            # i += 1
            # end of l1 for loop

        avg_of_avg_rmsd_for_c = 0.0
        if len(clusters[cluster_id]) > 0:
            avg_of_avg_rmsd_for_c = sum_of_avg_rmsd_for_c / len(clusters[cluster_id])

        rmsd_data_dict[cluster_id] = (avg_of_avg_rmsd_for_c, rmsd_data_list_dict) # sorted(rmsd_data_list_dict, key=lambda x: x[0][1]) #order by loop average

    f = open(rmsd_data_fname,"wb")
    pickle.dump(rmsd_data_dict, f)
    f.close()

    logger.info('Done')
    logger.info('Time taken: ' + str(round((time.time() - start_time), 3)) + ' seconds.\n')

    return rmsd_data_dict

def load_inter_family_alignment_and_rmsd_data(clusters, loop_list, corresponding_loop_data, input_fname_base, alignment_tool, partial_pdbx_dir, alignment_dir, graphs_and_pickles_dir):
    # graph_fname = os.path.join(graphs_and_pickles_dir, input_fname_base + '_' + alignment_tool + '.z.graph')

    if alignment_tool == 'TMalign':
        # alignment_data = load_inter_family_TM_alignment_data(input_fname_base, alignment_tool, alignment_dir, graphs_and_pickles_dir, graph_fname, clusters, corresponding_loop_data)
        alignment_data = load_inter_family_TM_alignment_data(input_fname_base, alignment_tool, alignment_dir, graphs_and_pickles_dir, clusters, corresponding_loop_data)
    else:
        # alignment_data = load_inter_family_alignment_data(input_fname_base, alignment_tool, alignment_dir, graphs_and_pickles_dir, graph_fname, clusters, corresponding_loop_data, previous_graph_file_reused)
        alignment_data = load_inter_family_alignment_data(input_fname_base, alignment_tool, alignment_dir, graphs_and_pickles_dir, clusters, corresponding_loop_data)

    if alignment_data == None:
        return None, None

    rmsd_data_dict = generate_inter_family_rmsd_data(input_fname_base, alignment_tool, partial_pdbx_dir, graphs_and_pickles_dir, alignment_data, clusters, loop_list)

    return alignment_data, rmsd_data_dict

def find_mean_rmsd_with_specific_cluster(cluster_id, clusters, pairwise_align_details):
    rmsd_list = []
    for (j, r2, rmsd, align_length) in pairwise_align_details:
        if r2 not in clusters[cluster_id]:
            continue
        rmsd_list.append(rmsd)
    return sum(rmsd_list) / float(len(rmsd_list))

def find_best_aligned_pair_of_specific_cluster(cluster_id, clusters, pairwise_align_details, align_len_threshold):
    # max align loop

    max_align_length = max_loop_index = 0
    max_length_loop = ''
    max_loop_rmsd = 1000.0

    # count = 0
    for (j, r2, rmsd, align_length) in pairwise_align_details:
        if r2 not in clusters[cluster_id]:
            continue
        # if r1 == '5J7L_DA:1711-1717_1741-1745':
        #     print(j, r2, rmsd, align_length)
        # count += 1
        if is_better_alignment_score((rmsd, align_length), (max_loop_rmsd, max_align_length), align_len_threshold, is_normalized_score):
            max_align_length = align_length
            max_length_loop = r2
            max_loop_index = j
            max_loop_rmsd = rmsd
    
    # print(count)
    return max_loop_index, max_length_loop, max_loop_rmsd, max_align_length

def extract_inter_family_current_rmsd_data_dict(rmsd_data_dict, cluster_id, loops1, loops2):

    for i in range(len(loops1)):
        loops1[i] = str(strToNode(loops1[i]))
    for i in range(len(loops2)):
        loops2[i] = str(strToNode(loops2[i]))

    new_rmsd_data_list_dict = {}
    # total_rmsd = 0.0
    rmsd_align_len_list = []
    _, rmsd_data_list_dict = rmsd_data_dict[cluster_id]
    not_found_count = 0
    for (i, r1) in rmsd_data_list_dict:
        if r1 in loops1:
            _, pairwise_align_details = rmsd_data_list_dict[(i, r1)]
            # total_rmsd_for_r1 = 0.0
            rmsd_align_len_list_for_r1 = []
            new_pairwise_align_details = []
            for (j, r2, rmsd, align_length) in pairwise_align_details:
                if r2 in loops2:
                    new_pairwise_align_details.append((j, r2, rmsd, align_length))
                    rmsd_align_len_list_for_r1.append((rmsd, align_length))
                    # total_rmsd_for_r1 += rmsd

            avg_rmsd_for_r1, total_align_len_for_r1 = get_weighted_avg_rmsd(rmsd_align_len_list_for_r1)
            # avg_rmsd_for_r1 = 0.0
            # if len(new_pairwise_align_details) > 0:
            #     avg_rmsd_for_r1 = total_rmsd_for_r1 / len(new_pairwise_align_details)
            new_rmsd_data_list_dict[(i, r1)] = (avg_rmsd_for_r1, sorted(new_pairwise_align_details, key=lambda x: x[2]))
            rmsd_align_len_list.append((avg_rmsd_for_r1, total_align_len_for_r1))
            # total_rmsd += avg_rmsd_for_r1
        else:
            # print r1 + ' not found'
            not_found_count += 1
            #sys.exit()

    # print str(not_fount_count) + ' loops not found.'
    avg_rmsd, total_align_len = get_weighted_avg_rmsd(rmsd_align_len_list)
    # avg_rmsd = 0.0
    # if len(new_rmsd_data_list_dict) > 0:
    #     avg_rmsd = total_rmsd / len(new_rmsd_data_list_dict)

    return (avg_rmsd, new_rmsd_data_list_dict)

def regenerate_families_from_superimposition_output(superimposition_output_dir):
    output_fname = os.path.join(superimposition_output_dir, 'subfamily_cluster.csv')
    if not os.path.exists(output_fname):
        print('Output subfamily file not found. Exiting.')
        sys.exit()
    fp = open(output_fname)
    lines = fp.readlines()
    fp.close()

    families = {}
    families_shortcoded = {}
    subfamilies = {}
    # subfamily_count = {}
    for line in lines:
        pieces = line.strip().split(',')
        subfamily_id = pieces[0]
        subpieces = subfamily_id.strip().split('-')
        family_id = '-'.join(subpieces[:-1])

        family_id_short = known_motif_shortcode[family_id.lower()]
        subfamily_id = family_id_short + subpieces[-1].lstrip('Sub')

        if family_id not in families:
            families[family_id] = []
        if family_id_short not in families_shortcoded:
            families_shortcoded[family_id_short] = []
        # if family_id not in subfamily_count:
        #     subfamily_count[family_id] = 0
        if family_id_short not in subfamilies:
            subfamilies[family_id_short] = {}
        if subfamily_id not in subfamilies[family_id_short]:
            subfamilies[family_id_short][subfamily_id] = []

        families[family_id] += pieces[1:]
        families_shortcoded[family_id_short] += pieces[1:]
        subfamilies[family_id_short][subfamily_id] += pieces[1:]
        # subfamily_count[family_id] += 1

    return families, families_shortcoded, subfamilies#, subfamily_count

def find_all_aligned_loops_of_specific_cluster(cluster_id, clusters, pairwise_align_details):
    item_list = []

    # count = 0
    for (j, r2, rmsd, align_length) in pairwise_align_details:
        if r2 not in clusters[cluster_id]:
            continue

        item_list.append((j, r2, rmsd, align_length))
        # if cluster_id == 'Sarcin-ricin':
            # print('rmsd')
        # print('rmsd', rmsd, align_length)
        # if align_length >= align_len_threshold and rmsd >= rmsd_threshold:
            # count += 1
            # print('filtered_rmsd', rmsd, align_length)
        
    
    # print(count)
    # return count
    return item_list

def find_avg_rmsd_of_all_filtered_pairs(filtered_pair_list):
    sum_rmsd = 0.0
    sum_aln_len = 0
    for (j, r2, rmsd, align_length) in filtered_pair_list:
        sum_rmsd += rmsd
        sum_aln_len += align_length

    return sum_rmsd/float(len(filtered_pair_list)), sum_aln_len/float(len(filtered_pair_list))

def generate_relative_graph_among_motif_families(families, loop_list, familywise_alignment_data, familywise_rmsd_data_dict, input_fname_base, alignment_tool, alignment_dir, graphs_and_pickles_dir, partial_pdbx_dir, output_dir, rmsd_threshold, participating_motif_instance_threshold):

    corresponding_loop_data = {}
    for c1_id in families:
        corresponding_loop_data[c1_id] = []
        for c2_id in families:
            if c1_id == c2_id:
                continue
            corresponding_loop_data[c1_id] += families[c2_id]

    alignment_data, rmsd_data_dict = load_inter_family_alignment_and_rmsd_data(families, loop_list, corresponding_loop_data, input_fname_base, alignment_tool, partial_pdbx_dir, alignment_dir, graphs_and_pickles_dir)

    # print('Checking alignment_data')
    # print(len(familywise_alignment_data['E-loop'][strToNode('4V88_A6:626-630_967-971')]))
    # print(len(familywise_alignment_data['Tandem-shear'][strToNode('5J7L_AA:1415-1418_1480-1483')]))

    # print('Checking rmsd_data')
    # ky = list(familywise_rmsd_data_dict['E-loop'][1].keys())[0]
    # print(ky)
    # print(len(familywise_rmsd_data_dict['E-loop'][1][ky][1]))
    # ky = list(familywise_rmsd_data_dict['Tandem-shear'][1].keys())[0]
    # print(ky)
    # print(len(familywise_rmsd_data_dict['Tandem-shear'][1][ky][1]))

    # print('#####')

    # print('Checking inter-family alignment_data')
    # print(len(alignment_data['E-loop'][strToNode('4V88_A6:626-630_967-971')]))
    # print(len(alignment_data['Tandem-shear'][strToNode('5J7L_AA:1415-1418_1480-1483')]))

    # print('Checking inter-family rmsd_data')
    # ky = list(rmsd_data_dict['E-loop'][1].keys())[0]
    # print(ky)
    # print(len(rmsd_data_dict['E-loop'][1][ky][1]))
    # ky = list(rmsd_data_dict['Tandem-shear'][1].keys())[0]
    # print(ky)
    # print(len(rmsd_data_dict['Tandem-shear'][1][ky][1]))

    # sys.exit()


    length_adjusted_rmsd_score_dict = rmsd_data_dict
    length_adjusted_familywise_rmsd_score_dict = familywise_rmsd_data_dict
    
    if is_length_adjusted_score:
        length_adjusted_rmsd_score_dict = generate_length_adjusted_rmsd_score(rmsd_data_dict)
        length_adjusted_familywise_rmsd_score_dict = generate_length_adjusted_rmsd_score(familywise_rmsd_data_dict)

    
    familywise_current_rmsd_data_dict = {}
    for cluster_id in families:
        loops = families[cluster_id]
        logger.info('Extracting familywise rmsd data dict for ' + cluster_id)
        familywise_current_rmsd_data_dict[cluster_id] = extract_current_rmsd_data_dict(length_adjusted_familywise_rmsd_score_dict, cluster_id, loops)
        # ky = list(familywise_current_rmsd_data_dict[cluster_id][1].keys())[0]
        # print(ky)
        # print(len(familywise_current_rmsd_data_dict[cluster_id][1][ky][1]))
        # _, familywise_cluster_pairwise_alignment_details = familywise_current_rmsd_data_dict[cluster_id]
        logger.info('Completed familywise extracting rmsd data dict for ' + cluster_id)
    # sys.exit()
    current_rmsd_data_dict = {}
    for cluster_id in families:
        loops = families[cluster_id]
        corresponding_loops = corresponding_loop_data[cluster_id]
        # loops = corresponding_loop_data[cluster_id]
        logger.info('Extracting rmsd data dict for ' + cluster_id)
        current_rmsd_data_dict[cluster_id] = extract_inter_family_current_rmsd_data_dict(length_adjusted_rmsd_score_dict, cluster_id, loops, corresponding_loops)
        # ky = list(current_rmsd_data_dict[cluster_id][1].keys())[0]
        # print(ky)
        # print(len(current_rmsd_data_dict[cluster_id][1][ky][1]))
        # _, cluster_pairwise_alignment_details = current_rmsd_data_dict[cluster_id]
        logger.info('Completed extracting rmsd data dict for ' + cluster_id)


    print('\n')
    inter_family_adj_mat = {}
    inter_family_adj_mat_avg = {}
    # inter_family_adj_data = {}

    participating_r1_list_dict = {}
    participating_r2_list_dict = {}
    for cluster1_id in families:
        # cluster1_id = 'Sarcin-ricin'
        # cluster1_id = 'E-loop'
        
        cluster1_id_short = cluster1_id
        if cluster1_id.lower() in known_motif_shortcode:
            cluster1_id_short = known_motif_shortcode[cluster1_id.lower()]

        if cluster1_id_short not in participating_r1_list_dict:
            participating_r1_list_dict[cluster1_id_short] = {}
        if cluster1_id_short not in participating_r2_list_dict:
            participating_r2_list_dict[cluster1_id_short] = {}

        # if cluster1_id_short not in inter_family_adj_data:
        #     inter_family_adj_data[cluster1_id_short] = {}

        _, cluster1_pairwise_alignment_details = current_rmsd_data_dict[cluster1_id]
        # _, cluster1_pairwise_alignment_details = familywise_current_rmsd_data_dict[cluster1_id]
        
        # align_len_threshold1 = generate_align_length_threshold(cluster1_pairwise_alignment_details)
        align_len_threshold1 = generate_align_length_threshold(familywise_current_rmsd_data_dict[cluster1_id][1])
        # print(cluster1_pairwise_alignment_details)
        # sys.exit()
        for cluster2_id in families:
            # print('cl1, cl2')
            # print(cluster1_id, cluster2_id)
            # cluster2_id = 'E-loop'
            # cluster2_id = 'Sarcin-ricin'
            # cluster2_id = 'Tetraloop-receptor'
            
            if cluster1_id == cluster2_id:
                # print('continueing')
                continue

            cluster2_id_short = cluster2_id
            if cluster2_id.lower() in known_motif_shortcode:
                cluster2_id_short = known_motif_shortcode[cluster2_id.lower()]

            if cluster2_id_short not in participating_r1_list_dict[cluster1_id_short]:
                participating_r1_list_dict[cluster1_id_short][cluster2_id_short] = []
            if cluster2_id_short not in participating_r2_list_dict[cluster1_id_short]:
                participating_r2_list_dict[cluster1_id_short][cluster2_id_short] = []
            # if cluster2_id_short not in inter_family_adj_data[cluster1_id_short]:
            #     inter_family_adj_data[cluster1_id_short][cluster2_id_short] = []

            _, cluster2_pairwise_alignment_details = current_rmsd_data_dict[cluster2_id]
            # _, cluster2_pairwise_alignment_details = familywise_current_rmsd_data_dict[cluster2_id]

            # align_len_threshold2 = generate_align_length_threshold(cluster2_pairwise_alignment_details)
            align_len_threshold2 = generate_align_length_threshold(familywise_current_rmsd_data_dict[cluster2_id][1])


            # if cluster1_id == 'E-loop' and cluster2_id == 'Tandem-shear':
            #     print('Sending RMSD data for E-loop threshold')
            #     ky = list(familywise_current_rmsd_data_dict[cluster1_id][1].keys())[0]
            #     print(ky)
            #     print(len(familywise_current_rmsd_data_dict[cluster1_id][1][ky][1]))
            #     print(familywise_current_rmsd_data_dict[cluster1_id][1][ky][0])
            #     # print(sorted(familywise_current_rmsd_data_dict[cluster1_id][1][ky][1], key = lambda x: x[1]))
            #     for it in sorted(familywise_current_rmsd_data_dict[cluster1_id][1][ky][1], key = lambda x: x[1]):
            #         print(it)
            #     # sys.exit()
            #     # print(familywise_current_rmsd_data_dict[cluster1_id][1])
            #     print('Calculated threshold: ')
            #     print(align_len_threshold1)
            #     print('Sending RMSD data for Tandem-shear threshold')
            #     ky = list(familywise_current_rmsd_data_dict[cluster2_id][1].keys())[0]
            #     print(ky)
            #     print(len(familywise_current_rmsd_data_dict[cluster2_id][1][ky][1]))
            #     print(familywise_current_rmsd_data_dict[cluster2_id][1][ky][0])
            #     # print(sorted(familywise_current_rmsd_data_dict[cluster1_id][1][ky][1], key = lambda x: x[1]))
            #     for it in sorted(familywise_current_rmsd_data_dict[cluster2_id][1][ky][1], key = lambda x: x[1]):
            #         print(it)
            #     print('Calculated threshold: ')
            #     print(align_len_threshold2)

            # if cluster1_id == 'Tandem-shear' and cluster2_id == 'E-loop':
            #     print('Sending RMSD data for Tandem-shear threshold')
            #     ky = list(familywise_current_rmsd_data_dict[cluster1_id][1].keys())[0]
            #     print(ky)
            #     print(len(familywise_current_rmsd_data_dict[cluster1_id][1][ky][1]))
            #     print(familywise_current_rmsd_data_dict[cluster1_id][1][ky][0])
            #     # print(sorted(familywise_current_rmsd_data_dict[cluster1_id][1][ky][1], key = lambda x: x[1]))
            #     for it in sorted(familywise_current_rmsd_data_dict[cluster1_id][1][ky][1], key = lambda x: x[1]):
            #         print(it)
            #     print('Calculated threshold: ')
            #     print(align_len_threshold1)
            #     print('Sending RMSD data for E-loop threshold')
            #     ky = list(familywise_current_rmsd_data_dict[cluster2_id][1].keys())[0]
            #     print(ky)
            #     print(len(familywise_current_rmsd_data_dict[cluster2_id][1][ky][1]))
            #     print(familywise_current_rmsd_data_dict[cluster2_id][1][ky][0])
            #     # print(sorted(familywise_current_rmsd_data_dict[cluster1_id][1][ky][1], key = lambda x: x[1]))
            #     for it in sorted(familywise_current_rmsd_data_dict[cluster2_id][1][ky][1], key = lambda x: x[1]):
            #         print(it)
            #     print('Calculated threshold: ')
            #     print(align_len_threshold2)


            # align_len_threshold = round(0.8 * min(align_len_threshold1, align_len_threshold2))
            align_len_threshold = max(align_len_threshold1, align_len_threshold2)
            # max value ensures that both motifs can have their own structure in the alignment

            # print(align_len_threshold1)
            # print(align_len_threshold2)
            # print(align_len_threshold)


            # print('Between ' + cluster1_id + ' and ' + cluster2_id)

            best_loop_rmsd = 1000.0
            best_align_length = 0
            best_r1 = ''
            best_r2 = ''
            mean_rmsd_list = []
            pair_count = 0
            total_count = 0
            # rmsd_threshold = rmsd_threshold_for_merging
            filtered_pair_list = []
            # participating_r1_list = []
            for (i, r1) in sorted(cluster1_pairwise_alignment_details, key=lambda x: cluster1_pairwise_alignment_details[x][0]):
                avg_rmsd, pairwise_align_details = cluster1_pairwise_alignment_details[(i, r1)]
                # print(avg_rmsd, pairwise_align_details)
                # sys.exit()
                if len(pairwise_align_details) > 0:
                    # print(align_len_threshold)
                    # sys.exit()

                    mean_rmsd_list.append(find_mean_rmsd_with_specific_cluster(cluster2_id, families, pairwise_align_details))
                    aligned_list = find_all_aligned_loops_of_specific_cluster(cluster2_id, families, pairwise_align_details)
                    total_count += len(aligned_list)
                    filtered_list = list(filter(lambda x: (x[2] <= rmsd_threshold and x[3] >= align_len_threshold), aligned_list))
                    filtered_pair_list += filtered_list
                    # pair_count += find_all_aligned_pairs_of_specific_cluster(cluster2_id, families, pairwise_align_details, align_len_threshold, 1.5)
                    pair_count += len(filtered_list)
                    if len(filtered_list) > 0 and r1 not in participating_r1_list_dict[cluster1_id_short][cluster2_id_short]:
                        participating_r1_list_dict[cluster1_id_short][cluster2_id_short].append(r1)

                    if len(filtered_list) > 0:
                        for (j, r2, rmsd, align_length) in filtered_list:
                            if r2 not in participating_r2_list_dict[cluster1_id_short][cluster2_id_short]:
                                participating_r2_list_dict[cluster1_id_short][cluster2_id_short].append(r2)

                    j, r2, max_loop_rmsd, max_align_length = find_best_aligned_pair_of_specific_cluster(cluster2_id, families, pairwise_align_details, align_len_threshold)
                    
                    if max_align_length >= align_len_threshold and max_loop_rmsd <= rmsd_threshold:
                    # inter_family_adj_data[cluster1_id_short][cluster2_id_short].append((r1, r2, max_loop_rmsd, max_align_length))
                        if is_better_alignment_score((max_loop_rmsd, max_align_length), (best_loop_rmsd, best_align_length), align_len_threshold, is_normalized_score):
                            best_loop_rmsd = max_loop_rmsd
                            best_align_length = max_align_length
                            best_r1 = r1
                            best_r2 = r2

            # print('pair_count', pair_count)
            # print('total_count', total_count)
            # print(cluster1_id, cluster2_id)
            # print(len(families[cluster1_id]), len(families[cluster2_id]))
            # print(best_align_length)
            # print(best_r1, best_r2, best_loop_rmsd, best_align_length)
            # pair_avg = 100.0 * pair_count / float(len(families[cluster1_id]) * len(families[cluster2_id]))
            pair_avg = 100.0 * pair_count / float(total_count)
            # cluster1_id_short = known_motif_shortcode[cluster1_id.lower()]
            # cluster2_id_short = known_motif_shortcode[cluster2_id.lower()]

            if cluster1_id_short not in inter_family_adj_mat:
                inter_family_adj_mat[cluster1_id_short] = {}
            if cluster1_id_short not in inter_family_adj_mat_avg:
                inter_family_adj_mat_avg[cluster1_id_short] = {}
            # if cluster2_id not in inter_family_adj_mat:
            #     inter_family_adj_mat[cluster2_id] = {}
            avg_mean_rmsd = sum(mean_rmsd_list) / float(len(mean_rmsd_list))
            inter_family_adj_mat_avg[cluster1_id_short][cluster2_id_short] = avg_mean_rmsd

            if best_r1 == '' or best_r2 == '':
                continue
            if input_index_type == 'pdb':
                best_r1 = convert_a_loop_from_FASTA_to_PDB(best_r1)
                best_r2 = convert_a_loop_from_FASTA_to_PDB(best_r2)

            # avg_mean_rmsd = sum(mean_rmsd_list) / float(len(mean_rmsd_list))
            avg_rmsd_of_all_filtered_pairs, avg_aln_len_of_all_filtered_pairs = find_avg_rmsd_of_all_filtered_pairs(filtered_pair_list)
            inter_family_adj_mat[cluster1_id_short][cluster2_id_short] = (best_r1, best_r2, best_loop_rmsd, best_align_length, avg_mean_rmsd, pair_avg, pair_count, avg_rmsd_of_all_filtered_pairs, avg_aln_len_of_all_filtered_pairs)
            # inter_family_adj_mat[cluster2_id][cluster1_id] = 

            print('Between ' + cluster1_id + ' and ' + cluster2_id, end=',')
            print(' the best aligned pair (r1, r2) is (' + best_r1 + ', ' + best_r2 + ') with alignment length ' + str(best_align_length) + ' and RMSD ' + str(round(best_loop_rmsd, 2)) + ' A.')
            # print(known_motif_shortcode[cluster1_id.lower()] + ',' + known_motif_shortcode[cluster2_id.lower()], end=',')
            # print(best_r1, end=',')
            # print(best_r2, end=',')
            # print(best_loop_rmsd, end=',')
            # print(best_align_length, end=',')
            # print(str(round(avg_mean_rmsd, 2)))

            # print('Between ' + cluster1_id + ' and ' + cluster2_id)
            # print(best_r1)
            # print(best_r2)
            # print(best_loop_rmsd)
            # print(best_align_length)
            # sys.exit()
            # break

    # print('participating_r1_list_dict')
    # print(participating_r1_list_dict)
    # for fam1 in participating_r1_list_dict:
    #     for fam2 in participating_r1_list_dict[fam1]:
    #         if len(participating_r1_list_dict[fam1][fam2]) > 0:
    #             print(fam1 + ' - ' + fam2)
    #             motif_list = participating_r1_list_dict[fam1][fam2]
    #             motif_list = list(map(lambda x: convert_a_loop_from_FASTA_to_PDB(x), motif_list))
    #             for r1 in motif_list:
    #                 org_type = ''
    #                 rna_type = ''
    #                 sub = find_subfamily(fam1, subfamilies, r1)
    #                 pdb_chain, regions = r1.strip().split(':')
    #                 if pdb_chain in pdb_organism_details:
    #                     org_type = pdb_organism_details[pdb_chain][3][:10] + '.' if len(pdb_organism_details[pdb_chain][3]) > 10 else pdb_organism_details[pdb_chain][3]
    #                     rna_type = pdb_organism_details[pdb_chain][0]
    #                     sub = find_subfamily(fam1, subfamilies, r1)
    #                 print(r1 + ',' + sub + ', ' + org_type + ', ' + rna_type)
    #             # print('\n'.join(motif_list))
    #             print('\n')

    # print('participating_r2_list_dict')
    # print(participating_r2_list_dict)
    # for fam1 in participating_r2_list_dict:
    #     for fam2 in participating_r2_list_dict[fam1]:
    #         if len(participating_r2_list_dict[fam1][fam2]) > 0:
    #             print(fam1 + ' - ' + fam2)
    #             motif_list = participating_r2_list_dict[fam1][fam2]
    #             motif_list = list(map(lambda x: convert_a_loop_from_FASTA_to_PDB(x), motif_list))
    #             print('\n'.join(motif_list))
    #             print('\n')

    # draw_figure(families_shortcoded, subfamilies, inter_family_adj_mat, participating_r1_list_dict)

    families_shortcoded = {}
    for family_id in families:
        family_id_short = family_id
        if family_id.lower() in known_motif_shortcode:
            family_id_short = known_motif_shortcode[family_id.lower()]

        families_shortcoded[family_id_short] = copy.deepcopy(families[family_id])

    draw_figure(families_shortcoded, alignment_tool, inter_family_adj_mat, participating_r1_list_dict, rmsd_threshold, participating_motif_instance_threshold, output_dir)
    # draw_figure_with_mean_rmsd(families_shortcoded, subfamilies, inter_family_adj_mat_avg)

def find_subfamily(cluster_id, subfamilies, r):
    for sub_fam in subfamilies[cluster_id]:
        if r in subfamilies[cluster_id][sub_fam]:
            return sub_fam

def find_parent(parent, i):
    if parent[i] == i:
        return i
    return find_parent(parent, parent[i])

def do_unite(parent, rank, x, y):
    x_parent = find_parent(parent, x)
    y_parent = find_parent(parent, y)

    if rank[x_parent] < rank[y_parent]:
        parent[x_parent] = y_parent
    elif rank[x_parent] > rank[y_parent]:
        parent[y_parent] = x_parent
    else:
        parent[y_parent] = x_parent
        rank[x_parent] += 1

def extract_edges(inter_family_adj_mat, ordered_family_list, edge_source='subfamily'):
    edge_list = []
    # parent = []
    # rank = []
    # selected_edges = []
    # other_edges = []

    # for i in range(len(ordered_family_list)):
    #     parent.append(i)
    #     rank.append(0)

    for i in range(len(ordered_family_list)):
        for j in range(i+1, len(ordered_family_list)):
            # print(i, j)
            if ordered_family_list[i] in inter_family_adj_mat and ordered_family_list[j] in inter_family_adj_mat[ordered_family_list[i]]:
                if edge_source == 'subfamily':
                    (best_r1, best_r2, best_loop_rmsd, best_align_length, mean_rmsd, pair_avg, pair_count, avg_rmsd_of_all_filtered_pairs, avg_aln_len_of_all_filtered_pairs) = inter_family_adj_mat[ordered_family_list[i]][ordered_family_list[j]]
                    edge_list.append((i, j, best_r1, best_r2, best_loop_rmsd, best_align_length, mean_rmsd, pair_avg, pair_count, avg_rmsd_of_all_filtered_pairs, avg_aln_len_of_all_filtered_pairs))
                else:
                    mean_rmsd = inter_family_adj_mat[ordered_family_list[i]][ordered_family_list[j]]
                    edge_list.append((i, j, mean_rmsd))

    return edge_list

def KruskalMST(inter_family_adj_mat, ordered_family_list, edge_source='subfamily'):
    print('inside MST')
    edge_list = []
    parent = []
    rank = []
    selected_edges = []
    other_edges = []

    for i in range(len(ordered_family_list)):
        parent.append(i)
        rank.append(0)

    for i in range(len(ordered_family_list)):
        for j in range(i+1, len(ordered_family_list)):
            # print(i, j)
            if ordered_family_list[i] in inter_family_adj_mat and ordered_family_list[j] in inter_family_adj_mat[ordered_family_list[i]]:
                if edge_source == 'subfamily':
                    (best_r1, best_r2, best_loop_rmsd, best_align_length, mean_rmsd, pair_avg, pair_count, avg_rmsd_of_all_filtered_pairs, avg_aln_len_of_all_filtered_pairs) = inter_family_adj_mat[ordered_family_list[i]][ordered_family_list[j]]
                    edge_list.append((i, j, best_r1, best_r2, best_loop_rmsd, best_align_length, mean_rmsd, pair_avg, pair_count, avg_rmsd_of_all_filtered_pairs, avg_aln_len_of_all_filtered_pairs))
                else:
                    mean_rmsd = inter_family_adj_mat[ordered_family_list[i]][ordered_family_list[j]]
                    edge_list.append((i, j, mean_rmsd))

    if edge_source == 'subfamily':
        for edge in sorted(edge_list, key=lambda x: (x[4], x[5])):
            # print(edge[4], edge[5])
            x = find_parent(parent, edge[0])
            y = find_parent(parent, edge[1])

            if x != y:
                selected_edges.append(edge)
                do_unite(parent, rank, x, y)
            else:
                other_edges.append(edge)
    else:
        for edge in sorted(edge_list, key=lambda x: (x[2])):
            # print(edge[4], edge[5])
            x = find_parent(parent, edge[0])
            y = find_parent(parent, edge[1])

            if x != y:
                selected_edges.append(edge)
                do_unite(parent, rank, x, y)
            else:
                other_edges.append(edge)

    print(selected_edges)
    return selected_edges, other_edges

def write_participating_motif_instances_to_file(family_id_short1, motif_list1, family_id_short2, motif_list2, output_dir):
    output_fname = os.path.join(output_dir, family_id_short1 + '_' + family_id_short2 + '.csv')
    
    if input_index_type == 'pdb':
        motif_list1 = list(map(lambda x: convert_a_loop_from_FASTA_to_PDB(x), motif_list1))
        motif_list2 = list(map(lambda x: convert_a_loop_from_FASTA_to_PDB(x), motif_list2))

    fp = open(output_fname, 'w')
    
    fp.write(get_motif_family_fullname(family_id_short1) + ',')
    fp.write(','.join(motif_list1))
    fp.write('\n')

    fp.write(get_motif_family_fullname(family_id_short2) + ',')
    fp.write(','.join(motif_list2))
    fp.write('\n')

    fp.close()

# def get_coordinates_in_circle(n):
#     return_list = []
#     for i in range(n):
#         theta = float(i)/n*2*math.pi
#         x = np.cos(theta)
#         y = np.sin(theta)
#         return_list.append((x,y))
#     return return_list

# def draw_figure(families, subfamilies, inter_family_adj_mat, participating_r1_list_dict):
def draw_figure(families, alignment_tool, inter_family_adj_mat, participating_r1_list_dict, rmsd_threshold, participating_motif_instance_threshold, output_dir):
    plt_fig = plt.figure(frameon = False)
    plt_fig = plt.figure()
    plt_fig.set_size_inches(9, 9)

    filename = os.path.join(output_dir, 'Graph_ScanX.png')
    if alignment_tool == 'TMalign':
        filename = os.path.join(output_dir, 'Graph_RNAalign.png')
    # print(subfamily_count)
    # G = nx.DiGraph(directed=True)
    G = nx.Graph()
    circular_positions = get_coordinates_in_circle(len(families), 0.0, 0.0, 2.0)

    # priority_order = ["GNRA", "GNAA", "GNGA", "KT", "rKT", "SR", "CL", "EL", "HT", "TS", "TR", "L1C", "TL", "RS"]
    priority_order = ["GNRA", "GNAA", "GNGA", "KT", "SR", "TL", "RS", "CL", "L1C", "EL", "HT", "rKT", "TR", "TS"]
    # colors, colors_tuple = get_sub_family_colors(len(families))
    colors, colors_tuple = get_family_colors(len(families))
    # print(colors)

    fixed_positions = {}
    fixed_nodes = []
    node_colors = []
    node_sizes = []

    edge_list = []
    edge_colors = []

    # fixed_positions['X'] = (0.0, 0.0)
    # fixed_nodes.append('X')
    # node_colors.append('salmon')

    ordered_family_list = []
    for family_id in priority_order:
        if family_id in families:
            ordered_family_list.append(family_id)


    # print(ordered_family_list)
    # print(family_id)
    # print(ordered_family_list.index(family_id))
    
    for i, family_id in enumerate(ordered_family_list):
        # node = known_motif_shortcode[family_id.lower()]
        node = family_id
        fixed_positions[node] = circular_positions[i]
        fixed_nodes.append(node)
        node_colors.append(colors[ordered_family_list.index(family_id)])
        node_sizes.append(1800)

    
    # print(fixed_positions)
    # print(fixed_nodes)

    G.add_nodes_from(fixed_nodes)

    # selected_edges, other_edges = KruskalMST(inter_family_adj_mat, ordered_family_list)
    # print('selected_edges')
    # print(selected_edges)
    # print('other_edges')
    # print(other_edges)
    # selected_edges = selected_edges+other_edges

    selected_edge_list = []
    selected_edge_list2 = []
    selected_edge_list3 = []
    edge_labels = {}
    all_edge_list = extract_edges(inter_family_adj_mat, ordered_family_list)
    # for edge in selected_edges:
    edge_details_fname = os.path.join(output_dir, alignment_tool + '_edge_details.tsv')
    fp = open(edge_details_fname, 'w')
    fp.write('MF1\tMF2\tMean RMSD\tMean alignment length\tParticipating instances from MF1\tParticipating instances from MF2\tEdge label\n')
    for edge in all_edge_list:
        (i, j, best_r1, best_r2, best_loop_rmsd, best_align_length, mean_rmsd, pair_avg, pair_count, avg_rmsd_of_all_filtered_pairs, avg_aln_len_of_all_filtered_pairs) = edge
        # sub1 = ''
        # sub2 = ''
        # for subfamily_id in subfamilies[ordered_family_list[i]]:
        #     if best_r1 in subfamilies[ordered_family_list[i]][subfamily_id]:
        #         sub1 = subfamily_id
        #         break
        # for subfamily_id in subfamilies[ordered_family_list[j]]:
        #     if best_r2 in subfamilies[ordered_family_list[j]][subfamily_id]:
        #         sub2 = subfamily_id
        #         break

        # if sub1 == '' or sub2 == '':
        #     print('check')
        #     print(best_r1)
        #     print(best_r2)
        #     sys.exit()
        # else:
        # edge_list.append((sub1, sub2, best_loop_rmsd))
        # if best_loop_rmsd < 1.0:
        fam1_participating_r1_p = 100.0 * len(participating_r1_list_dict[ordered_family_list[i]][ordered_family_list[j]]) / len(families[ordered_family_list[i]])
        fam2_participating_r1_p = 100.0 * len(participating_r1_list_dict[ordered_family_list[j]][ordered_family_list[i]]) / len(families[ordered_family_list[j]])

        max_participating_r1_p = min(fam1_participating_r1_p, fam2_participating_r1_p)

        # fp.write(str(round(avg_rmsd_of_all_filtered_pairs, 2)), str(round(avg_aln_len_of_all_filtered_pairs, 2)))
        edge_label = str(round(avg_rmsd_of_all_filtered_pairs, 2)) + ', ' + str(round(avg_aln_len_of_all_filtered_pairs, 2)) + ', ' + str(round(max_participating_r1_p, 1)) + '%'
        
        if avg_rmsd_of_all_filtered_pairs <= rmsd_threshold:
            fp.write(ordered_family_list[i] + '\t' + ordered_family_list[j] + '\t' + str(round(avg_rmsd_of_all_filtered_pairs, 2)) + '\t' + str(round(avg_aln_len_of_all_filtered_pairs, 2)) + '\t' + str(len(participating_r1_list_dict[ordered_family_list[i]][ordered_family_list[j]])) + ' / ' + str(len(families[ordered_family_list[i]])) + ' (' + str(round(fam1_participating_r1_p, 2)) + '%)' + '\t' + str(len(participating_r1_list_dict[ordered_family_list[j]][ordered_family_list[i]])) + ' / ' + str(len(families[ordered_family_list[j]])) + ' (' + str(round(fam2_participating_r1_p, 2)) + '%)' + '\t' + edge_label)
        # if best_loop_rmsd <= 1.5:
            # if pair_avg >= 5.0:
            if max_participating_r1_p >= participating_motif_instance_threshold:

                fp.write('\n')

                write_participating_motif_instances_to_file(ordered_family_list[i], participating_r1_list_dict[ordered_family_list[i]][ordered_family_list[j]], ordered_family_list[j], participating_r1_list_dict[ordered_family_list[j]][ordered_family_list[i]], output_dir)
                # selected_edge_list.append((sub1, sub2, best_loop_rmsd))
                selected_edge_list.append((ordered_family_list[i], ordered_family_list[j], 1.0))
                # edge_labels[(sub1, sub2)] =  str(round(best_loop_rmsd, 2)) + ',' + str(best_align_length) + ',' + str(round(pair_avg, 2)) + ',' + str(pair_count) + ',' + str(round(fam1_participating_r1_p, 1)) + ',' + str(round(fam2_participating_r1_p, 1)) + ',' + str(len(participating_r1_list_dict[ordered_family_list[i]][ordered_family_list[j]])) + ',' + str(len(participating_r1_list_dict[ordered_family_list[j]][ordered_family_list[i]]))
                # edge_labels[(sub1, sub2)] =  str(round(best_loop_rmsd, 2)) + ', ' + str(best_align_length) + ', ' + str(pair_count) + ', ' + str(round(max_participating_r1_p, 1)) + '%'
                # edge_labels[(sub1, sub2)] =  str(round(best_loop_rmsd, 2)) + ', ' + str(best_align_length) + ', ' + str(len(participating_r1_list_dict[ordered_family_list[i]][ordered_family_list[j]])) + ', ' + str(len(participating_r1_list_dict[ordered_family_list[j]][ordered_family_list[i]])) + ', ' + str(round(max_participating_r1_p, 1)) + '%'
                # edge_labels[(sub1, sub2)] =  str(round(best_loop_rmsd, 2)) + ', ' + str(best_align_length) + ', ' + str(round(max_participating_r1_p, 1)) + '%'
                # edge_labels[(sub1, sub2)] =  str(round(avg_rmsd_of_all_filtered_pairs, 2)) + ', ' + str(round(avg_aln_len_of_all_filtered_pairs, 2)) + ', ' + str(len(participating_r1_list_dict[ordered_family_list[i]][ordered_family_list[j]])) + ', ' + str(len(participating_r1_list_dict[ordered_family_list[j]][ordered_family_list[i]])) + ', ' + str(round(max_participating_r1_p, 1)) + '%'
                
                # edge_labels[(ordered_family_list[i], ordered_family_list[j])] =  str(round(avg_rmsd_of_all_filtered_pairs, 2)) + ', ' + str(round(avg_aln_len_of_all_filtered_pairs, 2)) + ', ' + str(len(participating_r1_list_dict[ordered_family_list[i]][ordered_family_list[j]])) + ', ' + str(len(participating_r1_list_dict[ordered_family_list[j]][ordered_family_list[i]])) + ', ' + str(round(fam1_participating_r1_p, 1)) + '%' + ', ' + str(round(fam2_participating_r1_p, 1)) + '%'


                # edge_labels[(ordered_family_list[i], ordered_family_list[j])] =  str(round(avg_rmsd_of_all_filtered_pairs, 2)) + ', ' + str(round(avg_aln_len_of_all_filtered_pairs, 2)) + ', ' + ordered_family_list[i] + ', ' + ordered_family_list[j] + ', ' + str(round(fam1_participating_r1_p, 1)) + '%' + ', ' + str(round(fam2_participating_r1_p, 1)) + '%'
                # edge_labels[(ordered_family_list[i], ordered_family_list[j])] =  str(round(avg_rmsd_of_all_filtered_pairs, 2)) + ', ' + str(round(avg_aln_len_of_all_filtered_pairs, 2)) + ', ' + str(round(max_participating_r1_p, 1)) + '%'
                edge_labels[(ordered_family_list[i], ordered_family_list[j])] = edge_label
            else:

                fp.write('\tBelow ' + str(round(participating_motif_instance_threshold, 1)) + '% threshold\n')
                
                # selected_edge_list3.append((sub1, sub2, best_loop_rmsd))
                selected_edge_list3.append((ordered_family_list[i], ordered_family_list[j], 1.0))
                # edge_labels[(sub1, sub2)] =  str(round(best_loop_rmsd, 2)) + ',' + str(best_align_length) + ',' + str(round(pair_avg, 2)) + ',' + str(pair_count) + ',' + str(round(fam1_participating_r1_p, 1)) + ',' + str(round(fam2_participating_r1_p, 1)) + ',' + str(len(participating_r1_list_dict[ordered_family_list[i]][ordered_family_list[j]])) + ',' + str(len(participating_r1_list_dict[ordered_family_list[j]][ordered_family_list[i]]))
                # edge_labels[(sub1, sub2)] =  str(round(best_loop_rmsd, 2)) + ', ' + str(best_align_length) + ', ' + str(pair_count) + ', ' + str(round(max_participating_r1_p, 1)) + '%'
                # edge_labels[(sub1, sub2)] =  str(round(best_loop_rmsd, 2)) + ', ' + str(best_align_length) + ', ' + str(len(participating_r1_list_dict[ordered_family_list[i]][ordered_family_list[j]])) + ', ' + str(len(participating_r1_list_dict[ordered_family_list[j]][ordered_family_list[i]])) + ', ' + str(round(max_participating_r1_p, 1)) + '%'
                # edge_labels[(sub1, sub2)] =  str(round(best_loop_rmsd, 2)) + ', ' + str(best_align_length) + ', ' + str(round(max_participating_r1_p, 1)) + '%'
                # edge_labels[(sub1, sub2)] =  str(round(avg_rmsd_of_all_filtered_pairs, 2)) + ', ' + str(round(avg_aln_len_of_all_filtered_pairs, 2)) + ', ' + str(len(participating_r1_list_dict[ordered_family_list[i]][ordered_family_list[j]])) + ', ' + str(len(participating_r1_list_dict[ordered_family_list[j]][ordered_family_list[i]])) + ', ' + str(round(max_participating_r1_p, 1)) + '%'
                
                # edge_labels[(ordered_family_list[i], ordered_family_list[j])] =  str(round(avg_rmsd_of_all_filtered_pairs, 2)) + ', ' + str(round(avg_aln_len_of_all_filtered_pairs, 2)) + ', ' + str(len(participating_r1_list_dict[ordered_family_list[i]][ordered_family_list[j]])) + ', ' + str(len(participating_r1_list_dict[ordered_family_list[j]][ordered_family_list[i]])) + ', ' + str(round(fam1_participating_r1_p, 1)) + '%' + ', ' + str(round(fam2_participating_r1_p, 1)) + '%'

                
                # edge_labels[(ordered_family_list[i], ordered_family_list[j])] =  str(round(avg_rmsd_of_all_filtered_pairs, 2)) + ', ' + str(round(avg_aln_len_of_all_filtered_pairs, 2)) + ', ' + ordered_family_list[i] + ', ' + ordered_family_list[j] + ', ' + str(round(fam1_participating_r1_p, 1)) + '%' + ', ' + str(round(fam2_participating_r1_p, 1)) + '%'
                # edge_labels[(ordered_family_list[i], ordered_family_list[j])] =  str(round(avg_rmsd_of_all_filtered_pairs, 2)) + ', ' + str(round(avg_aln_len_of_all_filtered_pairs, 2)) + ', ' + str(round(max_participating_r1_p, 1)) + '%'
                edge_labels[(ordered_family_list[i], ordered_family_list[j])] = edge_label
        else:
            # selected_edge_list2.append((sub1, sub2, best_loop_rmsd))
            selected_edge_list2.append((ordered_family_list[i], ordered_family_list[j], 1.0))
        # edge_colors.append('black')
            # edge_labels[(sub1, sub2)] =  str(round(best_loop_rmsd, 2)) + ', ' + str(best_align_length)

        # fp.write('\n')

    fp.close()

    label_options = {
    'font_size': 14, # 10
    # 'pos':nx.spring_layout(G, pos=fixed_positions, fixed=fixed_nodes),
    'pos': fixed_positions,
    'edge_labels': edge_labels,
    'label_pos': 0.65,
    }

    graph_options = {
    'node_color': node_colors,
    'edge_color': edge_colors,
    'edgelist': edge_list,
    # 'node_size': [len(v) * 300 for v in G.nodes()],
    'node_size': node_sizes,
    'node_shape': 'o', # s, o
    'font_size': 14, # 10
    'width': 1,
    # 'arrows': False,
    # 'arrowstyle': '->',
    # 'arrowsize': 12,
    # 'pos':nx.spring_layout(G, pos=fixed_positions),
    'pos':fixed_positions,
    'with_labels': True,    #node label
    'alpha': 0.75,
    # 'edge_labels': edge_labels,
    }
    
    nx.draw(G, **graph_options)

    # node_options = {
    # 'node_color': node_colors,
    # 'node_size': node_sizes,
    # 'node_shape': 'o', # s, o
    # 'font_size': 20, # 10
    # 'width': 1,
    # 'with_labels': True,    #node label
    # }
    # edge_options = {
    # 'width': 3,
    # 'edge_color': 'gray',
    # }

    nx.draw_networkx_edges(G, pos=fixed_positions, edgelist=selected_edge_list, width=3, edge_color='gray')
    # nx.draw_networkx_edges(G, pos=fixed_positions, style='dotted', edgelist=selected_edge_list2, width=3, edge_color='gray')
    nx.draw_networkx_edges(G, pos=fixed_positions, style='dotted', edgelist=selected_edge_list3, width=3, edge_color='darkgray')
    nx.draw_networkx_edge_labels(G, **label_options)

    # nx.draw_networkx_edges(G, pos=fixed_positions, edgelist=selected_edge_list, **edge_options)
    # nx.draw_networkx_edges(G, pos=fixed_positions, edgelist=selected_edge_list2, style='dotted', **edge_options)
    # nx.draw_networkx_nodes(G, pos=fixed_positions, **node_options)

    plt.savefig(filename, format="PNG")
    plt_fig.clf()

def get_family_colors(max_component_count = 100):
    # PyMol and pyplot common colors
    # (pyplot order) 'red','green','blue','brown','firebrick','salmon','darksalmon','chocolate','orange','wheat','olive','yellow','limegreen','aquamarine','teal','cyan','lightblue','skyblue','violet','purple','magenta','pink','lightpink'
    # (sorted) 'aquamarine','blue','brown','chocolate','cyan','darksalmon','firebrick','green','lightblue','lightpink','limegreen','magenta','olive','orange','pink','purple','red','salmon','skyblue','teal','violet','wheat','yellow'
    subfamily_colors_by_name = ['tomato', 'lightgreen', 'cornflowerblue', 'cyan', 'lightpink', 'wheat', 'orange', 'lightblue', 'violet', 'darksalmon', 'yellowgreen', 'teal', 'olive', 'yellow', 'purple', 'brown', 'magenta', 'red', 'green', 'blue']
    subfamily_colors_dict = {'red': [1.0, 0.0, 0.0],
                            'green': [0.0, 1.0, 0.0],
                            'blue': [0.0, 0.0, 1.0],
                            'tomato': [1.0, 0.39, 0.28], 
                            'lightgreen': [0.56, 0.93, 0.56], 
                            'cornflowerblue': [0.39, 0.58, 0.93], 
                            'cyan': [0.0, 1.0, 1.0], 
                            'brown': [0.65, 0.32, 0.17], 
                            'lightpink': [1.00, 0.75, 0.87],                             
                            'purple': [0.75, 0.00, 0.75],
                            'wheat': [0.99, 0.82, 0.65],
                            'teal': [0.00, 0.75, 0.75],
                            'orange': [1.0, 0.5, 0.0],
                            'lightblue': [0.75, 0.75, 1.0],                    
                            'violet': [1.0, 0.5, 1.0],
                            'magenta': [1.0, 0.0, 1.0], 
                            'darksalmon': [0.73, 0.55, 0.52],
                            'olive': [0.77, 0.70, 0.00],
                            'yellow': [1.0, 1.0, 0.0],
                            'yellowgreen': [0.60, 0.80, 0.20]}

    subfamily_colors_rgb = []
    for color_name in subfamily_colors_by_name:
        subfamily_colors_rgb.append(subfamily_colors_dict[color_name])

    random.seed(3)

    for i in range (max_component_count * 10):
        if len(subfamily_colors_rgb) >= max_component_count:
            break

        rand_ind1 = random.randint(0,len(subfamily_colors_rgb) - 1)
        rand_ind2 = random.randint(0,len(subfamily_colors_rgb) - 1)
        rand_ind3 = random.randint(0,len(subfamily_colors_rgb) - 1)

        # mix 3 random existing color to generate a new color
        new_color = [(x + y + z)/3.0 for x, y, z in zip(subfamily_colors_rgb[rand_ind1], subfamily_colors_rgb[rand_ind2], subfamily_colors_rgb[rand_ind3])]

        if not similar_to_existing_color(subfamily_colors_rgb, new_color):
            subfamily_colors_rgb.append(new_color)
            subfamily_colors_by_name.append(new_color)       

    # If need more class, assign random colors
    for i in range (max_component_count * 100):
        if len(subfamily_colors_rgb) >= max_component_count:
            break

        new_color = [random.random(), random.random(), random.random()]
        if not similar_to_existing_color(subfamily_colors_rgb, new_color):
            subfamily_colors_rgb.append(new_color)
            subfamily_colors_by_name.append(new_color)       

    # Generate tuple to make it compatible with pyplot (as list is needed for pymol)
    subfamily_colors_tuple = []
    for item in subfamily_colors_by_name:
        if type(item) == 'list' and len(item) == 3:
            subfamily_colors_tuple.append((item[0], item[1], item[2]))
        else:
            subfamily_colors_tuple.append(item)

    return subfamily_colors_by_name, subfamily_colors_tuple

def draw_figure_with_mean_rmsd(families, alignment_tool, subfamilies, inter_family_adj_mat_avg):
    plt_fig = plt.figure(frameon = False)
    plt_fig = plt.figure()
    plt_fig.set_size_inches(12, 12)

    filename = 'test_mean.png'
    # print(subfamily_count)
    # G = nx.DiGraph(directed=True)
    G = nx.Graph()
    circular_positions = get_coordinates_in_circle(len(families), 0.0, 0.0, 2.0)

    priority_order = ["GNRA", "GNAA", "GNGA", "KT", "rKT", "SR", "CL", "EL", "HT", "TS", "TR", "L1C", "TL", "RS"]
    colors, colors_tuple = get_family_colors(len(families))
    # print(colors)

    fixed_positions = {}
    fixed_nodes = []
    node_colors = []
    node_sizes = []

    edge_list = []
    edge_colors = []

    # fixed_positions['X'] = (0.0, 0.0)
    # fixed_nodes.append('X')
    # node_colors.append('salmon')

    ordered_family_list = []
    for family_id in priority_order:
        if family_id in families:
            ordered_family_list.append(family_id)


    print(ordered_family_list)
    print(family_id)
    print(ordered_family_list.index(family_id))
    
    for i, family_id in enumerate(ordered_family_list):
        # node = known_motif_shortcode[family_id.lower()]
        node = family_id
        fixed_positions[node] = circular_positions[i]
        fixed_nodes.append(node)
        node_colors.append(colors[ordered_family_list.index(family_id)])
        node_sizes.append(3500)

    
    print(fixed_positions)
    print(fixed_nodes)

    G.add_nodes_from(fixed_nodes)

    edge_source='family'
    selected_edges, other_edges = KruskalMST(inter_family_adj_mat_avg, ordered_family_list, edge_source)
    selected_edges += other_edges
    selected_edge_list = []
    selected_edge_list2 = []
    selected_edge_list3 = []
    selected_edge_list4 = []
    edge_labels = {}
    for edge in selected_edges:
        # (i, j, best_r1, best_r2, best_loop_rmsd, best_align_length, mean_rmsd, pair_avg, pair_count) = edge
        (i, j, mean_rmsd) = edge

        if mean_rmsd <= 2.0:
            selected_edge_list.append((ordered_family_list[i], ordered_family_list[j], mean_rmsd))
            edge_labels[(ordered_family_list[i], ordered_family_list[j])] = str(round(mean_rmsd, 2))
        elif mean_rmsd <= 2.5:
            selected_edge_list2.append((ordered_family_list[i], ordered_family_list[j], mean_rmsd))
            edge_labels[(ordered_family_list[i], ordered_family_list[j])] = str(round(mean_rmsd, 2))
        # edge_colors.append('black')
        # edge_labels[(ordered_family_list[i], ordered_family_list[j])] = str(round(mean_rmsd, 2))

    # for edge in other_edges:
    #     (i, j, best_r1, best_r2, best_loop_rmsd, best_align_length, mean_rmsd) = edge

    #     if mean_rmsd < 3.0:
    #         selected_edge_list3.append((ordered_family_list[i], ordered_family_list[j], mean_rmsd))
    #     else:
    #         selected_edge_list4.append((ordered_family_list[i], ordered_family_list[j], mean_rmsd))
    #     # edge_colors.append('black')
    #     edge_labels[(ordered_family_list[i], ordered_family_list[j])] = str(round(mean_rmsd, 2))


    # for a_node_in_circle in outgoing_nodes_from_circular_points:
    #     points = get_n_fixed_points_next_to_circle(len(outgoing_nodes_from_circular_points[a_node_in_circle]), outside_angle, fixed_positions[a_node_in_circle], fixed_positions[get_next_node(a_node_in_circle, circular_nodes)])
    #     for i, node in enumerate(outgoing_nodes_from_circular_points[a_node_in_circle]):
    #         fixed_positions[node] = points[i]
    #         fixed_nodes.append(node)

    label_options = {
    'font_size': 16, # 10
    # 'pos':nx.spring_layout(G, pos=fixed_positions, fixed=fixed_nodes),
    'pos': fixed_positions,
    'edge_labels': edge_labels,
    }

    graph_options = {
    'node_color': node_colors,
    'edge_color': edge_colors,
    'edgelist': edge_list,
    # 'node_size': [len(v) * 300 for v in G.nodes()],
    'node_size': node_sizes,
    'node_shape': 'o', # s, o
    'font_size': 20, # 10
    'width': 1,
    # 'arrows': False,
    # 'arrowstyle': '->',
    # 'arrowsize': 12,
    # 'pos':nx.spring_layout(G, pos=fixed_positions),
    'pos':fixed_positions,
    'with_labels': True,    #node label
    # 'edge_labels': edge_labels,
    }
    nx.draw_networkx_edge_labels(G, **label_options)
    nx.draw(G, **graph_options)


    nx.draw_networkx_edges(G, pos=fixed_positions, edgelist=selected_edge_list, width=3, edge_color='gray')
    nx.draw_networkx_edges(G, pos=fixed_positions, style='dotted', edgelist=selected_edge_list2, width=3, edge_color='gray')
    # nx.draw_networkx_edges(G, pos=fixed_positions, edgelist=selected_edge_list3, width=3, edge_color='lightgray')
    # nx.draw_networkx_edges(G, pos=fixed_positions, style='dotted', edgelist=selected_edge_list4, width=3, edge_color='lightgray')
    plt.savefig(filename, format="PNG")
    plt_fig.clf()

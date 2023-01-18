import sys
sys.path.append('src/scripts')
from utils import *

def get_family_or_subfamily_id(loop, annotated_dict):
    for ann_id in annotated_dict:
        loops = annotated_dict[ann_id]
        if str(strToNode(loop)) in loops:
            return ann_id

    return ''

def csv_to_list(lines):
    list_of_lists = []
    for line in lines:
        pieces = line.strip().split(',')
        list_of_lists.append(list(map(lambda x: x.strip(), pieces)))

    return list_of_lists

def sort_all_loop_segments(families):
    new_families = {}
    for fam_id in families:
        new_families[fam_id] = []
        loops = families[fam_id]
        for loop in loops:
            new_families[fam_id].append(str(strToNode(loop)))
    return new_families

def main():
    # input_file_name = 'output/IL-all/subfamily_cluster.csv'
    input_file_name = 'output/IL_TMalign/subfamily_cluster.csv'
    annotated_motif_family_file_name = 'data/supercluster_IL.in'
    annotated_motif_subfamily_file_name = 'data/subfamily_cluster_IL.csv'

    mixed_families = {}
    fp = open(input_file_name)
    lines = fp.readlines()
    loop_list = csv_to_list(lines)
    loop_cnt = 0
    for item in loop_list:
        mixed_families[item[0]] = item[1:]
        loop_cnt += len(item[1:])
    print(loop_cnt)
    fp.close()

    annotated_families = {}
    fp = open(annotated_motif_family_file_name)
    lines = fp.readlines()
    loop_list = csv_to_list(lines)
    loop_cnt = 0
    for item in loop_list:
        annotated_families[item[0]] = item[1:]
        loop_cnt += len(item[1:])
    print(loop_cnt)
    fp.close()

    annotated_subfamilies = {}
    fp = open(annotated_motif_subfamily_file_name)
    lines = fp.readlines()
    loop_list = csv_to_list(lines)
    loop_cnt = 0
    for item in loop_list:
        annotated_subfamilies[item[0]] = item[1:]
        loop_cnt += len(item[1:])
    print(loop_cnt)
    fp.close()

    mixed_families = sort_all_loop_segments(mixed_families)
    annotated_families = sort_all_loop_segments(annotated_families)
    annotated_subfamilies = sort_all_loop_segments(annotated_subfamilies)

    # sys.exit()

    # mixed_fam_stat = {}
    for mixed_fam_id in mixed_families:
        family_stat = {}
        subfamily_stat = {}
        loops = mixed_families[mixed_fam_id]
        total_loop_cnt = len(loops)
        for loop in loops:
            fam_id = get_family_or_subfamily_id(loop, annotated_families)
            subfam_id = get_family_or_subfamily_id(loop, annotated_subfamilies)

            if fam_id not in family_stat:
                family_stat[fam_id] = 0
            family_stat[fam_id] += 1

            if subfam_id not in subfamily_stat:
                subfamily_stat[subfam_id] = 0
            subfamily_stat[subfam_id] += 1

        # mixed_fam_stat[mixed_fam_id] = (family_stat, subfamily_stat)
        print('==========================================================')
        print(mixed_fam_id + ' (' + str(total_loop_cnt) + '):')
        print('family_stat:')
        for fam_id in family_stat:
            p = family_stat[fam_id] * 100 / total_loop_cnt
            print(fam_id + ': ' + str(family_stat[fam_id]) + ' (' + str(round(p, 2)) + '%)')
        print('\nsubfamily_stat:')
        for subfam_id in subfamily_stat:
            p = subfamily_stat[subfam_id] * 100 / total_loop_cnt
            print(subfam_id + ': ' + str(subfamily_stat[subfam_id]) + ' (' + str(round(p, 2)) + '%)')
        print('==========================================================\n\n')

        # sys.exit()

if __name__ == '__main__':
    main()
#!/bin/env python
# coding=utf8

u"""
Created by Zhang dan at 2021.05.18

All the class and function about command line parameters

Last modified at 2021.05.18

can only run at 4103
"""

import os
import re
import io
import gzip
import subprocess
import sys
import glob
import celescope
import click
import pandas as pd
from collections import defaultdict, Counter
from itertools import combinations, permutations, islice
from xopen import xopen
from celescope.tools.utils import format_number, seq_ranges, read_fasta, genDict
from celescope.tools.report import reporter
from celescope.tools.__init__ import __PATTERN_DICT__
from celescope.tools.Chemistry import Chemistry
from celescope.tools.utils import add_log

def ord2chr(q, offset=33):
    return chr(int(q) + offset)


# 生成错配字典
def generate_mis_seq(seq, n=1, bases='ACGTN'):
    # 以随机bases中的碱基替换seq中的n个位置，产生的错配字典
    # 返回字典，错配序列为key，
    # (正确序列，错配碱基数目，错配碱基位置，原始碱基，新碱基)组成的元组
    # 作为字典的值

    length = len(seq)
    assert length >= n, "err number should not be larger than sequence length!"
    res = {}
    seq_arr = list(seq)
    pos_group = list(combinations(range(0, length), n))
    bases_group = list(permutations(bases, n))

    for g in pos_group:
        for b in bases_group:
            seq_tmp = seq_arr[:]
            mis_num = n
            raw_tmp = []
            for i in range(n):
                raw_base = seq_tmp[g[i]]
                new_base = b[i]

                if raw_base == new_base:
                    mis_num -= 1

                raw_tmp.append(raw_base)
                seq_tmp[g[i]] = new_base

            if mis_num != 0:
                res[''.join(seq_tmp)] = (seq, mis_num, ','.join(
                    [str(i) for i in g]), ','.join(raw_tmp), ','.join(b))
    return(res)


def generate_seq_dict(seqlist, n=1):
    seq_dict = {}
    with open(seqlist, 'r') as fh:
        for seq in fh:
            seq = seq.strip()
            if seq == '':
                continue
            seq_dict[seq] = (seq, 0, -1, 'X', 'X')
            for k, v in generate_mis_seq(seq, n).items():
                # duplicate key
                if k in seq_dict:
                    generate_seq_dict.logger.warning('barcode %s, %s\n%s, %s' %
                                                     (v, k, seq_dict[k], k))
                else:
                    seq_dict[k] = v
    return seq_dict


def parse_pattern(pattern):
    # 解析接头结构，返回接头结构字典
    # key: 字母表示的接头, value: 碱基区间列表
    # eg.: C8L10C8L10C8U8T30
    # defaultdict(<type 'list'>:
    # {'C': [[0, 8], [18, 26], [36, 44]], 'U': [[44, 52]], 'L': [[8, 18], [26, 36]], 'T': [[52, 82]]})
    pattern_dict = defaultdict(list)
    p = re.compile(r'([CLUNT])(\d+)')
    tmp = p.findall(pattern)
    if not tmp:
        parse_pattern.logger.error(f'Invalid pattern: {pattern}')
        sys.exit()
    start = 0
    for item in tmp:
        end = start + int(item[1])
        pattern_dict[item[0]].append([start, end])
        start = end
    return pattern_dict


def get_scope_bc(bctype):
    import celescope
    root_path = os.path.dirname(celescope.__file__)
    linker_f = glob.glob(f'{root_path}/data/chemistry/{bctype}/linker*')[0]
    whitelist_f = f'{root_path}/data/chemistry/{bctype}/bclist'
    return linker_f, whitelist_f


def read_fastq(f):
    """
    Return tuples: (name, sequence, qualities).
    qualities is a string and it contains the unmodified, encoded qualities.
    """
    i = 3
    for i, line in enumerate(f):
        if i % 4 == 0:
            assert line.startswith('@'), ("Line {0} in FASTQ file is expected to start with '@', "
                                          "but found {1!r}".format(i + 1, line[:10]))
            name = line.strip()[1:]
        elif i % 4 == 1:
            sequence = line.strip()
        elif i % 4 == 2:
            line = line.strip()
            assert line.startswith('+'), ("Line {0} in FASTQ file is expected to start with '+', "
                                          "but found {1!r}".format(i + 1, line[:10]))
        elif i % 4 == 3:
            qualities = line.rstrip('\n\r')
            yield name, sequence, qualities
    if i % 4 != 3:
        raise Exception("FASTQ file ended prematurely")


def low_qual(quals, minQ='/', num=2):
    # print(ord('/')-33)           14
    return True if len([q for q in quals if q < minQ]) > num else False


def no_polyT(seq, strictT=0, minT=10):
    # strictT requires the first nth position to be T
    if seq[:strictT] != 'T' * strictT or seq.count('T') < minT:
        return True
    else:
        return False


def no_barcode(seq_arr, mis_dict, err_tolerance=1,barcode_corrected_num =0):
    tmp = [mis_dict[seq][0:2] if seq in mis_dict else (
        'X', 100) for seq in seq_arr]
    err = sum([t[1] for t in tmp])
    if err > err_tolerance:
        return True
    else:
        if err > 0:
            barcode_corrected_num += 1
            return ''.join([t[0] for t in tmp])
        else:
            return "correct"


def no_linker(seq, mis_dict, err_tolerance=2, linker_corrected_num = 0):
    tmp = mis_dict[seq][0:2] if seq in mis_dict else (
        'X', 100) 
    err = tmp[1]
    if err > err_tolerance:
        return True
    else:
        if err > 0:
            linker_corrected_num += 1
            return tmp[0]
        else:
            return "correct"


def check_seq(seq_file, pattern_dict, seq_abbr):
    length = 0
    for item in pattern_dict[seq_abbr]:
        start = item[0]
        end = item[1]
        length += end - start
    with open(seq_file, 'r') as fh:
        for seq in fh:
            seq = seq.strip()
            if seq == '':
                continue
            if len(seq) != length:
                raise Exception(
                    f'length of L in pattern ({length}) do not equal to length in {seq_file} ({len(seq)}) !')


@click.command(
    short_help="get new cell barcode including linker sequence from singleron singlecell NGS data"
)
@click.option(
    "-s","--sample", type=str, required=True,
    help="sample name."
)
@click.option(
    "-o","--outdir", type=click.Path(), required=True,
    help="Path to outdir directory."
)
@click.option(
    "-f1", "--fq1", type=str, required=True,
    help="Required, FASTQ R1 reads. Multiple FASTQ files are seperated by comma."
)
@click.option(
    "-f2", "--fq2", type=str, required=True,
    help="Required, FASTQ R2 reads. Multiple FASTQ files are seperated by comma."
)
def barcode(sample,outdir,fq1,fq2,nopolyT = False,noLinker = False,probe_file = None,allowNoPolyT = False,allowNoLinker = False,lowQual = 0, lowNum = 2,chemistry = "scopeV2.2.1"):
    # init
    fq1_list = fq1.split(",")
    fq2_list = fq2.split(",")

    # check dir
    if not os.path.exists(outdir):
        os.system('mkdir -p %s' % outdir)

    # get chemistry
    if chemistry == 'auto':
            ch = Chemistry(fq1)
            chemistry = ch.check_chemistry()
    else:
        chemistry = chemistry

    # get linker and whitelist
    bc_pattern = __PATTERN_DICT__[chemistry]
    if (bc_pattern):
        (linker, whitelist) = get_scope_bc(chemistry)
    else:
        bc_pattern = pattern
        linker = linker
        whitelist = whitelist
    if (not linker) or (not whitelist) or (not bc_pattern):
        raise Exception("invalid chemistry or [pattern,linker,whitelist]")
        
    # parse pattern to dict, C8L10C8L10C8U8
    # defaultdict(<type 'list'>, {'C': [[0, 8], [18, 26], [36, 44]], 'U':
    # [[44, 52]], 'L': [[8, 18], [26, 36]]})
    pattern_dict = parse_pattern(bc_pattern)

    # check linker
    check_seq(linker, pattern_dict, "L")

    bool_T = True if 'T' in pattern_dict else False
    bool_L = True if 'L' in pattern_dict else False

    C_len = sum([item[1] - item[0] for item in pattern_dict['C']])

    barcode_qual_Counter = Counter()
    umi_qual_Counter = Counter()
    C_U_base_Counter = Counter()
    lowQual = ord2chr(lowQual)

    # generate list with mismatch 1, substitute one base in raw sequence with
    # A,T,C,G
    barcode_dict = generate_seq_dict(whitelist, n=1)
    linker_dict = generate_seq_dict(linker, n=2)

    # prepare
    out_fq2 = outdir + '/' + sample + '_readid_bc.txt'
    fh3 = xopen(out_fq2, 'w')

    (total_num, clean_num, no_polyT_num, lowQual_num,
        no_linker_num, no_barcode_num) = (0, 0, 0, 0, 0, 0)
    Barcode_dict = defaultdict(int)

    if nopolyT:
        fh1_without_polyT = xopen(outdir + '/noPolyT_1.fq', 'w')
        fh2_without_polyT = xopen(outdir + '/noPolyT_2.fq', 'w')

    if noLinker:
        fh1_without_linker = xopen(outdir + '/noLinker_1.fq', 'w')
        fh2_without_linker = xopen(outdir + '/noLinker_2.fq', 'w')

    bool_probe = False
    if probe_file and probe_file != 'None':
        bool_probe = True
        count_dic = genDict(dim=3)
        valid_count_dic = genDict(dim=2)
        probe_dic = read_fasta(probe_file)
        reads_without_probe = 0

    # process
    fq_number = len(fq1_list)
    if fq_number != len(fq2_list):
        raise Exception('fastq1 and fastq2 do not have same file number!')
    for i in range(fq_number):
        fh1 = xopen(fq1_list[i])
        fh2 = xopen(fq2_list[i])
        g1 = read_fastq(fh1)
        g2 = read_fastq(fh2)

        while True:
            try:
                (header1, seq1, qual1) = next(g1)
                (header2, seq2, qual2) = next(g2)
            except BaseException:
                break
            if total_num > 0 and total_num % 1000000 == 0:
                print("finish!")
            total_num += 1

            # polyT filter
            if bool_T and (not allowNoPolyT):
                polyT = seq_ranges(seq1, pattern_dict['T'])
                if no_polyT(polyT):
                    no_polyT_num += 1
                    if nopolyT:
                        fh1_without_polyT.write(
                            '@%s\n%s\n+\n%s\n' % (header1, seq1, qual1))
                        fh2_without_polyT.write(
                            '@%s\n%s\n+\n%s\n' % (header2, seq2, qual2))
                    continue

            # lowQual filter
            C_U_quals_ascii = seq_ranges(
                qual1, pattern_dict['C'] + pattern_dict['U'])
            # C_U_quals_ord = [ord(q) - 33 for q in C_U_quals_ascii]
            if low_qual(C_U_quals_ascii, lowQual, lowNum):
                lowQual_num += 1
                continue

            # linker filter
            if bool_L and (not allowNoLinker):
                linker_arr = [seq_ranges(seq1, [i]) for i in pattern_dict['L']]
                raw_linker = ''.join(linker_arr)
                res_linker = no_linker(raw_linker, linker_dict)
                if res_linker is True:
                    no_linker_num += 1
                    continue
                elif res_linker == "correct":
                    linker = raw_linker
                else:
                    linker = res_linker

            # barcode filter
            barcode_arr = [seq_ranges(seq1, [i]) for i in pattern_dict['C']]
            raw_cb = ''.join(barcode_arr)
            res = no_barcode(barcode_arr, barcode_dict)

            if res is True:
                no_barcode_num += 1
                continue
            elif res == "correct":
                cb = raw_cb
            else:
                cb = res


            umi = seq_ranges(seq1, pattern_dict['U'])
            Barcode_dict[cb] += 1
            clean_num += 1
            read_name_probe = 'None'

            if bool_probe:
                # valid count
                valid_count_dic[cb][umi] += 1

                # output probe UMi and read count
                find_probe = False
                for probe_name in probe_dic:
                    probe_seq = probe_dic[probe_name]
                    probe_seq = probe_seq.upper()
                    if seq1.find(probe_seq) != -1:
                        count_dic[probe_name][cb][umi] += 1
                        read_name_probe = probe_name
                        find_probe = True
                        break

                if not find_probe:
                    reads_without_probe += 1

            barcode_qual_Counter.update(C_U_quals_ascii[:C_len])
            umi_qual_Counter.update(C_U_quals_ascii[C_len:])
            C_U_base_Counter.update(raw_cb + umi)

            # readid and cell barcode assignment
            cb1 = cb[0:8]
            cb2 = cb[8:16]
            cb3 = cb[16:24]
            lk1 = linker[0:16]
            lk2 = linker[16:32]
            lk3 = linker[32:33]
            new_cb = cb1 + lk1 + cb2 + lk2 + cb3 + lk3
            bc_id = cb1 + cb2 + cb3 + "_" + umi + "_" + read_name_probe + "_" + str(total_num)
            fh3.write(f'@{bc_id}\n{new_cb}\n')
    fh3.close()

if __name__ == '__main__':
    barcode()

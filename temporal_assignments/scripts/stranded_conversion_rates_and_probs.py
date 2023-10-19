import os
import argparse
import sys
import random
from joblib import Parallel, delayed
import pysam
from gtfparse import read_gtf
import pandas as pd
import itertools
from collections import defaultdict
import numpy as np
from tqdm import tqdm
from scipy import sparse
from scipy.stats import binom

def chunks(l, n):
    c = itertools.count()
    return [list(it) for _, it in itertools.groupby(l, lambda x: next(c)//n)]

def conversionCounting(bamfile, contig, start, end, strand, geneID, qcutoff, subsample, barcodes):
    bam = pysam.AlignmentFile(bamfile, 'rb')
    output = {geneID: {}}
    for read in bam.fetch(contig, start, end):

        if subsample < 1.0 and random.random() > subsample:
            continue
        
        if read.has_tag('GE'):
            gene = read.get_tag('GE')
        else:
            if read.has_tag('GI'):
                gene = read.get_tag('GI')
            else:
                continue

        if gene == geneID:
            cell = read.get_tag('BC')

            if barcodes and not cell in barcodes:
                continue
            
            if not cell in output[geneID]:
                output[geneID][cell] = {'+': {'coverage': {'A': 0, 'C': 0, 'G': 0, 'T': 0}, 'conversions': {'aC': 0, 'aG': 0, 'aT': 0, 'cA': 0, 'cG': 0, 'cT': 0, 'gA': 0, 'gC': 0, 'gT': 0, 'tA': 0, 'tC': 0, 'tG': 0}},
                                        '-': {'coverage': {'A': 0, 'C': 0, 'G': 0, 'T': 0}, 'conversions': {'aC': 0, 'aG': 0, 'aT': 0, 'cA': 0, 'cG': 0, 'cT': 0, 'gA': 0, 'gC': 0, 'gT': 0, 'tA': 0, 'tC': 0, 'tG': 0}}, 'pairs': []}
            pair_t = 0
            pair_tC = 0
            for pos in read.get_aligned_pairs(with_seq=True):
                base = pos[2]
                if base is not None:
                    if base.isupper() and base != 'N':
                        output[geneID][cell][strand]['coverage'][base] += 1
                        if (strand == '+' and base == 'T') or (strand == '-' and base == 'A'):
                            pair_t += 1
                    elif base.islower() and base != 'n':
                        if read.query_qualities[pos[0]] > qcutoff:
                            output[geneID][cell][strand]['coverage'][base.upper()] += 1
                            conv = (
                                ''.join([base, read.query_sequence[pos[0]]]))
                            output[geneID][cell][strand]['conversions'][conv] += 1
                            if (strand == '+' and conv == 'tC') or (strand == '-' and base == 'aG'):
                                pair_tC += 1
                                pair_t += 1
            output[geneID][cell]['pairs'].append((pair_t, pair_tC))
    return(output)


def conversionCounting_PE(bamfile, contig, start, end, strand, geneID, qcutoff, subsample=1.0, barcodes=False):
    bam = pysam.AlignmentFile(bamfile, 'rb')
    output = {geneID: {}}
    for reads in read_pair_generator(bam, contig, start, end):
        if subsample < 1.0 and random.random() > subsample:
            continue

        read1, read2 = reads
        cell = read1.get_tag('BC')
        if barcodes and not cell in barcodes:
            continue
        if not cell in output[geneID]:
            output[geneID][cell] = {'+': {'coverage': {'A': 0, 'C': 0, 'G': 0, 'T': 0}, 'conversions': {'aC': 0, 'aG': 0, 'aT': 0, 'cA': 0, 'cG': 0, 'cT': 0, 'gA': 0, 'gC': 0, 'gT': 0, 'tA': 0, 'tC': 0, 'tG': 0}},
                                    '-': {'coverage': {'A': 0, 'C': 0, 'G': 0, 'T': 0}, 'conversions': {'aC': 0, 'aG': 0, 'aT': 0, 'cA': 0, 'cG': 0, 'cT': 0, 'gA': 0, 'gC': 0, 'gT': 0, 'tA': 0, 'tC': 0, 'tG': 0}}, 'pairs': []}
        pair_t = 0
        pair_tC = 0
        for read in (read1, read2):
            for pos in read.get_aligned_pairs(with_seq=True):
                base = pos[2]
                if base is not None:
                    if base.isupper() and base != 'N':
                        output[geneID][cell][strand]['coverage'][base] += 1
                        if (strand == '+' and base == 'T') or (strand == '-' and base == 'A'):
                            pair_t += 1
                    elif base.islower() and base != 'n':
                        if read.query_qualities[pos[0]] > qcutoff:
                            output[geneID][cell][strand]['coverage'][base.upper()] += 1
                            conv = (
                                ''.join([base, read.query_sequence[pos[0]]]))
                            output[geneID][cell][strand]['conversions'][conv] += 1
                            if (strand == '+' and conv == 'tC') or (strand == '-' and base == 'aG'):
                                pair_tC += 1
                                pair_t += 1
        output[geneID][cell]['pairs'].append((pair_t, pair_tC))
    return(output)


def read_pair_generator(bam, contig, start, end):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(contig, start, end):
        if not read.has_tag('GE'):
            if not read.has_tag('GI'):
                continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


def conversionCountingGrouped(bamfile, contigList, startList, endList, strandList, geneIDList, l, seq_layout, prefix='conv_',
                              subsample=1.0, barcodes=False):
    
    listChunksContig, listChunksStart, listChunksEnd, listChunksStrand, listChunksGeneID = chunks(
        contigList, l), chunks(startList, l), chunks(endList, l), chunks(strandList, l), chunks(geneIDList, l)
    output = {}
    print(f"Counting conversions in {len(listChunksContig)} chunks")
    for j in tqdm(range(0, len(listChunksContig))):
        if seq_layout == "single_end":
            cr = Parallel(n_jobs=int(o.threads), verbose=3, backend='loky')(delayed(conversionCounting)(o.bam, seqname, start, end, strand, geneID, o.quality, subsample, barcodes)
                                                                            for seqname, start, end, strand, geneID in zip(listChunksContig[j], listChunksStart[j], listChunksEnd[j], listChunksStrand[j], listChunksGeneID[j]))
        else:
            cr = Parallel(n_jobs=int(o.threads), verbose=3, backend='loky')(delayed(conversionCounting_PE)(o.bam, seqname, start, end, strand, geneID, o.quality, subsample, barcodes)
                                                                            for seqname, start, end, strand, geneID in zip(listChunksContig[j], listChunksStart[j], listChunksEnd[j], listChunksStrand[j], listChunksGeneID[j]))
        for i in range(0, len(cr)):
            for geneID in cr[i].keys():
                for cell in cr[i][geneID].keys():
                    if not cell in output:
                        output[cell] = {'+': {'coverage': {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0},
                                              'conversions': {'aC': 0, 'aG': 0, 'aT': 0, 'aN': 0, 'cA': 0, 'cG': 0, 'cT': 0, 'cN': 0, 'gA': 0, 'gC': 0, 'gT': 0, 'gN': 0, 'tA': 0, 'tC': 0, 'tG': 0, 'tN': 0, 'nA': 0, 'nC': 0, 'nG': 0, 'nT': 0, 'nN': 0}},
                                        '-': {'coverage': {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0},
                                              'conversions': {'aC': 0, 'aG': 0, 'aT': 0, 'aN': 0, 'cA': 0, 'cG': 0, 'cT': 0, 'cN': 0, 'gA': 0, 'gC': 0, 'gT': 0, 'gN': 0, 'tA': 0, 'tC': 0, 'tG': 0, 'tN': 0, 'nA': 0, 'nC': 0, 'nG': 0, 'nT': 0, 'nN': 0}},
                                        'sparse_matrix': sparse.lil_matrix((200, 200), dtype=np.int32)}
                    for strand in ('-', '+'):
                        for stat in cr[i][geneID][cell][strand].keys():
                            for count in cr[i][geneID][cell][strand][stat].keys():
                                output[cell][strand][stat][count] += cr[i][geneID][cell][strand][stat][count]

                    # add read pair related conversion patterns
                    with open('./data/%sdetails.txt'%prefix,'a+') as yetmoreoutput:
                        yetmoreoutput.write("%s\t%s\t%s\t%s\n" % (cell, geneID, ",".join([str(val[1]) for val in cr[i][geneID][cell]['pairs']]),
                                                                  ",".join([str(val[0]) for val in cr[i][geneID][cell]['pairs']])))
                    
                    
                    for item in (cr[i][geneID][cell]['pairs']):
                        output[cell]['sparse_matrix'][item[1], item[0]] += 1

                        
        del cr
    return(output)


##############
# Calculate pc
##############

# Bisection-search for p_c
def estimateP_c(Akn, p_e, cell_id):
    l = p_e + 1e-7
    r = 1
    p_c0 = (l+r)/2.
    Akn0 = Akn.copy()  # makes a deep copy to allow in-place modifications
    p_c = p_c0
    Mkn = createMkn(Akn0, p_e)
    Akn = Akn0
    while r-l >= 10e-8:
        Akn = EstepAkn(Akn, Mkn, p_c)
        p_c_old = p_c
        p_c = MstepP_c(Akn)
        if p_c < p_c_old:
            r = p_c
        else:
            l = p_c
    return p_c, p_e, cell_id


def createMkn(Akn, p_e):  # Left out from Akn
    M = np.zeros(Akn.shape)
    for n in range(Akn.shape[1]):
        for k in range(Akn.shape[0]):
            Ekn = np.sum(Akn[(k+1):, n]) * binom.pmf(k, n, p_e)
            if Ekn > 0.01 * Akn[k, n]:
                M[k, n] = 1
    return M


def EstepAkn(Akn, Mkn, p_c):  # Alters Akn in place - modify initial step
    for k in range(Akn.shape[0]):
        for n in range(Akn.shape[1]):
            if Mkn[k, n] == 1:
                num = 0
                denom = 0
                for kp in range(Mkn.shape[0]):
                    if Mkn[kp, n] == 1:
                        num = num + binom.pmf(k, n, p_c) * Akn[kp, n]
                        denom = denom + binom.pmf(kp, n, p_c)
                Akn[k, n] = num / denom
    return Akn


def MstepP_c(Akn):
    num = 0
    denom = 0
    for k in range(Akn.shape[0]):
        for n in range(Akn.shape[1]):
            num = num + k * Akn[k, n]
            denom = denom + n * Akn[k, n]
    p_c = num / denom
    return p_c


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bam', required=True,
                        help='bam file to process.')
    parser.add_argument('-g', '--gtf', required=True,
                        help='gene annotations file.')
    parser.add_argument('-n', '--name', required=True,
                        help='output files will start with <name><...>')
    parser.add_argument('-t', '--threads', required=True,
                        help='CPU threads to use')
    parser.add_argument('-c', '--cutoff', required=False, default=100000, type = int,
                        help='Number of reads per cell required for processing cell stats.')
    parser.add_argument('-q', '--quality', required=False,
                        default=30, type=int,
                        help='Phred quality cutoff for mismatching bases')
    parser.add_argument('-s','--subsample', default=1.0, type=float,
                        help='Only process a fraction of reads defined by the subsample float')
    parser.add_argument('-w', '--barcode-file', default=False,help ='File with cell whitelisted barcodes to process')
    parser.add_argument('--pc_calc', default = 'pc', help='Shall conversion probabilites be calculated? Options: [pc/ no_pc]')
    parser.add_argument('--seq_layout', default = "paired_end", help='Sequencing layout. Options: [single_end/ paired_end]')
    o = parser.parse_args()

    if not os.path.exists("./data"):
        os.mkdir("./data")

    pc_calc = o.pc_calc
    seq_layout = o.seq_layout

    print(f"Mandatory flags set to:\npc_calc: {pc_calc}\nseq_layout: {seq_layout}\n***\n")


    if pc_calc not in ['pc', 'no_pc']:
        sys.exit("Error, pc_calc needs to be either pc or no_pc")
    
    if seq_layout not in ["single_end", "paired_end"]:
        sys.exit("Error, seq_layout needs to be single_end or paired_end")

    # parse gene annotations
    print("Reading GTF file...")
    # This returns output to stdout like this: INFO:root:Extracted GTF attributes...
    gtf = read_gtf(o.gtf)
    print("Done reading GTF file!")

    genes = gtf[gtf['feature'] == 'gene']
    starts = list(genes['start'])
    ends = list(genes['end'])
    seqnames = list(genes['seqname'])
    strands = list(genes['strand'])
    geneIDs = list(genes['gene_id'])

    print("Starting processing...\n***\n")

    if seq_layout == 'single_end':
        if o.subsample < 1.0:
            prefix = '%s.SE.%f' % (o.name,o.subsample)
        else:
            prefix = '%s.SE.' % o.name
    else:
        if o.subsample < 1.0:
            prefix = '%s.PE.%f' % (o.name,o.subsample)
        else:
            prefix = '%s.PE.' %	o.name

    if o.barcode_file:
        barcodes = []
        with open(o.barcode_file,'r') as inf:
            for line in inf:
                barcodes.append(line[:-1])
    else:
        barcodes=False


    print("Collecting conversion and gene coverage information...")

    output = conversionCountingGrouped(
        o.bam, seqnames, starts, ends, strands, geneIDs, 10000, seq_layout, prefix=prefix, subsample=o.subsample, barcodes=barcodes)
    print('Finished collecting conversion and coverage information.')

    cells, strands, convs, rates = [], [], [], []

    print("Calculating conversion rates...")
    unique_cells = set([])
    for cell in output.keys():
        if sum([sum(output[cell]['+']['coverage'].values()), sum(output[cell]['-']['coverage'].values())]) > o.cutoff:
            unique_cells.add(cell)
            for strand in ('+', '-'):
                output[cell][strand]['conversionRates'] = {}
                for conv in output[cell][strand]['conversions'].keys():
                    ref = conv[0].upper()
                    try:
                        rate = output[cell][strand]['conversions'][conv] / \
                            float(output[cell][strand]['coverage'][ref])
                    except ZeroDivisionError:
                        rate = 0
                    output[cell][strand]['conversionRates'][conv] = rate
                    cells.append(cell)
                    strands.append(strand)
                    convs.append(conv)
                    rates.append(rate)
    print('Finished calculating conversion rates.')

    print("Saving conversion rates to file...")
    df = pd.DataFrame([cells, strands, convs, rates])
    print()
    df = df.transpose()
    df.columns = ['cell', 'strand', 'conversionType', 'conversionRate']
    df.to_csv(f'./data/{prefix}conversionRates.csv')


    # compute Pc
    print('Starting to compute Pcs for cells with sufficient coverage.')
    if pc_calc == 'pc':

        pcs = Parallel(n_jobs=int(o.threads), verbose=3, backend='loky')(delayed(estimateP_c)(
            output[cell]['sparse_matrix'].toarray(), output[cell]['+']['conversionRates']['aG'], cell) for cell in list(unique_cells))
        # contains pc, pe, cellid for each processed cell.
        # This could be placed in a new output directory.
        with open('./data/%spc_pe.txt' % prefix, 'w') as pcof:
            pcof.write("cell\tpc\tpe\n")
            for item in pcs:
                pcof.write("%s\t%.5f\t%.5f\n" % (item[2], item[0], item[1]))
        print('Finished computing cell-specific Pcs')

    print('All done!')

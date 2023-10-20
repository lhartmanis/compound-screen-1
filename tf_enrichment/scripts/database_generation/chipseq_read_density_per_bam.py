import pandas, pysam
import argparse, os, sys, subprocess

def intarray(strrepr):
    return [int(v) for v in strrepr.strip(',').split(',')]

def error(message):
    print("Fatal error:", message, file=sys.stderr)
    exit(1)

def warn(message):
    print("Warning:", message, file=sys.stderr)

def chrformat(chrom):
    chrom = str(chrom)
    return chrom if chrom.startswith('chr') else 'chr'+chrom

if '__main__' == __name__:
    parser = argparse.ArgumentParser()
    parser.add_argument('bamfile')
    parser.add_argument('tablefile')
    parser.add_argument('name')
    parser.add_argument('--delete_bam', action='store_true')
    parser.add_argument('--RPKM', action='store_true')
    parser.add_argument('--nrows', type=int)
    parser.add_argument('--region', choices=['promoter', 'polyAsite', 'genebody', 'upstream', 'upstreamproximal', 'downstreamproximal', 'downstreamctrl'], default='promoter')
    o = parser.parse_args()
    
    txdata = pandas.read_csv('/home/danielr/from_crick2/meta/resources/Homo_sapiens.GRCh38.95.chr.genesymbols_genepred.txt', sep='\t', header=None, names=['name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames'], converters={'name':str, 'chrom':str, 'strand':str, 'txStart':int, 'txEnd':int, 'cdsStart':int, 'cdsEnd':int, 'exonCount':int, 'exonStarts':intarray, 'exonEnds':intarray, 'score':float, 'name2':str, 'cdsStartStat':str, 'cdsEndStat':str, 'exonFrames':str}, nrows=o.nrows)
    
    bamfile = o.bamfile
    tablefile = o.tablefile
    if not os.path.exists(bamfile + '.bai'):
        subprocess.check_call(['samtools', 'index', bamfile])
    
    gene_reads = dict()
    gene_locus_lengths = dict()
    samfile = pysam.AlignmentFile(bamfile, "rb")
    
    if o.RPKM:
        billionmappedreads = sum(info.mapped for info in samfile.get_index_statistics())/1e9
    
    for i, txrow in txdata[['name2', 'chrom', 'txStart', 'txEnd', 'strand']].iterrows():
        row = txrow
        if 'promoter' == o.region:
            if txrow['strand'] == '+':
                start = txrow['txStart'] - 1000
                end = txrow['txStart'] + 100
            else:
                start = txrow['txEnd'] - 100
                end = txrow['txEnd'] + 1000
        elif 'upstreamproximal' in o.region:
            if row['strand'] == '+':
                start = row['txStart'] - 5000
                end = row['txStart'] - 1000
            else:
                start = row['txEnd'] + 1000
                end = row['txEnd'] + 5000
        elif 'downstreamproximal' == o.region:
            if row['strand'] == '+':
                start = row['txStart'] + 100
                end = row['txStart'] + 5000
            else:
                start = row['txEnd'] - 5000
                end = row['txEnd'] - 100
        elif 'polyAsite' == o.region:
            if row['strand'] == '+':
                start = row['txEnd'] - 500
                end = row['txEnd'] + 500
            else:
                start = row['txStart'] - 500
                end = row['txStart'] + 500
        elif 'genebody' == o.region:
            if row['txEnd'] - row['txStart'] <= 5500:
                continue
            if row['strand'] == '+':
                start = row['txStart'] + 5000
                end = row['txEnd'] - 500
            else:
                start = row['txStart'] + 500
                end = row['txEnd'] - 5000
        elif 'upstream' == o.region:
            if row['strand'] == '+':
                start = row['txStart'] - 20000
                end = row['txStart'] - 5000
            else:
                start = row['txEnd'] + 5000
                end = row['txEnd'] + 20000
        elif 'downstreamctrl' == o.region:
            if row['strand'] == '+':
                start = row['txEnd'] + 5000
                end = row['txEnd'] + 20000
            else:
                start = row['txStart'] - 20000
                end = row['txStart'] - 5000
        chromosome = chrformat(txrow['chrom'])
        
        if chromosome in ('chrM', 'chrMT', 'chrY'): continue
        
        
        name2 = txrow['name2']
        gene_reads[name2] = set()
        gene_locus_lengths[name2] = end - start
        try:
            fetcher = samfile.fetch(chromosome, start, end)
        except ValueError:
            continue
        
        for read in fetcher:
            gene_reads[name2].add(read.query_name)
    if o.RPKM:
        pandas.Series({TF:len(reads)/gene_locus_lengths[TF]/billionmappedreads for TF, reads in gene_reads.items()}, name=o.name).to_frame().to_csv(tablefile)
    else:
        pandas.Series({TF:len(reads) for TF, reads in gene_reads.items()}, name=o.name).to_frame().to_csv(tablefile)
    
    samfile.close()
    
    
    if o.delete_bam:
        os.remove(bamfile+'.bai')
        os.remove(bamfile)
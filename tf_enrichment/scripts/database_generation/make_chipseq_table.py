import pandas, numpy, interlap
import argparse, os, collections, sys, random, time, gzip

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

def calc_length(boundaries):
    opened = 0
    length = 0
    for chrom, loc, boundtype in sorted(boundaries):
        if opened == 0:
            assert boundtype == 'begin'
            firststart = loc
            opened = 1
        elif boundtype == 'begin':
            opened += 1
        elif opened == 1:
            length += loc - firststart
            opened = 0
        else:
            opened -= 1
        assert opened >= 0
    assert opened == 0
    return length

if '__main__' == __name__:
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder_in', default='/home/danielr/from_crick2/work/features_of_genes/K562_chipseq/')
    parser.add_argument('--genepred_annotation', default='/home/danielr/from_crick2/meta/resources/Homo_sapiens.GRCh38.95.chr.genesymbols_genepred.txt')
    parser.add_argument('--table_out', default='/home/danielr/from_crick2/results/features_of_genes/human/K562_chipseq_to_GRCh38.95.chr.ensembl.txt')
    parser.add_argument('--out_metadata_prefilter')
    parser.add_argument('--out_metadata_filtered')
    parser.add_argument('--maxfiles', type=int)
    parser.add_argument('--region', choices=['promoter', 'proximal', 'polyAsite', 'genebody', 'upstream', 'upstreamproximal', 'downstreamproximal', 'downstreamctrl'], default=['promoter'], nargs='+')
    parser.add_argument('--suffix', default='')
    parser.add_argument('--PPKM', action='store_true')
    parser.add_argument('--min_neglog10_qValue', type=float)
    o = parser.parse_args()
    
    random.seed(0)
    
    # first: select the technically best experiment per transcription factor, this step takes 7 seconds
    
    metadata = pandas.read_table(os.path.join(o.folder_in, 'metadata.tsv'))
    metadata = metadata[metadata['Biosample genetic modifications methods'].isnull() & metadata['Biosample treatments'].isnull() & metadata['Biosample genetic modifications categories'].isnull()]
    metadata = metadata[metadata['Output type'].isin(('optimal IDR thresholded peaks', 'pseudoreplicated IDR thresholded peaks'))]
    
    metadata['num_non-compliant'] = metadata['Audit NOT_COMPLIANT'].str.split(',').str.len().fillna(0)
    metadata['num_warnings'] = metadata['Audit WARNING'].str.split(',').str.len().fillna(0)
    metadata['avg_fragsize'] = metadata['Library size range'].apply(lambda T: numpy.mean([float(v) for v in str(T).split('-')]))
    metadata['num_warnings_plus_100kdivSize'] = metadata['num_warnings'] + 100000/metadata['Size']
    num_choices = metadata['Experiment target'].value_counts()
    
    metadata = metadata.sort_values(by=['Output type', 'num_non-compliant', 'num_warnings_plus_100kdivSize'], ascending=[True, True, True])
    if o.out_metadata_prefilter:
        metadata.to_csv(o.out_metadata_prefilter, sep='\t', index=False)
    
    metadata = metadata.drop_duplicates(subset='Experiment target', keep='first')
    metadata['num_choices'] = list(num_choices.reindex(metadata['Experiment target']))
    metadata = metadata.sort_values(by='File accession')
    
    if o.maxfiles:
        metadata = metadata.head(o.maxfiles)
    
    if o.out_metadata_filtered:
        metadata.to_csv(o.out_metadata_filtered, sep='\t', index=False)
    
    
    # then, make an index of promoters, this step takes 20 seconds
    
    interlapindex = collections.defaultdict(interlap.InterLap)
    txdata = pandas.read_csv(o.genepred_annotation, sep='\t', header=None, names=['name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames'], converters={'name':str, 'chrom':str, 'strand':str, 'txStart':int, 'txEnd':int, 'cdsStart':int, 'cdsEnd':int, 'exonCount':int, 'exonStarts':intarray, 'exonEnds':intarray, 'score':float, 'name2':str, 'cdsStartStat':str, 'cdsEndStat':str, 'exonFrames':str})
    
    
    genedata = pandas.DataFrame(index=txdata['name2'].unique())
    genedata['tx_ids'] = txdata[['name2', 'name']].groupby('name2').agg({'name':lambda v: ','.join(v)})['name']
    
    
    genedata_default = pandas.Series(None, index=genedata.index)
    locus_boundaries = collections.defaultdict(list)
    
    # should the regions be allowed to overlap? no, have the user specify overlap by using several options in --regions instead
    if 'promoter' in o.region:
        for i, row in txdata[['name2', 'chrom', 'txStart', 'txEnd', 'strand']].iterrows():
            if row['strand'] == '+':
                start = row['txStart'] - 1000
                end = row['txStart'] + 100
            else:
                start = row['txEnd'] - 100
                end = row['txEnd'] + 1000
            interlapindex[chrformat(row['chrom'])].add((start, end, {'gene':row['name2']}))
            locus_boundaries[row['name2']].extend([(row['chrom'], start, 'begin'), (row['chrom'], end, 'end')])
            genedata_default[row['name2']] = 0 if o.PPKM else False
    if 'proximal' in o.region:
        for i, row in txdata[['name2', 'chrom', 'txStart', 'txEnd', 'strand']].iterrows():
            if row['strand'] == '+':
                start = row['txStart'] - 5000
                end = row['txStart'] - 1000
            else:
                start = row['txEnd'] + 1000
                end = row['txEnd'] + 5000
            interlapindex[chrformat(row['chrom'])].add((start, end, {'gene':row['name2']}))
            locus_boundaries[row['name2']].extend([(row['chrom'], start, 'begin'), (row['chrom'], end, 'end')])
            if row['strand'] == '+':
                start = row['txStart'] + 100
                end = row['txStart'] + 5000
            else:
                start = row['txEnd'] - 5000
                end = row['txEnd'] - 100
            interlapindex[chrformat(row['chrom'])].add((start, end, {'gene':row['name2']}))
            locus_boundaries[row['name2']].extend([(row['chrom'], start, 'begin'), (row['chrom'], end, 'end')])
            genedata_default[row['name2']] = 0 if o.PPKM else False
    if 'upstreamproximal' in o.region:
        for i, row in txdata[['name2', 'chrom', 'txStart', 'txEnd', 'strand']].iterrows():
            if row['strand'] == '+':
                start = row['txStart'] - 5000
                end = row['txStart'] - 1000
            else:
                start = row['txEnd'] + 1000
                end = row['txEnd'] + 5000
            interlapindex[chrformat(row['chrom'])].add((start, end, {'gene':row['name2']}))
            locus_boundaries[row['name2']].extend([(row['chrom'], start, 'begin'), (row['chrom'], end, 'end')])
            genedata_default[row['name2']] = 0 if o.PPKM else False
    if 'downstreamproximal' in o.region:
        for i, row in txdata[['name2', 'chrom', 'txStart', 'txEnd', 'strand']].iterrows():
            if row['strand'] == '+':
                start = row['txStart'] + 100
                end = row['txStart'] + 5000
            else:
                start = row['txEnd'] - 5000
                end = row['txEnd'] - 100
            interlapindex[chrformat(row['chrom'])].add((start, end, {'gene':row['name2']}))
            locus_boundaries[row['name2']].extend([(row['chrom'], start, 'begin'), (row['chrom'], end, 'end')])
            genedata_default[row['name2']] = 0 if o.PPKM else False
    if 'polyAsite' in o.region:
        for i, row in txdata[['name2', 'chrom', 'txStart', 'txEnd', 'strand']].iterrows():
            if row['strand'] == '+':
                start = row['txEnd'] - 500
                end = row['txEnd'] + 500
            else:
                start = row['txStart'] - 500
                end = row['txStart'] + 500
            interlapindex[chrformat(row['chrom'])].add((start, end, {'gene':row['name2']}))
            locus_boundaries[row['name2']].extend([(row['chrom'], start, 'begin'), (row['chrom'], end, 'end')])
            genedata_default[row['name2']] = 0 if o.PPKM else False
    if 'genebody' in o.region:
        for i, row in txdata[['name2', 'chrom', 'txStart', 'txEnd', 'strand']].iterrows():
            if row['txEnd'] - row['txStart'] <= 5500:
                continue
            if row['strand'] == '+':
                start = row['txStart'] + 5000
                end = row['txEnd'] - 500
            else:
                start = row['txStart'] + 500
                end = row['txEnd'] - 5000
            interlapindex[chrformat(row['chrom'])].add((start, end, {'gene':row['name2']}))
            locus_boundaries[row['name2']].extend([(row['chrom'], start, 'begin'), (row['chrom'], end, 'end')])
            genedata_default[row['name2']] = 0 if o.PPKM else False
    if 'upstream' in o.region:
        for i, row in txdata[['name2', 'chrom', 'txStart', 'txEnd', 'strand']].iterrows():
            if row['strand'] == '+':
                start = row['txStart'] - 20000
                end = row['txStart'] - 5000
            else:
                start = row['txEnd'] + 5000
                end = row['txEnd'] + 20000
            interlapindex[chrformat(row['chrom'])].add((start, end, {'gene':row['name2']}))
            locus_boundaries[row['name2']].extend([(row['chrom'], start, 'begin'), (row['chrom'], end, 'end')])
            genedata_default[row['name2']] = 0 if o.PPKM else False
    if 'downstreamctrl' in o.region:
        for i, row in txdata[['name2', 'chrom', 'txStart', 'txEnd', 'strand']].iterrows():
            if row['strand'] == '+':
                start = row['txEnd'] + 5000
                end = row['txEnd'] + 20000
            else:
                start = row['txStart'] - 20000
                end = row['txStart'] - 5000
            interlapindex[chrformat(row['chrom'])].add((start, end, {'gene':row['name2']}))
            locus_boundaries[row['name2']].extend([(row['chrom'], start, 'begin'), (row['chrom'], end, 'end')])
            genedata_default[row['name2']] = 0 if o.PPKM else False
    
    
    #print(genedata_default)
    assert o.region == ['genebody'] or not any(genedata_default.isna())
    
    
    # and then populate a table with results from each file, this 1 second per file, and there are 288 files
    
    
    peaksinfile = collections.defaultdict(int)
    for TF, accession in zip(metadata['Experiment target'], metadata['File accession']):
        TF = TF.split('-human')[0] + o.suffix
        #print("Starting at", accession, "for", TF, "at", time.asctime())
        filename = accession + '.bed.gz'
        genedata[TF] = genedata_default.copy()
        numtrue = 0
        with gzip.open(os.path.join(o.folder_in, filename), 'rt') as infh:
            for line in infh:
                p = line.rstrip().split('\t')
                middle = (int(p[1]) + int(p[2]))//2
                chromosome = chrformat(p[0])
                if o.min_neglog10_qValue is not None:
                    neglog10_qValue = float(p[8])
                    if neglog10_qValue < o.min_neglog10_qValue:
                        #print('skipped', p, 'due to low qValue', neglog10_qValue)
                        continue
                if chromosome in interlapindex:
                    peaksinfile[TF] += 1
                    genes = [match[2]['gene'] for match in interlapindex[chromosome].find((middle, middle))]
                    genes = [gene for gene in genes if not genedata.loc[gene, TF]]
                    if len(set(genes)) == 1:
                        if o.PPKM:
                            genedata.loc[genes[0], TF] += 1
                        else:
                            genedata.loc[genes[0], TF] = True
                        numtrue += 1
                    elif len(genes) > 0:
                        if o.PPKM:
                            genedata.loc[random.choice(genes), TF] += 1
                        else:
                            genedata.loc[random.choice(genes), TF] = True     # only set one True per binding site, to avoid double-counting during statistical testing
                        numtrue +=1
                elif "_random" not in chromosome and "chrUn" not in chromosome and chromosome != 'chrM':
                    warn("{} not in {},... in {} ({})".format(chromosome, ','.join(list(interlapindex)[:3]), accession, TF))
            if numtrue == 0:
                warn("No matches " + TF)
    
    print(genedata.shape)
    genedata = genedata.loc[:, genedata.any(axis=0)]  # remove all-False columns
    print(genedata.shape)
    
    
    if o.PPKM:
        tx_ids = genedata['tx_ids']
        del genedata['tx_ids']
        gene_locus_lengths = {gene:calc_length(boundaries) for gene, boundaries in locus_boundaries.items()}
        
        if set(genedata.index) - set(gene_locus_lengths.keys()):
            print(set(genedata.index) - set(gene_locus_lengths.keys()))
            genedata = genedata.loc[genedata.index.intersection(gene_locus_lengths.keys())]
        gene_locus_length_Series = pandas.Series({gene:gene_locus_lengths[gene] for gene in genedata.index})
        depth_Series = pandas.Series({TF:peaksinfile[TF] for TF in genedata.columns})
        genedata = genedata.divide(depth_Series)*1e9
        genedata = genedata.divide(gene_locus_length_Series, axis=0)
        
    
    genedata.to_csv(o.table_out, sep='\t')
    
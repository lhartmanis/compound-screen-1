import argparse, pandas

def intarray(strrepr):
    return [int(v) for v in strrepr.strip(',').split(',')]

def load_genepred(path):
    return pandas.read_csv(path, sep='\t', header=None, names=['name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames'], converters={'name':str, 'chrom':str, 'strand':str, 'txStart':int, 'txEnd':int, 'cdsStart':int, 'cdsEnd':int, 'exonCount':int, 'exonStarts':intarray, 'exonEnds':intarray, 'score':float, 'name2':str, 'cdsStartStat':str, 'cdsEndStat':str, 'exonFrames':str})

def commajoin(txdata, namecol):
    return txdata[['name2', namecol]].groupby('name2').agg({namecol:lambda v: ','.join(list(set(v)))})[namecol]

if '__main__' == __name__:
    parser = argparse.ArgumentParser()
    parser.add_argument('--genepred_annotation', default='/mnt/crick/rawdata/danielr/meta/resources/Homo_sapiens.GRCh38.95.chr.genesymbols_genepred.txt')
    parser.add_argument('--genepred_annotation_othername', default='/mnt/crick/rawdata/danielr/meta/resources/Homo_sapiens.GRCh38.95.chr.ensemblnames_genepred.txt.gz', nargs='?')
    parser.add_argument('--table_out', default='/mnt/crick/rawdata/danielr/results/features_of_genes/human/features_v2_GRCh38.95.chr.ensembl.txt')
    parser.add_argument('--txdata_debug')
    o = parser.parse_args()
    
    txdata = load_genepred(o.genepred_annotation)
    
    txdata['isCoding'] = txdata['cdsStart'] - txdata['cdsEnd'] != 0
    txdata['mRNAlength'] = txdata['exonEnds'].apply(sum) - txdata['exonStarts'].apply(sum)
    txdata['premRNAlength'] = txdata['txEnd'] - txdata['txStart']
    
    genedata = pandas.DataFrame(index=txdata['name2'].unique())
    genedata['tx_ids'] = commajoin(txdata, 'name')
    
    if o.genepred_annotation_othername is not None:
        ensname_txdata = load_genepred(o.genepred_annotation_othername)
        if list(ensname_txdata['name']) != list(txdata['name']): raise Exception
        txdata['gene_id'] = list(ensname_txdata['name2'])
        genedata['gene_ids'] = commajoin(txdata, 'gene_id')
        # I've checked, no ENSG id has several gene symbols
    
    genedata['hasIntron'] = txdata[['name2', 'exonCount']].groupby('name2').mean() > 1
    genedata['isCoding'] = txdata[['name2', 'isCoding']].groupby('name2').any()
    genedata['hasNonCodingIsoform'] = ~txdata[['name2', 'isCoding']].groupby('name2').all()
    genedata['mRNAlength'] = txdata[['name2', 'mRNAlength']].groupby('name2').mean()
    genedata['premRNAlength'] = txdata[['name2', 'premRNAlength']].groupby('name2').mean()
    
    
    genedata.to_csv(o.table_out, sep='\t')
    
    if o.txdata_debug:
        txdata.to_csv(o.txdata_debug, sep='\t')
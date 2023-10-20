import pandas, collections, argparse

def obo_accession_name_pairs(path):
    ID, name, skip = None, None, False
    with open(path) as infh:
        for line in infh:
            if line.startswith('id:'):
                ID = line.split(': ')[-1].strip()
            elif line.startswith('name:'):
                name = line.split(': ')[-1].strip()
            elif line.startswith('def: "OBSOLETE'):
                skip = True
            elif line.startswith('[Term]') or line.startswith('[Typedef]'):
                if ID is not None and not skip:
                    yield ID, name
                ID, name, skip = None, None, False

if '__main__' == __name__:
    parser = argparse.ArgumentParser()
    parser.add_argument('--obo', default='/home/danielr/from_crick2/work/features_of_genes/GO/go.obo')
    parser.add_argument('--goa', default='/home/danielr/from_crick2/work/features_of_genes/GO/goa_human_isoform.gaf.gz')
    parser.add_argument('-o', '--output', default='/home/danielr/from_crick2/results/features_of_genes/GO.tsv')
    o = parser.parse_args()
    obo_file = o.obo
    goa_file = o.goa
    
    GO_to_term = dict(obo_accession_name_pairs(obo_file))
    
    geneGO = pandas.read_csv(goa_file, sep='\t', comment='!', header=None, names=['DB', 'DB_ID', 'Symbol', 'Qualifier', 'GO_ID', 'DB_ref', 'Evidence', 'with_from', 'Aspect', 'DB_name', 'DB_syn', 'DB_type', 'Taxon', 'Date', 'AssignedBy', 'Extension', 'GeneProductForm'], usecols=['Symbol', 'GO_ID'])
    
    
    GO_to_symbols = collections.defaultdict(set)
    for symbol, GOid in zip(geneGO['Symbol'], geneGO['GO_ID']):
        GO_to_symbols[GOid].add(symbol)
    
    genedata = pandas.DataFrame(index=geneGO['Symbol'].unique())
    for GOid, symbols in GO_to_symbols.items():
        if GOid in GO_to_term:
            genedata.loc[symbols, GOid+': '+GO_to_term[GOid]] = True
    genedata = genedata.fillna(False)
    genedata = genedata.loc[:, genedata.sum()>=3]
    genedata.to_csv(o.output, sep='\t')
import argparse, pandas

if '__main__' == __name__:
    parser = argparse.ArgumentParser()
    parser.add_argument('other_table')
    parser.add_argument('bed_file', nargs='+')
    parser.add_argument('outfile')
    o = parser.parse_args()
    
    prev_table = pandas.read_table(o.other_table, index_col=0, usecols=[0,1])
    new_table = pandas.DataFrame(index=prev_table.index)
    
    
    for bedfile in o.bed_file:
        bed_df = pandas.read_table(bedfile, header=None, names=['chrom', 'start', 'end', 'gene_1', 'qual', 'strand'])
        genes = [gene.split('_')[0] for gene in bed_df['gene_1']]
        colname = bedfile.split('/')[-1].split('.bed')[0]
        new_table[colname] = pandas.Series({gene:True for gene in genes})
    new_table.fillna(False, inplace=True)
    new_table.to_csv(o.outfile, sep='\t')
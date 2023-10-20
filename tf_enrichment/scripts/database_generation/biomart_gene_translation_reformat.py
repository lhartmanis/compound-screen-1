import pandas, argparse

if '__main__' == __name__:
    parser = argparse.ArgumentParser()
    parser.add_argument('table_in')
    parser.add_argument('table_out')
    o = parser.parse_args()
    
    df = pandas.read_table(o.table_in).dropna()
    df.index = df['Gene name']
    df.index.name = None
    del df['Gene name']
    df.columns = ['gene_ids']
    df.to_csv(o.table_out, sep='\t')

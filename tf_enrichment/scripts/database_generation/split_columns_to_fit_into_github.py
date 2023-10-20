import pandas, argparse

if '__main__' == __name__:
    parser = argparse.ArgumentParser()
    parser.add_argument('splits', type=int)
    parser.add_argument('table_in')
    o = parser.parse_args()

    table = pandas.read_table(o.table_in, index_col=0)

    num_cols = 1+len(table.columns)//o.splits

    columns_printed = []
    base_name = o.table_in.rsplit('.txt', 1)[0]
    for split in range(o.splits):
        df = table[table.columns[split*num_cols:(1+split)*num_cols]]
        columns_printed.extend(df.columns)
        df.to_csv(base_name+'.part'+str(1+split)+'.txt', sep='\t')
    assert columns_printed == list(table.columns)
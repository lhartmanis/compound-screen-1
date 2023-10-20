import argparse, pandas

if '__main__' == __name__:
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--feature_table', required=True)
    parser.add_argument('-b', '--ctrl_table', required=True)
    parser.add_argument('-m', '--method', choices=['subtract', 'subtract_min0', 'add_divide'], required=True)
    parser.add_argument('-o', '--output_table', required=True)
    o = parser.parse_args()
    
    if '{}' in o.output_table:
        o.output_table = o.output_table.replace('{}', o.feature_table.split('/')[-1].rsplit('.',1)[0])
    
    fg = pandas.read_table(o.feature_table, index_col=0)
    bg = pandas.read_table(o.ctrl_table, index_col=0)
    
    if len(fg.columns) != len(bg.columns): raise Exception
    
    
    bg = bg.loc[fg.index.intersection(bg.index), :]
    bg.columns = fg.columns
    
    if o.method == 'subtract':
        out = fg - bg
    elif o.method == 'subtract_min0':
        out = fg - bg
        out[out < 0] = 0
    elif o.method == 'add_divide':
        avg_bg_per_sample = bg.mean()
        out = fg.add(avg_bg_per_sample, axis=1) / bg.add(avg_bg_per_sample, axis=1)
    
    out.to_csv(o.output_table, sep='\t')
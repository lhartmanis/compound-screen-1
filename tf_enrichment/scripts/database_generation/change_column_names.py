import pandas, argparse, re

def multiple_replace(string, rep_dict):
    # https://stackoverflow.com/questions/6116978/how-to-replace-multiple-substrings-of-a-string
    pattern = re.compile("|".join([re.escape(k) for k in sorted(rep_dict,key=len,reverse=True)]), flags=re.DOTALL)
    return pattern.sub(lambda x: rep_dict[x.group(0)], string)

if '__main__' == __name__:
    parser = argparse.ArgumentParser()
    parser.add_argument('table_tsv_in')
    parser.add_argument('table_tsv_out')
    parser.add_argument('-r', '--replace', nargs=2, action='append', metavar=('from', 'to'))
    o = parser.parse_args()
    
    df = pandas.read_table(o.table_tsv_in, index_col=0)
    replacements = dict(o.replace)
    df.columns = [multiple_replace(colname, replacements) for colname in df.columns]
    df.to_csv(o.table_tsv_out, sep='\t')
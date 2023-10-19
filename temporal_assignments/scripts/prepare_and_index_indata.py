import pandas as pd
import argparse
import subprocess
import os, re
import glob

if __name__ == "__main__":
    parser = argparse.ArgumentParser("This script takes the details file and conversion probabilities and returns a tabix indexed file ready for pi_g estimation")
    parser.add_argument('-d', '--work_dir', help = "Path to working directory. This folder shall include the details file & the Pc and Pe estimate file")
    args = parser.parse_args()

    work_dir = args.work_dir

    re_prob = re.compile(r'pc_pe')
    re_det  = re.compile(r'details')

    det_path  = list(filter(re_det.search, glob.glob(os.path.join(work_dir, '*'))))[0]
    prob_path = list(filter(re_prob.search, glob.glob(os.path.join(work_dir, '*'))))[0]

    print("Loading data...")
    details = pd.read_csv(det_path, sep = '\t', names = ['BC', 'gene_id', 'convs', 'content'])
    probs   = pd.read_csv(prob_path, sep = '\t')
    print("\tdone!")

    print("Merging files...")
    all_res = pd.merge(details, probs, left_on = 'BC', right_on = 'cell')
    all_res.drop('cell', axis = 1, inplace = True)
    all_res.loc[:, 'num_reads'] = all_res.convs.apply(lambda x: len(x.split(',')))
    print("\tdone!")

    bc_dict   = dict(zip(sorted(all_res.BC.unique()),range(all_res.BC.nunique())))
    gene_dict = dict(zip(sorted(all_res.gene_id.unique()), range(all_res.gene_id.nunique())))

    all_res.loc[:, 'BC_idx'] = all_res.BC.map(bc_dict)
    all_res.loc[:, 'gene_idx']   = all_res.gene_id.map(gene_dict)

    all_res = all_res.sort_values(by = ['BC', 'gene_id'])[['BC_idx', 'gene_idx', 'BC', 'gene_id', 'convs', 'content', 'pc', 'pe', 'num_reads']].reset_index(drop = True)
    all_res.drop('BC_idx', axis = 1, inplace = True)
    all_res.reset_index(inplace = True)
    all_res.rename({'index': '#data_idx'}, axis = 1, inplace = True)

    print("Saving results to file...")
    all_res.to_csv(os.path.join(work_dir, 'combined_results.txt'), sep = '\t', index = False)
    print("\tdone!")

    print("Compressing file...")
    subprocess.run(f"bgzip {os.path.join(work_dir, 'combined_results.txt')}", shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    print("\tdone!")

    print("Tabix indexing file...")
    subprocess.run(f"tabix {os.path.join(work_dir, 'combined_results.txt.gz')} -s 1 -b 2 -e 2 -c '#' -0 -f", shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    print("\tdone!")
    
    num_comps = all_res.shape[0]
    with open(os.path.join(work_dir, "total_comparisons.txt"), 'w') as fh:
        fh.write(str(num_comps))
    
    print("All done!")

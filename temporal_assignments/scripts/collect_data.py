import pandas as pd
import os, glob, re, shutil
import multiprocessing as mp
import argparse

def concatenate_results(work_folder):
    files_in_folder = glob.glob(os.path.join(work_folder, "*"))
    dfs = []
    
    for file in files_in_folder:
        dfs.append(pd.read_csv(file, sep = '\t', header = None))
    
    final_df = pd.concat(dfs)
    final_df.columns = ['cell-X-gene', 'conf5', 'average', 'conf95', 'numreads']
    final_df   = final_df.reset_index(drop = True)
    name_cols  = final_df.loc[:, 'cell-X-gene'].str.split('__', expand = True).rename({0: 'BC', 1: 'gene_id'}, axis = 1)
    final_df = pd.concat([name_cols, final_df.drop('cell-X-gene', axis = 1)], axis = 1)
    final_df.rename({"average": "pi_g"}, axis = 1, inplace = True)
    final_df.loc[:, "conf_width"] = final_df.loc[:, "conf95"] - final_df.loc[:, "conf5"]

    return final_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--base_path", help = "Path to working directory where pi_g inference results are stored")
    parser.add_argument("-p", '--num_procs', help = "Number of processes to use", type = int)
    parser.add_argument("-f", '--handle_folders', help = "How to handle intermediary folders [remove/ keep]. Default is keep", type = str, default = "keep")
    args = parser.parse_args()

    base_path = args.base_path
    num_procs = args.num_procs
    handle_folders = args.handle_folders

    print("Starting to collect output...")
    
    all_files = glob.glob(os.path.join(base_path, '*'))
    pattern   = re.compile("completed_pigs")
    all_pigs  = list(filter(pattern.search, all_files))
    all_pigs.append(os.path.join(base_path, "pi_g_results")) 

    if len(all_pigs) < num_procs:
        num_procs = len(all_pigs)

    print(f"Concatenating data from {len(all_pigs)} folders using {num_procs} parallel processes...")
    pool = mp.Pool(processes = num_procs)
    jobs = [pool.apply_async(concatenate_results, (work_folder, )) for work_folder in all_pigs]

    results = [job.get() for job in jobs]
    out_df = pd.concat(results)
    out_df = out_df.sort_values(by = ['BC', 'gene_id']).reset_index(drop = True)
    print("Saving results to file...")
    out_df.to_csv("data/pi_g_df.txt.gz", sep = '\t')

    if handle_folders == "remove":
        print("Removing intermediary folders...")
        for folder in all_pigs:
            shutil.rmtree(folder)
    else:
        print("Keeping intermediary files and folders")

    print("All done!")

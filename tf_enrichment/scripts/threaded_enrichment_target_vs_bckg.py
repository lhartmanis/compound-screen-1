import subprocess
import glob, os
import datetime
from time import time
import argparse
import multiprocessing as mp

def run_gsea(folder, target_name, background_name, sig_test, chip_type):
    regions = ['upstream', 'upstreamproximal', 'downstreamproximal' ,'promoter', 'genebody', 'polyAsite', 'proximal', 'downstreamctrl']
    script_path = "/mnt/davidson/danielr/scripts/features_of_genes/compare_gene_set_to_features.py"
    fdr_thresh = 0.1
    print(f"Working on {folder}")
    target_list     = os.path.join(folder, target_name)
    background_list = os.path.join(folder, background_name)

    comp_dict = {'target_vs_bckg': {'i': target_list,
                                    'b': background_list}}

    for comparison in comp_dict.keys():
        print(f"Working on {comparison}...")
        outfolder = os.path.join(folder, comparison)
        if not os.path.exists(outfolder):
            os.mkdir(outfolder)

        for region in regions:
            print(f"Working on {region}...")
            save_folder = os.path.join(outfolder, region)
            if not os.path.exists(save_folder):
                os.mkdir(save_folder)

            if chip_type == "filtered_peaks":
                chipfile = f"/compound-screen/tf_enrichment/chip_data/{region}.txt"
                
            else:
                print("Could not recognize ChIP file. Exiting...")
                sys.exit()

            if not os.path.exists(chipfile):
                continue

            volcano  = os.path.join(save_folder, f"{region}_volcano.pdf")
            pbars    = os.path.join(save_folder, f"{region}_pbars.pdf")
            results  = os.path.join(save_folder, f"{region}_chip_enrichment.txt")

            func_call = f"python3 {script_path} -i {comp_dict[comparison]['i']} -b {comp_dict[comparison]['b']} -o {results} -F {chipfile} -v {volcano} --Pbars_plot {pbars} --FDR_in_plots {fdr_thresh} --list_hits --quanttest {sig_test}"
            with open(os.path.join(save_folder, 'func_call.txt'), 'w') as fh:
                fh.write(func_call)
            subprocess.run(func_call, shell = True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infolder', type = str, help = 'Path to a folder containing GSEA data in one or more subfolders. It is a mess and I need to change it soon.')
    parser.add_argument('-p', '--processes', type = int, help = 'Number of concurrent processes to utilize')
    parser.add_argument('--target_name', type = str, help = 'Name of target gene list. Default is target.txt', default = "target.txt")
    parser.add_argument('--background_name', type = str, help = 'Name of background gene list. Default is background.txt', default = "background.txt")
    parser.add_argument('-c', '--chip_type', type = str, help = "Decide what ChIP-data to use. Implemented files are: [filtered_peaks]. Default is filtered_peaks", default = 'filtered_peaks')
    parser.add_argument('-s', '--sig_test', type = str, help = "Decide which significance test to use. Implemented tests are: [Ksamp_AndersonDarling/ Wilcoxon_ranksums]. Default is Wilcoxon ranksums", default = "Wilcoxon_ranksums")
    args = parser.parse_args()

    target_name = args.target_name
    background_name = args.background_name
    folder_path = args.infolder
    n_jobs    = args.processes
    sig_test  = args.sig_test
    chip_type = args.chip_type

    folders = glob.glob(os.path.join(folder_path, '*'))
    print(folders)

    tic = time()
    print(f"Starting at {datetime.datetime.now().strftime('%H:%M:%S')}")
    pool = mp.Pool(processes = n_jobs)
    jobs = [pool.apply_async(run_gsea, (folder, target_name, background_name, sig_test, chip_type, )) for folder in folders]
    results = [job.get() for job in jobs]

    print(f"Done at {datetime.datetime.now().strftime('%H:%M:%S')}")
    toc = time()
    print(f"Took {toc-tic} seconds")

import pandas, numpy, requests
import argparse, os, collections, sys, subprocess
from concurrent import futures

def intarray(strrepr):
    return [int(v) for v in strrepr.strip(',').split(',')]

def error(message):
    print("Fatal error:", message, file=sys.stderr)
    exit(1)

def warn(message):
    print("Warning:", message, file=sys.stderr)

def chrformat(chrom):
    chrom = str(chrom)
    return chrom if chrom.startswith('chr') else 'chr'+chrom

def run_one_bam(bamfile, tablefile, URL, TFname, try_short, region, useRPKM):
    print(TFname)
    os.chdir('/home/danielr/from_crick2/work/features_of_genes/K562_read_density')
    if not os.path.exists(bamfile):
        subprocess.check_call(['wget', URL])  # takes 80s, is a sorted bam file
        if not os.path.exists(bamfile):
            print(bamfile, URL)
            raise Exception
    cmd = ['python3', '/home/danielr/from_crick2/scripts/features_of_genes/database_generation/chipseq_read_density_per_bam.py', bamfile, tablefile, TFname, '--region', region]
    if try_short: cmd += ['--nrows', '10']
    else: cmd += []#['--delete_bam']
    if useRPKM: cmd += ['--RPKM']
    print(cmd)
    subprocess.check_call(cmd)
    return tablefile

def get_bam_info(accession):
    html = requests.get('https://www.encodeproject.org/files/%s/'%accession).text
    assembly = html[html.index('Assembly</dt><dd>')+17:].split('</dd>',1)[0]
    filesize = html[html.index('File size</dt><dd>')+18:].split('</dd>',1)[0]
    numB = float(filesize.split()[0])
    if filesize.endswith('GB'): numB*=1e9
    elif filesize.endswith('MB'): numB*=1e6
    elif filesize.endswith('kB'): numB*=1e3
    return assembly, numB, accession

if '__main__' == __name__:
    parser = argparse.ArgumentParser()
    parser.add_argument('region', choices=['promoter', 'polyAsite', 'genebody', 'upstream', 'upstreamproximal', 'downstreamproximal', 'downstreamctrl'])
    parser.add_argument('--try_short', action='store_true')
    parser.add_argument('--pop_bams', type=int, default=0)
    parser.add_argument('--RPKM', action='store_true')
    o = parser.parse_args()
    
    try_short = o.try_short
    region = o.region
    
    
    
    # first: select the technically best experiment per transcription factor, this step takes 7 seconds
    metadata = pandas.read_table('/home/danielr/from_crick2/work/features_of_genes/K562_chipseq/metadata.tsv')
    metadata = metadata[metadata['Biosample genetic modifications methods'].isnull() & metadata['Biosample treatments'].isnull() & metadata['Biosample genetic modifications categories'].isnull()]
    metadata = metadata[metadata['Output type'].isin(('optimal IDR thresholded peaks', 'pseudoreplicated IDR thresholded peaks'))]
    metadata['num_non-compliant'] = metadata['Audit NOT_COMPLIANT'].str.split(',').str.len().fillna(0)
    metadata['num_warnings'] = metadata['Audit WARNING'].str.split(',').str.len().fillna(0)
    metadata['avg_fragsize'] = metadata['Library size range'].apply(lambda T: numpy.mean([float(v) for v in str(T).split('-')]))
    metadata['num_warnings_plus_100kdivSize'] = metadata['num_warnings'] + 100000/metadata['Size']
    num_choices = metadata['Experiment target'].value_counts()
    metadata = metadata.sort_values(by=['Output type', 'num_non-compliant', 'num_warnings_plus_100kdivSize'], ascending=[True, True, True])
    metadata['rank'] = list(range(len(metadata)))
    metadata = metadata.drop_duplicates(subset='Experiment accession', keep='first')
    metadata['num_choices'] = list(num_choices.reindex(metadata['Experiment target']))
    metadata['ID'] = '/experiments/' + metadata['Experiment accession'] + '/'
    
    experimentreport = pandas.read_table('/home/danielr/from_crick2/work/features_of_genes/K562_read_density/experiment_report_2020_12_23_8h_12m.tsv', skiprows=1)    
    experimentreport = experimentreport.merge(metadata, suffixes=[None, '_metadata'], left_on='ID', right_on='ID', validate='1:1', how='inner').sort_values(by='rank', na_position='last') # inner merge instead of left merge to remove treated or genetically modified samples
    experimentreport = experimentreport.drop_duplicates('Target of assay', keep='first').reset_index()
    
    # second, load a list of bam files
    bam_URLs = dict()
    with open('/home/danielr/from_crick2/work/features_of_genes/K562_read_density/files-13.txt', 'rt') as infh:
        for line in infh:
            if line.startswith('http'):
                bam_URLs[line.split('/')[-1].split('.bam')[0]] = line.strip()
    
    processpool = futures.ProcessPoolExecutor(10)
    
    # create a file per sample
    
    jobs = []
    for rowindex, TFrow in experimentreport.iterrows():
        if try_short and rowindex > 5: break
        file_accessions = [p.strip('/').split('/')[-1] for p in TFrow["Files"].split(',')]
        bam_accessions = [acc for acc in file_accessions if acc in bam_URLs]
        bam_infos = [info for info in sorted(get_bam_info(a) for a in bam_accessions) if '38' in info[0]]
        for pop_i in range(o.pop_bams):
            if os.path.exists('/home/danielr/from_crick2/work/features_of_genes/K562_read_density/'+bam_infos[-1][-1]+'.bam'): break
            if len(bam_infos) > 1:
                bam_infos.pop()
        print(rowindex)
        if len(bam_infos) >= 1:
            accession = bam_infos[-1][-1]
            experimentreport.loc[rowindex, 'chosenbam_sizeB'] = bam_infos[-1][1]
            experimentreport.loc[rowindex, 'chosenbam_accession'] = accession
            experimentreport.loc[rowindex, 'chosenbam_assembly'] = bam_infos[-1][0]
        else:
            warn("No BAM file for {}".format(gene))
            continue
        gene = TFrow['Target of assay']
        bamfile = '/home/danielr/from_crick2/work/features_of_genes/K562_read_density/{}.bam'.format(accession)
        tablefile = '/home/danielr/from_crick2/work/features_of_genes/K562_read_density/{}{}.csv'.format(accession, region)
        if os.path.exists(tablefile):
            jobs.append(processpool.submit(str, tablefile))
            continue
        
        try:
            URL = bam_URLs[accession]
        except KeyError:
            warn("No BAM file for {} ({})".format(accession, gene))
            continue
        print("Sending", gene)
        jobs.append(processpool.submit(run_one_bam, bamfile, tablefile, URL, gene + '_K562_'+region+('_D' if o.RPKM else '_R'), try_short, region, o.RPKM))
    
    experimentreport.to_csv('experimentreport_after_merging.tsv', sep='\t')
    
    # collect output to one file
    series_collection = dict()
    for job in jobs:
        tablefile = job.result()
        ser = pandas.read_csv(tablefile, index_col=0, squeeze=True)
        print(ser.head())
        #ser = ser.iloc[:, 0]
        series_collection[ser.name] = ser
    pandas.DataFrame(series_collection).to_csv('/home/danielr/from_crick2/work/features_of_genes/K562_read_density/try_truncated.txt' if try_short else '/home/danielr/from_crick2/results/features_of_genes/normreads_proximalsplit_K562/K562_chipseq_'+region+'_RPKMs_to_GRCh38.95.chr.ensembl.txt' if o.RPKM else '/home/danielr/from_crick2/results/features_of_genes/readcount_proximalsplit_K562_v2/K562_chipseq_'+region+'_reads_to_GRCh38.95.chr.ensembl.txt', sep='\t')
    
    # took 6 hours to run, excluding downloads
import pandas, numpy
from scipy import interpolate

# based on stattest_nonzero_fraction_and_values_v7_csv.py

def loadlist(path):
    with open(path, 'r') as infh:
        return [l.strip() for l in infh]

def binary(value):
    return int(value > 0)

def adjustP(column):
    Pval = column.dropna()
    if len(Pval)==0: return None
    return pandas.Series(multipletests(Pval, method="fdr_tsbh")[1], index=Pval.index)

def rank(sortby, handleties=0):
	""" return rank for each value in sortby, in the same order """
	zippedin = list(zip(sortby, list(range(len(sortby)))))
	if handleties:
		zippedin.sort()
		ranks = list(range(len(zippedin)))
		lastsameindex = -1
		lastsamevalue = None
		for ii in range(len(zippedin)):
			if zippedin[ii][0] != lastsamevalue:
				if lastsameindex != -1:
					targetrank = sum(ranks[lastsameindex:ii])/float(ii-lastsameindex)
					for jj in range(lastsameindex, ii):
						ranks[jj] = targetrank
				lastsameindex = ii
				lastsamevalue = zippedin[ii][0]
		if lastsameindex != -1:
			ii = len(zippedin)
			targetrank = sum(ranks[lastsameindex:ii])/float(ii-lastsameindex)
			for jj in range(lastsameindex, ii):
				ranks[jj] = targetrank
	else:
		import random
		random.shuffle(zippedin)
		zippedin.sort(key=lambda o:o[0])
		ranks = list(range(len(zippedin)))
	zippedout = list(zip([z[1] for z in zippedin], ranks))
	zippedout.sort()
	return [z[1] for z in zippedout]

def do_some_permutations(func, c, lena, n):
    import random
    arr = []
    while len(arr) < n:
        random.shuffle(c)
        arr.append(abs(func(c[:lena]) - func(c[lena:])))
    return arr

def expanding_permutation_test(in_func, a, b, controls=[100, 1000, 10000], proc=1):
    import random, math
    func = lambda V: math.log10(in_func(V))
    arr = [abs(func(a) - func(b))]
    if str(arr[0])=='nan': return float('nan'), float('nan')
    c = a.append(b, ignore_index=True).values
    lena = len(a)
    for n_ctrl in controls:
        if n_ctrl in controls[3:] and proc > 1:
            pool = futures.ProcessPoolExecutor(proc)
            jobs = []
            ctrls_per_process = (n_ctrl - len(arr))//proc
            ctrls_sent = 0
            for i in range(proc):
                jobs.append(pool.submit(do_some_permutations, func, c, lena, ctrls_per_process))
                ctrls_sent += ctrls_per_process
            while len(arr)+ctrls_sent <= n_ctrl:
                random.shuffle(c)
                arr.append(abs(func(c[:lena]) - func(c[lena:])))
            for job in jobs:
                arr += job.result()
        else:
            while len(arr) <= n_ctrl:
                random.shuffle(c)
                arr.append(abs(func(c[:lena]) - func(c[lena:])))
        noNA = [v for v in arr if str(v) != 'nan']
        r = rank(noNA)
        if n_ctrl-r[0] >= 20 or r[0] < len(noNA)/2:
            break
    
    if len(noNA) < 2:
        return float('nan'), float('nan')
    
    Pval = 1-float(r[0])/(len(r)+1)
    
    conf_width_values = numpy.percentile(noNA[1:], [2.5, 97.5])
    
    if conf_width_values[1]==conf_width_values[0]:
        return Pval, float('nan')
    else:
        return Pval, 10**(conf_width_values[1]-conf_width_values[0])

def kdtree_query_and_interpolation(kdtree, point, param_lookup, interpolation_points):
    # returns k_on, k_off, k_syn
    if kdtree is None: raise ValueError
    import warnings
    warnings.filterwarnings("ignore")
    try:
        v = interpolation_points(point)
        if 'nan' in str(v): raise Exception 
        return v
    except:
        d, closest_lookup_index = kdtree.query(point)
        return param_lookup[closest_lookup_index]

def lookup_and_ML(lookupML_support, point, value_series, runML):
    from scipy.optimize import minimize
    kdtree, param_lookup, interpolation_points, prob_fun, time, deg, proxy_cols = lookupML_support
    point = tuple(v for v,p in zip(point, ['freq_adj0', 'size_adj0', 'size_cov']) if p in proxy_cols)
    looked_up = kdtree_query_and_interpolation(kdtree, point, param_lookup, interpolation_points)
    if not runML: return looked_up
    x0 = tuple(k/deg for k in looked_up)
    bnds = tuple((k/100, k*100) for k in x0)
    for prob_function in prob_fun:
        fun = make_loglikelihood_fun(time*deg, value_series, prob_function)
        res = minimize(fun, x0, method='L-BFGS-B', bounds=bnds)
        if res.success:
            return tuple(k*deg for k in res.x)
    else:
        return float('nan'), float('nan'), float('nan')

def make_c_function(prec, path='/home/antonl/programs/hyp_new_git/hyp1f1_arb/hyp1f1_arb_vector.so'):
    import ctypes
    c_lib = ctypes.CDLL(path)
    c_lib.calculate_probs.argtypes =  [ctypes.POINTER(ctypes.c_double), ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_ulong, ctypes.c_ulong]
    c_lib.calculate_probs.restype = None
    prec_ulong = ctypes.c_ulong(prec)
    def prepared_function(kon_double, koff_double, ksyn_double, t_double, nmax):

        kon_c_double = ctypes.c_double(kon_double)
        koff_c_double = ctypes.c_double(koff_double)
        ksyn_c_double = ctypes.c_double(ksyn_double)
        t_c_double = ctypes.c_double(t_double)
        nmax_ulong = ctypes.c_ulong(nmax)

        res = (ctypes.c_double*(nmax+1))()

        c_lib.calculate_probs(res, kon_c_double, koff_c_double, ksyn_c_double, t_c_double, nmax_ulong, prec_ulong)
        
        res_list = list(res)
        
        return res_list
    return prepared_function

def make_loglikelihood_fun(t, vals, prob_fun):
    from collections import Counter
    vals_dict = Counter(vals)
    nmax = int(max(vals))
    def loglikelihood(x):
        kon = x[0]
        koff = x[1]
        ksyn = x[2]
        
        prob_list = prob_fun(kon,koff,ksyn,t,nmax)
        ll_dict = dict(enumerate(numpy.log(numpy.clip(prob_list, 1e-16,1))))
        
        return -numpy.sum([n_obs*ll_dict[count] for count,n_obs in vals_dict.items()])
    return loglikelihood

def calc_nums(lookupML_support, value_series, params, runML):
    freq_adj0 = value_series.apply(binary).mean()
    if freq_adj0 == 0:
        return 0, 0, float('nan'), float('nan'), float('nan'), float('nan'), float('nan'), float('nan')
    size_adj0 = value_series.mean() / freq_adj0
    size_cov = numpy.std([v for v in value_series if v > 0])/size_adj0
    if freq_adj0 == 1:
        return len(value_series), 1, size_adj0, size_cov, float('inf'), float('nan'), float('nan'), float('nan')
    kon, koff, ksyn = lookup_and_ML(lookupML_support, (freq_adj0, size_adj0, size_cov), value_series, runML)[:3]
    ret = [value_series.apply(binary).sum(), freq_adj0, size_adj0, size_cov]
    for param in params:
        if param == 'kon': ret.append(kon)
        elif param == 'koff': ret.append(koff)
        elif param == 'ksyn': ret.append(ksyn)
        elif param == 'burstsize': ret.append(ksyn/koff)
        elif param == 'burstfrequency': ret.append(1/(1/kon + 1/koff))
        elif param == 'expressionrate':ret.append(kon*ksyn/(kon+koff))
        elif param == 'meanoccupancy': ret.append(kon/(kon+koff))
        else: raise ValueError
    return ret

def get_param(lookupML_support, param, runML, value_series):
    freq_adj0 = numpy.mean([binary(v) for v in value_series])
    if freq_adj0 == 0: return float('nan')
    if freq_adj0 == 1: return float('nan')
    size_adj0 = value_series.mean() / freq_adj0
    size_cov = numpy.std([v for v in value_series if v > 0])/size_adj0
    kon, koff, ksyn = lookup_and_ML(lookupML_support, (freq_adj0, size_adj0, size_cov), value_series, runML)[:3]
    if param == 'kon': return kon
    if param == 'koff': return koff
    if param == 'ksyn': return ksyn
    if param == 'burstsize': return ksyn/koff
    if param == 'burstfrequency': return 1/(1/kon + 1/koff)
    if param == 'expressionrate': return kon*ksyn/(kon+koff)
    if param == 'meanoccupancy': return kon/(kon+koff)
    raise ValueError

def one_gene(gene, df1_gene, df2_gene, freq_size_adj0_kdtree, param_lookup, permutations, interpolation_points, proc, time, deg, hyp1f1_arb_vector_path, MLprecision, proxy_cols, params_out):
    import warnings, functools
    #warnings.filterwarnings("ignore")
    
    if freq_size_adj0_kdtree is None:
        freq_size_adj0_kdtree, param_lookup, interpolation_points = build_KD_tree(param_lookup, True, proxy_cols)
    
    prob_fun = [make_c_function(prec, hyp1f1_arb_vector_path) for prec in MLprecision]
    lookupML_support = freq_size_adj0_kdtree, param_lookup, interpolation_points, prob_fun, time, deg, proxy_cols
    
    row = pandas.Series(name=gene)
    res = calc_nums(lookupML_support, df1_gene, params_out, False)
    row['G1_binary_sum'], row["G1_binary_avg"], row["G1_nonzero_avg"], row["G1_nonzero_cov"] = res[:4]
    for i, param in enumerate(params_out): row["G1_"+param] = res[4+i]
    res = calc_nums(lookupML_support, df2_gene, params_out, False)
    row['G2_binary_sum'], row["G2_binary_avg"], row["G2_nonzero_avg"], row["G2_nonzero_cov"] = res[:4]
    for i, param in enumerate(params_out): row["G2_"+param] = res[4+i]
    res = calc_nums(lookupML_support, df1_gene, params_out, True)
    for i, param in enumerate(params_out): row["G1_ML_"+param] = res[4+i]
    res = calc_nums(lookupML_support, df2_gene, params_out, True)
    for i, param in enumerate(params_out): row["G2_ML_"+param] = res[4+i]
    
    for param in params_out:
        if row["G1_binary_avg"] not in (0, 1) and row["G2_binary_avg"] not in (0, 1) and row["G1_nonzero_avg"] != 1 and row["G2_nonzero_avg"] != 1:
            fetch_func = functools.partial(get_param, lookupML_support, param, False)
            # 95%confratio = permutations cause random ratio betwenn 2.5 and 75.5 percentiles that's this big, ie this is the limit for being able to find a p=0.05
            row["G1vsG2_"+param+"_P"], row["G1vsG2_"+param+"_95%confratio"] = expanding_permutation_test(fetch_func, df1_gene, df2_gene, permutations, proc)
        else:
            row["G1vsG2_"+param+"_P"], row["G1vsG2_"+param+"_95%confratio"] = None, None
    
    return row

def hash_dataframe(params_to_zeros_table):
    return int(pandas.util.hash_pandas_object(params_to_zeros_table).sum())

def build_KD_tree(params_to_zeros, from_one_gene, proxy_cols):
    params_to_zeros_table = pandas.read_table(params_to_zeros).dropna().drop_duplicates(subset=proxy_cols, keep=False)
    freq_size_adj0_kdtree = spatial.KDTree(params_to_zeros_table[proxy_cols])
    param_lookup = list(zip(params_to_zeros_table[open_chance_name], params_to_zeros_table[close_chance_name], params_to_zeros_table[transcribe_chance_name], *(params_to_zeros_table[col] for col in proxy_cols)))
    
    interpolation_points = interpolate.LinearNDInterpolator([param_lookup_entry[3:] for param_lookup_entry in param_lookup], [param_lookup_entry[:3] for param_lookup_entry in param_lookup])
    
    return freq_size_adj0_kdtree, param_lookup, interpolation_points

def strictlyDecreasing(l):
    # adapted from https://stackoverflow.com/questions/3755136/pythonic-way-to-check-if-a-list-is-sorted-or-not/4404056#4404056
    for i, el in enumerate(l[1:]):
        if el >= l[i]:
            return False
    return True
def strictlyIncreasing(l):
    for i, el in enumerate(l[1:]):
        if el <= l[i]:
            return False
    return True

if '__main__' == __name__:
    import tqdm
    import collections, argparse
    from statsmodels.stats.multitest import multipletests
    from scipy import spatial
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--params_to_zeros", required=True)
    parser.add_argument('--max_genes', type=int)
    parser.add_argument('--bootstrap', '--permutations', dest='permutations', type=int, nargs='+', default=[100, 1000, 10000, 100000])
    parser.add_argument('--proc', type=int, default=30)
    parser.add_argument("csv1", metavar='csv_G1')
    parser.add_argument("csv2", metavar='csv_G2')
    parser.add_argument("table_out")
    parser.add_argument('--reuse_lookups', action='store_true')
    parser.add_argument('--proxy_cols', nargs='+', choices=['freq_adj0', 'size_adj0', 'size_cov'], default=['freq_adj0', 'size_adj0'])
    parser.add_argument('--params_out', nargs='+', choices=['kon', 'koff', 'ksyn', 'burstsize', 'burstfrequency', 'expressionrate', 'meanoccupancy'], default=['burstfrequency', 'burstsize'])
    parser.add_argument('-t', '--time', type=float, required=True)
    parser.add_argument('-d', '--degradationrate', required=True)
    parser.add_argument('--prec', '--MLprecision', type=int, default=[10000], nargs='+')
    parser.add_argument('--genes', nargs='+')
    o = parser.parse_args()
    load_lookups_per_gene = not o.reuse_lookups
    o.proc = [o.proc, 1]
    
    if not strictlyIncreasing(o.permutations): raise Exception
    if not strictlyIncreasing(o.prec): raise Exception
    
    import os
    scriptfolder = os.path.dirname(os.path.realpath(__file__))
    hyp1f1_arb_vector_path = os.path.join(scriptfolder, 'hyp1f1_arb_vector.so')
    if not os.path.exists(hyp1f1_arb_vector_path): raise Exception('Cannot find hyp1f1_arb_vector.so in ' + scriptfolder)
    
    global open_chance_name, close_chance_name, transcribe_chance_name
    open_chance_name, close_chance_name, transcribe_chance_name = 'open_chance', 'close_chance', 'transcribe_chance'
    
    # load data
    df1 = pandas.read_csv(o.csv1, index_col=0, nrows=o.max_genes).T
    df2 = pandas.read_csv(o.csv2, index_col=0, nrows=o.max_genes).T
    
    if not df1.columns.equals(df2.columns): raise Exception("Not the same row labels (genes) in the two input files")
    
    if load_lookups_per_gene:
        freq_size_adj0_kdtree = None
        param_lookup = o.params_to_zeros
        interpolation_points = None
    else:
        freq_size_adj0_kdtree, param_lookup, interpolation_points = build_KD_tree(o.params_to_zeros, False, o.proxy_cols)
    
    try:
        degradationrate = float(o.degradationrate)
    except ValueError:
        dfD = pandas.read_csv(o.degradationrate, index_col=0)
    else:
        dfD = None
    
    # test on frequency of non-zero
    from concurrent import futures
    pool = futures.ProcessPoolExecutor(o.proc[0])
    jobs = []
    for gene in df1:
        if o.genes is not None and gene not in o.genes: continue
        jobs.append(pool.submit(one_gene, gene, df1[gene], df2[gene], freq_size_adj0_kdtree, param_lookup, o.permutations, interpolation_points, o.proc[1], o.time, degradationrate if dfD is None else dfD.loc[gene, 'degradationrate'], hyp1f1_arb_vector_path, o.prec, o.proxy_cols, o.params_out))
    table_out = pandas.DataFrame([job.result() for job in tqdm.tqdm(jobs)])
    
    if all(len(table_out["G1vsG2_"+param+"_P"].dropna())==0 for param in o.params_out):
        from warnings import warn
        warn("No P-values")
    else:
        for param in o.params_out:
            Pcol = "G1vsG2_"+param+"_P"
            if o.genes is None:
                table_out.insert(table_out.columns.get_loc(Pcol)+1, Pcol.replace("_P", "_Padj"),  adjustP(table_out[Pcol]))
            else:
                table_out.insert(table_out.columns.get_loc(Pcol)+1, Pcol.replace("_P", "_Padj"), None)
    table_out.to_csv(o.table_out, sep='\t')
    


def addD(v,D):
    return v*(1-D)+D
def beyondD(v,D):
    return (v-D)/(1-D)
def calc_Aexp_het(extraA,B,C,D):
    # B = observed risk10
    # C = observed risk01
    # D = observed risk00
    # extraA = modelled extra risk11 beyond heterogeneous model sum of B and C under a background of D
    
    if extraA == 0:
        return addD(beyondD(B,D) + beyondD(C,D) - beyondD(B,D)*beyondD(C,D),D)
    else:
        xAexp = beyondD(B,D) + beyondD(C,D) - beyondD(B,D)*beyondD(C,D)
        return addD(xAexp + extraA - extraA*xAexp, D)

def calc_Aexp_same(extraA,B,C,D):
    Aexp = (B+C)/2
    if extraA == 0:
        return Aexp
    else:
        xAexp = beyondD(Aexp,D)
        return addD(xAexp + extraA - extraA*xAexp, D)

def calc_Aexp_R10(extraA,B,C,D):
    Aexp = B
    if extraA == 0:
        return Aexp
    else:
        xAexp = beyondD(Aexp,D)
        return addD(xAexp + extraA - extraA*xAexp, D)

def calc_Aexp_R01(extraA,B,C,D):
    Aexp = C
    if extraA == 0:
        return Aexp
    else:
        xAexp = beyondD(Aexp,D)
        return addD(xAexp + extraA - extraA*xAexp, D)

def calc_Aexp_mult(factorA,B,C,D):
    Aexp = B * C / D
    if factorA == 1:
        return Aexp
    else:
        return factorA * Aexp

def calc_risk(sim):
    A = sim.loc[(sim['TF1']==True) & (sim['TF2']==True), 'geneset'].mean()
    B = sim.loc[(sim['TF1']==True) & (sim['TF2']==False), 'geneset'].mean()
    C = sim.loc[(sim['TF1']==False) & (sim['TF2']==True), 'geneset'].mean()
    D = sim.loc[(sim['TF1']==False) & (sim['TF2']==False), 'geneset'].mean()
    return A,B,C,D

def generic_calc_extra_risk(A,B,C,D, Aexp_func, ratio=False):
    Aexp, Aobs = Aexp_func(1 if ratio else 0, B, C, D), A
    if ratio:
        return Aobs/Aexp if Aexp > 0 else float('inf')
    else:
        return (beyondD(Aobs,D)-beyondD(Aexp,D))/(1-beyondD(Aexp,D)) if Aexp < 1 else Aobs-Aexp

def calc_extra_risk_het(A,B,C,D):
    return generic_calc_extra_risk(A,B,C,D,calc_Aexp_het,False)

def calc_extra_risk_same(A,B,C,D):
    return generic_calc_extra_risk(A,B,C,D,calc_Aexp_same,False)

def calc_extra_risk_mult(A,B,C,D):
    return generic_calc_extra_risk(A,B,C,D,calc_Aexp_mult,True)

def calc_extra_risk_R10(A,B,C,D):
    return generic_calc_extra_risk(A,B,C,D,calc_Aexp_R10,False)

def calc_extra_risk_R01(A,B,C,D):
    return generic_calc_extra_risk(A,B,C,D,calc_Aexp_R01,False)

def parse_contig_table_repr(contig_table_repr):
    return [int(s.strip('[] ')) for s in contig_table_repr.split(',')]

def load_features(feature_file):
    return pandas.read_csv(feature_file, index_col=0, sep='\t')

def load_genes(paths, genedata, id_type):
    genes = list()
    for path in paths:
        if not os.path.exists(path):
            genes.append(path)
        else:
            with (gzip.open if path.endswith('.gz') else open)(path, 'rt') as infh:
                genes.extend([line.strip() for line in infh if line.strip()])
    try:
        return translate_genes(genes, genedata, id_type)
    except GeneTranslationError:
        return set()

def translate_genes(genes, genedata, id_type="guess"):
    genes = [str(gene) for gene in genes]
    if id_type == 'symbol' or (id_type == 'guess' and len(genedata.index.intersection(genes)) >= 1):
        return list(genes)
    elif id_type == 'mixed':
        translations = dict()
        for colname in reversed(genedata.columns):
            translations.update({ID:sym for sym, IDstr in genedata[colname].items() for ID in str(IDstr).split(',')})
        translations.update({sym:sym for sym in genedata.index})
        translated = []
        for gene in genes:
            if gene in translations:
                translated.append(translations[gene])
            else:
                translated.append(None)
                warn('Could not identify '+gene)
        return redundant_to_None(translated)
    else:
        for colname in genedata.columns:
            tested_cols = []
            if (colname == id_type) or (id_type == 'guess' and colname.endswith('_ids')):
                translation = {part:key for key, val in genedata[colname].items() for part in str(val).split(',')}
                translated = redundant_to_None(translation.get(ID, None) for ID in genes)
                if (colname == id_type) or (id_type == 'guess' and not all(t is None for t in translated)):
                    return translated
                tested_cols.append(colname)
        raise GeneTranslationError(tested_cols)

def expand_folder_list_into_file_list(paths:list, recursive:bool) -> list:
    expanded = True
    while expanded:
        expanded = False
        paths_new = []
        for path in paths:
            if os.path.isdir(path):
                expanded = recursive
                paths_new.extend(os.path.join(path, filename) for filename in os.listdir(path))
            else:
                paths_new.append(path)
        paths = paths_new
    return paths


def error(message):
    print("Fatal error:", message, file=sys.stderr)
    exit(1)

def warn(message):
    print("Warning:", message, file=sys.stderr)


'''
metric testing code:

def ratio_logitdiff_v():
	ValS = pandas.Series(numpy.exp(numpy.random.normal(size=1000)))
	ratio = ValS[:100].median()/ValS[100:].median()
	ValR = ValS.rank(method='min', pct=True)
	logitdiff = special.logit(ValR[:100].mean())-special.logit(ValR[100:].mean())
	stat, P = stats.ranksums(ValS[:100], ValS[100:])
	return ratio, logitdiff, stat, ValS[:100], ValS[100:]
	

>>> ratio, logitdiff, stat = zip(*[ratio_logitdiff_v()[:3] for j in range(1000)])
>>> pylab.plot(numpy.log(ratio), stat, 'k.')
[<matplotlib.lines.Line2D object at 0x11d8ddf28>]
>>> pylab.show()
>>> stats.spearmanr(ratio, stat)
SpearmanrResult(correlation=0.8669380985674183, pvalue=4.579020955008183e-304)
>>> pylab.clf()
>>> pylab.plot(logitdiff, stat, 'k.')
[<matplotlib.lines.Line2D object at 0x11d80eb70>]
>>> pylab.show()
>>> pylab.clf()
>>> pylab.plot(ratio, logitdiff, 'k.')
[<matplotlib.lines.Line2D object at 0x11d8522b0>]
>>> pylab.show()
>>> stats.spearmanr(logitdiff, stat)
SpearmanrResult(correlation=0.9999999249999172, pvalue=0.0)
'''

def flip_values(column):
    if column.apply(lambda v: isinstance(v, bool)).all():
        return ~column
    elif column.apply(lambda v: isinstance(v, int) or isinstance(v, float)).all():
        return -column
    else:
        warn("A feature column has values that were neither Boolean, integer nor floating point, will only flip any Boolean values")
        return column.replace({True:False, False:True})

def sign(v):
    return 1 if v>=0 else -1

def id_type_parse(string):
    if string in ('symbol', 'guess', 'mixed') or string.endswith('_ids'): return string
    raise argparse.ArgumentTypeError("valid values are: guess, mixed, symbol, tx_ids, gene_ids or anything ending on _ids")

def nr_list(seq):
    # https://stackoverflow.com/questions/480214/how-do-you-remove-duplicates-from-a-list-whilst-preserving-order
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

def redundant_to_None(seq):
    seen = set()
    seen_add = seen.add
    return [None if (x in seen or seen_add(x)) else x for x in seq]

class GeneTranslationError(KeyError):
    def __init__(self, tested_cols):
        super().__init__()
        self.tested_cols = tested_cols

def get_pvalue(distribution, test_against_factor=False, twotailed=True):
    import numpy
    test_against = 1 if test_against_factor else 0
    p_1tail = numpy.mean([r<=test_against for r in distribution]+[0,1])
    if twotailed:
        return 2*p_1tail if p_1tail <= 0.5 else 2*(1-p_1tail)
    else:
        return p_1tail

def magnitudes(last_one, maximum):
    last_one = min(last_one, maximum)
    total = 0
    while total < maximum:
        yield min(last_one, maximum) - total, min(last_one, maximum)
        total = last_one
        last_one = last_one*10

def calc_interaction(interactiondata, o, n_per_quarter):
    risks = calc_risk(interactiondata)
    expectedR11_add = calc_Aexp_het(0, *risks[1:])
    if min(expectedR11_add, 1-expectedR11_add) * n_per_quarter[0] < o.interaction_settings[3]:
        return None
    if min(risks[1:3]) < risks[3]: # protective factors in the absense of the other will violate formula assumptions
        return None
    import numpy
    
    extraA_add_from_bootstrap = []
    extraA_R10_from_bootstrap = []
    extraA_R01_from_bootstrap = []
    factorA_mult_from_bootstrap = []
    for additional_bootstrap_repl_num, total_bootstrap_repl_num in magnitudes(100, int(o.interaction_settings[2])):
        for j in range(additional_bootstrap_repl_num):
            bootstrap_repl_risks = calc_risk(interactiondata.sample(frac=1, replace=True)) # this line probably takes some time
            v = calc_extra_risk_het(*bootstrap_repl_risks)
            if -1 <= v <= 1: extraA_add_from_bootstrap.append(v)
            v = calc_extra_risk_R01(*bootstrap_repl_risks)
            if -1 <= v <= 1: extraA_R01_from_bootstrap.append(v)
            v = calc_extra_risk_R10(*bootstrap_repl_risks)
            if -1 <= v <= 1: extraA_R10_from_bootstrap.append(v)
            v = calc_extra_risk_mult(*bootstrap_repl_risks)
            if str(v) != 'nan': factorA_mult_from_bootstrap.append(v)
        
        p_values_so_far = [get_pvalue(extraA_R10_from_bootstrap, False, True), get_pvalue(extraA_R01_from_bootstrap, False, True), get_pvalue(extraA_add_from_bootstrap, False, True), get_pvalue(factorA_mult_from_bootstrap, True, True)]
        if min(pval*total_bootstrap_repl_num*factor for pval,factor in zip(p_values_so_far, o.interaction_bootstrap_ratios)) > 20:
            break
    
    newrow = dict()
    try:
        newrow.update({name+'_eq10':val for name,val in zip(['median'], numpy.quantile(extraA_R10_from_bootstrap, [0.5]))})
        newrow.update({name+'_eq01':val for name,val in zip(['median'], numpy.quantile(extraA_R01_from_bootstrap, [0.5]))})
        newrow.update({name+'_add':val for name,val in zip(['CI95low', 'CI90low', 'median', 'CI90high', 'CI95high'], numpy.quantile(extraA_add_from_bootstrap, [0.025, 0.05, 0.5, 0.95, 0.975]))})
        newrow.update({name+'_mult':val for name,val in zip(['median'], numpy.quantile(factorA_mult_from_bootstrap, [0.5]))})
    except IndexError:
        return None
    newrow.update(dict(zip(['n11', 'n10', 'n01', 'n00'], n_per_quarter)))
    newrow.update(dict(zip(['R11', 'R10', 'R01', 'R00'], risks)))
    newrow['R11expected_add'] = expectedR11_add
    newrow['R11expected_mult'] = calc_Aexp_mult(1, *risks[1:])
    newrow['scaled_add0_mult1'] = (newrow['R11']-newrow['R11expected_add'])/( newrow['R11expected_mult']- newrow['R11expected_add'])
    newrow['RERI'] = 1 + (newrow['R11'] - newrow['R10'] - newrow['R01'])/newrow['R00']
    newrow['Attributable_proportion_AP'] = newrow['RERI'] * newrow['R00'] / newrow['R11']
    newrow['Synergy_index_S'] = (newrow['R11']/newrow['R00'] - 1)/(newrow['R10']/newrow['R00'] + newrow['R01']/newrow['R00'])
    newrow['p_nominal_eq10_2tailed'] = get_pvalue(extraA_R10_from_bootstrap, False, True)
    newrow['p_nominal_eq01_2tailed'] = get_pvalue(extraA_R01_from_bootstrap, False, True)
    newrow['p_nominal_add_2tailed'] = get_pvalue(extraA_add_from_bootstrap, False, True)
    newrow['p_nominal_mult_2tailed'] = get_pvalue(factorA_mult_from_bootstrap, True, True)
    newrow['sort_by_tmp'] = get_pvalue(extraA_add_from_bootstrap, False, False)
    return newrow

def calc_interaction_status(o, row):
    def check_p(nullhypothesis):
        return row['padj_%s_2tailed'%nullhypothesis] > 0.05
    def check_up(nullhypothesis, mid):
        return row['median_'+nullhypothesis] > mid
    def is_lowest_p(nullhypothesis):
        return row['p_nominal_%s_2tailed'%nullhypothesis] <= 2/o.interaction_settings[2]
    possibly_mult = check_p('mult')
    possibly_add = check_p('add')
    possibly_eq10 = check_p('eq10')
    possibly_eq01 = check_p('eq01')
    above_add = check_up('add', 0)
    above_mult = check_up('mult', 1)
    above_R10 = check_up('eq10', 0)
    above_R01 = check_up('eq01', 0)
    
    if any(is_lowest_p(nullhypothesis) and check_p(nullhypothesis) for nullhypothesis in ('eq10', 'eq01', 'add', 'mult')):
        return 0, 'too few bootstrap replicates'
    elif possibly_add:
        if possibly_mult or possibly_eq10 or possibly_eq01:
            return 11, 'multiple model fit: insufficient data'
        else:
            return 4, 'additive fit: independent contributions'
    elif above_add:
        if possibly_mult:
            return 2, 'multiplicative fit: synergy'
        elif above_mult:
            return 1, 'positive deviation from multiplicative: synergy'
        else:
            return 3, 'between additive and multiplicative: synergy'
    else:
        if possibly_eq01 and possibly_eq10:
            return 9, 'R11 matches R10 or R01: single contribution'#'R11 matches R10 or R01: redundancy'
        elif possibly_eq01:
            return 6, 'R11 matches R01: single contribution'#'R11 matches R01: redundancy'
        elif possibly_eq10:
            return 7, 'R11 matches R10: single contribution'#'R11 matches R10: redundancy'
        elif above_R10 and above_R01:
            return 5, 'below additive and above R10 and R01: partial redundancy'
        elif above_R10 or above_R01:
            return 8, 'between R10 and R01: redundancy'
        else:
            return 10, 'below R10 and R01: anti-synergy'

def add_padj(interaction_out, nullhypothesis, o):
    interaction_out['padj_%s_2tailed'%nullhypothesis] = multitest.multipletests(interaction_out['p_nominal_%s_2tailed'%nullhypothesis], is_sorted=False, method='fdr_bh')[1]
    if interaction_out['p_nominal_%s_2tailed'%nullhypothesis].min() < 10/o.interaction_settings[2] and interaction_out['padj_%s_2tailed'%nullhypothesis].min() > 0.05:
        warn("Some results affected by too few bootstrap replicates (for %s)"%nullhypothesis)

def test_quantities(A, B, o):
    if o.quanttest == 'KolmogorovSmirnov':
        return stats.ks_2samp(A,B)
    elif o.quanttest == 'Ksamp_AndersonDarling':
        import warnings
        warnings.filterwarnings("ignore", message="p-value capped: true value larger than 0.25")
        warnings.filterwarnings("ignore", message="p-value floored: true value smaller than 0.001")
        stat, crit, P = stats.anderson_ksamp([A, B])
        return stat, P
    elif o.quanttest == 'MannWhitneyU':
        return stats.mannwhitneyu(A, B, use_continuity=False, alternative='two-sided')
    elif o.quanttest == 'Wilcoxon_ranksums':
        return stats.ranksums(A,B)
    else:
        raise ValueError(o.quanttest + ' option not found in test_quantities function')

'''
def biserial(cat, cont):
    try:
        from pypair import association
    except:
        try:
            has_warned = biserial.has_warned
        except AttributeError:
            warn("Could not find the pypair module, using point-biserial correlation instead of rank-biserial correlation")
            biserial.has_warned = True
        from scipy import stats
        return stats.pointbiserialr(cat.astype(bool), cont.astype(float))[0]
    else:
        return association.binary_continuous(cat.astype(bool), cont.astype(float), 'rank_biserial')
'''
def biserial(cat, cont):
    from scipy import stats
    return stats.pointbiserialr(cat.astype(bool), cont.astype(float))[0]

def one_minus_correlation_customizable(V1, V2, allow_mixed, default_value=1):
    # V1 and V2 are equal-length vectors, either all-Boolean or all-numerical
    table = pandas.DataFrame({'V1':V1, 'V2':V2}).dropna()
    if len(table) == 0:
        return default_value
    V1_categorical = table.V1.apply(lambda v: isinstance(v, bool)).all()
    V2_categorical = table.V2.apply(lambda v: isinstance(v, bool)).all()
    
    if V1_categorical and V2_categorical:
        from sklearn import metrics
        return 1-metrics.matthews_corrcoef(table.V1.astype(bool), table.V2.astype(bool))
    elif (not V1_categorical) and (not V2_categorical):
        from scipy import stats
        return 1-stats.spearmanr(table.V1.astype(float), table.V2.astype(float))[0]
    elif not allow_mixed:
        return default_value
    elif V1_categorical:
        return 1-biserial(table.V1, table.V2)
    else:
        return 1-biserial(table.V2, table.V1)

def one_minus_correlation(V1, V2):
    return one_minus_correlation_customizable(V1, V2, True)

def one_minus_correlation_same_type(V1, V2):
    return one_minus_correlation_customizable(V1, V2, False)


def feature_clustering_panel(featuretable, plot_values=False):
    # featuretable is a subset of genedata, with features as columns and genes as rows
    import seaborn
    from scipy import cluster, spatial
    condensed_dists = spatial.distance.pdist(featuretable.T, metric=one_minus_correlation)
    try:
        linked_cols = cluster.hierarchy.linkage(condensed_dists, method='complete')
    except ValueError:
        warn('(feature_clustering_plot) Linkage failed for distance matrix '+ repr(condensed_dists))
        return
    if plot_values:
        vals = featuretable.astype(float).divide(featuretable.max().apply(lambda m: max(1, m)))
        seaborn.clustermap(vals, col_linkage=linked_cols, row_cluster=False)
        return vals
    else:
        square_dists = pandas.DataFrame(1-spatial.distance.squareform(condensed_dists), columns=featuretable.columns, index=featuretable.columns)
        seaborn.clustermap(square_dists, row_linkage=linked_cols, col_linkage=linked_cols, center=0, vmax=1, vmin=-1, yticklabels=True, xticklabels=True)
        return square_dists

def average_feature_values(generow):
    generow = generow.dropna()
    if len(generow) == 0:
        return None
    midval = generow.median()
    if generow.apply(lambda v: isinstance(v, bool)).all():
        return bool(math.ceil(midval))
    else:
        return midval

def hcluster(featuretable, o):
    from scipy import cluster, spatial
    condensed_dists = spatial.distance.pdist(featuretable.T, metric=one_minus_correlation)
    try:
        linked_cols = cluster.hierarchy.linkage(condensed_dists, method='complete')
    except ValueError:
        return []
    return zip(featuretable.columns, cluster.hierarchy.fcluster(linked_cols, 1-o.cluster_cutoff))

if '__main__' == __name__:
    
    default_addition_database = '/home/danielr/from_crick2/results/features_of_genes/human/added_gene_lists.tsv'
    
    import argparse, pandas, gzip, sys, os, math, random, tqdm
    from scipy import stats, special
    from statsmodels.stats import multitest
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '-fg', '--geneset_file', nargs='+', required=True, help='Plain text file with a gene identifier on each line, as test set', action='append')
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-b', '-bg', '--background_geneset_file', nargs='+', help='Plain text file with a gene identifier on each line, as background set', action='append')
    parser.add_argument('-F', '-f', '--gene_feature_files', '--gene_features', default=['/mnt/davidson/danielr/results/features_of_genes/human/'], nargs='+', dest='gene_features')
    parser.add_argument('--limit_to_features', nargs='+')
    parser.add_argument('-t', '--id_type', metavar="{symbol, gene_ids, tx_ids, guess, mixed, *_ids}", default='mixed', type=id_type_parse)
    parser.add_argument('-s', '--sortP', action='store_true', help='Sort by P value')
    parser.add_argument('-r', '--ranked', action='store_true')
    parser.add_argument('-a', '--add_to_database', nargs='?', const=default_addition_database)
    parser.add_argument('-n', '--geneset_name')
    parser.add_argument('-u', '--update_database_entry', action='store_true')
    parser.add_argument('-e', '--expression_file', action='append', help='tab-separated plain text file with genes in the first column and a header row')
    parser.add_argument('-E', '--export_background', '--export_expression_matching')
    parser.add_argument('-P', '--params_expr_matching', nargs=4, type=float, metavar=('max_rank_difference', 'max_unmatched_genes_min', 'max_unmatched_genes_fraction', 'max_iterations'), default=[100, 2, 0.1, float('inf')])
    parser.add_argument('--rm_expr_unmatched', action='store_true')
    parser.add_argument('-S', '--background_subsampling', '--expression_match_subsampling', type=float, metavar='fraction_of_genes')
    parser.add_argument('-v', '--volcano_plot', '--vulcano_plot', dest='vulcano_plot')
    parser.add_argument('-B', '--Pbars_plot')
    parser.add_argument('-C', '--feature_clustering_plotprefix')
    parser.add_argument('-c', '--cluster_cutoff', type=float, metavar='correlation', default=None, nargs='?', const=0.7)
    parser.add_argument('--clustering_set', choices=['input', 'all'], default='all')
    parser.add_argument('--FDR_in_plots', type=float, default=0.05)
    parser.add_argument('--list_hits', action='store_true')
    parser.add_argument('--filter_on_TF', nargs='+')
    parser.add_argument('--filter_on_TF_combine', choices=['intersection', 'union'], default='union')
    parser.add_argument('--filter_on_TF_matching', choices=['exact', 'word', 'partial', 'startswith', 'regularexpression'], default='word')
    parser.add_argument('--random_geneset', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--interaction_output')
    parser.add_argument('--interaction_settings', nargs=6, metavar=('padj_cutoff_single_factor', 'min_n_per_group', 'bootstrap_replicates', 'min_expected_true11', 'max_factors', 'num_processes'), type=float, default=[0.05, 50, 10000, 10, float('inf'), 20], help='default: 0.05 50 10000 10 inf 20')
    parser.add_argument('--interaction_bootstrap_ratios', metavar=('eq10', 'eq01', 'add', 'mult'), default=[0.01, 0.01, 1, 1], nargs=4, type=float)
    parser.add_argument('--flip_feature_values', action='store_true')
    parser.add_argument('--flip_depleted_for_interaction', action='store_true')
    parser.add_argument('--quanttest', choices=['Wilcoxon_ranksums', 'MannWhitneyU', 'KolmogorovSmirnov', 'Ksamp_AndersonDarling'], default='Wilcoxon_ranksums')
    parser.add_argument('--combine_by_clustering', type=float, metavar='correlation', nargs='?', const=0.7, help='slow, can take over an hour')
    parser.add_argument('--export_feature_table', help=argparse.SUPPRESS)
    o = parser.parse_args()
    
    o.interaction_bootstrap_ratios = [v/max(o.interaction_bootstrap_ratios) for v in o.interaction_bootstrap_ratios] # should contain 1, and might have lower values including zero
    
    if o.add_to_database is None and (o.geneset_name is not None or o.update_database_entry):
        o.add_to_database = default_addition_database
    
    if o.ranked and (o.expression_file is not None):
        warn("The expression matching for background does not do anything for ranked input (--ranked and -e in arguments)")
    
    if o.ranked and (o.background_geneset_file is not None):
        warn("Tests with ranked input are not against background (--ranked and -b in arguments)")
    
    if o.ranked:
        if pandas.Index([5, 8, 9, 1, 2]).intersection(pandas.Index([9, 1, 8])).equals(pandas.Index([8, 9, 1])):
            pandas_newer_intersection_order = True
        elif pandas.Index([5, 8, 9, 1, 2]).intersection(pandas.Index([2, 1, 9, 10])).equals(pandas.Index([2, 1, 9])):
            pandas_newer_intersection_order = False
            warn("Switching code to support an outdated version of the python module pandas")
        else:
            error("Your installed version of the python module pandas is not supported")
    else:
        pandas_newer_intersection_order = True
    
    gene_features_files = expand_folder_list_into_file_list(o.gene_features, recursive=True)
    genedata = load_features(gene_features_files[0])
    for feature_file in gene_features_files[1:]:
        nextdata = load_features(feature_file)
        nextdata = nextdata[nextdata.columns.difference(genedata.columns)]
        genedata = genedata.merge(nextdata, how='outer', left_index=True, right_index=True)
    
    if o.limit_to_features:
        genedata = genedata[[col for col in genedata.columns if col.endswith('_ids')] + o.limit_to_features]
    
    if len(genedata) == 0:
        error("Empty feature database")
    
    if len(genedata.dropna())==0:
        warn("No gene identifier is found in every column in gene feature databases, could be sign of identifier mismatch between the files from -F")
    
    if o.filter_on_TF is not None:
        use_cols = set()
        for searchterm in o.filter_on_TF:
            if o.filter_on_TF_matching == 'exact':
                if searchterm in genedata.columns:
                    use_cols.add(searchterm)
                else:
                    warn('(--filter_on_TF) Did not find '+searchterm)
            else:
                if o.filter_on_TF_matching == 'word':
                    found_col = set(col for col in genedata.columns if searchterm in col.replace('_', ' ').split())
                elif o.filter_on_TF_matching == 'partial':
                    found_col = set(col for col in genedata.columns if searchterm in col)
                elif o.filter_on_TF_matching == 'startswith':
                    found_col = set(col for col in genedata.columns if col.startswith(searchterm))
                elif o.filter_on_TF_matching == 'regularexpression':
                    import re
                    prog = re.compile(searchterm)
                    found_col = set(col for col in genedata.columns if prog.match(col))
                if found_col:
                    use_cols.update(found_col)
                else:
                    warn('(--filter_on_TF) Did not match '+searchterm)
        if not use_cols:
            error('(--filter_on_TF) Did not match any columns in the gene feature database')
        if o.filter_on_TF_combine == 'intersection':
            gene_filter = genedata[use_cols].all(axis=1)
        elif o.filter_on_TF_combine == 'union':
            gene_filter = genedata[use_cols].any(axis=1)
        if gene_filter.all():
            warn('(--filter_on_TF) No gene filtering being done')
        elif not gene_filter.any():
            error("(--filter_on_TF) Removed all genes")
        cols_remaining = genedata.columns.difference(use_cols)
        if cols_remaining.empty:
            error("(--filter_on_TF) Removed all feature columns")
        genedata = genedata.loc[gene_filter, cols_remaining]
        
    
    genes = None
    for geneset_files in o.geneset_file:
        genes_l = load_genes(geneset_files, genedata, o.id_type)
        genes_l = pandas.Index(genes_l)
        genes = genes_l if genes is None else genes.intersection(genes_l)
    
    if o.background_geneset_file is not None:
        bg_genes = None
        for background_geneset_files in o.background_geneset_file:
            bg_genes_l = load_genes(background_geneset_files, genedata, o.id_type)
            if any('\t' in gene for gene in bg_genes_l if gene is not None):
                warn('There are tabs in the background file, note that -b will treat whole rows as gene identifiers, it is -e that takes a tab-delimited file')
            bg_genes_l = pandas.Index(bg_genes_l)
            bg_genes = bg_genes_l if bg_genes is None else bg_genes.intersection(bg_genes_l)
    
    if o.add_to_database is not None:
        if o.geneset_name is None and len(o.geneset_file) > 1:
            error("You need to set -n/--geneset_name to use -a/--add_to_database with multiple input")
        geneset_name = o.geneset_name if o.geneset_name is not None else o.geneset_file[0]
        
        if os.path.exists(o.add_to_database):
            existing = load_features(o.add_to_database)
        else:
            existing = pandas.DataFrame()
        
        if not o.update_database_entry and geneset_name in existing.columns:
            warn("Will not add {} to {} as the column already exists".format(geneset_name, o.add_to_database))
        else:
            if o.update_database_entry:
                del existing[geneset_name]
            
            if o.ranked:
                newcol = pandas.Series(range(len(genes)), index=genes)
            else:
                newcol = pandas.Series(False, index=bg_genes.union(genes) if o.background_geneset_file is not None else genedata.index)
                newcol.loc[genes.intersection(newcol.index) if pandas_newer_intersection_order else newcol.index.intersection(genes)] = True
            newcol.name = geneset_name
            newcol = newcol.to_frame()
            
            if not existing.empty:
                new_database = existing.merge(newcol, how='outer', left_index=True, right_index=True)
            else:
                new_database = newcol
            new_database.to_csv(o.add_to_database, sep='\t')
    
    genes = genes.intersection(genedata.index) if pandas_newer_intersection_order else genedata.index.intersection(genes)
    if len(genes)==0:
        error("Input gene set did not overlap genes in gene feature database")
    
    if o.background_geneset_file is not None:
        bg_genes = bg_genes.intersection(genedata.index)
        if len(bg_genes)==0:
            error("The background gene set did not overlap genes in gene feature database")
        
        if not genes.isin(bg_genes).all():
            warn("Not all genes in -i are in the background list -b, adding them to the background")
            bg_genes = bg_genes.union(genes)
        
        genedata = genedata.loc[bg_genes, :]
    
    if o.expression_file is not None:
        bg_genes2 = None
        
        for expression_file in o.expression_file:
            
            expression = pandas.read_table(expression_file, index_col=0)
            if len(expression.columns) == 0:
                error("No columns in the expression file beyond the index, is it really tab-separated?")
            expression = expression.mean(axis = 1, skipna = True)
            try:
                expression.index = translate_genes(expression.index, genedata, o.id_type)
            except GeneTranslationError as err:
                print(err.tested_cols)
                error("The row indexes in the expression table are of an unrecognized gene identifier type, not found in: %s, if this doesn't include the identifier type you might be ale to add them using -F and a relevant file"%(", ".join(["symbols"]+err.tested_cols)))
            expression = expression.reset_index().dropna().set_index('index') # dropping genes that couldn't be translated to gene symbols
            expression.index.name = None
            expression = expression[0]
            
            
            gene_to_rank = expression.sort_values().rank(method="average")
            max_rank_difference = int(o.params_expr_matching[0])
            max_unmatched_genes = max(int(o.params_expr_matching[1]), len(genes) * o.params_expr_matching[2])
            max_iterations_expr_matching = o.params_expr_matching[3]
            
            matches = []
            unmatched = set()
            missing = set()
            search_sign = -1
            new_matches = list(genes)
            iterations_done = 0
            import collections; matching_log = collections.defaultdict(list)
            while len(unmatched) < max_unmatched_genes or iterations_done <= 1:
                if len(new_matches) == 0: break
                if iterations_done >= max_iterations_expr_matching: break
                iterations_done += 1
                matches.extend(new_matches)
                new_matches = []
                for gene in genes:
                    if gene in unmatched or gene in missing: continue
                    try:
                        geneloc = gene_to_rank.index.get_loc(gene)
                    except KeyError:
                        missing.add(gene)
                        continue
                    if geneloc == 0: search_sign = 1
                    elif geneloc+1 == len(gene_to_rank): search_sign = -1
                    label = gene_to_rank.index[geneloc+search_sign]
                    search_sign = -search_sign
                    search_dist = 1
                    
                    while label is None or label in genes or (abs(gene_to_rank[label] - gene_to_rank[gene]) > max_rank_difference):
                        if geneloc == 0: search_sign = 1
                        elif geneloc+search_dist >= len(gene_to_rank): search_sign = -1
                        label = gene_to_rank.index[geneloc+search_sign*search_dist]
                        if search_sign > 0 or geneloc+search_dist >= len(gene_to_rank): search_dist += 1
                        search_sign = -search_sign
                        if search_dist > 50 or (search_dist > 3 and label not in genes): break
                    else:
                        new_matches.append(label)
                        gene_to_rank.drop(label, inplace=True)
                        matching_log[gene].append(label)
                        continue
                    unmatched.add(gene)
            
            bg_genes2_l = pandas.Index(matches)
            bg_genes2_l = bg_genes2_l.intersection(genedata.index)
            if len(bg_genes2_l)==0:
                error("The background gene set, derived from the expression file, did not overlap genes in gene feature database or in the given background set")
            if len(bg_genes2_l) == len(genes):
                error("Failed to generate an expression-matched background gene set")
            
            if o.rm_expr_unmatched:
                genes = genes.difference(missing)
                bg_genes2_l = bg_genes2_l.difference(missing)
            
            bg_genes2 = bg_genes2_l if bg_genes2 is None else bg_genes2.intersection(bg_genes2_l)
        
        if len(bg_genes2)==0:
            error("The combined background gene set, derived from the expression files, did not overlap genes in gene feature database or in the given background set")
        if len(bg_genes2) == len(genes):
            error("Failed to generate a combined expression-matched background gene set, equals the positive set")
        
        genedata = genedata.loc[bg_genes2, :]
        
        
    if o.background_subsampling:
        # take 100*o.expression_match_subsampling % of the expression matched genes (which are bg_genes2 minus genes) as a new smaller negative set
        expression_matched_genes = genedata.index.difference(genes)
        sampled_genes = random.sample(list(expression_matched_genes), round(len(expression_matched_genes)*o.background_subsampling))
        bg_genes3 = pandas.Index(sampled_genes).union(genes)
        genedata = genedata.loc[bg_genes3, :]
    
    if o.export_background is not None:
        with open(o.export_background, 'wt') as outfh:
            if o.expression_file is not None and len(o.expression_file)==1:
                for gene in genes:
                    print(gene, file=outfh)
                    for bggene in matching_log[gene]:
                        if bggene in genedata.index:
                            print(bggene, file=outfh)
            else:
                for gene in genedata.index:
                    print(gene, file=outfh)
    
    if o.random_geneset:
        genes = pandas.Index(random.sample(list(genedata.index), len(genes)))
    
    if o.combine_by_clustering is not None:
        from scipy import cluster, spatial
        import collections
        featuretable = genedata[[col for col in genedata.columns if not col.endswith('_ids')]]
        if o.clustering_set != 'all':
            featuretable = featuretable.loc[genes, :]
        condensed_dists = spatial.distance.pdist(featuretable.T, metric=one_minus_correlation_same_type)
        linked_cols = cluster.hierarchy.linkage(condensed_dists, method='complete')
        genedata = pandas.DataFrame(index=genedata.index)
        clusternum_to_features = collections.defaultdict(list)
        for feature, clusternum in zip(featuretable.columns, cluster.hierarchy.fcluster(linked_cols, 1-o.combine_by_clustering)):
            clusternum_to_features[clusternum].append(feature)
        for featurenum, features in clusternum_to_features.items():
            genedata['+'.join(features)] = featuretable[features].apply(average_feature_values, axis=1)
        genedata = genedata.dropna(how='all')
    
    if o.flip_feature_values:
        for colname in genedata.columns:
            if colname.endswith('_ids'): continue
            genedata[colname] = flip_values(genedata[colname])
    flippedI = {f:False for f in genedata.columns}
    def name_factor(f):
        return '!'+f if o.flip_feature_values ^ flippedI[f] else f
    
    if o.export_feature_table is not None:
        genedata.to_csv(o.export_feature_table, sep='\t')
    
    results = pandas.DataFrame()
    
    for colname in tqdm.tqdm(genedata.columns, desc='Factor enrichment'):
        if colname.endswith('_ids'): continue
        extracolumns = dict()
        
        column = genedata[colname].dropna()
        
        if len(column.unique()) == 1: continue
        
        pos_set = column[genes.intersection(column.index) if pandas_newer_intersection_order else column.index.intersection(genes)]
        if o.ranked:
            columnranked = pandas.Series([(v+0.5)/len(pos_set) for v in range(len(pos_set))], index=pos_set.index)
        else:
            neg_set = column[column.index.difference(genes)]
        
        if column.apply(lambda v: isinstance(v, bool)).all():
            # column only has True/False values
            
            if o.ranked:
                pos_set_ranks = columnranked[column]
                neg_set_ranks = columnranked[~column]
                logitdiff = special.logit(pos_set_ranks.mean()) - special.logit(neg_set_ranks.mean())
                ratio = float('nan')
                stat, P = test_quantities(pos_set_ranks, neg_set_ranks, o)
                
                if o.list_hits:
                    extracolumns['hitlist'] = ','.join([gene for gene in (genes.intersection(column.index) if pandas_newer_intersection_order else column.index.intersection(genes)) if column[gene]])
            else:
                cont_table = [[list(pos_set).count(True), list(pos_set).count(False)], [list(neg_set).count(True), list(neg_set).count(False)]]
                if cont_table[0][0] == 0 and cont_table[1][0] == 0:
                    logitdiff = 0
                    ratio = 1
                    P = 1
                elif sum(cont_table[0]) == 0:
                    warn("No genes for testing "+colname)
                    continue
                elif sum(cont_table[1]) == 0:
                    warn("No background (excluding the input set) for testing "+colname)
                    continue
                else:
                    if all(v>=10 for A in cont_table for v in A) and all(v>=10 for A in stats.contingency.expected_freq(cont_table) for v in A):
                        chi2, P, dof, expected = stats.chi2_contingency(cont_table, True)
                        try: oddsratio = cont_table[0][0]*cont_table[1][1]/cont_table[0][1]/cont_table[1][0]
                        except ZeroDivisionError: oddsratio = float('inf')
                    else:
                        oddsratio, P = stats.fisher_exact(cont_table, 'two-sided')
                    logitdiff = special.logit(cont_table[0][0]/sum(cont_table[0])) - special.logit(cont_table[1][0]/sum(cont_table[1]))
                    ratio = oddsratio
                    
                    if o.list_hits:
                        try:
                            extracolumns['relativerisk'] = cont_table[0][0]*sum(cont_table[1])/sum(cont_table[0])/cont_table[1][0]
                        except ZeroDivisionError:
                            extracolumns['relativerisk'] = None
                    if o.list_hits or o.interaction_output is not None:
                        extracolumns['contig_table'] = repr(cont_table)
                    if o.list_hits:
                        extracolumns['hitlist'] = ','.join([gene for gene in (genes.intersection(column.index) if pandas_newer_intersection_order else column.index.intersection(genes)) if column[gene]])
        elif column.apply(lambda v: isinstance(v, int) or isinstance(v, float)).all():
            # column only has numerical values
            
            if o.list_hits:
                extracolumns['pos_set'] = str(pos_set.to_dict())
            
            if o.ranked:
                stat, P = stats.spearmanr(columnranked, pos_set)
                ratio = '=correlation'
                logitdiff = stat
            else:
                stat, P = test_quantities(pos_set, neg_set, o)
                
                if o.list_hits:
                    extracolumns['neg_set'] = str(neg_set.to_dict())
                
                columnranked = column.rank(method='min', pct=True)
                pos_set_ranks = columnranked[column.index.intersection(genes)]
                neg_set_ranks = columnranked[column.index.difference(genes)]
                logitdiff = special.logit(pos_set_ranks.mean()) - special.logit(neg_set_ranks.mean()) # correlates well and linearly with log(ratio), and gets the same scale as log(ratio) [when using lognormal distribution as exp(normal())], i.e. 0.3 as logitdiff means ~0.3 as log(ratio); correlates perfectly linearly with stat from stats.ranksums, i.e. consult this value to see the direction. So: as correct as the test metric, and nearly as interpretable as the metric since the value is so similar to ln(ratio)
                if (column >= 0).all():
                    ratio = pos_set.median()/neg_set.median()
                    if (ratio > 1  and logitdiff < 0) or (ratio < 1  and logitdiff > 0):
                        ratio = float('nan')
                else:
                    ratio = float('nan')
            
        else:
            warn("The variable type for column "+ colname+ " is neither boolean nor numerical")
            continue
        
        results.loc[colname, 'P'] = P
        results.loc[colname, 'logitdiff'] = logitdiff
        results.loc[colname, 'ratio'] = ratio
        
        for key, val in extracolumns.items():
            results.loc[colname, key] = val
        
    
    if o.sortP:
        results = results.sort_values(by='P')
    
    nonNaN_P = results['P'].dropna()
    results['Padj_tsBH'] = pandas.Series(multitest.multipletests(nonNaN_P, is_sorted=False, method='fdr_tsbh')[1], index=nonNaN_P.index)
    results_out = results.copy()
    if o.cluster_cutoff:
        results_out['cluster'] = None
        features_to_show = results.index[results['Padj_tsBH'] <= o.FDR_in_plots]
        for feature, cluster in hcluster((genedata.loc[:, features_to_show] if o.clustering_set == 'all' else genedata.loc[genes, features_to_show]), o):
            results_out.loc[feature, 'cluster'] = cluster
    results_out.index = map(name_factor, results_out.index)
    results_out.to_csv(o.output, sep='\t')
    
    if o.feature_clustering_plotprefix is not None:
        import matplotlib
        matplotlib.use('Agg')
        from matplotlib import pyplot
        
        features_to_show = results.index[results['Padj_tsBH'] <= o.FDR_in_plots]
        if len(features_to_show) < 2:
            warn('Too few significant results to plot as cluster map')
        else:
            
            #genedata.loc[genes, features_to_show].to_csv(o.feature_clustering_plotprefix+'dump.csv')
            pyplot.clf()
            feature_clustering_panel(genedata.loc[genes, features_to_show])
            pyplot.gcf().suptitle('Input genes, correlation matrix')
            pyplot.savefig(o.feature_clustering_plotprefix+'1.pdf')
            tab.to_csv(o.feature_clustering_plotprefix+'1.csv')
            
            pyplot.clf()
            tab = feature_clustering_panel(genedata.loc[:, features_to_show])
            pyplot.gcf().suptitle('Input+background genes, correlation matrix')
            pyplot.savefig(o.feature_clustering_plotprefix+'2.pdf')
            tab.to_csv(o.feature_clustering_plotprefix+'2.csv')
            
            pyplot.clf()
            feature_clustering_panel(genedata.loc[genes, features_to_show], True)
            pyplot.gcf().suptitle('Input genes, data matrix')
            pyplot.savefig(o.feature_clustering_plotprefix+'3.pdf')
            pyplot.clf()
            feature_clustering_panel(genedata.loc[:, features_to_show], True)
            pyplot.gcf().suptitle('Input+background genes, data matrix')
            pyplot.savefig(o.feature_clustering_plotprefix+'4.pdf')
        
    
    if o.vulcano_plot is not None:
        import matplotlib
        matplotlib.use('Agg')
        from matplotlib import pyplot
        import seaborn, numpy, random
        random.seed(0)
        df = results_out
        if any(math.isinf(v) for v in df['ratio']):
            warn("(volcano plot) Some ratios are infinite, these might trigger a warning about posx and posy and that these data points are not plotted")
        if any(v == 0 for v in df['ratio']):
            warn("(volcano plot) Some ratios are zero, these might trigger a warning about division by zero and log2, and posx and posy, and that these data points are not plotted")
        df['log2_ratio'] = numpy.log2(df['ratio'])
        df['-log10_P'] = -numpy.log10(df['P'])
        seaborn.scatterplot(data=df, x='log2_ratio', y='-log10_P')
        highest_significant_P = df.loc[df['Padj_tsBH']<=o.FDR_in_plots, '-log10_P'].min()
        pyplot.axhline(highest_significant_P, linestyle=':')
        text_elements = []
        dfPsort = df.sort_values(by='-log10_P', ascending=False)
        for s, x, y in zip(dfPsort.index, dfPsort['log2_ratio'], dfPsort['-log10_P']):
            if y >= highest_significant_P:
                col = '#' + ('%02x'%random.randint(0,80))*3
                text_elements.append(pyplot.text(x, y, s, fontname='Arial', fontsize=6, color=col))
        '''try:
            import adjustText
        except ImportError:
            warn("Could not import the module adjustText")
        else:
            adjustText.adjust_text(text_elements[:20])'''
        pyplot.savefig(o.vulcano_plot)
    
    if o.Pbars_plot is not None:
        import matplotlib
        matplotlib.use('Agg')
        from matplotlib import pyplot
        import seaborn, numpy
        df = results_out[results_out['Padj_tsBH']<=o.FDR_in_plots].copy()
        if df.empty:
            warn("No bar plot created since there are no results at {}% FDR".format(100*o.FDR_in_plots))
        else:
            df['Test_done'] = df.index
            df['-log10_P'] = -numpy.log10(df['P'])
            df['log2_ratio'] = numpy.log2(df['ratio'])
            highest_significant_P = df['-log10_P'].min()
            df = df.head(20)
            pyplot.clf()
            zeros = [0 for i in range(len(df))]
            try: pathcollection = pyplot.scatter(zeros, zeros, c=df['log2_ratio'], cmap='flare')
            except ValueError: pathcollection = pyplot.scatter(zeros, zeros, c=df['log2_ratio'], cmap='rocket_r')
            pyplot.colorbar().set_label('log2_ratio', rotation=270)
            pathcollection.remove()
            for palette in ('flare', 'rocket_r', None):
                try: seaborn.barplot(data=df, x='-log10_P', y='Test_done', hue='log2_ratio', palette=palette, dodge=False)
                except ValueError:
                    if palette is None:
                        print(palette, df.head())
                        raise
                else:
                    break
            pyplot.gca().get_legend().set_visible(False)
            highest_significant_P = df.loc[df['Padj_tsBH']<=o.FDR_in_plots, '-log10_P'].min()
            pyplot.axvline(highest_significant_P, linestyle=':')
            pyplot.tight_layout()
            pyplot.savefig(o.Pbars_plot)
    
    del results_out
    
    if o.interaction_output is not None:
        # pick factors that fit all assumptions: enriched, binary, and big enough values to bootstrap
        results = results.sort_values(by='P').dropna()
        if o.flip_depleted_for_interaction:
            significant_single_factors = results.index[(results['Padj_tsBH']<=o.interaction_settings[0])]
            flippedI.update({f:(results.loc[f, 'logitdiff']<0) for f in significant_single_factors})
            for f, flip in flippedI.items():
                if flip:
                    genedata[f] = flip_values(genedata[f])
        else:
            significant_single_factors = results.index[(results['Padj_tsBH']<=o.interaction_settings[0]) & (results['logitdiff']>0)]
        
        #print(significant_single_factors, [(min(parse_contig_table_repr(results.loc[TF, 'contig_table'])) >= o.interaction_settings[1], genedata[TF].dropna().apply(lambda v: isinstance(v, bool)).all()) for TF in significant_single_factors])
        
        significant_single_factors = [TF for TF in significant_single_factors if min(parse_contig_table_repr(results.loc[TF, 'contig_table'])) >= o.interaction_settings[1] and genedata[TF].dropna().apply(lambda v: isinstance(v, bool)).all()]
        
        
        if o.interaction_settings[4] < len(significant_single_factors):
            significant_single_factors = significant_single_factors[:int(o.interaction_settings[4])]
        
        
        import itertools, numpy, functools
        interaction_outT = pandas.DataFrame()
        multiprocess = o.interaction_settings[5] > 1.01
        
        if multiprocess:
            from concurrent import futures
            pool = futures.ProcessPoolExecutor(int(o.interaction_settings[5]))
            newrow_futures = []
        
        for factor1, factor2 in tqdm.tqdm(list(itertools.combinations(significant_single_factors, 2)), desc='Interaction preparation' if multiprocess else 'Interaction calculation'):
            interactiondata = genedata[[factor1, factor2]].dropna().copy()
            interactiondata.columns = ['TF1', 'TF2']
            
            n_per_quarter = [((interactiondata['TF1']==v1) & (interactiondata['TF2']==v2)).sum() for v1 in (True, False) for v2 in (True,False)]
            
            if min(n_per_quarter) < o.interaction_settings[1]: continue
            
            interactiondata['geneset'] = False
            interactiondata.loc[interactiondata.index.intersection(genes), 'geneset'] = True
            
            if multiprocess:
                newrow_futures.append((pool.submit(calc_interaction, interactiondata, o, n_per_quarter), factor1, factor2))
            else:
                newrow = calc_interaction(interactiondata, o, n_per_quarter)
                if newrow is not None:
                    interaction_outT[name_factor(factor1) + '_x_' + name_factor(factor2)] = pandas.Series(newrow)
        
        if multiprocess:
            for fut, factor1, factor2 in tqdm.tqdm(newrow_futures, desc='Interaction calculation'):
                newrow = fut.result()
                if newrow is not None:
                    interaction_outT[name_factor(factor1) + '_x_' + name_factor(factor2)] = pandas.Series(newrow)
        
        interaction_out = interaction_outT.T
        
        if len(interaction_out) == 0:
            error("No interactions could be tested")
        
        
        add_padj(interaction_out, 'eq10', o)
        add_padj(interaction_out, 'eq01', o)
        add_padj(interaction_out, 'add', o)
        add_padj(interaction_out, 'mult', o)
        interaction_out[['interaction_sort_order', 'interaction_type']] = interaction_out.apply(functools.partial(calc_interaction_status, o), axis=1).apply(pandas.Series)
        
        if o.sortP:
            interaction_out = interaction_out.sort_values(by=['interaction_sort_order', 'sort_by_tmp'])
        del interaction_out['sort_by_tmp']
        
        interaction_out.to_csv(o.interaction_output, sep='\t')
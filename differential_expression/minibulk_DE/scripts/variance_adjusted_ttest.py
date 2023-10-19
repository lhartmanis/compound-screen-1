import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
from matplotlib.backends.backend_pdf import PdfPages

from sklearn.linear_model import LinearRegression
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy import stats
import statsmodels.api as sm

def train_regressors(exp_data, meta, outpath, rnatype = 'New', print_conds = False, kdeplot = False, kdecmap = 'Blues', replicate_column = 'compound_dose', max_data_points = 50_000):

    print(f"Training regressors. Using {replicate_column} as replicates and restricting to {max_data_points} data points")

    replicates = meta.loc[:, replicate_column].unique()
    can_use  = 0
    cant_use = 0

    namelist = []
    meanlist = []
    std_list = []
    cv2_list = []

    data_counter = 0

    for rep in replicates:
        while data_counter < max_data_points:
            current_meta = meta[meta.loc[:, replicate_column] == rep]
            comp_times   = current_meta.compound_time.unique()
            for comp_time in comp_times:
                lab_times = current_meta.query(f"compound_time == '{comp_time}'").label_time.unique()
                for lab_time in lab_times:
                    bcs = current_meta.query(f"compound_time == '{comp_time}' & label_time == '{lab_time}'").BC.values
                    # This check is needed if we're starting from an already subsetted dataframe
                    bcs = [bc for bc in bcs if bc in exp_data.columns.values]
                    if len(bcs) < 3:
                        cant_use += 1
                        continue
                    else:
                        can_use += 1

                    # Keep genes that have at least three not NA replicates
                    data = exp_data.loc[:, bcs][exp_data.loc[:, bcs].notna().sum(1) >= 3]

                    treatment_mean = data.mean(1)
                    treatment_std  = data.std(1)
                    treatment_cv2  = (treatment_std / treatment_mean) ** 2

                    genes = data.index.values
                    names = [f"{gene}-{comp_time}{rep}_{lab_time}4sU" for gene in genes]

                    namelist.append(names)
                    meanlist.append(treatment_mean.values)
                    std_list.append(treatment_std.values)
                    cv2_list.append(treatment_cv2.values)
                    data_counter = len([i for sublist in cv2_list for i in sublist])

    namelist = [i for sublist in namelist for i in sublist]
    meanlist = [i for sublist in meanlist for i in sublist]
    std_list = [i for sublist in std_list for i in sublist]
    cv2_list = [i for sublist in cv2_list for i in sublist]

    out_df = pd.DataFrame({'name': namelist,
                           'mean_val': meanlist,
                           'stdev': std_list,
                           'cv2': cv2_list})

    out_df.dropna(inplace = True)
    print(f"Shape of training data: {out_df.shape}")

    if print_conds:
        print(f"Conditions used to train regressor: {can_use}")
        print(f"Conditions not used: {cant_use}")

    # Removing zeros and outliers
    out_df = out_df[(out_df.cv2 != 0) & (out_df.mean_val != 0)]
    out_df.loc[:, 'log10_cv2_zscore'] = stats.zscore(np.log10(out_df.cv2))
    out_df = out_df[np.abs(out_df.log10_cv2_zscore) < 3]

    xvals = np.log10(out_df.mean_val).values.reshape(-1, 1)
    yvals = np.log10(out_df.cv2).values.reshape(-1, 1)

    linreg = LinearRegression().fit(xvals, yvals)
    regx = np.linspace(xvals.min(), xvals.max(), 1000)
    regy = linreg.predict(regx.reshape(-1, 1))

    lowess_res = lowess(yvals.ravel(), xvals.ravel(), return_sorted = False, xvals = regx)
    lowess_dict = {'y_values': yvals.ravel(),
                   'x_values': xvals.ravel()}

    fig, ax = plt.subplots(figsize = (10, 5))
    # Plot data values
    ax.scatter(xvals, yvals, s = 10, alpha = 0.3)

    # Plot densities
    if kdeplot:
        sns.kdeplot(x = xvals.ravel(),
                    y = yvals.ravel(),
                    cmap = kdecmap)

    # Plot linear regression fit
    ax.plot(regx, regy, color = 'k', linestyle = '-', linewidth = 2, alpha = 0.7,
            label = 'Linear regression')

    # Plot lowess regression fit
    ax.plot(regx, lowess_res, color = 'firebrick', linestyle = '-', linewidth = 2, alpha = 0.8,
            label = 'Loess regression')

    ax.tick_params(labelsize = 12)
    ax.set_ylabel("log$_{10}$ CV$^2$", fontsize = 12)
    ax.set_xlabel("log$_{10}$ mean rpkm", fontsize = 12)

    ax.legend(fontsize = 14)
    fig.suptitle(f"{rnatype} RNA", fontsize = 14)

    sns.despine()
    plt.savefig(os.path.join(outpath, 'regression_training.png'))

    return out_df, linreg, lowess_dict

def run_linear_predict(regressor, log10_mean_vals):
    linear_res = regressor.predict(log10_mean_vals.reshape(-1, 1))
    return linear_res.ravel()

def run_lowess_predict(lowess_dict, log10_mean_vals):
    lowess_res = lowess(lowess_dict['y_values'],
                        lowess_dict['x_values'],
                        return_sorted = False,
                        xvals = log10_mean_vals)
    return lowess_res

def plot_regressions_test_data(indata, log10_cv2_column, meanexp_column, ax, linear_predictor, lowess_dict, title = '', kdecmap = 'Blues_r', scatter_color = 'C0', linreg_color = 'k', loess_color = 'darkred'):

    plot_df = indata[(indata.loc[:, log10_cv2_column] != 0) & (indata.loc[:, meanexp_column] != 0)]
    plot_df.loc[:, 'log10_cv2_zscore'] = stats.zscore(plot_df.loc[:, log10_cv2_column])
    plot_df = plot_df[np.abs(plot_df.log10_cv2_zscore) < 3]

    ax.scatter(x = np.log10(plot_df.loc[:, meanexp_column]),
               y = plot_df.loc[:, log10_cv2_column], alpha = 0.3, s = 10, c = scatter_color)

    sns.kdeplot(x = np.log10(plot_df.loc[:, meanexp_column]),
                y = plot_df.loc[:, log10_cv2_column], ax = ax, cmap = kdecmap)

    pred_range = np.linspace(np.log10(plot_df.loc[:, meanexp_column]).min(),
                             np.log10(plot_df.loc[:, meanexp_column]).max())


    ax.plot(pred_range, run_linear_predict(linear_predictor, pred_range),
            label = 'Linear regression', linewidth = 2, c = linreg_color)
    ax.plot(pred_range, run_lowess_predict(lowess_dict, pred_range),
            label = 'Lowess regression', linewidth = 2, c = loess_color)

    ax.tick_params(labelsize = 12)
    ax.set_xlabel('log$_{10}$ mean RPKM', fontsize = 14)
    ax.set_ylabel('log$_{10}$ CV$^{2}$', fontsize = 14)
    ax.set_title(title, fontsize = 14)
    ax.legend(fontsize = 12)

# This could be two functions!
def add_predicted_cv2(indata, linear_predictor, lowess_dict):
    work_data = indata.copy()
    work_data.loc[:, 'low_linear_log10_cv2']  = run_linear_predict(linear_predictor, np.log10(work_data.low_mean.values))
    work_data.loc[:, 'high_linear_log10_cv2'] = run_linear_predict(linear_predictor, np.log10(work_data.high_mean.values))

    work_data.loc[:, 'low_lowess_log10_cv2']   = run_lowess_predict(lowess_dict, np.log10(work_data.low_mean.values))
    work_data.loc[:, 'high_lowess_log10_cv2']  = run_lowess_predict(lowess_dict, np.log10(work_data.high_mean.values))

    return work_data

# Make this into a function that gets called with indata and names of variance and mean columns
def compare_and_append_stdevs(indata):
    work_data = indata.copy()

    # Calculate standard deviation from CV² for low and high dose treatments
    low_stdev_linear_list  = np.sqrt(10 ** (work_data.low_linear_log10_cv2)  * work_data.loc[:, 'low_mean']**2)
    high_stdev_linear_list = np.sqrt(10 ** (work_data.high_linear_log10_cv2) * work_data.loc[:, 'high_mean']**2)

    low_stdev_lowess_list  = np.sqrt(10 ** (work_data.low_lowess_log10_cv2)  * work_data.loc[:, 'low_mean']**2)
    high_stdev_lowess_list = np.sqrt(10 ** (work_data.high_lowess_log10_cv2) * work_data.loc[:, 'high_mean']**2)

    # Instantiate new columns with standard deviations computed from data
    work_data.loc[:, 'low_stdev_to_use_linear']  = work_data.loc[:, 'low_stdev']
    work_data.loc[:, 'high_stdev_to_use_linear'] = work_data.loc[:, 'high_stdev']

    work_data.loc[:, 'low_stdev_to_use_lowess']  = work_data.loc[:, 'low_stdev']
    work_data.loc[:, 'high_stdev_to_use_lowess'] = work_data.loc[:, 'high_stdev']

    # Compare observed and predicted variances. Save indeces where predicted is higher than observed
    change_idx_low_linear  = np.where(work_data.low_linear_log10_cv2  > work_data.low_log10_cv2)[0]
    change_idx_high_linear = np.where(work_data.high_linear_log10_cv2 > work_data.high_log10_cv2)[0]

    change_idx_low_lowess  = np.where(work_data.low_lowess_log10_cv2  > work_data.low_log10_cv2)[0]
    change_idx_high_lowess = np.where(work_data.high_lowess_log10_cv2 > work_data.high_log10_cv2)[0]

    # Change to predicted CV² where the predicted variance is larger than observed
    work_data.loc[:, 'low_stdev_to_use_linear'].iloc[change_idx_low_linear] = low_stdev_linear_list[change_idx_low_linear]
    work_data.loc[:, 'high_stdev_to_use_linear'].iloc[change_idx_high_linear] = high_stdev_linear_list[change_idx_high_linear]

    work_data.loc[:, 'low_stdev_to_use_lowess'].iloc[change_idx_low_lowess] = low_stdev_lowess_list[change_idx_low_lowess]
    work_data.loc[:, 'high_stdev_to_use_lowess'].iloc[change_idx_high_lowess] = high_stdev_lowess_list[change_idx_high_lowess]


    return work_data

def make_res_dict(indata, meta, comparison_dict, linear_predictor, lowess_dict, subsample_g1, subsample_g2 ,target_column = 'compound_name', pcount = 0.001, label_time = '1hr', compound_time = '1hr'):

    name_dict = {}
    for key in comparison_dict.keys():
        g1names = '-'.join([str(i) for i in comparison_dict[key]['group1']])
        g2names = '-'.join([str(i) for i in comparison_dict[key]['group2']])
        comparison = f"{g1names}_vs_{g2names}"
        name_dict[key] = comparison

    analysis_meta = meta.query(f'label_time == "{label_time}" & compound_time == "{compound_time}"')
    res_dict = {}

    for idx, comp in enumerate(list(comparison_dict.keys())):
        g1 = analysis_meta[analysis_meta.loc[:, target_column].isin(comparison_dict[comp]['group1'])].BC.values
        g2 = analysis_meta[analysis_meta.loc[:, target_column].isin(comparison_dict[comp]['group2'])].BC.values

        if subsample_g1:
            g1 = np.random.choice(g1, replace = False, size = subsample_g1)
        if subsample_g2:
            g2 = np.random.choice(g2, replace = False, size = subsample_g2)

        if (len(g1) < 1) or (len(g2) < 1):
            print(comp)
            print(f"Can't use comparison {comparison_dict[comp]['group1']} vs {comparison_dict[comp]['group2']}, too few replicates")
            continue

        name_list = []

        low_group_list  = []
        high_group_list = []

        low_size_list   = []
        high_size_list  = []

        low_mean_list   = []
        high_mean_list  = []

        low_stdev_list  = []
        high_stdev_list = []

        low_sum_list    = []
        high_sum_list   = []

        low_log10_cv2_list  = []
        high_log10_cv2_list = []

        low_bcs  = []
        high_bcs = []

        for gname, gdata in indata.iterrows():
            low_data  = gdata.loc[g1].dropna()
            high_data = gdata.loc[g2].dropna()
            if (len(low_data) == 0) or (len(high_data) == 0):
                continue
            else:
                name_list.append(gname)

                low_bcs.append(g1)
                high_bcs.append(g2)

                low_group_list.append(low_data.values)
                high_group_list.append(high_data.values)

                low_size_list.append(len(low_data))
                high_size_list.append(len(high_data))

                low_mean_list.append(low_data.mean())
                high_mean_list.append(high_data.mean())

                low_stdev_list.append(low_data.std())
                high_stdev_list.append(high_data.std())

                low_log10_cv2  = np.log10(((low_data + pcount).std()  / (low_data + pcount).mean()) ** 2)
                high_log10_cv2 = np.log10(((high_data + pcount).std() / (high_data + pcount).mean()) ** 2)

                low_log10_cv2_list.append(low_log10_cv2)
                high_log10_cv2_list.append(high_log10_cv2)

                low_sum_list.append(low_data.sum())
                high_sum_list.append(high_data.sum())

        rd = {'gene_id':name_list,
              'low_vals': low_group_list, 'high_vals': high_group_list,
              'low_bcs': low_bcs, 'high_bcs': high_bcs,
              'low_size': low_size_list, 'high_size': high_size_list,
              'low_mean': low_mean_list, 'high_mean': high_mean_list,
              'low_stdev': low_stdev_list, 'high_stdev': high_stdev_list,
              'low_log10_cv2': low_log10_cv2_list, 'high_log10_cv2': high_log10_cv2_list,
              'low_sum': low_sum_list, 'high_sum': high_sum_list}

        results = pd.DataFrame(rd)

        results = results[(results.low_size > 1) & (results.high_size > 1)]
        results = results[(results.low_sum > 0)  & (results.high_sum > 0)]

        results = pd.concat([results.iloc[:, :5],
                             results.iloc[:, 5:].replace({np.inf: np.nan, -np.inf: np.nan})], axis = 1)

        results.dropna(inplace = True)
        results.reset_index(drop = True, inplace = True)

        if results.shape[0] == 0:
            continue

        print(f'Predicting dispersion based on mean RPKMs for comparison {idx+1}/{len(list(comparison_dict.keys()))}...')
        results = add_predicted_cv2(results, linear_predictor, lowess_dict)
        results = compare_and_append_stdevs(results)

        res_dict[comp] = results

    return res_dict, name_dict

def calculate_t_statistic(group1, group2, equal_variance, stdev_group1 = None, stdev_group2 = None):
    g1s = group1.size
    g2s = group2.size

    g1mean = group1.mean()
    g2mean = group2.mean()

    if not stdev_group1:
        print("Calculating standard deviation for group1 from data...")
        g1var = group1.var()
    else:
        g1var = stdev_group1 ** 2

    if not stdev_group2:
        print("Calculating standard_deviation for group2 from data...")
        g2var = group2.var()
    else:
        g2var = stdev_group2 ** 2

    if equal_variance:
        pooled_stdev = np.sqrt(((g1s - 1) * g1var + (g2s - 1) * g2var) / (g1s + g2s - 2))
        t = (g1mean - g2mean) / (pooled_stdev * (np.sqrt(g1s**-1 + g2s**-1)))

    else:
        denom = np.sqrt((g1var / g1s) + (g2var / g2s))
        t  = (g1mean - g2mean) / denom

    return t

def equal_var_df(group1, group2):
    df = group1.size + group2.size - 2
    return df

def welch_satterthwaite_df(group1, group2):
    df = (group1.var() / group1.size + group2.var() / group2.size)**2 / ((group1.var() / group1.size)**2 / (group1.size - 1) + (group2.var() / group2.size)**2 / (group2.size - 1))
    return df

def perform_ttest(group1, group2, std_group1, std_group2, equal_variance = 'infer'):
    if equal_variance == 'infer':
        var_ratio = std_group1**2 / std_group2**2
        if (var_ratio < 2) & (var_ratio > 0.5):
            equal_variance = True
        else:
            equal_variance = False

    t = calculate_t_statistic(group1, group2, equal_variance, std_group1, std_group2)

    if equal_variance:
        df = equal_var_df(group1, group2)

    else:
        df = welch_satterthwaite_df(group1, group2)

    p_value = stats.t.sf(np.abs(t), df = df) * 2
    return t, p_value

def run_adjusted_ttest(indata, var_correction_type = 'both', equal_variance = 'infer', plot = False, returnfig = True, plot_title = '', fdr_method = "fdr_bh"):
    if var_correction_type not in ['both', 'linear', 'lowess']:
        raise ValueError("Correction type needs to be one of 'both', 'linear' or 'lowess'")

    print("Calculating p-values...")

    linear_res = indata.apply(lambda x: perform_ttest(x.low_vals, x.high_vals,
                                                      x.low_stdev_to_use_linear,
                                                      x.high_stdev_to_use_linear,
                                                      equal_variance = equal_variance),
                                                      axis = 1)

    linear_df = pd.DataFrame({'tstat': [i[0] for i in linear_res],
                              'pval':  [i[1] for i in linear_res]})

    lowess_res = indata.apply(lambda x: perform_ttest(x.low_vals, x.high_vals,
                                                      x.low_stdev_to_use_lowess,
                                                      x.high_stdev_to_use_lowess,
                                                      equal_variance = equal_variance),
                                                      axis = 1)

    lowess_df = pd.DataFrame({'tstat': [i[0] for i in lowess_res],
                              'pval':  [i[1] for i in lowess_res]})

    print("Performing multiple test correction...")
    adjusted_pvals_linear = sm.stats.multipletests(linear_df.loc[:,'pval'].values,
                                                   is_sorted = False,
                                                   method = fdr_method)[1]

    adjusted_pvals_lowess = sm.stats.multipletests(lowess_df.loc[:,'pval'].values,
                                                   is_sorted = False,
                                                   method = fdr_method)[1]

    return_dict = {'linear_corrected': {'tstat': linear_df.loc[:, 'tstat'].values,
                                        'unadjusted_pvals': linear_df.loc[:, 'pval'].values,
                                        'adjusted_pvals': adjusted_pvals_linear},

                   'lowess_corrected': {'tstat': lowess_df.loc[:, 'tstat'].values,
                                        'unadjusted_pvals': lowess_df.loc[:, 'pval'].values,
                                        'adjusted_pvals': adjusted_pvals_lowess}}
    if plot:
        print("Plotting...")
        fig, axes = plt.subplots(2, 2, figsize = (10, 10), sharex = 'all')
        ax1, ax2 = axes[0]
        ax3, ax4 = axes[1]

        kwargs_unadj = {'bins': 50, 'alpha' : 0.8, 'ec' : 'k', 'color': 'coral'}
        kwargs_adj   = {'bins': 50, 'alpha' : 0.8, 'ec' : 'k', 'color': 'mediumpurple'}

        ax1.hist(return_dict['linear_corrected']['unadjusted_pvals'], **kwargs_unadj)
        ax2.hist(return_dict['linear_corrected']['adjusted_pvals'], **kwargs_adj)
        ax3.hist(return_dict['lowess_corrected']['unadjusted_pvals'], **kwargs_unadj)
        ax4.hist(return_dict['lowess_corrected']['adjusted_pvals'], **kwargs_adj)
        fig.suptitle(plot_title, fontsize = 15)

        for idx, ax in enumerate(axes.ravel()):
            ax.tick_params(labelsize = 12)
            if idx == 0:
                ax.set_title('Unadjusted p-values', fontsize = 15)
                ax.set_ylabel('Stdev lin. reg. corrected\nOccurrence', fontsize = 14)
            if idx == 1:
                ax.set_title('Adjusted p-values', fontsize = 15)
            if idx == 2:
                ax.set_ylabel('Stdev lowess corrected\nOccurence', fontsize = 14)
                ax.set_xlabel('P-value', fontsize = 14)
            if idx == 3:
                ax.set_xlabel("P-value", fontsize = 14)

        sns.despine()

    print("Done performing DE testing!")
    
    if returnfig:
        return return_dict, fig
    else:
        return return_dict

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--indata_path', required = True, help = 'Path to expression data')
    parser.add_argument('-m', '--meta_path', required = True, help = 'Path to meta data file')
    parser.add_argument('-o', '--outfolder', required = True, help = 'Path to outfolder')
    parser.add_argument('-r', '--rna_type', required = True, help = '[New/ Total/ Old]')
    parser.add_argument('-e', '--equal_variance', default = 'infer', help = 'Assume equal variances in t-test?\nDefault = infer')
    parser.add_argument('-l', '--label_time', default = '1hr', help = 'Label time to analyze.\nDefault = 1hr')
    parser.add_argument('-c', '--compound_time', default = '1hr', help = 'Compound treatment time to analyze.\nDefault = 1hr')
    parser.add_argument('-p', '--pseudocount', type = float, default = 0.001, help = 'What pseudocount to use.\nDefault = 0.001')
    parser.add_argument('-v', '--var_correction_type', default = 'both', help = "How to correct variance? [both/ linear/ lowess]. Lowess correction is slow with many data points owing to the need to locally fit the regression function.\nDefault = both")
    parser.add_argument('-t', '--test_dict_path', help = "Path to dictionary specifying the tests to run.")
    parser.add_argument('--max_data_points', type = int, default = 50000, help = "How many datapoints to use when training the regressors? More points will make the lowess regression slow.\nDefault = 50000")
    parser.add_argument('--replicate_column', type = str, help = "What column in the meta file will be used to specify replicates for regression training and DE testing?")
    parser.add_argument('--subsample_g1', type = int, help = "Subsample this many conditions from group 1", default = None)
    parser.add_argument('--subsample_g2', type = int, help = "Subsample this many conditions from group 2", default = None)
    parser.add_argument('--conds_to_keep', type = str, help = "Path to file containing barcodes to keep for further analysis")

    args = parser.parse_args()

    indata_path = args.indata_path
    meta_path = args.meta_path
    outpath = args.outfolder
    rna_type = args.rna_type
    equal_variance = args.equal_variance
    label_time = args.label_time
    compound_time = args.compound_time
    pseudocount = args.pseudocount
    var_correction_type = args.var_correction_type
    max_data_points = args.max_data_points
    replicate_column = args.replicate_column
    test_dict_path = args.test_dict_path
    subsample_g1 = args.subsample_g1
    subsample_g2 = args.subsample_g2
    conds_to_keep_path = args.conds_to_keep

    if not os.path.exists(outpath):
        os.mkdir(outpath)

    with open(test_dict_path, 'r') as fh:
        comparison_dict = json.load(fh)

    print(f"Differential expression tests will be run on the following comparisons:\n{comparison_dict.keys()}")

    print(f"Files will be saved to {outpath}\n")
    print("Reading data...\n")
    exp_data = pd.read_csv(indata_path, sep = '\t', index_col = 0)
    meta = pd.read_csv(meta_path, sep = '\t')

    if conds_to_keep_path:
        with open(conds_to_keep_path, 'r') as fh:
            conds_to_keep = fh.read().splitlines()
        meta = meta[meta.BC.isin(conds_to_keep)]

    common_conds = np.intersect1d(exp_data.columns.values, meta.BC.values)

    exp_data = exp_data.loc[:, common_conds]
    meta = meta[meta.BC.isin(common_conds)]

    print("Fitting mean expression / variance regressors...\n")
    out_df, linear_predictor, lowess_dict = train_regressors(exp_data, meta, outpath, rnatype = rna_type, print_conds = False, kdeplot = True, kdecmap = 'Blues_r', replicate_column = replicate_column, max_data_points = max_data_points)

    print("Estimating variance based on mean expression...\n")
    res_dict, name_dict = make_res_dict(exp_data, meta, comparison_dict, linear_predictor, lowess_dict, subsample_g1, subsample_g2, target_column = replicate_column, pcount = pseudocount, label_time = label_time, compound_time = compound_time)

    # Save regression plots to file
    # This becomes quite slow, A better alternative is to save all plots as png:s
    # and then knit those together into multiple PDF pages
    print("Saving regression plots...")
    with PdfPages(os.path.join(outpath, 'regression_results.pdf')) as pdf:
        for idx, key in enumerate(res_dict.keys()):
            print(f"Saving regression plot {idx + 1}/{len(res_dict.keys())}")
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (15, 5))
            plot_regressions_test_data(res_dict[key], 'low_log10_cv2', 'low_mean', ax1, linear_predictor, lowess_dict, title = f'{rna_type} RNA\n{key} low treatment group')
            plot_regressions_test_data(res_dict[key], 'high_log10_cv2', 'high_mean', ax2, linear_predictor, lowess_dict, title = f'{rna_type} RNA\n{key} high treatment group')
            sns.despine()
            pdf.savefig(figure = fig)

    # Run t-tests
    print("Running tests...\n")

    all_breaks = {}
    unadjusted_pvals_dict = {}
    adjusted_pvals_dict = {}

    with PdfPages(os.path.join(outpath, 'pvalue_distributions.pdf')) as pdf:
        for breakpoint in res_dict.keys():
            temp_dict, fig = run_adjusted_ttest(res_dict[breakpoint], var_correction_type = var_correction_type,
                                                plot = True, plot_title = breakpoint, returnfig = True)

            all_breaks[breakpoint] = temp_dict
            pdf.savefig(figure = fig)

    for idx, breakpoint in enumerate(all_breaks.keys()):
        data_part = res_dict[breakpoint].loc[:, ['gene_id', 'low_vals', 'high_vals', 'low_bcs', 'high_bcs']]
        linear_part = pd.DataFrame(all_breaks[breakpoint]['linear_corrected']).drop('tstat', axis = 1).rename({'unadjusted_pvals': 'unadjusted_pvals_linear_corr', 'adjusted_pvals': 'adjusted_pvals_linear_corr'}, axis = 1)
        lowess_part = pd.DataFrame(all_breaks[breakpoint]['lowess_corrected']).drop('tstat', axis = 1).rename({'unadjusted_pvals': 'unadjusted_pvals_lowess_corr', 'adjusted_pvals': 'adjusted_pvals_lowess_corr'}, axis = 1)

        results_df = pd.concat([data_part, linear_part, lowess_part], axis = 1)

        results_df.loc[:, 'mean_ratio'] = results_df.apply(lambda x: x.high_vals.mean() / x.low_vals.mean(), axis = 1)
        results_df.loc[:, 'mean_diff']  = results_df.apply(lambda x: x.high_vals.mean() - x.low_vals.mean(), axis = 1)
        results_df.loc[:, 'diff_over_untreat']  = results_df.apply(lambda x: (x.high_vals.mean() - x.low_vals.mean()) / x.low_vals.mean(), axis = 1)

        print(f"Saving results {idx+1}/{len(all_breaks.keys())}, {breakpoint}")
        results_df.to_csv(os.path.join(outpath, f"{breakpoint}_DE_results.txt"), sep = '\t')

    # Two last plots
    fig, axes = plt.subplots(2, 2, figsize = (15, 15))
    ax1, ax2 = axes[0]
    ax3, ax4 = axes[1]

    kwargs = {'bins':50, 'alpha':0.7, 'ec':'k'}

    for breakpoint in list(all_breaks.keys())[:]:
        try:
            sns.kdeplot(all_breaks[breakpoint]['linear_corrected']['unadjusted_pvals'], label = breakpoint, ax = ax1, shade = True)
        except np.linalg.LinAlgError:
            print("No significant genes found with unadjusted p-values, leaving plot area empty...")
        
        try:
            sns.kdeplot(all_breaks[breakpoint]['linear_corrected']['adjusted_pvals'],   label = breakpoint, ax = ax2, shade = True)
        except np.linalg.LinAlgError:
            print("No significant genes found with adjusted p-values, leaving plot area empty...")
        
        try:
            sns.kdeplot(all_breaks[breakpoint]['lowess_corrected']['unadjusted_pvals'], label = breakpoint, ax = ax3, shade = True)
        except np.linalg.LinAlgError:
            print("No significant genes found with unadjusted p-values, leaving plot area empty...")

        try:
            sns.kdeplot(all_breaks[breakpoint]['lowess_corrected']['adjusted_pvals'],   label = breakpoint, ax = ax4, shade = True)
        except np.linalg.LinAlgError:
            print("No significant genes found with adjusted p-values, leaving plot area empty...")

    [ax.tick_params(labelsize = 12) for ax in axes.ravel()]

    for idx, ax in enumerate(axes.ravel()):
        if idx == 0:
            ax.set_title('Unadjusted p-values', fontsize = 15)
            ax.set_ylabel('Stdev lin. reg. corrected\nDensity', fontsize = 14)
        if idx == 1:
            ax.set_title('Adjusted p-values', fontsize = 15)
            ax.set_ylabel("")
        if idx == 2:
            ax.set_ylabel('Stdev lowess corrected\nDensity', fontsize = 14)
            ax.set_xlabel('P-value', fontsize = 14)
        if idx == 3:
            ax.set_xlabel("P-value", fontsize = 14)
            ax.set_ylabel("")

    sns.despine()
    ax4.legend(loc = 'upper left', fontsize = 12)
    plt.savefig(os.path.join(outpath, 'pvalue_density.pdf'))

    fig, axes = plt.subplots(2, 2, figsize = (15, 15))
    ax1, ax2 = axes[0]
    ax3, ax4 = axes[1]

    kwargs = {'bins':50, 'alpha':0.5, 'ec':'k'}

    for breakpoint in list(all_breaks.keys()):
        ax1.hist(all_breaks[breakpoint]['linear_corrected']['unadjusted_pvals'], label = breakpoint, **kwargs)
        ax2.hist(all_breaks[breakpoint]['linear_corrected']['adjusted_pvals'],   label = breakpoint, **kwargs)
        ax3.hist(all_breaks[breakpoint]['lowess_corrected']['unadjusted_pvals'], label = breakpoint, **kwargs)
        ax4.hist(all_breaks[breakpoint]['lowess_corrected']['adjusted_pvals'],   label = breakpoint, **kwargs)


    [ax.tick_params(labelsize = 12) for ax in axes.ravel()]

    for idx, ax in enumerate(axes.ravel()):
        ax.tick_params(labelsize = 12)
        if idx == 0:
            ax.set_title('Unadjusted p-values', fontsize = 15)
            ax.set_ylabel('Stdev lin. reg. corrected\nOccurrence', fontsize = 14)
        if idx == 1:
            ax.set_title('Adjusted p-values', fontsize = 15)
            ax.set_ylabel("")
        if idx == 2:
            ax.set_ylabel('Stdev lowess corrected\nOccurrence', fontsize = 14)
            ax.set_xlabel('P-value', fontsize = 14)
        if idx == 3:
            ax.set_xlabel("P-value", fontsize = 14)
            ax.set_ylabel("")

    sns.despine()
    ax4.legend(loc = 'upper left', fontsize = 12)
    plt.savefig(os.path.join(outpath, 'pvalue_hist.pdf'))

    print("\nAll done!")

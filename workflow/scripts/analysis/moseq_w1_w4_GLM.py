import os, sys
import h5py
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.formula.api import glm
from sklearn.model_selection import train_test_split

# import util functions from utils module
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(os.getcwd())
sys.path.append(parent_dir)

from utils.population import activity_at_phase


def get_GLM_and_prediction(syl_ratio_mx, pop_at_phase, test_size=0.33, glm_min_pval=0.99):

    # separate train / test
    X_train, X_test, y_train, y_test = train_test_split(syl_ratio_mx, pop_at_phase, test_size=test_size)
    
    # train glm to get contributions of each syllable
    data = np.column_stack([y_train, X_train])
    columns = ['state'] + ["x%d" % x for x in range(X_train.shape[1])]
    AM_df = pd.DataFrame(data, columns=columns)

    model = glm('state ~ ' + ' + '.join(columns[1:]), data=AM_df).fit()
    #model.summary()
    
    glm_coeffs = dict([(i, coef) for i, coef in enumerate(model.params[1:]) if model.pvalues[1:][i] < glm_min_pval])
    #if len(glm_coeffs) == 0:
    #    return 0, 0, model.params, model.pvalues
    target_fit = np.zeros(len(y_test))
    for idx, coef in glm_coeffs.items():
        target_fit += coef * X_test[:, idx]

    corr, pval = stats.pearsonr(target_fit, y_test)
    
    return corr, pval, model.params, model.pvalues


iter_count = 100  # how much to shuffle and train / test split

with h5py.File(snakemake.input[0], 'r') as f:
    tl = np.array(f['processed']['timeline'])
    tgt_mx = np.array(f['processed']['target_matrix'])
    events = np.array(f['processed']['sound_events'])
with h5py.File(snakemake.input[1], 'r') as f:
    syl_ratio_mx = np.array(f['syl_ratio_mx'])
    idxs_srm_tl  = np.array(f['idxs_srm_tl'])
    
for j, phase in enumerate([1, 2, 3, 4]):
    w_pca = activity_at_phase(os.path.dirname(snakemake.input[0]), phase, do_pca=True)
    w_int = np.interp(tl[idxs_srm_tl][:, 0], events[:, 0], w_pca)
    
    # original
    corr, pval, params, pvalues = get_GLM_and_prediction(syl_ratio_mx, w_int)
    
    # diff train / split
    corr_mx_chun = np.zeros([iter_count, 2])
    for k in range(iter_count):
        corr_s, pval_s, _, _ = get_GLM_and_prediction(syl_ratio_mx, w_int)
        corr_mx_chun[k] = (corr_s, pval_s)
        
    # shuffled
    corr_mx_shuf = np.zeros([iter_count, 2])  # coeff, pval for each shuffle
    for k in range(iter_count):
        syl_ratio_mx_s = syl_ratio_mx.copy()
        np.random.shuffle(syl_ratio_mx_s)
        corr_s, pval_s, _, _ = get_GLM_and_prediction(syl_ratio_mx_s, w_int)
        corr_mx_shuf[k] = (corr_s, pval_s)
    
    grp_name = 'W' + str(j+1)
    with h5py.File(snakemake.output[0], 'a') as f:
        if not grp_name in f:
            f.create_group(grp_name)
        tgt_grp = f[grp_name]

        for val in ['corr_glm_fit_orig', 'glm_fit_params', 'glm_fit_pvalues', 'corr_glm_fit_shuf', 'corr_glm_fit_chun']:
            if val in tgt_grp:
                del tgt_grp[val]

        tgt_grp.create_dataset('corr_glm_fit_orig', data=np.array([corr, pval]))
        tgt_grp.create_dataset('glm_fit_params',    data=params)
        tgt_grp.create_dataset('glm_fit_pvalues',   data=pvalues)
        tgt_grp.create_dataset('corr_glm_fit_shuf', data=corr_mx_shuf)
        tgt_grp.create_dataset('corr_glm_fit_chun', data=corr_mx_chun)
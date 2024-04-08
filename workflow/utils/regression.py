import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from statsmodels.formula.api import glm
from scipy import stats


def get_GLM_and_prediction(X, y, test_size=0.33, glm_min_pval=0.95):

    # separate train / test
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size)
    
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
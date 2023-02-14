##########
# Matthew Muller
# 11/24/2022
#
# Stats Functions to help 
##########

##########
# Library Imports
import sys
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import scipy.stats as st

from sklearn.metrics import roc_auc_score, roc_curve, auc

##########
# Set/Append Working directory
sys.path.append('.')

##########
# Import Functions

##########
# Code Below

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), st.sem(a)
    h = se * st.t.ppf((1 + confidence) / 2., n-1)
    low, high = st.t.interval(confidence, len(data)-1, loc=np.mean(data), scale=st.sem(data))
    return m, low, high

def bootstrap_auc_confidence(y_pred, y_true, ci =  0.95,
                  n_bootstraps = 1000, rng_seed = 100,
                  plot_histogram = False):
    '''
    Summary: Binary target implementation of auc bootstrapping for determining a confidence interval

    y_pred (np.array) : numpy array of predicted values, usually given from the predict_proba method.
    y_true (np.array) : numpy array of true label values (1 is the presumed target)
    n_bootstraps (int) : integer value of the number of bootstraps
    rng_seed (str) : random seed to set a random state, which 
    plot_histogram (bool) : Plot a histogram of AUC values (default = False)

    output (list) : working_roc_auc, confidence_lower, confidence_upper
    '''

    working_roc_auc = roc_auc_score(y_true, y_pred)
    bootstrapped_scores = []
    rng = np.random.RandomState(rng_seed)
    for i in range(n_bootstraps):
        # bootstrap by sampling with replacement on the prediction indices
        indices = rng.randint(0, len(y_pred), len(y_pred))
        if len(np.unique(y_true[indices])) < 2:
            # We need at least one positive and one negative sample for ROC AUC
            # to be defined: reject the sample
            continue

        score = roc_auc_score(y_true[indices], y_pred[indices])
        bootstrapped_scores.append(score)
    
    if plot_histogram:
        import matplotlib.pyplot as plt
        plt.hist(bootstrapped_scores, bins=50)
        plt.title('Histogram of the bootstrapped ROC AUC scores')
        plt.show()

    sorted_scores = np.array(bootstrapped_scores)
    sorted_scores.sort()

    working_roc_auc, confidence_lower, confidence_upper = mean_confidence_interval(sorted_scores, confidence=ci)

    print(f"Confidence interval for the score: {working_roc_auc:0.3f} [{confidence_lower:0.3f} - {confidence_upper:0.3f}]")
    return working_roc_auc, confidence_lower, confidence_upper
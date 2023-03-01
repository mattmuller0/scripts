##########
# Matthew Muller
# 11/24/2022
#
# Plotting Functions to help 
##########

##########
# Library Imports
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import sys

from sklearn.decomposition import PCA, KernelPCA, NMF
from sklearn.metrics import RocCurveDisplay, confusion_matrix, ConfusionMatrixDisplay
from sklearn.metrics import roc_auc_score, roc_curve, auc
from sklearn.model_selection import StratifiedKFold, KFold

from sklearn.base import clone
import scipy.stats as st
import statsmodels.stats.api as sms
from scipy.stats import kruskal, zscore

##########
# Set/Append Working directory
sys.path.append('/Users/muller/Documents/RugglesLab')


##########
# Import Functions
from Pytools.stats import mean_confidence_interval

##########
# Code Below

def plot_scree(
    pca : PCA(), components = 50, 
    figsize = (8,6), fontsize = 20, 
    save_path = None
    ):
    '''
    pca (sklearn.decomposition.pca()) : sklearn pca instance fitted to data
    components (int) : number of components to plot on y-axis
    figsize (tuple) : size of graph
    fontsize (int) : title fontsize
    '''
    plt.figure(figsize=figsize)
    plt.plot(np.arange(pca.n_components_)[:components], 
             pca.explained_variance_ratio_[:components],
             'o-')
    plt.title('Scree Plot', fontsize = fontsize)
    plt.xlabel('Principal Component')
    plt.ylabel('Proportion of Variance Explained')
    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()

def plot_pca(
    data,  x = 0, y = 1, labels = None, 
    figsize=(8,6), 
    save_path=None
    ):
    '''
    data (np.matrix) : numpy matrix of sklearn pca reduced data
    x (int) : integer value of first PC (default 0)
    y (int) : integer value of second PC (default 1)
    save_path (str) : string pointing where to save image
    figsize (tuple) : size of graph
    '''
    plt.figure(figsize=figsize)
    sns.scatterplot(pd.DataFrame(data), x=x, y=y, hue=labels)
    plt.title(f'PCA: PC {x} versus PC {y}', fontsize = 20)
    plt.xlabel(f'PC {x}')
    plt.ylabel(f'PC {y}')
    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()

def plot_nmf(
    data,  x = 0, y = 1, labels = None, 
    figsize = (8,6),
    save_path = None
    ):
    '''
    data (np.matrix) : numpy matrix of sklearn nmf reduced data
    x (int) : integer value of first column (default 0)
    y (int) : integer value of second column (default 1)
    save_path (str) : string pointing where to save image
    figsize (tuple) : size of graph
    '''
    plt.figure(figsize=figsize)
    sns.scatterplot(pd.DataFrame(data), x=x, y=y, hue=labels)
    plt.title(f'NMF:  Column {x} versus Column {y}', fontsize = 20)
    plt.xlabel(f'Column {x}')
    plt.ylabel(f'Column {y}')
    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()

def plot_tsne(
    data,  x = 0, y = 1, labels = None,
    figsize = (8,6),
    save_path = None
    ):
    '''
    data (np.matrix) : numpy matrix of sklearn tsne reduced data
    x (int) : integer value of first column (default 0)
    y (int) : integer value of second column (default 1)
    save_path (str) : string pointing where to save image
    figsize (tuple) : size of graph
    '''
    plt.figure(figsize=figsize)
    sns.scatterplot(pd.DataFrame(data), x=x, y=y, hue=labels)
    plt.title(f't-SNE {x} versus t-SNE {y}', fontsize = 20)
    plt.xlabel(f't-SNE {x}')
    plt.ylabel(f't-SNE {y}')
    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()

def plot_roc_curve(
    y_true, score, 
    save_path = None
    ):
    '''
    y_true (np.array) : numpy array of the testing labels
    score (np.array) : numpy array of the prediction values (using model.predict_proba)
    save_path (str) : string pointing where to save image
    '''
    fpr, tpr, _ = roc_curve(y_true, score)
    roc_auc = auc(fpr, tpr)
    roc_display = RocCurveDisplay(fpr=fpr, tpr=tpr,
                                  roc_auc=roc_auc).plot()
    plt.title('ROC Curve', fontsize = 14)
    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()

def plot_confusion_matrix(
    y_true, score, 
    save_path=None, labels=None):
    '''
    y_true (np.array) : numpy array of the testing labels
    score (np.array) : predictions of model (using model.predict)
    save_path (str) : string pointing where to save image
    labels (list) : list of labels used for classes
    '''
    cm = confusion_matrix(y_true, score)
    cm_display = ConfusionMatrixDisplay(cm, display_labels=labels).plot()
    plt.title('Confusion Matrix', fontsize = 14)
    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()

def plot_training_roc_curve_ci(
    model,
    X, y,
    cv_splits = 6,
    fill = True,
    title = None,
    save_path = None,
    ):
    '''
    model (joblib) : model to be copied for determining training confidence
    X (np.array or pd.DataFrame) : np array of features used
    y (np.array) : list of labels used for classes
    save_path (str) : string pointing where to save image
    '''
    if not isinstance(X, np.ndarray):
        X = X.to_numpy()

    cv = StratifiedKFold(n_splits=cv_splits)
    model = clone(model)

    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)

    fig, ax = plt.subplots(figsize=(6, 6))
    for fold, (train, test) in enumerate(cv.split(X, y)):
        model.fit(X[train], y[train])
        viz = RocCurveDisplay.from_estimator(
            model,
            X[test],
            y[test],
            name=f"ROC fold {fold}",
            alpha=0.3,
            lw=1,
            ax=ax,
        )
        plt.cla() # This removes each individual interation
        interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
        interp_tpr[0] = 0.0
        tprs.append(interp_tpr)
        aucs.append(viz.roc_auc)
    ax.plot([0, 1], [0, 1], "k--", label="chance level (AUC = 0.5)")
    
    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = mean_confidence_interval(aucs)[0] - mean_confidence_interval(aucs)[1]
    ax.plot(
        mean_fpr,
        mean_tpr,
        color="b",
        label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
        lw=2,
        alpha=0.8,
    )

    # This is some array magic to make mean_confidence_interval play well along a numpy axis=0 selection
    tprs = np.array(tprs)
    ci_tpr = [mean_confidence_interval(tprs[:,idx].flatten()) for idx in range(tprs.shape[1])]
    ci_tpr = np.array(ci_tpr)[:,0] - np.array(ci_tpr)[:,1] # mean - (mean + ci)

    tprs_upper = np.minimum(mean_tpr + ci_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - ci_tpr, 0)
    if fill:
        ax.fill_between(
            mean_fpr,
            tprs_lower,
            tprs_upper,
            color="grey",
            alpha=0.2,
            label="95% Confidence Interval",
        )
    if not title:
        title="Mean ROC curve with 95% Confidence Interval"
    ax.set(
        xlim=[0, 1],
        ylim=[0, 1],
        xlabel="False Positive Rate",
        ylabel="True Positive Rate",
        title=title,
    )
    ax.axis("square")
    ax.grid(False)
    ax.legend(loc="lower right")
    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()

def plot_roc_curve_ci(
    model,
    X, y,
    cv_splits = 5,
    stratified = True,
    title = None,
    random_state = 100,
    save_path = None,
    ):
    '''
    model (joblib) : model to be copied for determining training confidence
    X (np.array or pd.DataFrame) : np array of features used
    y (np.array) : list of labels used for classes
    save_path (str) : string pointing where to save image
    '''
    if isinstance(X, np.ndarray):    
        if stratified:
            cv = StratifiedKFold(n_splits=cv_splits)

            tprs = []
            aucs = []
            mean_fpr = np.linspace(0, 1, 100)

            fig, ax = plt.subplots(figsize=(6, 6))
            for fold, (train, test) in enumerate(cv.split(X, y)):
                viz = RocCurveDisplay.from_estimator(
                    model,
                    X[train],
                    y[train],
                    name=f"ROC fold {fold}",
                    alpha=0.3,
                    lw=1,
                    ax=ax,
                )
                plt.cla() # This removes each individual interation
                interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
                interp_tpr[0] = 0.0
                tprs.append(interp_tpr)
                aucs.append(viz.roc_auc)
            ax.plot([0, 1], [0, 1], "k--", label="chance level (AUC = 0.5)")
            
            mean_tpr = np.mean(tprs, axis=0)
            mean_tpr[-1] = 1.0
            mean_auc = auc(mean_fpr, mean_tpr)
            std_auc = mean_confidence_interval(aucs)[0] - mean_confidence_interval(aucs)[1]

            ax.plot(
                mean_fpr,
                mean_tpr,
                color="b",
                label=f"ROC (AUC = {mean_auc:0.2f} $\pm$ {std_auc:0.2f})",
                lw=2,
                alpha=0.8,
            )

            # This is some array magic to make mean_confidence_interval play well along a numpy axis=0 selection
            tprs = np.array(tprs)
            ci_tpr = [mean_confidence_interval(tprs[:,idx].flatten()) for idx in range(tprs.shape[1])]
            ci_tpr = np.array(ci_tpr)[:,0] - np.array(ci_tpr)[:,1] # mean - (mean + ci)

            tprs_upper = np.minimum(mean_tpr + ci_tpr, 1)
            tprs_lower = np.maximum(mean_tpr - ci_tpr, 0)
            ax.fill_between(
                mean_fpr,
                tprs_lower,
                tprs_upper,
                color="grey",
                alpha=0.2,
                label="95% Confidence Interval",
            )

            if not title:
                title="Mean ROC curve with 95% Confidence Interval"
            ax.set(
                xlim=[0, 1],
                ylim=[0, 1],
                xlabel="False Positive Rate",
                ylabel="True Positive Rate",
                title="Mean ROC curve with 95% Confidence Interval",
            )
            ax.axis("square")
            ax.grid(False)
            ax.legend(loc="lower right")
            if save_path:
                plt.savefig(save_path)
            else:
                plt.show()

        if not stratified:
            tprs = []
            fprs = []
            aucs = []
            mean_fpr = np.linspace(0, 1, cv_splits)
            working_roc = roc_auc_score(y, model.predict_proba(X)[:,1])
            rng = np.random.RandomState(random_state)
            fig, ax = plt.subplots(figsize=(6, 6))

            for i in range(cv_splits):
                # bootstrap by sampling with replacement on the prediction indices
                indices = rng.randint(0, len(y), len(y))
                fpr, tpr, _ = roc_curve(y[indices], model.predict_proba(X[indices])[:,1])

                interp_tpr = np.interp(mean_fpr, fpr, tpr)
                interp_tpr[0] = 0.0
                tprs.append(interp_tpr)
                aucs.append( roc_auc_score(y[indices], model.predict_proba(X[indices])[:,1]) )
            ax.plot([0, 1], [0, 1], "k--", label="chance level (AUC = 0.5)")
            
            mean_tpr = np.mean(tprs, axis=0)
            mean_tpr[-1] = 1.0
            mean_auc = auc(mean_fpr, mean_tpr)
            std_auc = mean_confidence_interval(aucs)[0] - mean_confidence_interval(aucs)[1]
            
            # print(mean_fpr)
            # print(mean_tpr)
            
            ax.plot(
                mean_fpr,
                mean_tpr,
                color="b",
                label=f"ROC (AUC = {mean_auc:0.2f} $\pm$ {0.02:0.2f})",
                lw=2,
                alpha=0.8,
            )

            # This is some array magic to make mean_confidence_interval play well along a numpy axis=0 selection

            tprs = np.array(tprs)
            tprs_lower = []
            tprs_upper = []
            for idx in range(tprs.shape[1]):
                tpr_row = tprs[:,idx].flatten()
                ci_low, ci_up = st.t.interval(0.95, len(tpr_row)-1, loc=np.mean(tpr_row), scale=st.sem(tpr_row))
                # print(ci_low, ci_up, mean_tpr[idx])
                tprs_lower.append(ci_low)
                tprs_upper.append(ci_up)
            
            ax.fill_between(
                mean_fpr,
                tprs_lower,
                tprs_upper,
                color="grey",
                alpha=0.2,
                label="95% Confidence Interval",
            )

            if not title:
                title="Mean ROC curve with 95% Confidence Interval"
            ax.set(
                xlim=[0, 1],
                ylim=[0, 1],
                xlabel="False Positive Rate",
                ylabel="True Positive Rate",
                title="Mean ROC curve with 95% Confidence Interval",
            )
            ax.axis("square")
            ax.grid(False)
            ax.legend(loc="lower right")
            if save_path:
                plt.savefig(save_path)
            else:
                plt.show()

    if isinstance(X, pd.DataFrame):    
        if stratified:
            cv = StratifiedKFold(n_splits=cv_splits)

            tprs = []
            aucs = []
            mean_fpr = np.linspace(0, 1, 100)

            fig, ax = plt.subplots(figsize=(6, 6))
            for fold, (train, test) in enumerate(cv.split(X, y)):
                viz = RocCurveDisplay.from_estimator(
                    model,
                    X.iloc[train],
                    y[train],
                    name=f"ROC fold {fold}",
                    alpha=0.3,
                    lw=1,
                    ax=ax,
                )
                plt.cla() # This removes each individual interation
                interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
                interp_tpr[0] = 0.0
                tprs.append(interp_tpr)
                aucs.append(viz.roc_auc)
            ax.plot([0, 1], [0, 1], "k--", label="chance level (AUC = 0.5)")
            
            mean_tpr = np.mean(tprs, axis=0)
            mean_tpr[-1] = 1.0
            mean_auc = auc(mean_fpr, mean_tpr)
            std_auc = mean_confidence_interval(aucs)[0] - mean_confidence_interval(aucs)[1]

            ax.plot(
                mean_fpr,
                mean_tpr,
                color="b",
                label=f"ROC (AUC = {mean_auc:0.2f} $\pm$ {std_auc:0.2f})",
                lw=2,
                alpha=0.8,
            )

            # This is some array magic to make mean_confidence_interval play well along a numpy axis=0 selection
            tprs = np.array(tprs)
            ci_tpr = [mean_confidence_interval(tprs[:,idx].flatten()) for idx in range(tprs.shape[1])]
            ci_tpr = np.array(ci_tpr)[:,0] - np.array(ci_tpr)[:,1] # mean - (mean + ci)

            tprs_upper = np.minimum(mean_tpr + ci_tpr, 1)
            tprs_lower = np.maximum(mean_tpr - ci_tpr, 0)
            ax.fill_between(
                mean_fpr,
                tprs_lower,
                tprs_upper,
                color="grey",
                alpha=0.2,
                label="95% Confidence Interval",
            )

            if not title:
                title="Mean ROC curve with 95% Confidence Interval"
            ax.set(
                xlim=[0, 1],
                ylim=[0, 1],
                xlabel="False Positive Rate",
                ylabel="True Positive Rate",
                title="Mean ROC curve with 95% Confidence Interval",
            )
            ax.axis("square")
            ax.grid(False)
            ax.legend(loc="lower right")
            if save_path:
                plt.savefig(save_path)
            else:
                plt.show()

        if not stratified:
            tprs = []
            fprs = []
            aucs = []
            mean_fpr = np.linspace(0, 1, cv_splits)
            working_roc = roc_auc_score(y, model.predict_proba(X)[:,1])
            rng = np.random.RandomState(random_state)
            fig, ax = plt.subplots(figsize=(6, 6))

            for i in range(cv_splits):
                # bootstrap by sampling with replacement on the prediction indices
                indices = rng.randint(0, len(y), len(y))
                fpr, tpr, _ = roc_curve(y[indices], model.predict_proba(X.iloc[indices])[:,1])

                interp_tpr = np.interp(mean_fpr, fpr, tpr)
                interp_tpr[0] = 0.0
                tprs.append(interp_tpr)
                aucs.append( roc_auc_score(y[indices], model.predict_proba(X.iloc[indices])[:,1]) )
            ax.plot([0, 1], [0, 1], "k--", label="chance level (AUC = 0.5)")
            
            mean_tpr = np.mean(tprs, axis=0)
            mean_tpr[-1] = 1.0
            mean_auc = auc(mean_fpr, mean_tpr)
            std_auc = mean_confidence_interval(aucs)[0] - mean_confidence_interval(aucs)[1]
            
            # print(mean_fpr)
            # print(mean_tpr)
            
            ax.plot(
                mean_fpr,
                mean_tpr,
                color="b",
                label=f"ROC (AUC = {mean_auc:0.2f} $\pm$ {0.02:0.2f})",
                lw=2,
                alpha=0.8,
            )

            # This is some array magic to make mean_confidence_interval play well along a numpy axis=0 selection

            tprs = np.array(tprs)
            tprs_lower = []
            tprs_upper = []
            for idx in range(tprs.shape[1]):
                tpr_row = tprs[:,idx].flatten()
                ci_low, ci_up = st.t.interval(0.95, len(tpr_row)-1, loc=np.mean(tpr_row), scale=st.sem(tpr_row))
                # print(ci_low, ci_up, mean_tpr[idx])
                tprs_lower.append(ci_low)
                tprs_upper.append(ci_up)
            
            ax.fill_between(
                mean_fpr,
                tprs_lower,
                tprs_upper,
                color="grey",
                alpha=0.2,
                label="95% Confidence Interval",
            )

            if not title:
                title="Mean ROC curve with 95% Confidence Interval"
            ax.set(
                xlim=[0, 1],
                ylim=[0, 1],
                xlabel="False Positive Rate",
                ylabel="True Positive Rate",
                title="Mean ROC curve with 95% Confidence Interval",
            )
            ax.axis("square")
            ax.grid(False)
            ax.legend(loc="lower right")
            if save_path:
                plt.savefig(save_path)
            else:
                plt.show()


def plot_training_probas(
    model,
    X, y,
    cv_splits = 6,
    plot = 'boxplot',
    title = None,
    save_path = None,
    ):
    '''
    model (joblib) : model to be copied for determining training confidence
    X (np.array or pd.DataFrame) : np array of features used
    y (np.array) : list of labels used for classes
    save_path (str) : string pointing where to save image
    '''
    if not isinstance(X, np.ndarray):
        X = X.to_numpy()

    cv = StratifiedKFold(n_splits=cv_splits)
    model = clone(model)

    sns.set_theme(style="whitegrid", palette="colorblind")

    scores = []
    labels = []
    fig, ax = plt.subplots(figsize=(6, 6))
    for fold, (train, test) in enumerate(cv.split(X, y)):
        model.fit(X[train], y[train])

        labels += list(y[test])
        scores += list(model.predict_proba(X[test])[:,1])

    df = pd.DataFrame({'preds' : scores}).apply(zscore)
    df['labels'] = labels
    df['labels'] = df['labels'].map({0:"Normal", 1:"Hyper"})

    if plot == 'boxplot':
        sns.boxplot(
            df, x='labels', y='preds', 
            order=["Normal", "Hyper"], boxprops={'alpha': 0.75}
            ).set(ylabel="Prediction Score", xlabel=None, title=title)
    if plot == 'violinplot':
        sns.violinplot(
            df, x='labels', y='preds', 
            order=["Normal", "Hyper"], boxprops={'alpha': 0.75}
            ).set(ylabel="Prediction Score", xlabel=None, title=title)

    kruskal_wallis = kruskal(df.loc[df['labels'] == "Normal"]['preds'], df.loc[df['labels'] == "Hyper"]['preds'])
    if kruskal_wallis[1] < 0.05:
        # statistical annotation
        x1, x2 = 0, 1
        y, h = df['preds'].max()+0.02, 0.01
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.2, c='k')
        plt.text((x1+x2)*.5, y+h, f"p = {kruskal_wallis[1]:.4f}",
                 ha='center', va='bottom', color='k')
        print(f'p-value is {kruskal_wallis[1]}')
    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()

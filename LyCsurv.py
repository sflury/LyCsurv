import numpy as np
import pandas as pd
import sys
from lifelines import CoxPHFitter,WeibullAFTFitter
'''
Name:
    interp_ph

Purpose:
    Interpolate over the survival function for each input observation to predict
    the appropriate escape fraction.

Arugments:
    :dat (*pandas.DataFrame*): pandas DataFrame of observables used to train model
    :part (**): Cox partial proportional hazard values for survival function
    :base (**): predicted response corresponding to full range of observations

Returns:
    :predict (*np.ndarray*): Nx3 array of response values predicted by the Cox PH
            model. First row is the lower uncertainty. Middle row is the median.
            Last row is the upper uncertainty.
'''
def interp_ph(dat,part,base):
    # storage array of predicted values
    f = np.zeros((len(dat),4))
    # for each object
    for i in dat.index:
        # calculate the survival function
        surv = np.exp(-1.0*base*part[i])['baseline cumulative hazard'].to_numpy()
        # predicted "expected" values and related quantiles
        f[i,:3] = np.interp([0.1587,0.5,0.8413],surv[::-1],1-base.index[::-1])
        if np.max(surv) < 0.5:    f[i,3] = -1.
        if np.max(surv) < 0.1587: f[i,0] = 0.
        if np.min(surv) > 0.5:    f[i,3] = 1.
        if np.min(surv) > 0.8413: f[i,2] = 1.
    # convert percentiles into uncertainties
    f[:,[0,2]] = abs(np.diff(f[:,:3],axis=1))
    return f
'''
Name:
    cox_ph

Purpose:
    Perform a Cox proportional hazards regression on the input reference data

Arguments:
    :dat (*pandas.DataFrame*): pandas DataFrame containing columns named
            according to conventions in `params.lis` with values corresponding
            to the galaxy sample which the user desires to fit

Keyword Arguments:
    :resp (*str*): string indicating the desired response variable. Options are
            'f_esc(LyC)', 'f_esc(LyA)', 'f(LyC)', and 'f(LyA)'. Default is
            'f_esc(LyC)'.
    :verbose (*bool*): boolean indicating whether to print out details of the
            Cox proportional hazards regression. Default is `False`.

Returns:
    :cph_fit (*numpy.ndarray*): n by 4 array containing the median, lower and
            upper uncertainties corresponding to 0.16 and 0.84 quantiles, and
            an indicator of whether the survival function is always below (-1)
            or above (+1) the median of the predicted distribution
'''
def cox_ph(dat,resp='f_esc(LyC)',verbose=False):
    # read in data
    trn = pd.read_csv('./tab/lzlcs.csv')
    # censor array based on LyC detection
    if 'LyC' in resp:
        trn['censors'] = trn['P(>N|B)'] < 0.02275#0.00135#
    elif 'LyA' in resp:
        trn['censors'] = trn['f_esc(LyA)'] > 0.0
        trn['f_esc(LyA)'][trn['f_esc(LyA)']<=0] = abs(trn['f_esc(LyA)'][trn['f_esc(LyA)']<=0])
    c = np.where(trn['censors'].values)[0]
    u = [i for i in trn.index if i not in c]
    # predictor variables
    pred = [line.strip() for line in open('./tab/params.lis').readlines()
            if not line.startswith('#')]
    trn = trn.dropna(subset=pred).reset_index()
    # add in weights
    trn['weights'] = trn[resp+' err'].values**-2
    # switch to left censoring
    trn[resp] = 1-trn[resp]
    # Cox proportional hazard fit
    cph = CoxPHFitter()
    cph.fit(trn[pred+[resp,'censors']],resp,'censors',robust=True)#,weights_col='weights')
    if verbose:
        cph.print_summary()
        print('---\nx̄')
        print(trn[pred].mean(axis=0))
    # results, including quantiles on fesc predictions
    base = cph.baseline_cumulative_hazard_
    return interp_ph(dat[pred],cph.predict_partial_hazard(dat[pred]),base)
'''
Name:
    aft

Purpose:
    Perform parametric survival regression using the Accelerated Failure method
    assuming a generic Weibull distribution.

Arguments:
    :dat (*pandas.DataFrame*): pandas DataFrame containing columns named
            according to conventions in `params.lis` with values corresponding
            to the galaxy sample which the user desires to fit

Keyword Arguments:
    :resp (*str*): string indicating the desired response variable. Options are
            'f_esc(LyC)', 'f_esc(LyA)', 'f(LyC)', and 'f(LyA)'. Default is
            'f_esc(LyC)'.
    :verbose (*bool*): boolean indicating whether to print out details of the
            accelerated failure regression. Default is `False`.
    :intercept (*bool*): boolean indicating whether to include an intercept in
            the parametric hazard function

Returns:
    :aft_fit (*np.ndarray*): n by 3 array containing the median, lower and
            upper uncertainties corresponding to 0.16 and 0.84 quantiles
'''
def aft(dat,resp='f_esc(LyC)',verbose=False,intercept=True):
    if intercept:
        dat['Intercept'] = np.ones(len(dat))
    # read in data
    trn = pd.read_csv('./tab/lzlcs.csv')
    # censor array based on LyC detection
    if 'LyC' in resp:
        trn['censors'] = trn['P(>N|B)']>0.002275
    elif 'LyA' in resp:
        trn['censors'] = trn['f_esc(LyA)'] > 0.0
        trn['f_esc(LyA)'][trn['f_esc(LyA)']<=0] = abs(trn['f_esc(LyA)'][trn['f_esc(LyA)']<=0])
    c = np.where(trn['censors'].values)[0]
    u = [i for i in trn.index if i not in c]
    # predictor variables
    pred = [line.strip() for line in open('./tab/params.lis').readlines()
            if not line.startswith('#')]
    # add in weights
    trn['weights'] = trn[resp+' err'].values**-2
    # add in fit
    aft = WeibullAFTFitter()
    aft.fit(trn[pred+[resp,'censors']],duration_col=resp,event_col='censors',\
            fit_intercept=intercept)
    if verbose:
        aft.print_summary()
    #
    if intercept:
        pred += ['Intercept']
    #return aft.params_['lambda_'].values @ dat[pred].values.T
    f = np.zeros((len(dat),3))
    for i,q in enumerate([0.8413,0.5,0.1587]):
        f[:,i] = aft.predict_percentile(dat,p=q)
    f[f>1] = 1
    f[f<0] = 0
    f[:,[0,2]] = abs(np.diff(f,axis=1))
    return f

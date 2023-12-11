import numpy as np
import pandas as pd
import sys
from lifelines import CoxPHFitter,WeibullAFTFitter
# response variable (must call using quotes because of dashes and underscores)
resp = sys.argv[1]
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
    predict = np.zeros((len(dat),4))
    # for each object
    for i in dat.index:
        # calculate the survival function
        surv = np.exp(-1.0*base*part[i])['baseline cumulative hazard'].to_numpy()
        # predicted "expected" values and related quantiles
        predict[i,:3] = np.interp([0.8413,0.5,0.1587],surv[::-1],base.index[::-1])
        if np.max(surv) < 0.5: predict[i,3] = -1.
        if np.max(surv) < 0.1587: predict[i,0] = 0.
        if np.min(surv) > 0.5: predict[i,3] = 1.
        if np.min(surv) > 0.8413: predict[i,2] = 1.
    # convert percentiles into uncertainties
    predict[:,[0,2]] = abs(np.diff(predict[:,:3],axis=1))
    return predict
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
    :verbose (*bool*): boolean indicating whether to print out details of the
            Cox proportional hazards regression. Default is `False`.

Returns:
    :cph_fit (*numpy.ndarray*): n by 4 array containing the median, lower and
            upper uncertainties corresponding to 0.16 and 0.84 quantiles, and
            an indicator of whether the survival function is always below (-1)
            or above (+1) the median of the predicted distribution
'''
def cox_ph(dat,verbose=False):
    # read in data
    trn = pd.read_csv('./tab/lzlcs_surv.csv')
    # censor array based on LyC detection
    if 'LyC' in resp:
        trn['censors'] = trn['P(>N|B)'] < 0.02275#0.00135#
    elif 'LyA' in resp:
        trn['censors'] = trn['f_esc(LyA)'] > 0.0
        trn['f_esc(LyA)'][trn['f_esc(LyA)']<=0] = abs(trn['f_esc(LyA)'][trn['f_esc(LyA)']<=0])
    c = np.where(trn['censors'].values)[0]
    u = [i for i in trn.index if i not in c]
    # predictor variables
    pred = [line.strip() for line in open('./tab/pars.lis').readlines()
            if not line.startswith('#')]
    # add in weights
    trn['weights'] = trn[resp+' err'].values**-2
    # Cox proportional hazard fit
    cph = CoxPHFitter()
    cph.fit(trn[pred+[resp,'censors']],resp,'censors',robust=True)#,weights_col='weights')
    if verbose:
        cph.print_summary()
    # results, including quantiles on fesc predictions
    base = cph.baseline_cumulative_hazard_
    return interp_ph(dat[pred],cph.predict_partial_hazard(dat[pred]),base)
'''
Name:
    aft

Purpose:
    Perform parametric survival regression using the Accelerated Failure method
    assuming a generic Weibull distribution.
'''

dat[resp+'-cph']     = cph_interp[:,1]
dat[resp+'-cph.lo']  = cph_interp[:,0]
dat[resp+'-cph.up']  = cph_interp[:,2]
dat[resp+'-cph.lim'] = cph_interp[:,3]

cols = ['Object'] + [key for key in trn.keys() if resp in key]
dat[cols].to_csv(f'./out/dat_pred_{resp}.csv',index=False,float_format='%.6f')

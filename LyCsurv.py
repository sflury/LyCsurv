import numpy as np
import pandas as pd
import sys
from lifelines import CoxPHFitter,WeibullAFTFitter

def InterpPH(dat,part,base):
    '''
    Name:
        InterpPH

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

def ModAssess(trn,mod,cens,p,concord='harrell'):
    '''
    Name:
        ModAssess

    Purpose:
        Assess the quality of a survival model by testing it against the training
        data set (in this case, the LzLCS).

    Arguments:
        :trn (*np.ndarray*): Nx1 array of observed values
        :mod (*np.ndarray*): Nx1 array of predicted values
        :cens (*np.ndarray*): Nx1 array of censors as booleans

    Keyword Arguments:
        :concord (*str*): string indicating the method to use for concordance
            calculation. Currently supports Harrell+ 1996 and Uno+ 2011 methods.
            Default is \'harrell\'.

    Returns:
        :R2  (*float*): R^2 metric of the residuals
        :R2a (*float*): adjusted R^2 metric of the residuals
        :RMS (*float*): root-mean-square of the residuals
        :C   (*float*): concordance index
    '''
    trnLog = np.log10(trn)
    modLog = np.log10(mod)
    n = len(trn)
    R2   = 1 - np.sum(np.square(trnLog-modLog)) / \
               np.sum(np.square(trnLog-trnLog.mean()))
    R2a  = 1 - (1-R2)*(n-1)/(n-p-1)
    RMS  = np.sqrt(np.sum(np.square(trnLog-modLog))/n)
    if concord.lower() == 'harrell':
        from sksurv.metrics import concordance_index_censored
        C = concordance_index_censored(cens,trnLog,modLog)[0]
    elif concord.lower() == 'uno':
        from sksurv.metrics import concordance_index_ipcw
        C = concordance_index_ipcw(trnLog,modLog)[0]
    #elif concord.lower() == 'auc':
    #    from sksurv.metrics import cumulative_dynamic_auc
    return R2,R2a,RMS,C

def CoxPH(dat,resp='f_esc(LyC)',verbose=False,StatsVerbose=False):
    '''
    Name:
        CoxPH

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
        :StatsVerbose (*bool*): boolean indicating whether to perform and return
                statistical assessments of the model using the training set.
                Default is `False`.

    Returns:
        :cph_fit (*numpy.ndarray*): n by 4 array containing the median, lower and
                upper uncertainties corresponding to 0.16 and 0.84 quantiles, and
                an indicator of whether the survival function is always below (-1)
                or above (+1) the median of the predicted distribution
        :ModAssess (*tuple*): (_optional_) 4x1 tuple of model assessments
                containing the R^2, adjusted R^2, RMS, and concordance index.
                Only returned if `StatsVerbose` set to `True`.
    '''
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
    if StatsVerbose:
        mod = InterpPH(dat[pred],cph.predict_partial_hazard(dat[pred]),base)
        return mod, ModAssess(1-trn[resp],1-mod[:,1],len(pred),trn['censors'])
    else:
        return InterpPH(dat[pred],cph.predict_partial_hazard(dat[pred]),base)

def AFT(dat,resp='f_esc(LyC)',verbose=False,intercept=True,StatsVerbose=False):
    '''
    Name:
        AFT

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
        :StatsVerbose (*bool*): boolean indicating whether to perform and return
                statistical assessments of the model using the training set.
                Default is `False`.

    Returns:
        :aft_fit (*np.ndarray*): n by 3 array containing the median, lower and
                upper uncertainties corresponding to 0.16 and 0.84 quantiles
        :ModAssess (*tuple*): (_optional_) 4x1 tuple of model assessments
                containing the R^2, adjusted R^2, RMS, and concordance index.
                Only returned if `StatsVerbose` set to `True`.
    '''
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
    if StatsVerbose:
        return f, ModAssess(trn[resp],f[:,1],len(pred),trn['censors'])
    else:
        return f

class Train(object):
    '''
    Name:
        Train

    Purpose:
        Train a specified survival model on reference data and assess the results.

    Keyword Arguments:
        :resp (*str*): string indicating the desired response variable. Options are
                'f_esc(LyC)', 'f_esc(LyA)', 'f(LyC)', and 'f(LyA)'. Default is
                'f_esc(LyC)'.
        :method (*str*): string indicating the method of survival analysis to be
                used in the training run: 'CoxPH' or 'AFT'. Default is 'CoxPH'.

    Attributes:
        :train (*numpy.ndarray*): 88x2 array containing the observed and predicted
                response variable
        :stats (*tuple*): (_optional_) 4x1 tuple of model assessments
                containing the R^2, adjusted R^2, RMS, and concordance index.
                Only returned if `StatsVerbose` set to `True`.
        :resp (*str*): response variable (corresponds to `resp` input)
        :meth (*str*): method used for training (corresponds to `method` input)
    '''
    def __init__(self,resp='f_esc(LyC)',method='CoxPH',intercept=True,verbose=False):
        self.resp = resp
        self.meth = method
        # read in data
        trn = pd.read_csv('./tab/lzlcs.csv')
        # censor array based on LyC detection
        if 'LyC' in resp:
            trn['censors'] = trn['P(>N|B)'] < 0.02275#0.00135#
        elif 'LyA' in resp:
            trn['censors'] = trn['f_esc(LyA)'] > 0.0
            trn['f_esc(LyA)'][trn['f_esc(LyA)']<=0] = abs(trn['f_esc(LyA)'][trn['f_esc(LyA)']<=0])
        self.c = np.where(trn['censors'].values)[0]
        self.u = [i for i in trn.index if i not in self.c]
        # predictor variables
        pred = [line.strip() for line in open('./tab/params.lis').readlines()
                if not line.startswith('#')]
        trn = trn.dropna(subset=pred).reset_index()
        # add in weights
        trn['weights'] = trn[resp+' err'].values**-2
        # Cox proportional hazard fit
        if method == 'CoxPH':
            # switch to left censoring
            trn[resp] = 1-trn[resp]
            cph = CoxPHFitter()
            cph.fit(trn[pred+[resp,'censors']],resp,'censors',robust=True)#,weights_col='weights')
            if verbose:
                cph.print_summary()
                print('---\nx̄')
                print(trn[pred].mean(axis=0))
            # results, including quantiles on fesc predictions
            base = cph.baseline_cumulative_hazard_
            mod = InterpPH(trn[pred],cph.predict_partial_hazard(trn[pred]),base)[:,1]
            self.stats = ModAssess(trn[resp],1-mod,trn['censors'],len(pred))
            trn[resp] = 1-trn[resp]
        elif method == 'AFT':
            aft = WeibullAFTFitter()
            aft.fit(trn[pred+[resp,'censors']],duration_col=resp,event_col='censors',\
                    fit_intercept=intercept)
            if verbose:
                aft.print_summary()
            #
            if intercept:
                pred += ['Intercept']
            #return aft.params_['lambda_'].values @ dat[pred].values.T
            mod = aft.predict_percentile(trn,p=0.5)
            self.stats = ModAssess(trn[resp],mod,trn['censors'],len(pred))

        else:
            print('Model method not recognized. Supported methods are either\n'+\
                  '\'CoxPH\' for Cox proportional hazards or \'AFT\' for \n'+\
                  'accelerated failure time.')
        self.train = pd.DataFrame({f'{resp} obs':trn[resp],f'{resp} pred':mod})

    def pprint(self):
        for s,l in zip(self.stats,['R^2','adj R^2','RMS','concord']):
            print(f'{l: >10s}  :  {s:5.3f}')

    def plot(self):
        import matplotlib.pyplot as plt
        line = [self.train.iloc[:,0].min()//1,self.train.iloc[:,0].max()//1+1]
        plt.loglog(line,line,color='xkcd:cerulean')
        plt.scatter(self.train.iloc[self.u,0],self.train.iloc[self.u,1],\
                        marker='$\u21A4$',edgecolor='xkcd:peach',s=100,zorder=4)
        plt.scatter(self.train.iloc[self.c,0],self.train.iloc[self.c,1],\
                        c='xkcd:peach',edgecolor='black',zorder=4)
        plt.xlabel(self.train.keys()[0])
        plt.ylabel(self.train.keys()[1])
        if 'f_esc' in self.resp:
            plt.xticks([0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1.0],\
                       ['1e-3','2e-3','5e-3',0.01,0.02,0.05,0.1,0.2,0.5,1.0])
            plt.yticks([0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1.0],\
                       ['1e-3','2e-3','5e-3',0.01,0.02,0.05,0.1,0.2,0.5,1.0])
            plt.xlim(0.0075,1)
            plt.ylim(0.0075,1)
            plt.gca().set_aspect('equal')
        plt.subplots_adjust(bottom=0.15,left=0.18,right=0.95,top=0.95)
        plt.show()

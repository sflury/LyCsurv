# LyCsurv

Based on Jaskot et al. (2024a,b), the `LyCsurv` code predicts Lyman continuum (LyC) escape fractions (fesc) from given input variables using the Low-redshift Lyman Continuum Survey (LzLCS, Flury et al 2022a) combined with results from the literature. The authors request that any use of this code cite Jaskot et al. (2024a). A BibTeX reference is provided at the end of this document.

## Practical Details

The `LyCsurv` scripts depend on several software packages which are readily installed using `pip` from the command line. These packages include `numpy`, `matplotlib`, `lifelines`, and `scikit-survival`. For example
```
$ pip install lifelines
```

Calling the `LyCsurv.CoxPH` or `LyCsurv.AFT` functions will return an Nx3 array of the predicted fesc and its lower and upper uncertainties in the first, zeroth, and second rows, respectively, for each of N galaxies. In the worked examples below, N=1.

All files and scripts necessary to run `LyCsurv` are contained in this repository, including combined LzLCS+archival data (`./tab/lzlcs.csv`), a file to control input parameters (`./tab/params.lis`), and the source code (`LyCsurv.py`). To run `LyCsurv`, simply download the zipped directory or clone the directory using `git`. To run locally, navigate into the `LyCsurv` (or `LyCsurv-main`) directory and simply use an `ipython` terminal or write a new script based on the examples below.

Variables for input are controlled via a table `./tab/params.lis`. To exclude a variable from the fitting routine, simply include a \# at the start of the line just as in python comments. To include a variable, simply delete the \#. In the examples below, only the `O32` and `beta1550` variables are included. The params file is the only file the user should alter in any way.

Input data should come in the form of a dictionary or `pandas` dataframe (the latter is recommended). All columns to be used with `LyCsurv` must be named such that they match the names in the `params.lis` table.

For details, API documentation for `LyCsurv` is available [here](https://github.com/sflury/LyCsurv/wiki/API).

### Example Usage - Cox Proportional Hazards
``` python
>>> import pandas as pd
>>> from LyCsurv import *
>>> # pandas DataFrame
>>> # should contain user’s input data table
>>> # with measurements for the required input variables
>>> dat = pd.DataFrame({'O32':[1.2],'beta1550':[-2.4]})
>>> # Cox proportional hazards fit
>>> cph_fit = CoxPH(dat)
>>> print(f'proportional hazards fit fesc(LyC) : {cph_fit[:,1][0]:.3f}'+\
...       f'-{cph_fit[:,0][0]:.3f}+{cph_fit[:,2][0]:.3f}')
```

which prints the $f_{esc}^{\rm LyC}$ predicted by the Cox PH model to the command line as

```
proportional hazards fit fesc(LyC) : 0.147-0.106+0.359
```

Here, the output from `LyCsurv.cox_ph` as printed indicates a predicted $f_{esc}^{\rm LyC} = 0.146^{+0.359}_{-0.106}$.

### Example Usage - Accelerated Failure Time
``` python
>>> import pandas as pd
>>> from LyCsurv import *
>>> # pandas DataFrame
>>> # should contain user’s input data table
>>> # with measurements for the required input variables
>>> dat = pd.DataFrame({'O32':[1.2],'beta1550':[-2.4]})
>>> # accelerated failure time with Weibull distribution
>>> aft_fit = AFT(dat)
>>> print(f'accelerated failures fit fesc(LyC) : {aft_fit[:,1][0]:.3f}'+\
...       f'-{aft_fit[:,0][0]:.3f}+{aft_fit[:,2][0]:.3f}')
```

which prints the $f_{esc}^{\rm LyC}$ predicted by the AFT model to the command line as

``` 
accelerated failures fit fesc(LyC) : 0.369-0.312+0.631
```

Here, the output from `LyCsurv.AFT` as printed indicates a predicted $f_{esc}^{\rm LyC} = 0.369^{+0.631}_{-0.312}$.

### Example Usage - Statistical Assessments
``` python
>>> import pandas as pd
>>> from LyCsurv import *
>>> # assess quality of trained model by
>>> # instantiating LyCsurv.Train
>>> train = Train(method='AFT')
>>> train.pprint()
>>> train.plot()
```

which prints to the terminal command line using the `pprint` method

```
     R^2  :  0.482
 adj R^2  :  0.464
     RMS  :  0.612
 concord  :  0.209
```

and displays this figure using the `plot` method

!['example of a training model using UV slope and O32 ratio'](train_examp.png)


## The `params.lis` Input File

Included in the `./tab/` directory of `LyCsurv` is a table containing all of the possible inputs for predictor variables named according to the column headings in the reference training set of LzLCS+ galaxies contained in `/tab/lzlcs.csv`. To toggle a variable "off", add a \# before the variable name. To toggle a variable "on", simply omit the \#.

## Motivation

Survival analysis as implemented here and in the underlying `lifelines` software accounts for censoring (limits) in the response variable (the "y" variable) and allows for predictions incorporating multiple predictor variables (the "x" variables). Accounting for upper limits is critical for inferring LyC escape fractions due to the large number of non-detections in the sample. Allowing for a multivariate predictor can also account for the scatter seen in LyC escape diagnostics due to variations in galaxy properties (Flury et al 2022b). The methods outlined here are semi-parametric (Cox model of proportional hazards) and parametric (accelerated failure time model) in their treatment of the so-called hazard function $\lambda$. Each serves a unique purpose in the sort of LyC escape fraction predicted. We leave the selection of method to the user but provide details below.

## Some Technical Stats Details

The hazard $\lambda$ is definied as the ratio of the probability $P$ of response (in this case, LyC escape fraction $f_{esc}^{\rm LyC}$) to the survival $S=1-{\rm CDF}$ of the response. Under this definition, $\lambda$ and $P$ are used to determine $S$ from a set of calibrator measurements used to perform regression. Regression predicts not the response but rather the survival of the response, which requires an additional step of determining the value of $f_{esc}^{\rm LyC}$ which corresponds to $S=0.1587,~0.5,~0.8413$.

### Cox Proportional Hazards

Cox proportional hazards provides a semi-parametric prediction of the hazard function $\lambda$. The parametric component is the partial hazard function, a portion of the hazard inferred from a set of predictors **x** which are linearly combined. The non-parametric component is the baseline hazard $\lambda_0$, a part of the hazard function based solely on given responses rank-ordered by their values and accounting for upper limits. Together, the baseline and partial hazard predict $S(f_{esc}^{\rm LyC}|\mathbf{x})$ such that

$$ S(1-f_{esc}^{\rm LyC}) = \lambda_0 \exp\left( \beta \cdot ( \mathbf{x} - \bar{\mathbf{x}} )^\prime  \right) .$$

A fundamental nuance of the Cox model is that the $f_{esc}^{\rm LyC}$ predicted is the *typical line-of-sight* $f_{esc}^{\rm LyC}$ for that set of input variables and gives its typical variation within the LzLCS. If the input variables are sensitive to LOS properties (e.g., absorption line depth), the predicted fesc will be LOS-specific too. If only global galaxy properties are used AND if we assume the LzLCS sufficiently probes all random orientations, then, yes, the predicted fesc = the global fesc. (see Jaskot et al. 2024 for discussion).

### Accelerated Failure

Like the Cox proportional hazards model, the accelerated failure time (AFT) approach accounts for upper limits in the response and allows for a multivariate predictor. However, unlike the Cox model, AFT uses a fully parametric description of the hazard function so that

$$ S(f_{esc}^{\rm LyC}) = \left(\frac{f_{esc}^{\rm LyC}}{\lambda}\right)^\rho$$

where

$$ \lambda = \exp(\beta \cdot \mathbf{x}^\prime) .$$

Under this formalism, the `lifelines` AFT implementation assumes a Weibull distribution.

Unlike the Cox model, the AFT $f_{esc}^{\rm LyC}$ is predominantly line-of-sight and thus typically exceeds the Cox $f_{esc}^{\rm LyC}$ as a result of anisotropic LyC escape (see Jaskot et al. 2024 for discussion).

## BibTeX References
``` bibtex
@software{2024zndo..11392442F,
       author = {{Flury}, Sophia and {Jaskot}, Anne and {Silveyra}, Anneliese},
        title = "{LyCsurv}",
         year = 2024,
        month = may,
          eid = {10.5281/zenodo.11392442},
          doi = {10.5281/zenodo.11392442},
      version = {v0.1.0},
    publisher = {Zenodo},
          url = {https://github.com/sflury/LyCsurv},
      license = {GPL-3.0-or-later},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2024zndo..11392442F},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
@ARTICLE{Jaskot2024a,
       author = {{Jaskot}, Anne E. and {Silveyra}, Anneliese C. and {Plantinga}, Anna and {Flury}, Sophia R. and {Hayes}, Matthew and {Chisholm}, John and {Heckman}, Timothy and {Pentericci}, Laura and {Schaerer}, Daniel and {Trebitsch}, Maxime and {Verhamme}, Anne and {Carr}, Cody and {Ferguson}, Henry C. and {Ji}, Zhiyuan and {Giavalisco}, Mauro and {Henry}, Alaina and {Marques-Chaves}, Rui and {{\"O}stlin}, G{\"o}ran and {Saldana-Lopez}, Alberto and {Scarlata}, Claudia and {Worseck}, G{\'a}bor and {Xu}, Xinfeng},
        title = "{Multivariate Predictors of Lyman Continuum Escape. I. A Survival Analysis of the Low-redshift Lyman Continuum Survey}",
      journal = {\apj},
     keywords = {Astrostatistics, Reionization, High-redshift galaxies, Starburst galaxies, Interstellar medium, Ultraviolet astronomy, Radiative transfer, 1882, 1383, 734, 1570, 847, 1736, 1335, Astrophysics - Astrophysics of Galaxies},
         year = 2024,
        month = sep,
       volume = {972},
       number = {1},
          eid = {92},
        pages = {92},
          doi = {10.3847/1538-4357/ad58b9},
archivePrefix = {arXiv},
       eprint = {2406.10171},
 primaryClass = {astro-ph.GA},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2024ApJ...972...92J},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
@ARTICLE{Jaskot2024b,
       author = {{Jaskot}, Anne E. and {Silveyra}, Anneliese C. and {Plantinga}, Anna and {Flury}, Sophia R. and {Hayes}, Matthew and {Chisholm}, John and {Heckman}, Timothy and {Pentericci}, Laura and {Schaerer}, Daniel and {Trebitsch}, Maxime and {Verhamme}, Anne and {Carr}, Cody and {Ferguson}, Henry C. and {Ji}, Zhiyuan and {Giavalisco}, Mauro and {Henry}, Alaina and {Marques-Chaves}, Rui and {{\"O}stlin}, G{\"o}ran and {Saldana-Lopez}, Alberto and {Scarlata}, Claudia and {Worseck}, G{\'a}bor and {Xu}, Xinfeng},
        title = "{Multivariate Predictors of Lyman Continuum Escape. II. Predicting Lyman Continuum Escape Fractions for High-redshift Galaxies}",
      journal = {\apj},
     keywords = {Astrostatistics, Reionization, High-redshift galaxies, Starburst galaxies, Interstellar medium, Ultraviolet astronomy, Radiative transfer, 1882, 1383, 734, 1570, 847, 1736, 1335, Astrophysics - Astrophysics of Galaxies},
         year = 2024,
        month = oct,
       volume = {973},
       number = {2},
          eid = {111},
        pages = {111},
          doi = {10.3847/1538-4357/ad5557},
archivePrefix = {arXiv},
       eprint = {2406.10179},
 primaryClass = {astro-ph.GA},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2024ApJ...973..111J},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

## Zenodo DOI

[![DOI](https://zenodo.org/badge/730334029.svg)](https://zenodo.org/doi/10.5281/zenodo.11392442)


## Licensing
<a href="https://github.com/sflury/LyCsurv">LyCsurv</a> © 2024 by <a href="https://sflury.github.io">Sophia Flury</a>, Anne Jaskot, and Anneliese Silvera is licensed under <a href="https://creativecommons.org/licenses/by-nc-nd/4.0/">Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International</a>

<img src="https://mirrors.creativecommons.org/presskit/icons/cc.svg" style="max-width: 1em;max-height:1em;margin-left: .2em;"><img src="https://mirrors.creativecommons.org/presskit/icons/by.svg" style="max-width: 1em;max-height:1em;margin-left: .2em;"><img src="https://mirrors.creativecommons.org/presskit/icons/nc.svg" style="max-width: 1em;max-height:1em;margin-left: .2em;"><img src="https://mirrors.creativecommons.org/presskit/icons/nd.svg" style="max-width: 1em;max-height:1em;margin-left: .2em;">

This license enables reusers to copy and distribute LyCsurv in any medium or format in unadapted form only, for noncommercial purposes only, and only so long as attribution is given to the creator. CC BY-NC-ND 4.0 includes the following elements:

BY: credit must be given to the creator.
NC: Only noncommercial uses of the work are permitted.
ND: No derivatives or adaptations of the work are permitted.

You should have received a copy of the CC BY-NC-ND 4.0 along with this program. If not, see <https://creativecommons.org/licenses/by-nc-nd/4.0/>.

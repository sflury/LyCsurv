# LyCsurv

This code is developed based on Jaskot et al. (2024) to predict Lyman continuum (LyC) escape fractions (fesc) from given predictors using the concatenated Low-redshift Lyman Continuum Survey (LzLCS, Flury et al 2022a). While this code is provided publicly, the authors request that any use thereof be cited in any publications in which this code is used. A BibTeX reference is provided at the end of this document.

Survival analysis as implemented here and in the underlying `lifelines` software accounts for censoring (limits) in the response variable (the "y" variable) and allows for the predictor variable (the "x" variables) to be multivariate. Accounting for upper limits is critical for inferring LyC escape fractions due to the large number of non-detections in the sample. Allowing for a multivariate predictor can also account for the scatter seen in LyC escape diagnostics due to variations in galaxy properties (Flury et al 2022b). The methods outlined here are semi-parametric (Cox model of proportional hazards) and parametric (accelerated failure time model) in their treatment of the so-called hazard function $\lambda$. Each serves a unique purpose in the sort of LyC escape fraction predicted. We leave the selection of method to the user but provide details below.

The hazard $\lambda$ is definied as the ratio of the probability $P$ of response (in this case, LyC escape fraction $f_{esc}^{\rm LyC}$) to the survival $S=1-{\rm CDF}$ of the response. Under this definition, $\lambda$ and $P$ are used to determine $S$ from a set of calibrator measurements used to perform regression. Regression predicts not the response but rather the survival of the response, which requires an additional step of determining the value of $f_{esc}^{\rm LyC}$ which corresponds to $S=0.1587,~0.5,~0.8413$.

Variables for input are controlled via a table `./tab/params.lis`. To exclude a variable from the fitting routine, simply include a \# at the start of the line just as in python comments. To include a variable, simply delete the \#. In the examples below, only the `O32` and `beta1550` variables are included.

## Cox Proportional Hazards

Cox proportional hazards provides a semi-parametric prediction of the hazard function $\lambda$. The parametric component is the partial hazard function, a portion of the hazard inferred from a set of predictors **x** which are linearly combined. The non-parametric component is the baseline hazard $\lambda_0$, a part of the hazard function based solely on given responses rank-ordered by their values and accounting for upper limits. Together, the baseline and partial hazard predict $S(f_{esc}^{\rm LyC}|**x**})$ such that

$$ S(f_{esc}^{\rm LyC}) = \lambda_0 \exp\left( \boldsymbol\beta \cdot ( \mathbf{x} - \bar{\mathbf{x}} )^\prime  \right) .$$

LyCsurv solves for the $f_{esc}^{\rm LyC}$ which satisfies the median and $1\sigma$ confidence intervals of $S$ given **x** observed for a user-provided sample of galaxy measurements.

A fundamental nuance of the Cox model is that the $f_{esc}^{\rm LyC}$ predicted corresponds to the *global* escape fraction rather than the line-of-sight (see Jaskot et al. 2024 for discussion).

### Example Usage
``` python
import pandas as pd
from LyCsurv import *
# pandas DataFrame of data
dat = pd.DataFrame({'O32':[1.2],'beta1550':[-2.4]})
# Cox proportional hazards fit
cph_fit = cox_ph(dat)
print(f'proportional hazards fit fesc(LyC) : {cph_fit[:,1][0]:.3f}'+\
      f'-{cph_fit[:,0][0]:.3f}+{cph_fit[:,2][0]:.3f}')
```

## Accelerated Failure

Like the Cox proportional hazards model, the accelerated failure time (AFT) approach accounts for upper limits in the response and allows for a multivariate predictor. However, unlike the Cox model, AFT uses a fully parametric description of the hazard function so that

$$ S(f_{esc}^{\rm LyC}) = \left(\frac{f_{esc}^{\rm LyC}}{\lambda}\right)^\rho$$

where

$$ \lambda = \exp(\boldsymbol\beta \cdot \mathbf{x}^\prime) .$$

LyCsurv solves for the $f_{esc}^{\rm LyC}$ which satisfies the median and $1\sigma$ confidence intervals of $S$ given **x** observed for a user-provided sample of galaxy measurements.

Unlike the Cox model, the AFT $f_{esc}^{\rm LyC}$ is predominantly line-of-sight and thus typically exceeds the Cox $f_{esc}^{\rm LyC}$ as a result of anisotropic LyC escape (see Jaskot et al. 2024 for discussion).

### Example Usage
``` python
import pandas as pd
from LyCsurv import *
# pandas DataFrame of data
dat = pd.DataFrame({'O32':[1.2],'beta1550':[-2.4]})
# accelerated failure time with Weibull distribution
aft_fit = aft(dat)
print(f'accelerated failures fit fesc(LyC) : {aft_fit[:,1][0]:.3f}'+\
      f'-{aft_fit[:,0][0]:.3f}+{aft_fit[:,2][0]:.3f}')
```

## BibTeX Citation
``` bibtex
@ARTICLE{2024ApJ...J,
       author = {{Jaskot}, Anne E., et al.},
        title = "{}",
      journal = {\apj},
         year = 2024,
        month = ,
       volume = {},
       number = {},
        pages = {},
          doi = {} }
}
```

## Licensing

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

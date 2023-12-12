# LyCsurv

This code is developed based on Jaskot et al. (2024) to predict Lyman continuum (LyC) escape fractions (fesc) from given predictors using the concatenated Low-redshift Lyman Continuum Survey (LzLCS, Flury et al 2022a). Two options for predictions are included: semi-parametric Cox proportional hazards models and parametric accelerated failure models. Variables for input are controlled via a table `params.lis`.

Survival analysis as implemented here and in the underlying `lifelines` software accounts for censoring (limits) in the response variable (the "ordinate" or "y" variable) and allows for the predictor variable (the "abscissa" or "x" variables) to be multivariate. Accounting for upper limits is critical for inferring LyC escape fractions due to the large number of non-detections in the sample. Allowing for a multivariate predictor can also account for the scatter seen in LyC escape diagnostics due to variations in galaxy properties (Flury et al 2022b). The methods outlined here are semi-parametric (Cox model of proportional hazards) and parametric (accelerated failure time model). Each serves a unique purpose in the sort of LyC escape fraction predicted. We leave the selection of method to the user but provide details below.

## Cox Proportional Hazards
The Cox proportional hazards method provides a semi-parametric prediction of the hazard function $\lambda$ -- the ratio of probability $P$ of response (in this case, LyC escape fraction $f_{esc}^{\rm LyC}$) to the survival $S$ of an event or measurement. The parametric component is the partial hazard function, a portion of the hazard inferred from a set of predictors **x** which are linearly combined. The non-parametric component is the baseline hazard $\lambda_0$, a part of the hazard function based solely on given responses rank-ordered by their values and accounting for upper limits. Together, the baseline and partial hazard predict $S(f_{esc}^{\rm LyC}|**x**})$ such that
$$ S(f_{esc}^{\rm LyC}) = \lambda_0 \exp\left[ **\beta** \cdot ( **x** - **\bar{x}** )^\prime  \right] $$
. LyCsurv solves for the $f_{esc}^{\rm LyC}$ which satisfies $S(f_{esc}^{\rm LyC})=0.1587,0.5,0.8413$ given **x** observed for a user-provided sample of galaxy measurements.

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
Like the Cox proportional hazards method, the accelerated failure time (AFT) approach accounts for upper limits in the response and allows for a multivariate predictor. However, unlike the Cox model, AFT uses a fully parametric description of the hazard function so that
$$ S(f_{esc}^{\rm LyC}) = \left(\frac{f_{esc}^{\rm LyC}}{\lambda}\right)^\rho$$
where
$$ \lambda = \exp(**\beta**\cdot**x**^\prime) $$
. As before, LyCsurv solves for the $f_{esc}^{\rm LyC}$ which satisfies $S(f_{esc}^{\rm LyC})=0.1587,0.5,0.8413$ given **x** observed for a user-provided sample of galaxy measurements.
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

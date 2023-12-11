# LyCsurv

This code is developed based on Jaskot et al. (2024) to predict Lyman continuum (LyC) escape fractions (fesc) from given predictors using the concatenated Low-redshift Lyman Continuum Survey (LzLCS, Flury et al 2022). Two options for predictions are included: semi-parametric Cox proportional hazards models and parametric accelerated failure models. Variables for input are controlled via a table `params.lis`.

## Cox Proportional Hazards

### Example Usage
``` python
from LyCsurv import *

```

## Accelerated Failure

### Example Usage
``` python
from LyCsurv import *
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

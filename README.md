ASSISTant
=========

`ASSISTant` is a package to assist in performing **A**daptive
**S**ubgroup **S**election **I**n **S**equential **T**rials.

Install this package (after taking care of dependencies shown) the
usual way in R once it gets on CRAN or via:

```{r}
require("mvtnorm")
require("R6")
library(devtools)
install_github("bnaras/ASSISTant")
```

`ASSISTant` provides a three-stage adaptive group sequential clinical
trial design with provision for selecting a subgroup where the
treatment may be effective.  The package provides facilities for
design, exploration and analysis of such trials. A complete
implementation of the initial design for the DEFUSE-3 trial is
also provided as a vignette.

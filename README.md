<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- use knitr::knit(input="README.Rmd", output = "README.md") to generate README.md -->




[![Build Status](https://travis-ci.org/cells2numbers/migrationminer.svg?branch=master)](https://travis-ci.org/cells2numbers/migrationminer) 
[![codecov](https://codecov.io/gh/cells2numbers/migrationminer/branch/master/graph/badge.svg)](https://codecov.io/gh/cells2numbers/migrationminer)

# migrationminer
migrationminer is a R package to analyse and profile in vitro cell tracking and migration data. It is belongs to the [cytominer-verse](https://github.com/cytomining/) used for morphological profiling and allows to create temporal or dynamic profiles. 

It works well together with CellProiler and the CellProfiler tracking module. Beside that it is quite easy to parse other formats as long as the tracking data is available in tidy format. That means, that the data is available as csv where each row represents a cell with at least four columns

 * x coordinate  
 * y coordinate 
 * time point 
 * track label

## Installation 
migrationminer is not yet available from CRAN but you can install directly from github. First, you need to install the devtools package. You can do this from CRAN. Invoke R and then type

```{}
install.packages("devtools")
```

Load the package using 

```{}
library(devtools)
```
Now install migrationminer from github using
```{}
install_github("cells2numbers/migrationminer")
```

You should see something like 
```
Downloading GitHub repo cells2numbers/migrationminer@master
from URL https://api.github.com/repos/cells2numbers/migrationminer/zipball/master
Installing migrationminer
'/usr/local/Cellar/r/3.4.1_2/lib/R/bin/R' --no-site-file --no-environ --no-save --no-restore --quiet CMD INSTALL  \
  '/private/var/folders/qq/v5mwx8g15655gt708s86vys8yzcpll/T/RtmpUyrzzK/devtoolse7f945c7b360/cells2numbers-migrationminer-8a95c5d'  \
  --library='/usr/local/lib/R/3.4/site-library' --install-tests 

* installing *source* package ‘migrationminer’ ...
** R
** tests
** preparing package for lazy loading
** help
*** installing help indices
** building package indices
** testing if installed package can be loaded
* DONE (migrationminer)
Reloading installed migrationminer
```

If this does not work for you or you find any bugs, please file an issue! 

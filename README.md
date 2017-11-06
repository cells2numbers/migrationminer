[![Build Status](https://travis-ci.org/cells2numbers/neutrominer.svg?branch=master)](https://travis-ci.org/cells2numbers/neutrominer) 
[![Coverage Status](https://img.shields.io/codecov/c/github/cells2numbers/neutrominer/master.svg)](https://codecov.io/github/cells2numbers/neutrominer?branch=master)

# neutrominer
Neutrominer is a R package to analyse and profile in vitro cell tracking and migration data. It is belongs to the [cytominer-verse](https://github.com/cytomining/) used for morphological profiling and allows to create temporal or dynamic profiles. 

# Installation 
Neutrominer is not yet available from CRAN but you can install directly from github. First, you need to install the devtools package. You can do this from CRAN. Invoke R and then type

```
install.packages("devtools")
```

Load the package using 

```
library(devtools)
```
Now install neutrominer from github using
```
install_github("cells2numbers/neutrominer")
```

You should see something like 
```
Downloading GitHub repo cells2numbers/neutrominer@master
from URL https://api.github.com/repos/cells2numbers/neutrominer/zipball/master
Installing neutrominer
'/usr/local/Cellar/r/3.4.1_2/lib/R/bin/R' --no-site-file --no-environ --no-save --no-restore --quiet CMD INSTALL  \
  '/private/var/folders/qq/v5mwx8g15655gt708s86vys8yzcpll/T/RtmpUyrzzK/devtoolse7f945c7b360/cells2numbers-neutrominer-8a95c5d'  \
  --library='/usr/local/lib/R/3.4/site-library' --install-tests 

* installing *source* package ‘neutrominer’ ...
** R
** tests
** preparing package for lazy loading
** help
*** installing help indices
** building package indices
** testing if installed package can be loaded
* DONE (neutrominer)
Reloading installed neutrominer
```

If this does not work for you or you find any issue, please file an issue! 

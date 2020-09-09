# VisualVisium
R functions for visualizing Visium data. 

Installation: simply clone the repo and run from R (notice `/path/to/repo/` points to the parent directory, not the `VisualVisium` directory itself): 

```
library('devtools')
setwd("/path/to/repo/")
install('VisualVisium')
```

This helper package contains two functions: 

`cntmtr_tpl` loads the count matrix and the tissue positions list -- probably the two most important data objects produced by the Visium array -- into R objects. 

`visium_bands` groups Visium spots into bands that align with the perimeter of the tissue. 

Head over to the [demo](VisualVisium_demo/VisualVisium_demo.md) for usage examples.

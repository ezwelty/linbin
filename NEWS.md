# linbin 0.1.3

* Declare `rmarkdown` as a soft dependency (see https://github.com/yihui/knitr/issues/1864).
* Update external links to resolve redirects or to follow the CRAN canonical form.

# linbin 0.1.2

* Fixes a bug that caused results from `sample_events()` to be returned out of order in certain cases when `by = ` was used in a sampling function.

# linbin 0.1.1

* Functions from core packages are now called explicitly, i.e. `package::function()`, with the exception of `graphics` functions which are all imported to the namespace.
* Package citation and documentation points to the journal article [Multiscale analysis of river networks using the R package linbin](https://doi.org/10.1080/02755947.2015.1044764).

# linbin 0.1.0

* The linbin R package is released into the world! Further changes will be tracked here.

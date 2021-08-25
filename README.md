# DeeDee 

The package DeeDee provides a user-friendly framework for the visual comparison of two or more results from Differential Expression Analysis workflows. 

It contains eight exported functions: one for the data preparation, six for the actual graphical output-producing comparison functions, and the last one opens an interactive Shiny App. The App comprises the functionality of all functions, combining them with an easy to use graphical user interface.

### Installation 
```
library("remotes")
remotes::install_github("lea-rothoerl/DeeDee", 
                        dependencies = TRUE, 
                        build_vignettes = TRUE)
```

### Example

Save your DEA results as .RDS files or convert them DeeDee tables by running `deedee_prepare()` and save those to .RDS, .txt or .xlsx files. Open the DeeDee App by running 
```
deedee_app()
```
and use the interactive graphical user interface to analyze your data.


### License

MIT &copy; Lea Rothörl

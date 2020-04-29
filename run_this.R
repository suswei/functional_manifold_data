knitr::knit("index.Rmd")
rmarkdown::pandoc_convert("index.md", to = "latex", output = "./writing/index.tex",citeproc=FALSE)

# if index.Rmd is only changed in terms of writing, set cache=True to avoid long running times producing certain figures
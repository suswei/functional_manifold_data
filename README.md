# About

The tex files are written through index.Rmd. 

Any changes to the writing should be performed in index.Rmd. Afterwards, do "run_this.R" to produce a corresponding index.tex.

The wrapper latex files call index.tex.

Warning: bad idea to direclty modify index.tex as it gets overwritten when index.Rmd is knit.
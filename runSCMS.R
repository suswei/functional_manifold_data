library(reticulate)
scms = import_from_path("scms",path='.')
py_run_file("demo.py")

blah <- import("demo",convert=TRUE)
scmspath = blah$denoised
plot(scmspath[,1],scmspath[,2])

obj = EuclideanExamples("bananas",100,0)
data = obj$data
plot(data)
denoised = scms$scms(data,0.2)
points(denoised,col="red")

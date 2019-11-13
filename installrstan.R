data<-data.frame(a = rnorm(100,3,2), b = rnorm(100,1,2))
testmodel<-brm(a~b, data = data)

remove.packages("rstan")
if (file.exists(".RData"))(file.remove(".RData"))

#restart r

Sys.setenv(MAKEFLAGS = "-j4") # use four cores
remotes::install_github("stan-dev/rstan", ref = "develop", subdir = "rstan/rstan", build_opts = "")

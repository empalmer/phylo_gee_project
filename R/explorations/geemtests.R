install.packages("geeM")
library(geeM)

data("seizure", package="geepack")
seiz.l <- reshape(seizure,
                  varying=list(c("base","y1", "y2", "y3", "y4")),
                  v.names="y", times=0:4, direction="long")
seiz.l <- seiz.l[order(seiz.l$id, seiz.l$time),]
seiz.l$t <- ifelse(seiz.l$time == 0, 8, 2)
seiz.l$x <- ifelse(seiz.l$time == 0, 0, 1)

seiz <- geem(y~ x + trt + x:trt+ offset(log(t)), id=id,data = seiz.l, 
             family = poisson, corstr = "exchangeable")
seiz

mat <- matrix(c(1,0,2,3,0,
                0,1,3,3,0,
                2,3,1,2,0,
                3,3,2,1,4,
                0,0,0,4,1), nrow = 5)
seiz <- geem(y~ x + trt + x:trt+ offset(log(t)), id=id,data = seiz.l, 
             family = poisson, corstr = "userdefined", corr.mat = mat)
seiz

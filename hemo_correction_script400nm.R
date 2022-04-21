args = commandArgs(trailingOnly=TRUE)

data_file = args[1]
param_file = args[2]
output_file = args[3]
path = args[4]
#Rscript --vanilla hemo_correction.R  run18Residue.xlsx parameters.xlsx Output_fit.txt ~/data/R_test/fMRI_CION

setwd(path)
#install.packages("readxl")
library("readxl")

ehbo = 331462
ehbr = 261958
x405 = 0.0158 ###Harry's latest excitaiton X405 value  

n = 65

data <- read_xlsx(data_file, col_names = FALSE)
params <- read_xlsx(param_file, col_names = TRUE)

names(data) <- as.character(params$wavelength)

####### set baseline time points
base_time_points <- 10
time_points <- length(data[[1]]) - base_time_points

base_data <- apply(data[1:base_time_points, ], 2, mean)

data <- data[-(1:base_time_points),]


c1 <- c()
c2 <- c()
c3 <- c()
c4 <- c()
for(t in 1:length(data[[1]])){
  x = ehbo*x405+params$`HbO(e)`*params$X
  y = ehbr*x405+params$`HbR(e)`*params$X
  signal = data
  base = base_data
  z = as.numeric(log(signal[t,]/base, base=exp(1)))
  x2 <- sum(x**2)
  y2 <- sum(y**2)
  xy <- t(x)%*%y
  xz <- t(x)%*%z
  yz <- t(y)%*%z
  xI <- sum(x)
  yI <- sum(y)
  zI <- sum(z)
  A <- matrix(c(x2, xy, -xI, xy, y2, -yI, -xI, -yI, n**2), nrow = 3, byrow = TRUE)
  b <- c(-xz, -yz, zI )
  c <- solve(A, b)
  c1 <- c(c1, c[1])
  c2 <- c(c2, c[2])
  c3 <- c(c3, c[3])
  c4 <- c(c4, sum((z+c[1]*x+c[2]*y-c[3]*rep(1, n))**2))
}


data_fit = data.frame(c1,c2,c3,c4)
names(data_fit) = c("dChbo", "dChbr", "lnC(t)/C(0)", "goodness")
write.table(data_fit, file = output_file, row.names=FALSE)



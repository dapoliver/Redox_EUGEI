library(pmsampsize)

# CHR-T v HC
pmsampsize(type = "b", cstatistic = 0.938, parameters = 2, prevalence = 0.46) # miRNA
pmsampsize(type = "b", cstatistic = 0.943, parameters = 10, prevalence = 0.46) # miRNA+demo

# CHR-T v CHR-NT
pmsampsize(type = "b", cstatistic = 0.988, parameters = 5, prevalence = 0.19) # miRNA
pmsampsize(type = "b", cstatistic = 0.986, parameters = 7, prevalence = 0.19) # miRNA+demo

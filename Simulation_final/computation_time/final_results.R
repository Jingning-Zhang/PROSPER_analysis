#
#################################################
### runtime
#
#time_prscsx2 <- numeric()
#for (repl in 1:10){
#  time_prscsx2[repl] <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/computation_time/prscsx/results/",repl, "/two_runtime.rds"))/60
#}
#time_prscsx5 <- numeric()
#for (repl in 1:10){
#  time_prscsx5[repl] <- readRDS(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/computation_time/prscsx/results/",repl, "/five_runtime.rds"))/60
#}
#
#library(stringr)
#
#time_prsepr2 <- numeric()
#for (repl in 1:10){
#  tmp <- readLines(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/computation_time/prsepr_2eth.o1216861.",repl))
#  start_time <- strptime(str_split(str_split(tmp[1], "2022 ")[[1]][2], " ")[[1]][1], "%H:%M:%S")
#  end_time <- strptime(str_split(str_split(tmp[5], "2022 ")[[1]][2], " ")[[1]][1], "%H:%M:%S")
#  time_prsepr2[repl] <- as.numeric(end_time-start_time, units="mins")
#}
#
#time_prsepr5 <- numeric()
#for (repl in 1:10){
#  tmp <- readLines(paste0("/dcs04/nilanjan/data/jzhang2/MEPRS/codes/Simulation_final/computation_time/prsepr_5eth.o1216855.",repl))
#  start_time <- strptime(str_split(str_split(tmp[1], "2022 ")[[1]][2], " ")[[1]][1], "%H:%M:%S")
#  end_time <- strptime(str_split(str_split(tmp[5], "2022 ")[[1]][2], " ")[[1]][1], "%H:%M:%S")
#  time_prsepr5[repl] <- as.numeric(end_time-start_time, units="mins")
#}
#
#res_time <- data.frame(time_prsepr2,time_prsepr5,time_prscsx2,time_prscsx5)
#apply(res_time,2,mean)


################################################

#prscsx_2eth: qacct -j 1262471
#prscsx_5eth: qacct -j 1262470
#prsepr_2eth: qacct -j 1262473
#prsepr_5eth: qacct -j 1262472


#1262470 0.59846 prscsx_5et jzhang2      qw    11/30/2022 00:22:36                                    1 1-10:1
#1262471 0.50895 prscsx_2et jzhang2      qw    11/30/2022 00:22:51                                    1 1-10:1
#1262472 0.50469 prsepr_5et jzhang2      qw    11/30/2022 00:23:09                                    1 1-10:1
#1262473 0.50318 prsepr_2et jzhang2      qw    11/30/2022 00:23:12                                    1 1-10:1
#

### time (secs)
time_prscsx2 <- c(6461,6462,6461,6482,6529,6580,6654,6661,6872,7511)
time_prscsx5 <- c(33172,34521,34634,34983,36645,36649,36678,36713,36737,36771)
time_prsepr2 <- c(180,169,161,164,187,188,190,190,171,177)
time_prsepr5 <- c(500,547,407,424,391,394,358,364,345,353)

### memory (GB)
mem_prscsx2 <- c(795.797,715.168,824.406,824.414,815.406,824.414,815.395,807.781,665.145,665.363)/1000
mem_prscsx5 <- c(824.438,841.367,841.371,839.484,839.500,836.965,841.363,824.551,836.957,841.363)/1000
mem_prsepr2 <- c(2.227,2.243,2.244,2.241,2.218,2.220,2.284,2.237,2.239,2.247)
mem_prsepr5 <- c(2.363,2.330,2.338,2.338,2.338,2.344,2.341,2.338,2.367,2.369)

res_time <- data.frame(time_prscsx2,time_prscsx5,time_prsepr2,time_prsepr5)
apply(res_time,2,mean)/60
res_mem <- data.frame(mem_prscsx2,mem_prscsx5,mem_prsepr2,mem_prsepr5)
apply(res_mem,2,mean)

#all_chrom_prsepr_5eth.sh: qacct -j 1260601
#ru_wallclock 2597s
#maxvmem      33.811GB

#all_chrom_prsepr_2eth.sh: qacct -j 1260623
#ru_wallclock 1190s
#maxvmem      24.500GB


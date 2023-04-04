
Note that the format of codes in the original analysis here and the final package (https://github.com/Jingning-Zhang/PROSPER) 
is slightly different. The latter one is much well-annotated and cleaned so that people can use it easily.

1. The directory *lassosum2* is used to perform single-ancestry analysis, which is also the first step in PROSPER.

2. The directory *multi-lassosum* is used to clean the data into a standard format which can then be inputted to PROSPER.

3. The directory *multi-lassosum_train_from_single_eth* is used to perform PROSPER. Some optimal tuning parameter values (as 
described in https://github.com/Jingning-Zhang/PROSPER) are borrowed from the single-ancestry analysis (lassosum2), so the 
directory name has a subscript of "train_from_single_eth" but it is actually what the PROSPER did.

4. The directory *computation_time* is used to compare the computational time for PRS-CSx and PROSPER using simulated phenotype
on the chromosome 22. 


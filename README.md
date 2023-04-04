
# Codes for data analysis in PROSPER paper

This is a repo containing codes used for simulation and data analyses in the PROSPER paper. For the command line tool of PROSPER, please go to  https://github.com/Jingning-Zhang/PROSPER.


**PROSPER** is a command line tool based on R programming language. It is a powerful multi-ancestry PRS method that jointly models GWAS summary statistics from multiple populations by penalized regression and ensemble approach to improve the performance in minority populations.

Preprint manuscript will be put online very soon. Please contact Jingning Zhang (jzhan218@jhu.edu) for citation.


## Directory description

1. The directory of *Simulation_final* contains the codes used to perform simulation analysis in this manuscript.
2. The directory of *data_preparation* contains the codes used to clean real data, including UK Biobank, GLGC, and All of Us. The data cleaning and processing of 23andMe data has been described in a previous paper from [Zhang et al.](https://www.biorxiv.org/content/10.1101/2022.03.24.485519v5.abstract)
3. The directory of *real_data_analysis* contains the codes used to perform real data analysis in this manuscript. The training data are from 23andMe, GLGC, and All of Us. The tuning and testing data are from 23andMe (for 23andMe training) and UK Biobank (for GLGC and All of Us training). 
4. The directory of *results_in_paper* contains the codes used to plot figures and create numerical results. 
5. The directory of *pacakge* contains the codes used to generate the final command line tool. This includes cleaning of scripts, creating example input data, running the software, and finally producing example output results. The final command line tool is provided in https://github.com/Jingning-Zhang/PROSPER.

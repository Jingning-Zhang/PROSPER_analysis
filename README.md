
# Codes for data analysis in PROSPER paper

This is a repo containing codes used for simulation and data analyses in the PROSPER paper. 

For the command line tool of PROSPER, please go to https://github.com/Jingning-Zhang/PROSPER.

**PROSPER** is a command line tool based on R programming language. It is a powerful multi-ancestry PRS method that jointly models GWAS summary statistics from multiple populations by penalized regression and ensemble approach to improve the performance in minority populations.

Preprint manuscript is released at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10055041/. The peer-reviewed publication will be available soon. Please contact Jingning Zhang (jingningzhang238@gmail.com) if you have any questions regarding the citation.


## Directory Description

1. The directory of *Simulation_final* contains the codes used to perform simulation analysis in this manuscript.
2. The directory of *data_preparation* contains the codes used to clean real data, including UK Biobank, 23andMe, GLGC, and All of Us. 
3. The directory of *real_data_analysis* contains the codes used to perform real data analysis in this manuscript. The training data are from 23andMe, GLGC, and All of Us. The tuning and testing data are from 23andMe (for 23andMe training) and UK Biobank (for GLGC and All of Us training). 
4. The directory of *results_in_paper* contains the codes used to plot figures and create numerical results. 
5. The directory of *pacakge* contains the codes used to generate the final command line tool. This includes cleaning of scripts, creating example input data, running the software, and finally producing example output results. The final command line tool is provided in https://github.com/Jingning-Zhang/PROSPER.
6. The directory of *revision1* and *revision2* contains the codes for the revisions of our manuscript for publication.

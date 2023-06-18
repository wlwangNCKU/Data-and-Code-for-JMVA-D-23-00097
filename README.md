# Data-and-Code-for-JMVA-D-23-00097
Supplementary Material for "On moments of truncated multivariate normal/independent distributions" by Tsung-I Lin and Wan-Lun Wang

######## Author responsible for the code ########
For questions, comments or remarks about the code please contact responsible author, Wan-Lun Wang (wangwl@gs.ncku.edu.tw).

######## Configurations ########
The code was written/evaluated in R with the following software versions:
R version 4.1.1 (2021-08-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19044)

Matrix products: default

locale:
[1] LC_COLLATE=Chinese (Traditional)_Taiwan.950  LC_CTYPE=Chinese (Traditional)_Taiwan.950    LC_MONETARY=Chinese (Traditional)_Taiwan.950
[4] LC_NUMERIC=C                                 LC_TIME=Chinese (Traditional)_Taiwan.950    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
[1] compiler_4.1.1

######## Descriptions of the codes ######### 
Please extract the file "Data-and-Code.zip" to the "current working directory" of the R package.
The getwd() function shall determine an absolute pathname of the "current working directory".

Before running the codes 'fig1.r', 'fig2.r', 'simulation.r', and 'example.r', one needs to install the following R packages 

    install.packages("mvtnorm")  Version: 1.1-3
    install.packages("tmvtnorm") Version: 1.5
    install.packages("cubature") Version： 2.0.4.5
    install.packages("Bessel")   Version： 0.6-0

R codes for the implementation of our methodology are provided.

### Subfolder: ./function ###
./function
	contains the program (function) of
                 (1) 'mni.r' for calculating the probability density function (pdf) and cumulative distribution function (cdf) of multivariate normal/independent (MNI) distributions;
                 (2) 'tmni.r' for calculating the pdf and cdf of the truncated MNI distributions; and
		 (3) 'TMNImoment.r' for evaluting the first two moments of truncated MNI distributions.

### Subfolder: ./code ###
./code
       contains (1) 'fig1.r' main script for drawing the density curves of normal/independent and standard normal distributions;
                (2) 'fig2.r' main script for drawing the contours plots of bivariate normal/independent and bivariate normal distributions;

# Note for Section 2 - Preliminary:
To draw the density curves and contour lines as shown in Figures 1 and 2, please source the 'mni.r', 'tmni.r' and 'TMNImoment.r' scripts in subfolder "./function/",
and then run the 'fig1.r' and 'fig2.r' scripts in subfolder "./code/".
The results have been stored in "./results/" and the densities of the bivariate normal/independen distributions has been stored in './data/dMNI.Rdata'.

                (3) 'simulation.r' main script for re-generating part of intermediate results for simualtion (note: The cases of right-, left- and doubly-truncation should be done separately.);
                (4) 'SIMtab2tabS1.r' main script for re-producting the theoretical mean vectors and variance-covariance matrices of the MNI distributions with various truncations;
                (5) 'SIMfig.r' main script for re-producting the trace plots of empirical mean vectors and variance-covariance parameters for five TMNI distributions;

# Note for Section 5 - Simulation:
R code 'simulation.r' generates the intermediate results of Table 2 and Figures 3-4 in the manuscript and Figures S.1-S.4 in the supplementary materials.
Because the code takes a huge amount of time to run, we record these intermediate results so that one can use 
the R codes 'SIMtab2tabS1.r' and 'SIMfig.r' to obtain the final results based on files stored in "./results/right", "./results/left", and "./results/doubly" subfolders.

To reproduce the results presented in Table 2, run the script 'SIMtab2tabS1.r' inf the subfolder "./code/".
To reproduce the results presented in Figures 3-4 and Figures S.1-S.4, just load 'mMNI.RData' files in the "./data/" and the '.txt' files in the subfolders "./results/doubly", "./results/right/", and "./results/left/", 
and then run the script 'SIMfig.r' in the subfolder "./code/". 

The resulting Table 2 and Table S.1 has been stored in the "./results/" subfolder, and 
Figures 3-4, Figures S.1-S.2, and Figures S.3-S.4 have been stored in the "./results/doubly", "./results/right" and "./results/left" subfolders, respectively.

               (6) 'example.r' main script for re-producting the multivariate tail conditional expectation (MTCE) for the MNI distributions.
# Note for Section 6 - Application:
To draw the MTCE measures evaluated at specific quantiles, please source the 'mni.r', 'tmni.r' and 'TMNImoment.r' scripts in subfolder "./function/",
and then run the 'example.r' scripts in subfolder "./code/".

### Subfolder: ./data ###
./data
      contains
      (1) dMNI.RData: the densities of the bivariate normal/independen distributions with true values of parameters specified in Section 2;
      (2) mMNI.RData: the theoretical moments of the truncated multivariate normal/indpendent distributions with true values of parameters specified in Section 5;
      (3) example.Rdata: the numerical results of the examples presented in Section 6.

### Subfolder: ./results ###
./results
      contains 
        (1) fig1.eps: pdfs for normal, Student's t, slash, contaminated noraml, variance gamma and double exponential distributions under various settings of truncations;
        (2) fig2.eps: contour plots of bivaraite normal, Student's t, slash, contaminated normal, variance gamma and double exponential under various settings of truncations;
        (3) fig5.eps: MTCE allocated to the individual risk in 5-variate MVN and five MNI distributions evaluated at two specified quantities.
        (4) Table2.csv: theoretical moments of five cases of MNI distributions with double truncation;
        (5) TableS1.csv: theoretical moments of five cases of MNI distributions with right and left truncations.

./results/doubly; ./results/right; ./results/left
      contain intermediately numerical results (empirical sample mean and sample variances and covariances) 
      from the simulation studies under three settings of truncations over different sample sizes.
      
### Additional Remark ###
# Note: 
One can directly run each "source(.)" described in 'master.r' file in the seperate R session to obtain the results.

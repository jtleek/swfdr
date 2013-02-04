
README for swfdr folder
-------------------------

Copyright (C) 2011 Jeffrey T. Leek (http://www.biostat.jhsph.edu/~jleek/contact.html) and Leah R. Jager (jager@usna.edu)

Note: These functions were written on a Mac and may have difficulties when
         read on Windows machines. 


### getPvalues.R


This file contains the code to scrape the P-values from pubmed (either run it first, or use the already
calculated pvalueData.rda)



### calculateSwfdr.R


This file contains the function to estimate the science-wise false discovery rate


### journalAnalysis.R


This file contains the code to reproduce the quantities and figures in the Jager/Leek paper


### journalAnalysisHelp.R


This file contains a helper function necessary for journalAnalysis.R


### pvalueData.rda


The pre-computed p-value data used for the Jager/Leek paper in .rda format. 


### simulation.R


A simulation study comparing our estimates to the truth when the assumptions hold and when they are badly violated. 


## To reproduce the results in the Jager/Leek paper


Either: (1) Run getPvalues.R to obtain the p-value data, (2) Run journalAnalysis.R
Or: (2) Run journalAnalysis.R using the pre-computed pvalueData.rda

To get the simulated results you should run sensitivity.R (supplementary sensitivity analysis) and simulation.R (main text sensivity analysis)

Note, because of the bootstrapping calculations, these functions may take a while (think order hours) to run. 








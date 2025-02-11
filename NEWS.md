# nephro 1.5
* Date: 2025-02-10 
* Authors: Cristian Pattaro <cristian.pattaro@eurac.edu>, Janina Herold <Janina.Herold@klinik.uni-regensburg.de>

* New functions
  + Urinary albumin-to-creatinine and protein-to-creatinine ratios: UACR() and UPCR(), respectively
  + Urea-to-creatinine and BUN-to-creatinine ratios: UCR() and BCR(), respectively
  + Convertion functions Convert.mgdL.umolL() and Convert.mgdL.mmolL() for measurement unit conversion

* Updates
  + function name CKDEpi_RF.creat() -> changed into CKDEpi2021.creat()
  + function name CKDEpi_RF.creat.cys() -> changed into CKDEpi2021.creat.cys()

* Documentation fixes
  + In the manual, measurement unit for Cystatin C in the crea-cys CKDEpi has been corrected; many thanks to Isaac V. Cohen for spotting it!


# nephro 1.4
* Date: 2023-10-10 
* Authors: Janina Herold <Janina.Herold@klinik.uni-regensburg.de>, Thomas Winkler, Cristian.Pattaro <cristian.pattaro@eurac.edu>

* New functions
  + European Kidney Foundation Consortium (EKFC) formulas
  + FAS equations 
  + Bedside Schwartz equation (courtesy of Andrew Srisuwananukorn)

*Updates
  + function name CKDEpi.creat.rf() -> changed into CKDEpi_RF.creat()
  + function name CKDEpi.creat.cys.rf() -> changed into CKDEpi_RF.creat.cys()


# nephro 1.3
* Date: 2022-05-01 
* Author: Ryosuke Fujii <ryosuke.fujii@eurac.edu>

* New functions
  + CKDEpi.creat.rf() function: race-free CKD-EPI equation for serum creatinine
  + CKDEpi.creat.cys.rf() function: race-free CKD-EPI equation for serum creatinine and cystatin C
 
* Documentation fixes
  + In the manual, measurement unit for Cystatin C has been corrected


# nephro 1.2
* Date: 2017-05-05 
* Author: Cristian Pattaro <cristian.pattaro@eurac.edu>

* Bug fixes
  + Serious bug fixed in the CG(): previous wrong formula: eGFR <- (140-age) * wt / 72 * creatinine; current fixed formula: 
eGFR <- (140-age) * wt / (72 * creatinine); previous CG() results not valid - many thanks to Joao Sabino for uncovering the mistake!


# nephro 1.1
* Date: 2015-01-31 
* Author: Cristian Pattaro <cristian.pattaro@eurac.edu>

* Bug fixes  
  + all functions: a serious bug on NA values has been fixed - many thanks to Max Plisckhe (Medical University of Vienna) who uncovered that!
  + CG(): new Cockroft and Gault creatinine clearance estimation function added


# nephro 1.0
* Date: 2013-01-30 
* Author: Cristian Pattaro <cristian.pattaro@eurac.edu>

* initial version released to CRAN
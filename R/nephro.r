MDRD4 <- function(creatinine, sex, age, ethnicity, method="IDMS")
   { 
   if (!is.null(creatinine) & !is.null(sex) & !is.null(age) & !is.null(ethnicity))
      {
      creatinine <- as.numeric(creatinine)
      ethnicity <- as.numeric(ethnicity)      
      sex <- as.numeric(sex)
      age <- as.numeric(age)
      n <- length(creatinine)
      
      if (length(sex) == n & length(age) == n & length(ethnicity) == n)
         {
         eGFR <- creatinine^(-1.154)*age^(-0.203)
         eGFR[sex==0] <- eGFR[sex==0] * 0.742
         eGFR[ethnicity==1] <- eGFR[ethnicity==1] * 1.212
         eGFR * ifelse(method == 'IDMS', 175, 186)
         } else
         stop ("Different number of observations between variables") 
      } else
      stop ("Some variables are not defined") 
   }

MDRD6 <- function(creatinine, sex, age, albumin, BUN, ethnicity, method="IDMS")
   { 
   if (!is.null(creatinine) & !is.null(sex) & !is.null(age) & !is.null(albumin) & !is.null(BUN) & !is.null(ethnicity))
      {
      creatinine <- as.numeric(creatinine)
      ethnicity <- as.numeric(ethnicity)      
      albumin <- as.numeric(albumin)
      BUN <- as.numeric(BUN)
      sex <- as.numeric(sex)
      age <- as.numeric(age)

      n <- length(creatinine)
      
      if (length(sex) == n & length(age) == n & length(ethnicity) == n & length(albumin) == n & length(BUN) == n )
         {
         eGFR <- creatinine^(-0.999)*age^(-0.176)*BUN^(-0.17)*albumin^0.318
         eGFR[sex==0] <- eGFR[sex==0] * 0.762
         eGFR[ethnicity==1] <- eGFR[ethnicity==1] * 1.180
         eGFR * ifelse(method == 'IDMS', 161.5, 170)
         } else
         stop ("Different number of observations between variables") 
      } else
      stop ("Some variables are not defined") 
   }

Virga <- function(creatinine, sex, age, wt)
   {
   if (!is.null(creatinine) & !is.null(sex) & !is.null(age) & !is.null(wt))
      {
      creatinine <- as.numeric(creatinine)
      sex <- as.numeric(sex)
      age <- as.numeric(age)
      wt <- as.numeric(wt)
      n <- length(creatinine)
      
      if (length(sex) == n & length(age) == n & length(wt) == n)
         {
         eGFR <- numeric(n)
         eGFR[sex==1] <- (69.4 - 0.59*age[sex==1] + 0.79*wt[sex==1]) / creatinine[sex==1] - 3.0
         eGFR[sex==0] <- (57.3 - 0.37*age[sex==0] + 0.51*wt[sex==0]) / creatinine[sex==0] - 2.9
         eGFR
         } else
         stop ("Different number of observations between variables") 
      } else
      stop ("Some variables are not defined") 
   }

CKDEpi.creat <- function(creatinine, sex, age, ethnicity)
   { 
   if (!is.null(creatinine) & !is.null(sex) & !is.null(age) & !is.null(ethnicity))
      {
      creatinine <- as.numeric(creatinine)
      ethnicity <- as.numeric(ethnicity)      
      sex <- as.numeric(sex)
      age <- as.numeric(age)

      n <- length(creatinine)
      
      if (length(sex) == n & length(age) == n & length(ethnicity) == n)
         {
         k <- a <- numeric(n)
         k[sex==0] <- 0.7
         k[sex==1] <- 0.9
         k[is.na(sex)] <- NA
         a[sex==0] <- -0.329
         a[sex==1] <- -0.411
         a[is.na(sex)] <- NA
         one <- rep(1,n)
         eGFR <- apply(cbind(creatinine/k,one),1,min,na.rm=T)^a * apply(cbind(creatinine/k,one),1,max,na.rm=T)^-1.209 * 0.993^age
         eGFR[sex==0] <- eGFR[sex==0] * 1.018
         eGFR[ethnicity==1] <- eGFR[ethnicity==1] * 1.159
         141 * eGFR
         } else
         stop ("Different number of observations between variables") 
      } else
      stop ("Some variables are not defined") 
   }

CKDEpi.cys <- function(cystatin, sex, age)
   { 
   if (!is.null(cystatin) & !is.null(sex) & !is.null(age))
      {
      cystatin <- as.numeric(cystatin)
      sex <- as.numeric(sex)
      age <- as.numeric(age)

      n <- length(cystatin)
      
      if (length(sex) == n & length(age) == n)
         {
         k_sex <- rep(1,n)
         k_sex[sex==0] <- 0.932
         k_sex[is.na(sex)] <- NA
         one <- rep(1,n)
         133 * apply(cbind(cystatin/0.8,one),1,min,na.rm=T)^-0.499 * apply(cbind(cystatin/0.8,one),1,max,na.rm=T)^-1.328 * 0.996^age * k_sex
         } else
         stop ("Different number of observations between variables") 
      } else
      stop ("Some variables are not defined") 
   }

CKDEpi.creat.cys <- function(creatinine, cystatin, sex, age, ethnicity)
   {
   if (!is.null(cystatin) & !is.null(cystatin) & !is.null(sex) & !is.null(age) & !is.null(ethnicity))
      {
      creatinine <- as.numeric(creatinine)
      cystatin <- as.numeric(cystatin)
      sex <- as.numeric(sex)
      age <- as.numeric(age)
      ethnicity <- as.numeric(ethnicity)
      n <- length(creatinine)
      
      if (length(sex) == n & length(age) == n & length(ethnicity) == n & length(cystatin) == n)
         {
         k_sex <- k_ethnicity <- rep(1,n)
         k_sex[sex==0] <- 0.969
         k_sex[is.na(sex)] <- NA
         k_ethnicity[ethnicity==1] <- 1.08
         k_ethnicity[is.na(ethnicity)] <- NA
         k <- rep(0.9,n)
         k[sex==0] <- 0.7
         k[is.na(sex)] <- NA
         a <- rep(-0.207,n)
         a[sex==0] <- -0.248
         a[is.na(sex)] <- NA     
         one <- rep(1,n)
         CR <- cbind(creatinine/k,one)
         CY <- cbind(cystatin/0.8,one)                 
         135 * apply(CR,1,min,na.rm=T)^a * apply(CR,1,max,na.rm=T)^-0.601 * apply(CY,1,min,na.rm=T)^-0.375 * apply(CY,1,max,na.rm=T)^(-0.711) * 0.995^age * k_sex * k_ethnicity
         } else
         stop ("Different number of observations between variables") 
      } else
      stop ("Some variables are not defined") 
   }   
  
Stevens.cys1 <- function(cystatin)
   { 
   if (!is.null(cystatin))
      {
      cystatin <- as.numeric(cystatin)
      76.7*cystatin^(-1.19)
      } else
      stop ("Variable not defined") 
   }

Stevens.cys2 <- function(cystatin, sex, age, ethnicity)
   { 
   if (!is.null(cystatin) & !is.null(sex) & !is.null(age) & !is.null(ethnicity))
      {
      cystatin <- as.numeric(cystatin)
      ethnicity <- as.numeric(ethnicity)
      sex <- as.numeric(sex)
      age <- as.numeric(age)
      n <- length(cystatin)
      if (length(sex) == n & length(age) == n & length(ethnicity) == n)
         {
         eGFR <- cystatin^(-1.17) * age^(-0.13)
         eGFR[sex == 0] <- eGFR[sex == 0] * 0.91
         eGFR[ethnicity == 1] <- eGFR[ethnicity == 1] * 1.06
         127.7 *eGFR
         } else
         stop ("Different number of observations between variables") 
      } else
      stop ("Some variables are not defined") 
   }

Stevens.creat.cys <- function(creatinine, cystatin, sex, age, ethnicity)
   { 
   if (!is.null(cystatin) & !is.null(sex) & !is.null(age) & !is.null(ethnicity))
      {
      creatinine <- as.numeric(creatinine)
      cystatin <- as.numeric(cystatin)
      ethnicity <- as.numeric(ethnicity)
      sex <- as.numeric(sex)
      age <- as.numeric(age)
      n <- length(cystatin)
      if (length(sex) == n & length(age) == n & length(ethnicity) == n & length(cystatin) == n)
         {
         eGFR <- creatinine^(-0.65) * cystatin^(-0.57) * age^(-0.20)
         eGFR[sex == 0] <- eGFR[sex == 0] * 0.82
         eGFR[ethnicity == 1] <- eGFR[ethnicity == 1] * 1.11
         177.6 * eGFR    
         } else
         stop ("Different number of observations between variables") 
      } else
      stop ("Some variables are not defined") 
   }

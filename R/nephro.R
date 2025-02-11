# =============================================================================
#
# TITLE : nephro - functions to deal with nephrological outcomes
# AUTHOR: Cristian Pattaro, Ryosuke Fujii, Janina Herold
# DATE  : Jan 30, 2025

# =============================================================================

MDRD4 <- function(creatinine, sex, age, ethnicity, method="IDMS")
{ 
  if (!is.null(creatinine) & !is.null(sex) & !is.null(age) & !is.null(ethnicity))
  {
    if (is.null(method)) method <- "other"
    
    creatinine <- as.numeric(creatinine)
    ethnicity <- as.numeric(ethnicity)      
    sex <- as.numeric(sex)
    age <- as.numeric(age)
    n <- length(creatinine)
    
    if (length(sex) == n & length(age) == n & length(ethnicity) == n)
    {
      # Identify missing data and store the index
      idx <- c(1:n)[is.na(creatinine) | is.na(sex) | is.na(age) | is.na(ethnicity)]
      
      # Replace missing data with fake data to avoid problems with formulas
      creatinine[is.na(creatinine)] <- 10
      sex[is.na(sex)] <- 10
      age[is.na(age)] <- 10
      ethnicity[is.na(ethnicity)] <- 10
      
      # MDRD4 formula 
      eGFR <- creatinine^(-1.154)*age^(-0.203)
      eGFR[sex==0] <- eGFR[sex==0] * 0.742
      eGFR[ethnicity==1] <- eGFR[ethnicity==1] * 1.212
      
      # Restore missing data at the indexed positions
      eGFR[idx] <- NA
      
      # Output
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
    if (is.null(method)) method <- "other"
    
    creatinine <- as.numeric(creatinine)
    ethnicity <- as.numeric(ethnicity)      
    albumin <- as.numeric(albumin)
    BUN <- as.numeric(BUN)
    sex <- as.numeric(sex)
    age <- as.numeric(age)
    
    n <- length(creatinine)
    
    if (length(sex) == n & length(age) == n & length(ethnicity) == n & length(albumin) == n & length(BUN) == n )
    {
      # Identify missing data and store the index
      idx <- c(1:n)[is.na(creatinine) | is.na(sex) | is.na(age) | is.na(ethnicity) | is.na(BUN) | is.na(albumin)]
      
      # Replace missing data with fake data to avoid problems with formulas
      creatinine[is.na(creatinine)] <- 10
      sex[is.na(sex)] <- 10
      age[is.na(age)] <- 10
      ethnicity[is.na(ethnicity)] <- 10
      albumin[is.na(albumin)] <- 10
      BUN[is.na(BUN)] <- 10
      
      # MDRD6 formula 
      eGFR <- creatinine^(-0.999)*age^(-0.176)*BUN^(-0.17)*albumin^0.318
      eGFR[sex==0] <- eGFR[sex==0] * 0.762
      eGFR[ethnicity==1] <- eGFR[ethnicity==1] * 1.180
      
      # Restore missing data at the indexed positions
      eGFR[idx] <- NA
      
      # Output
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
      # Identify missing data and store the index
      idx <- c(1:n)[is.na(creatinine) | is.na(sex) | is.na(age) | is.na(wt)]
      
      # Replace missing data with fake data to avoid problems with formulas
      creatinine[is.na(creatinine)] <- 10
      sex[is.na(sex)] <- 10
      age[is.na(age)] <- 10
      wt[is.na(wt)] <- 10
      
      # Virgas' formula         
      eGFR <- numeric(n)
      eGFR[sex==1] <- (69.4 - 0.59*age[sex==1] + 0.79*wt[sex==1]) / creatinine[sex==1] - 3.0
      eGFR[sex==0] <- (57.3 - 0.37*age[sex==0] + 0.51*wt[sex==0]) / creatinine[sex==0] - 2.9
      
      # Restore missing data at the indexed positions
      eGFR[idx] <- NA
      
      # Output
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
      # Identify missing data and store the index
      idx <- c(1:n)[is.na(creatinine) | is.na(sex) | is.na(age) | is.na(ethnicity)]
      
      # Replace missing data with fake data to avoid problems with formulas
      creatinine[is.na(creatinine)] <- 10
      sex[is.na(sex)] <- 10
      age[is.na(age)] <- 10
      ethnicity[is.na(ethnicity)] <- 10
      
      # CKD-Epi equation
      k <- a <- numeric(n)
      k[sex==0] <- 0.7
      k[sex==1] <- 0.9
      a[sex==0] <- -0.329
      a[sex==1] <- -0.411
      one <- rep(1,n)
      eGFR <- apply(cbind(creatinine/k,one),1,min,na.rm=T)^a * apply(cbind(creatinine/k,one),1,max,na.rm=T)^-1.209 * 0.993^age
      eGFR[sex==0] <- eGFR[sex==0] * 1.018
      eGFR[ethnicity==1] <- eGFR[ethnicity==1] * 1.159
      
      # Restore missing data at the indexed positions
      eGFR[idx] <- NA
      
      # Output
      141 * eGFR
      
    } else
      stop ("Different number of observations between variables") 
  } else
    stop ("Some variables are not defined") 
}

CKDEpi2021.creat <- function (creatinine, sex, age) 
{
  if (!is.null(creatinine) & !is.null(sex) & !is.null(age))
  {
    creatinine <- as.numeric(creatinine)
    sex <- as.numeric(sex)
    age <- as.numeric(age)
    n <- length(creatinine)
    
    if (length(sex) == n & length(age) == n)
    {
      # Identify missing data and store the index
      idx <- c(1:n)[is.na(creatinine) | is.na(sex) | is.na(age)]
      
      # Replace missing data with fake data to avoid problems with formulas
      creatinine[is.na(creatinine)] <- 10
      sex[is.na(sex)] <- 10
      age[is.na(age)] <- 10
      
      # CKD-Epi equation
      k <- a <- numeric(n)
      k[sex == 0] <- 0.7
      k[sex == 1] <- 0.9
      a[sex == 0] <- -0.241
      a[sex == 1] <- -0.302
      one <- rep(1, n)
      eGFR <- apply(cbind(creatinine/k, one), 1, min, na.rm = T)^a * apply(cbind(creatinine/k, one), 1, max, na.rm = T)^-1.200 * 0.9938^age
      eGFR[sex == 0] <- eGFR[sex == 0] * 1.012
      
      # Restore missing data at the indexed positions
      eGFR[idx] <- NA
      
      # Output
      142 * eGFR
      
    } else
      stop("Different number of observations between variables")
  } else 
    stop("Some variables are not defined")
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
      # Identify missing data and store the index
      idx <- c(1:n)[is.na(cystatin) | is.na(sex) | is.na(age)]
      
      # Replace missing data with fake data to avoid problems with formulas
      cystatin[is.na(cystatin)] <- 10
      sex[is.na(sex)] <- 10
      age[is.na(age)] <- 10
      
      # CKD-Epi equation
      k_sex <- rep(1,n)
      k_sex[sex==0] <- 0.932
      one <- rep(1,n)
      eGFR <- 133 * apply(cbind(cystatin/0.8,one),1,min,na.rm=T)^-0.499 * apply(cbind(cystatin/0.8,one),1,max,na.rm=T)^-1.328 * 0.996^age * k_sex
      
      # Restore missing data at the indexed positions
      eGFR[idx] <- NA
      
      # Output
      eGFR
      
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
      # Identify missing data and store the index
      idx <- c(1:n)[is.na(creatinine) | is.na(cystatin) | is.na(sex) | is.na(age) | is.na(ethnicity)]
      
      # Replace missing data with fake data to avoid problems with formulas
      creatinine[is.na(creatinine)] <- 10
      cystatin[is.na(cystatin)] <- 10
      sex[is.na(sex)] <- 10
      age[is.na(age)] <- 10
      ethnicity[is.na(ethnicity)] <- 10
      
      # CKD-Epi equation
      k_sex <- k_ethnicity <- rep(1,n)
      k_sex[sex==0] <- 0.969
      k_ethnicity[ethnicity==1] <- 1.08
      k <- rep(0.9,n)
      k[sex==0] <- 0.7
      a <- rep(-0.207,n)
      a[sex==0] <- -0.248
      one <- rep(1,n)
      CR <- cbind(creatinine/k,one)
      CY <- cbind(cystatin/0.8,one)                 
      eGFR <- 135 * apply(CR,1,min,na.rm=T)^a * apply(CR,1,max,na.rm=T)^-0.601 * apply(CY,1,min,na.rm=T)^-0.375 * apply(CY,1,max,na.rm=T)^(-0.711) * 0.995^age * k_sex * k_ethnicity
      
      # Restore missing data at the indexed positions
      eGFR[idx] <- NA
      
      # Output
      eGFR
      
    } else
      stop ("Different number of observations between variables") 
  } else
    stop ("Some variables are not defined") 
}   

CKDEpi2021.creat.cys <- function(creatinine, cystatin, sex, age)
{
  if (!is.null(cystatin) & !is.null(cystatin) & !is.null(sex) & !is.null(age))
  {
    creatinine <- as.numeric(creatinine)
    cystatin <- as.numeric(cystatin)
    sex <- as.numeric(sex)
    age <- as.numeric(age)
    n <- length(creatinine)
    
    if (length(sex) == n & length(age) == n & length(creatinine) == n & length(cystatin) == n)
    {
      # Identify missing data and store the index
      idx <- c(1:n)[is.na(creatinine) | is.na(cystatin) | is.na(sex) | is.na(age)]
      
      # Replace missing data with fake data to avoid problems with formulas
      creatinine[is.na(creatinine)] <- 10
      cystatin[is.na(cystatin)] <- 10
      sex[is.na(sex)] <- 10
      age[is.na(age)] <- 10
      
      # CKD-Epi equation
      a <- b <- c <- numeric(n)
      a[sex == 0] <- 0.7
      a[sex == 1] <- 0.9
      b[sex == 0] <- -0.219
      b[sex == 1] <- -0.144
      c <- 0.8
      one <- rep(1, n)               
      eGFR <- apply(cbind(creatinine/a, one), 1, min, na.rm = T)^b * apply(cbind(creatinine/a, one), 1, max, na.rm = T)^-0.544 * apply(cbind(cystatin/c, one), 1, min, na.rm = T)^-0.323 * apply(cbind(cystatin/c, one), 1, max, na.rm = T)^-0.778 * 0.9961^age
      eGFR[sex == 0] <- eGFR[sex == 0] * 0.963
      
      # Restore missing data at the indexed positions
      eGFR[idx] <- NA
      
      # Output
      135 * eGFR
      
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
      # Identify missing data and store the index
      idx <- c(1:n)[is.na(cystatin) | is.na(sex) | is.na(age) | is.na(ethnicity)]
      
      # Replace missing data with fake data to avoid problems with formulas
      cystatin[is.na(cystatin)] <- 10
      sex[is.na(sex)] <- 10
      age[is.na(age)] <- 10
      ethnicity[is.na(ethnicity)] <- 10
      
      # Stevens' formula
      eGFR <- cystatin^(-1.17) * age^(-0.13)
      eGFR[sex == 0] <- eGFR[sex == 0] * 0.91
      eGFR[ethnicity == 1] <- eGFR[ethnicity == 1] * 1.06
      
      # Restore missing data at the indexed positions
      eGFR[idx] <- NA
      
      # Output
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
    n <- length(creatinine)
    
    if (length(sex) == n & length(age) == n & length(ethnicity) == n & length(cystatin) == n)
    {
      # Identify missing data and store the index
      idx <- c(1:n)[is.na(creatinine) | is.na(cystatin) | is.na(sex) | is.na(age) | is.na(ethnicity)]
      
      # Replace missing data with fake data to avoid problems with formulas
      creatinine[is.na(creatinine)] <- 10
      cystatin[is.na(cystatin)] <- 10
      sex[is.na(sex)] <- 10
      age[is.na(age)] <- 10
      ethnicity[is.na(ethnicity)] <- 10
      
      # Stevens' formula
      eGFR <- creatinine^(-0.65) * cystatin^(-0.57) * age^(-0.20)
      eGFR[sex == 0] <- eGFR[sex == 0] * 0.82
      eGFR[ethnicity == 1] <- eGFR[ethnicity == 1] * 1.11
      
      # Restore missing data at the indexed positions
      eGFR[idx] <- NA
      
      # Output
      177.6 * eGFR
      
    } else
      stop ("Different number of observations between variables") 
  } else
    stop ("Some variables are not defined")
}

CG <- function(creatinine, sex, age, wt)
{ 
  if (!is.null(creatinine) & !is.null(sex) & !is.null(age) & !is.null(wt))
  {
    creatinine <- as.numeric(creatinine)
    wt  <- as.numeric(wt)      
    sex <- as.numeric(sex)
    age <- as.numeric(age)
    n <- length(creatinine)
    
    if (length(sex) == n & length(age) == n & length(wt) == n)
    {
      # Identify missing data and store the index
      idx <- c(1:n)[is.na(creatinine) | is.na(sex) | is.na(age) | is.na(wt)]
      
      # Replace missing data with fake data to avoid problems with formulas
      creatinine[is.na(creatinine)] <- 10
      sex[is.na(sex)] <- 10
      age[is.na(age)] <- 10
      wt[is.na(wt)] <- 10
      
      # Cockroft-Gault formula 
      eGFR <- (140-age) * wt / (72 * creatinine)
      eGFR[sex==0] <- eGFR[sex==0] * 0.85
      
      # Restore missing data at the indexed positions
      eGFR[idx] <- NA
      
      # Output
      eGFR
      
    } else
      stop ("Different number of observations between variables") 
  } else
    stop ("Some variables are not defined") 
}

FAS.creat <- function(creatinine, sex, age)
{ 
  if (!is.null(creatinine) & !is.null(sex) & !is.null(age))
  {
    creatinine <- as.numeric(creatinine)     
    sex <- as.numeric(sex)
    age <- as.numeric(age)
    n <- length(creatinine)
    
    if (length(sex) == n & length(age) == n)
    {
      # Identify missing data and store the index
      idx <- c(1:n)[is.na(creatinine) | is.na(sex) | is.na(age)]
      
      # Replace missing data with fake data to avoid problems with formulas
      creatinine[is.na(creatinine)] <- 10
      sex[is.na(sex)] <- 10
      age[is.na(age)] <- 10
      
      
      Qcrea <- a <- age_term<-numeric(n)
      Qcrea[sex==0] <- 0.7
      Qcrea[sex==1] <- 0.9
      
      age_term[age>40]  <- 0.988^(age[age>40] - 40)
      age_term[age<=40] <- 1
      a<--1
      
      eGFR <- cbind(creatinine/Qcrea)^a * age_term
      eGFR[idx] <- NA
      
      # Output
      107.3 * eGFR
      
    } else
      stop ("Different number of observations between variables") 
  } else
    stop ("Some variables are not defined")
}

FAS.cys <- function(cystatin, sex, age)
{ 
  if (!is.null(cystatin) & !is.null(sex) & !is.null(age))
  {
    cystatin <- as.numeric(cystatin)     
    sex <- as.numeric(sex)
    age <- as.numeric(age)
    
    n <- length(cystatin)
    
    if (length(sex) == n & length(age) == n)
    {
      # Identify missing data and store the index
      idx <- c(1:n)[is.na(cystatin) | is.na(sex) | is.na(age)]
      
      # Replace missing data with fake data to avoid problems with formulas
      cystatin[is.na(cystatin)] <- 10
      sex[is.na(sex)] <- 10
      age[is.na(age)] <- 10
      
      a <- Qcys <- age_term <- numeric(n)
      
      Qcys[age<=70 & sex==0] <- 0.82
      Qcys[age<=70 & sex==1] <- 0.82
      Qcys[age>70  & sex==0] <- 0.863 + 0.01704*(age[age>70 & sex==0] - 70)
      Qcys[age>70  & sex==1] <- 0.863 + 0.01704*(age[age>70 & sex==1] - 70)
      age_term[age>40]  <- 0.988^(age[age>40] - 40)
      age_term[age<=40] <- 1
      a <- -1
      
      eGFR <- cbind(cystatin/Qcys)^a * age_term
      eGFR[idx] <- NA
      
      # Output
      107.3 * eGFR
      
    } else
      stop ("Different number of observations between variables") 
  } else
    stop ("Some variables are not defined") 
}

FAS.creat.cys <- function(creatinine, cystatin, sex, age)
{
  if (!is.null(creatinine) & !is.null(cystatin) & !is.null(sex) & !is.null(age))
  {
    creatinine <- as.numeric(creatinine)
    cystatin <- as.numeric(cystatin)
    sex <- as.numeric(sex)
    age <- as.numeric(age)
    n <- length(creatinine)
    
    if (length(sex) == n & length(age) == n & length(creatinine) == n & length(cystatin) == n)
    {
      # Identify missing data and store the index
      idx <- c(1:n)[is.na(cystatin) | is.na(sex) | is.na(age)|is.na(creatinine)]
      
      # Replace missing data with fake data to avoid problems with formulas
      creatinine[is.na(creatinine)] <- 10
      cystatin[is.na(cystatin)] <- 10
      sex[is.na(sex)] <- 10
      age[is.na(age)] <- 10
      
      Qcrea <- Qcys <- a <- age_term <- b <- numeric(n)
      
      Qcys[age<=70 & sex==0] <- 0.82
      Qcys[age<=70 & sex==1] <- 0.82
      Qcys[age>70  & sex==0] <- 0.863 + 0.01704*(age[age>70 & sex==0] - 70)
      Qcys[age>70  & sex==1] <- 0.863 + 0.01704*(age[age>70 & sex==1] - 70)
      Qcrea[sex==0] <- 0.7
      Qcrea[sex==1] <- 0.9
      
      age_term[age>40]  <- 0.988^(age[age>40] - 40)
      age_term[age<=40] <- 1
      a <- 0.5
      b <- -1
      
      eGFR <- cbind(a*(creatinine/Qcrea) + (1-a)*(cystatin/Qcys))^b * age_term
      eGFR[idx] <- NA
      
      #Output
      107.3 * eGFR
      
    } else
      stop ("Different number of observations between variables") 
  } else
    stop ("Some variables are not defined") 
}

EKFC.creat <- function(creatinine, sex, age)
{ 
  if (!is.null(creatinine) & !is.null(sex) & !is.null(age))
  {
    creatinine <- as.numeric(creatinine)     
    sex <- as.numeric(sex)
    age <- as.numeric(age)
    
    n <- length(creatinine)
    
    if (length(sex) == n & length(age) == n)
    {
      # Identify missing data and store the index
      idx <- c(1:n)[is.na(creatinine) | is.na(sex) | is.na(age)]
      
      # Replace missing data with fake data to avoid problems with formulas
      creatinine[is.na(creatinine)] <- 10
      sex[is.na(sex)] <- 10
      age[is.na(age)] <- 10
      
      Qcrea <- a <- age_term <- Qcrea_mgDl <- numeric(n)
      
      Qcrea_mgDl[age>25 & sex==0] <- 0.7
      Qcrea_mgDl[age>25 & sex==1] <- 0.9
      Qcrea[age>2 & age<=25 & sex==0] <- exp(3.080 + 0.177*age[age>2 & age<=25 & sex==0] - 0.223*log(age[age>2 & age<=25 & sex==0]) - 0.00596*age[age>2 & age<=25 & sex==0]^2 + 0.0000686*age[age>2 & age<=25 & sex==0]^3)
      Qcrea[age>2 & age<=25 & sex==1] <- exp(3.2   + 0.259*age[age>2 & age<=25 & sex==1] - 0.543*log(age[age>2 & age<=25 & sex==1]) - 0.00763*age[age>2 & age<=25 & sex==1]^2 + 0.0000790*age[age>2 & age<=25 & sex==1]^3)
      
      #SCr and Q in umol/L (to convert to mg/dL, divide by 88.4
      Qcrea_mgDl[age>2 & age<=25 & sex==0] <- Qcrea[age>2 & age<=25 & sex==0]/88.4
      Qcrea_mgDl[age>2 & age<=25 & sex==1] <- Qcrea[age>2 & age<=25 & sex==1]/88.4
      
      CreaQ <- creatinine/Qcrea_mgDl
      a[CreaQ<1]  <- -0.322
      a[CreaQ>=1] <- -1.132
      age_term[age>40] <- 0.990^(age[age>40] - 40)
      age_term[age<=40] <- 1
      
      eGFR <- cbind(creatinine/Qcrea_mgDl)^a*age_term
      eGFR[idx] <- NA
      
      #Output
      107.3 * eGFR
      
    } else
      stop ("Different number of observations between variables") 
  } else
    stop ("Some variables are not defined") 
}

EKFC.cys <- function(cystatin, sex, age)
{ 
  if (!is.null(cystatin) & !is.null(sex) & !is.null(age))
  {
    cystatin <- as.numeric(cystatin)     
    sex <- as.numeric(sex)
    age <- as.numeric(age)
    
    n <- length(cystatin)
    
    if (length(sex) == n & length(age) == n)
    {
      # Identify missing data and store the index
      idx <- c(1:n)[is.na(cystatin) | is.na(sex) | is.na(age)]
      
      # Replace missing data with fake data to avoid problems with formulas
      cystatin[is.na(cystatin)] <- 10
      sex[is.na(sex)] <- 10
      age[is.na(age)] <- 10
      
      age_term <- Q <- a <- numeric(n)
      Q[age>50 & sex==0] <- 0.79 + 0.005*(age[age>50 & sex==0] - 50)
      Q[age>50 & sex==1] <- 0.86 + 0.005*(age[age>50 & sex==1] - 50)
      Q[age>18 & age<=50 & sex==0] <- 0.79
      Q[age>18 & age<=50 & sex==1] <- 0.86
      Q[age<=18 & sex==1] <- 0.84101056 + 0.02936824*age[age<=18 & sex==1] - 0.01008673*age[age<=18 & sex==1]^2 + 0.00111224*age[age<=18 & sex==1]^3 + 0.00003552*age[age<=18 & sex==1]^4
      Q[age<=18 & sex==0] <- 0.96723605 + 0.07852349*age[age<=18 & sex==0] - 0.01562922*age[age<=18 & sex==0]^2 + 0.00108505*age[age<=18 & sex==0]^3 + 0.00002382*age[age<=18 & sex==0]^4
      CysQ <- cystatin/Q
      a[CysQ<1]  <- -0.322
      a[CysQ>=1] <- -1.132
      age_term[age>40]  <- 0.990^(age[age>40] - 40)
      age_term[age<=40] <- 1
      
      #for age > 40 years  
      eGFR <- cbind(cystatin/Q)^a * age_term
      
      eGFR[idx] <- NA
      107.3 * eGFR
      
    } else
      stop ("Different number of observations between variables") 
  } else
    stop ("Some variables are not defined") 
}

EKFC_SF.cys <- function(cystatin,age)
{ 
  if (!is.null(cystatin) & !is.null(age))
  {
    cystatin <- as.numeric(cystatin)     
    age <- as.numeric(age)
    n <- length(cystatin)
    
    if (length(age) == n)
    {
      # Identify missing data and store the index
      idx <- c(1:n)[is.na(cystatin) |  is.na(age)]
      
      # Replace missing data with fake data to avoid problems with formulas
      cystatin[is.na(cystatin)] <- 10
      age[is.na(age)] <- 10
      
      age_term<-numeric(n)
      Q <- a <- numeric(n)
      
      Q[age>50]           <- 0.83       + 0.005*(age[age>50] - 50)
      Q[age>18 & age<=50] <- 0.83
      Q[age<=18]          <- 0.96723605 + 0.07852349*age[age<=18] - 0.01562922*age[age<=18]^2 + 0.00108505*age[age<=18]^3 + 0.00002382*age[age<=18]^4
      
      CysQ <- cystatin/Q
      a[CysQ<1]  <- -0.322
      a[CysQ>=1] <- -1.132
      
      age_term[age>40]  <- 0.990^(age[age>40] - 40)
      age_term[age<=40] <- 1
      
      #for age > 40 years  
      eGFR <- cbind(cystatin/Q)^a * age_term
      eGFR[idx] <- NA
      
      #Outpout
      107.3 * eGFR
      
    } else
      stop ("Different number of observations between variables") 
  } else
    stop ("Some variables are not defined") 
}

Schwartz.Bedside <- function(creatinine, ht, age)
{
  if (!is.null(creatinine) & !is.null(ht) & !is.null(age))
  {
    creatinine <- as.numeric(creatinine)     
    ht <- as.numeric(ht)
    age <- as.numeric(age)
    n <- length(creatinine)
    
    if (length(ht) == n)
    {
      if ( (min(age)<1) | (max(age)>17)) cat("\nWarning: there are age values <1 or >17 years; for those children, eGFR values might be invalid\n") 
      
      # Identify missing data and store the index
      idx <- c(1:n)[is.na(creatinine) |  is.na(ht)]
      
      # Replace missing data with fake data to avoid problems with formulas
      creatinine[is.na(creatinine)] <- 10
      ht[is.na(ht)] <- 10
      
      # Bedside Schwartz formula
      eGFR <- 0.413 * ht / creatinine
      
      # Restore missing data at the indexed positions
      eGFR[idx] <- NA
      
      # Output
      eGFR
      
    } else
      stop ("Different number of observations between variables") 
  } else
    stop ("Some variables are not defined") 
}

UCR <- function(urea=NULL, creatinine=NULL)
{ 
  # urea-to-creatinine ratio in SI units
  # urea in mmol/L
  # creatinine in umol/L
  if (!is.null(creatinine) & !is.null(urea))
  {      
    creatinine <- as.numeric(creatinine)
    urea <- as.numeric(urea)
    n <- length(creatinine)
    
    if (length(urea) == n)
    {
      creatinine[creatinine==0] <- NA
      return (urea/creatinine/1000)
    } else
      stop ("Different number of observations between variables") 
  } else
    stop ("Some variables are not defined") 
}

BCR <- function(BUN=NULL, creatinine=NULL)
{ 
  # BUN-to-creatinine ratio in non-SI units
  # BUN in mg/dL
  # creatinine in mg/dL
  if (!is.null(creatinine) & !is.null(BUN))
  {      
    creatinine <- as.numeric(creatinine)
    BUN <- as.numeric(BUN)
    n <- length(creatinine)
    
    if (length(BUN) == n)
    {
      creatinine[creatinine==0] <- NA
      return (BUN/creatinine)
    } else
      stop ("Different number of observations between variables") 
  } else
    stop ("Some variables are not defined") 
}

UACR <- function(u.albumin=NULL, u.creatinine=NULL)
{ 
  # albumin in mg/dL
  # creatinine in mg/dL
  if (!is.null(u.creatinine) & !is.null(u.albumin))
  {      
    u.creatinine <- as.numeric(u.creatinine)/1000
    u.albumin <- as.numeric(u.albumin)
    n <- length(u.creatinine)
    
    if (length(u.albumin) == n)
    {
      u.creatinine[u.creatinine==0] <- NA
      return (u.albumin/u.creatinine)
    } else
      stop ("Different number of observations between variables") 
  } else
    stop ("Some variables are not defined") 
}

UPCR <- function(u.protein=NULL, u.creatinine=NULL)
{ 
  # protein in mg/dL
  # creatinine in mg/dL
  if (!is.null(u.creatinine) & !is.null(u.protein))
  {      
    u.creatinine <- as.numeric(u.creatinine)
    u.protein <- as.numeric(u.protein)
    n <- length(u.creatinine)
    
    if (length(u.protein) == n)
    {
      u.creatinine[u.creatinine==0] <- NA
      return (u.protein/u.creatinine)
    } else
      stop ("Different number of observations between variables") 
  } else
    stop ("Some variables are not defined") 
}

Convert.mgdL.umolL <- function(mg.dL = NULL, umol.L = NULL)
{
  if (!is.null(mg.dL) & is.null(umol.L))
  {
    umol.L <- mg.dL * 88.42
    cat("\n mg/dL converted to umol/L (conversion factor: 88.42)\n")
    return(umol.L)
  } else {
    if (is.null(mg.dL) & !is.null(umol.L))
    {
      mg.dL <- umol.L / 88.42
      cat("\n umol/L converted to mg/dL (conversion factor: 88.42)\n")
      return(mg.dL)
    } else stop ("No conversion was performed because either no measures to convert were provided or conversion unnecessary because both measures were provided.")    
  }
}

Convert.mgdL.mmolL <- function(mg.dL = NULL, mmol.L = NULL)
{
  if (!is.null(mg.dL) & is.null(mmol.L))
  {
    mmol.L <- mg.dL * 0.357
    cat("\n mg/dL converted to mmol/L (conversion factor: 0.357)\n")
    return(mmol.L)
  } else {
    if (is.null(mg.dL) & !is.null(mmol.L)) 
    {
      mg.dL <- mmol.L / 0.357
      cat("\n mmol/L converted to mg/dL (conversion factor: 0.357)\n")
      return(mg.dL)
    } else stop ("No conversion was performed because either no measures to convert were provided or conversion unnecessary because both measures were provided.")    
  }
}


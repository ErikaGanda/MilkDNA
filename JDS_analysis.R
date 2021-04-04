library(tidyverse)
library(multcomp)
library(car)

set.seed(50113)


# load data ---------------------------------------------------------------


AllData <- read.table("Inputs/CleanDNAprepData1.18.19.txt", sep="\t",  
                      fill = TRUE, header=TRUE)




# First summarize missing data patterns -----------------------------------

dat <- AllData %>% 
  
  # Select variables of interest
  dplyr::select(Assay, VariableKit, VariableSampleType, SpikeSet, LogCopiespermLofMilk) %>% 
  
  # Just take Innoculated Milk 
  filter(VariableSampleType == "InoculatedMilk") %>% 
  
  # Drop EZNA Food DNA because that seems to be the major problem with 
  # "heteroskedasticity"
  # Just show it but don't analyze. 
  filter(VariableKit != "EZFood") %>% 
  
  # split by Assay into list for parallel analysis 
  split(.$Assay) 




# Now summarise missing data patters 
dat %>% 
  
  # Create new missing binary variable
  map(~mutate(.x, missing = is.na(LogCopiespermLofMilk))) %>% 
  
  # Summarise by VariableKit and SpikeSet 
  map(~.x %>% 
        group_by(VariableKit, SpikeSet) %>% 
        summarise(n = n(), 
                  missing = sum(missing),
                  prop_missing = missing/n))

# Well there is only one missing value and its for CORDENA on the Second SpikeSet, 
# So just say that and then be done with it. 

# Although based on my findings in the next section I am guessing this is not totally correct... 
# Why no PViralDNA for SpikeSet 1 in The Listeria Assay?

# Linear model ------------------------------------------------------------

# Just take complete cases
dat <- dat %>% map(~filter(.x, complete.cases(.x))) %>% 
  
  # Just a bit of data munging to convert VariableKit to factor
  map(~mutate(.x, VariableKit = factor(VariableKit)))
    
  
# Now linear model with interactions
lms <- dat %>% 
    map(~lm(LogCopiespermLofMilk ~ VariableKit*SpikeSet - 1,
            data=.x))


# Identify worst offending outliers (those with cooks distance greater than 0.5)
outliers <- lms %>% 
  map(cooks.distance) %>% 
  map(~which(.x>0.5))
outliers
# only two found


# helper function
drop_if_non_zero <- function(x,y){
  if (length(y) > 0) x <- x[-y,]
  return(x)
}

lms <- dat %>% 
  
  # Drop outliers
  map2(outliers, ~drop_if_non_zero(.x,.y))  %>% 
  
  # Refit the model
  map(~lm(LogCopiespermLofMilk ~ VariableKit*SpikeSet - 1,
          data=.x))

# Check for identifiabilty 
dat %>% 
  map(~.x[,c("VariableKit", "SpikeSet")] %>% table())

# We are going to have to Drop PviralDNA on hte Listeria Monocytogenes assay. 
comparisons <-  c("Mastitis - COREDNA", 
                  "Pfood - COREDNA", 
                  "PSoilP - COREDNA", 
                  "PviralDNA - COREDNA",
                  "ZymoDNA - COREDNA",
                  "Pfood - Mastitis", 
                  "PSoilP - Mastitis", 
                  "PviralDNA - Mastitis",
                  "ZymoDNA - Mastitis", 
                  "PSoilP - Pfood",
                  "PviralDNA - Pfood",
                  "ZymoDNA - Pfood",    
                  "PviralDNA - PSoilP",
                  "ZymoDNA - PSoilP", 
                  "ZymoDNA - PviralDNA" ) %>% 
  paste(., "== 0")



# Now do pairwise comparisons 
comps <- lms[-3] %>% 
  
  # Use Tukey's all-pair comparison
  map2(comparisons[-3], ~glht(.x, mcp(VariableKit=.y) ))

# Handle case 3 specially due to missingness
compl <- comparisons
tmp <- grep("viral", comparisons)
compl <- compl[-tmp]
rm(tmp)
#hacking multcomp to make it work - this is a bug in their error checking
tmp_vcov <- vcov(lms[[3]])
tmp_vcov[is.na(tmp_vcov)] <- 1e-8 # this value doesn't matter
tmp_coef <- coef(lms[[3]]) 
tmp_coef[is.na(tmp_coef)] <- 1e-8 # this value doesn't matter

comps[[names(dat)[3]]] <- glht(lms[[3]], mcp(VariableKit = compl), vcov.=tmp_vcov, coef.=tmp_coef)


# Now summarise results
comps <- comps %>% 
  map(summary)

# Comps is then the key thing to use and investigate for reports. 
comps

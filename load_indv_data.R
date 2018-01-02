###------------------------------------------------------------------------
# Use this script to load data as individual proxies from .csv files, 
# instead of the complete dataset used in the analysis at once.
###------------------------------------------------------------------------

## Specify lakes and time period of analysis
lakes <- c('BB15','DU15','MA15','NW15','SH15','TL15','UR15') 
studyPeriod <- c(1550,2015)

## Specify directory on your machine containing data (ends with 'nowitna/data') 
dataDir <- "/Users/tylerhoecker/GitHub/nowitna/data/"

## Chronological data
## Age-depth models were built in Bacon v2.2 (http://www.chrono.qub.ac.uk/blaauw/bacon.html), 
## using the parameters described in the text (Table 2).  
dates <- read_csv(paste0(dataDir,'/ageData_sample_inventory.csv'))

ageModel.df <- map(paste0(dataDir,'ageModel_',lakes,'.csv'), read_csv) %>% 
  `names<-` (lakes) %>%
  bind_rows(.id = 'lake')

## Charcoal data
## Processed charcoal data, which have undergone peak analysis and interpolation 
## in CharAnalysis (https://github.com/phiguera/CharAnalysis), using the 
## parameters described in the text (Table 2).
char.ldf <- map(paste0(dataDir,'charResults_',lakes,'.csv'),read_csv)   

char.df <- char.ldf %>%
  bind_rows() %>%
  mutate(ageCE = 1950 - ageTop,
         peakYr = ifelse(peaksFinal == 1,ageCE,NA),
         peakInsigYr = ifelse(peaksInsig == 1,ageCE,NA),
         peakYr1 = ifelse(peaks1 == 1,ageCE,NA),
         peakYr2 = ifelse(peaks2 == 1,ageCE,NA),
         peakYr3 = ifelse(peaks3 == 1,ageCE,NA)) %>%
  filter(ageCE >= studyPeriod[1]) %>%
  group_by(lake) %>%
  mutate(length = max(ageCE) - min(ageCE)) 

## Load raw charcoal count data; derive and standardize charcoal accumulation 
## rates (# cm2 yr-1, "CHAR").
# Function for standardizing non-zero charcoal accumlations rates
transFunc <- function(x) {
  x = ifelse(x > 0, x, NA)
  logX = log(x)
  zX = (logX - mean(logX, na.rm = T)) / sd(logX, na.rm = T)
  expX = exp(zX)
  expX[is.na(expX)] = 0
  return(expX)
}
# Calculate and standardize CHAR
charData <- map(paste0(dataDir,'charData_',lakes,'.csv'),read_csv) %>%
  `names<-` (lakes) %>%
  bind_rows(.id = 'lake') %>% 
  group_by(lake) %>% 
  mutate(sedAcc = (cmTop - cmBot) / (ageTop - ageBot),
         rawChar = (charCount / charVol) * sedAcc,
         char = transFunc(rawChar)) %>% 
  rowwise() %>% 
  mutate(age = round(mean(c(ageTop,ageBot)))) %>% 
  select(char, age, lake)

## Build composite biomass burning recrod from individual records (no peak 
## analysis or interpolation, raw charcoal count data). 
## Execute the method for building a composite biomass burning record used 
## in Kelly et al. 2013 (doi: 10.1073/pnas.1305069110). The method estimates 
## the parameters of a zero-inflated log-normal (ZIL) distribution of pooled 
## charcoal counts in continuous moving windows of a user-defined width. 
## In this analysis, 5-year and 50-year window widths are used 
## (2.5 and 25 half-kernel widths, respectively).

## This portion of the analysis is relegated to a separate R script for 
## clarity of the workflow. Window widths and other parameters can be manipulated within the 'analysis_composite.r' script.
source(file = file.path("analysis_composite.R"))

## Observed (historic) fire event data.
## Using a GIS, we identified fire perimeters from 1950-2014 that overlaped 
## more than 50% of a 1 km buffer around each site. 
## Fire perimeter data from: https://afsmaps.blm.gov/imf/imf.jsp?site=firehistory 
observed.df <- read_csv(paste0(dataDir,'/observedfireData.csv'))
obs.1km <- filter(observed.df,distance == 1.0) %>%
  rename(obsCE = ageCE)

## Tree demography data.
## The bulk of these data are pith date estimates from tree cross sections developed by Paul Duffy (2006, Ph.D. Dissertation, University of Alaska-Fairbanks). A small proportion, those immediatly adjacent to lake cores, were developed by Meghan Foard and Philip Higuera.
tree.df <- read_csv(paste0(dataDir,'treeData.csv'))

tree.df %>%
  mutate(ageCE = plyr::round_any(pith,5)) %>%
  group_by(sp) %>%
  summarise(count = n()) 

## Climate data.
## Modeled temperature data were obtained from SNAP: http://ckan.snap.uaf.edu/dataset/historical-monthly-temperature-1-km-cru-ts. Tree-ring reconstructed temperature were developed by Gregory Wiles et al. (2014) and are available here: .
goa.df <- read_csv(paste0(dataDir,'climateData_GOA.csv'))
cru.df <- read_csv(paste0(dataDir,'climateData_CRU.csv'))
# Extract growing season temperatures from CRU data
cru.gs.df <- cru.df %>%
  filter(month %in% c(4:9)) %>%
  group_by(yearCE,variable) %>%
  summarise(mean = mean(value))

# Standardize and bin climate data to allow for direct comparison.
# Modify tree dataframe 
tree.cor.df <- tree.df %>%
  mutate(ageCE = plyr::round_any(pith,5)) %>%
  group_by(sp, ageCE) %>%
  summarise(count = n()) %>%
  tidyr::spread(key = sp, value = count) %>%
  mutate(bepa = ifelse(is.na(bepa),0,bepa),
         lala = ifelse(is.na(lala),0,lala),
         pigl = ifelse(is.na(pigl),0,pigl),
         pima = ifelse(is.na(pima),0,pima),
         potr = ifelse(is.na(potr),0,potr)) %>%
  mutate(tree.count = rowSums(.[c('bepa','lala','pigl','pima','potr')])) 
  
# Use GOA mean from 1550-2010 to standardize both climate datasets. 
goa.1900 <- goa.df %>%
  filter(yearCE >= 1900) 
goa.mean <- mean(goa.1900$temp)

# Bin data into universal 5-year means
goa.df$yearBins <- cut(goa.df$yearCE,include.lowest = T,right = F,
                       breaks = seq(min(goa.df$yearCE),max(goa.df$yearCE),5), 
                       labels = seq(min(goa.df$yearCE),2005,5))
goa.binned <- goa.df %>%
  group_by(yearBins) %>%
  summarise(bin.temp = mean(temp)) %>%
  mutate(zscore = (bin.temp - goa.mean)) %>%
  mutate(sign = ifelse(zscore >= 0,'positive','negative'))


cru.gs.df$yearBins <- cut(cru.gs.df$yearCE,include.lowest = T,right = F,
                          breaks = seq(1900,2010,5), 
                          labels = seq(1900,2005,5))
cru.binned <- cru.gs.df %>%
  group_by(yearBins,variable) %>%
  summarise(bin.mean = mean(mean)) %>%
  group_by(variable) %>%
  mutate(zscore = (bin.mean - mean(bin.mean))) %>%
  mutate(sign = ifelse(zscore >= 0,'positive','negative'))

# Created combined dataframe of standardized, binned proxies for 1550-1895
goa.1550 <- goa.binned %>%
  mutate(ageCE = as.numeric(as.character(.$yearBins))) %>%
  select(-sign,-yearBins) %>%
  filter(ageCE < 1900)

combined.1550_1895.df <- composite.df %>%
  select(ageCE, lowMean, highMean) %>%
  inner_join(., pctBurned, by = 'ageCE') %>%
  inner_join(., goa.1550, by = 'ageCE') %>%
  left_join(., tree.cor.df, by = 'ageCE') %>%
  select(lowMean, highMean, win.pct, bin.temp, tree.count)

# Created combined dataframe of standardized, binned proxies for 1900-2005
goa.1900 <- goa.binned %>%
  mutate(ageCE = as.numeric(as.character(.$yearBins))) %>%
  select(-sign,-yearBins) 

cru.combined <- cru.binned %>%
  mutate(ageCE = as.numeric(as.character(yearBins))) %>%
  select(ageCE,variable,bin.mean) %>%
  spread(key = variable, value = bin.mean)

combined.1900_2010.df <- composite.df %>%
  inner_join(., pctBurned, by = 'ageCE') %>%
  inner_join(., cru.combined, by = 'ageCE') %>%
  inner_join(., goa.1900, by = 'ageCE') %>%
  left_join(., tree.cor.df, by = 'ageCE') %>%
  select(lowMean, highMean, win.pct, cru.precip = precip, cru.temp = temp, goa.temp = bin.temp, tree.count)




  
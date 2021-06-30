# library(devtools)
# install_github('appling/unitted')
# install_github('USGS-R/streamMetabolizer') 
# install_github('USGS-R/mda.streams') 

library(mda.streams) # usgs package used to build dataset for powell center
library(dataRetrieval)
library(tidyverse)

# site list file from Powell Center data on Sciencebase
pw_sites <- read.table(file = 'C:/Users/Alice Carter/git/oxygen_timeseries/data/site_data.tsv',
           header = T, sep = '\t')

# all sites currently in database:
# list_sites()

# download data from NWIS for each of the sites in the powell center list
whatNWISsites(stateCd = "AZ", 
              parameterCd = c('00300'), 
              siteType = c('ST', 'ST-CA', 'ST-DCH', 'ST-TS', 'SP'),
              hasDataTypeCd = 'iv')

# iterate through each state because USGS doesn't let you query entire database
states <- c(state.abb, 'DC', 'PR')
sites <- data.frame()

for(state in states){
  print(state)
  ss <- try(whatNWISsites(stateCd = state,
                      parameterCd = '00300',
                      siteType = c('ST', 'ST-CA', 'ST-DCH', 'ST-TS', 'SP'),
                      hasDataTypeCd = 'iv'))
  if(inherits(ss, 'try-error')) next
  sites <- bind_rows(sites, ss)

  
}
all_sites <- sites

# for each site, check to see if it has at least 100 days of DO obs

site_info <- whatNWISdata(siteNumber = sites$site_no, 
             service = 'uv', 
             parameterCd = '00300') %>% 
  select(site_no, alt_va, dec_coord_datum_cd, 
                     begin_date, end_date, count_nu) %>%
  filter(count_nu >= 100)

sites <- right_join(sites, site_info, by = "site_no") 
sites <- sites[!duplicated(sites[,1:6]),]


n <- nrow(sites)
skip <- c()
d <- data.frame()

for(i in 1:n){
  print(i)
  dat <- try(readNWISuv(siteNumbers = sites$site_no[i], 
                        parameterCd = c('00010', '00060', '00025', '00300', '00301')))
  # temperature_C, discharge_cfs, baro_mmHg, DO_mgL, DO_persat
  if(inherits(dat, 'try-error')) {
    skip = c(skip, i)
    next}
  d <- bind_rows(d, dat)
  if (i%%10 == 0) {
    saveRDS(d,'C:/Users/Alice Carter/git/oxygen_timeseries/data/nwis_DO_data.rds' )
    print(paste('skipped numbers:', skip))
  }
}

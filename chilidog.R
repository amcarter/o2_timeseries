library(StreamPULSE)

phil_to_mike_format = function(dset, mod_dset, arrange=TRUE){

    nwis_ind = substr(dset$Site_ID, 1, 4) == 'nwis'
    dset$Site_ID[nwis_ind] = gsub('_', '-', dset$Site_ID[nwis_ind])
    dset$site = unname(sapply(dset$Site_ID, function(x) strsplit(x, '_')[[1]][2]))
    dset$site[is.na(dset$site)] = dset$Site_ID[is.na(dset$site)]
    dset$Site_ID = NULL

    #join region column from model outputs (which for streampulse == US state)
    dset = dset %>%
        left_join(mod_dset, 'site') %>%
        mutate(sitecode=paste(region, site, sep='_')) %>%
        filter(! is.na(sitecode)) %>%
        select(-year) %>%
        distinct()

    if(arrange){
        dset = arrange(dset, sitecode)
    }

    return(dset)
}

metr = as_tibble(readRDS('data/site_metrics.rds'))
metr = readRDS('data/metrics_bound.rds') %>%
    select(-one_of('Name', 'Source', 'Lat', 'Lon', 'COMID', 'VPU')) %>%
    as_tibble() %>%
    left_join(metr, by='Site_ID')

# these two years are problematic according to Mike. We're not sure why
dropsiteyears = data.frame(region=c('MD', 'MD', 'AZ', 'AZ'),
                           site=c('BARN', 'POBR', 'OC', 'MV'),
                           sitecode=c('MD_BARN', 'MD_POBR', 'AZ_OC', 'AZ_MV'),
                           year=c('2016', '2016', '2017', '2018'), stringsAsFactors=FALSE)
dropsiteyears_p = apply(dropsiteyears[, c('region', 'site', 'year')], 1, paste, collapse='_')



mods = query_available_results('all')[[1]] %>%
    as_tibble() %>%
    mutate_all(as.character)
modx = apply(mods[, c('region', 'site', 'year')], 1, paste, collapse='_')
mods = filter(mods, ! modx %in% dropsiteyears_p)

metr = phil_to_mike_format(metr, mods) %>%
    select(-site, -region)

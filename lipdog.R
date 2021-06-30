lips_plot = function(readfile, datlist, diagnostics, sitedata, quant_filt=NULL,
    standalone, outfile, filter_label=TRUE, ...){

    #fixed ylims to 5, -5 for standalone; didn't change standalone=F

    #quant_filt example: 'width_calc > 0.25'
    if(standalone){
        filt = readRDS(readfile)
    } else {
        filt = datlist
    }

    var_quant_filt = NULL
    if(! is.null(quant_filt)){
        quant_comp = strsplit(quant_filt, ' ')[[1]]
        qf = quantile(sitedata[, quant_comp[1], drop=TRUE], na.rm=TRUE,
            probs=as.numeric(quant_comp[3]))
        filt_sites = sitedata %>%
            filter_(paste(quant_comp[1], quant_comp[2], qf)) %>%
            pull(sitecode)
        filt = filt[names(filt) %in% filt_sites]

        var_quant_filt = paste0(quant_comp[1], ' ', quant_comp[2], ' ',
            as.numeric(quant_comp[3]) * 100, '%')
    }

    # nsites_included = sum(names(filt) %in%
    #     sitedata$sitecode[! is.na(sitedata$gpp_C_amp)])
    nsites_included = sum(sapply(filt, nrow) > 0)

    # fnz = names(filt)
    # for(i in 1:length(filt)){
    #     filt[[i]]$sitename = fnz[i]
    # }
    #
    # smry = consolidate_list_withnames(filt) %>%
    #     as_tibble() %>%
    #     group_by(sitename, DOY) %>%
    #     summarize_all(list(median=~median(., na.rm=TRUE),
    #         mean=~mean(., na.rm=TRUE),
    #         quant25=~quantile(., na.rm=TRUE)[2],
    #         quant75=~quantile(., na.rm=TRUE)[4]))
    #
    # cam = smry %>%
    #     group_by(sitename) %>%
    #     summarize(cumul_annual_mean=sum(GPP_C_filled_mean, na.rm=TRUE))


    smry = consolidate_list(filt) %>%
        as_tibble() %>%
        group_by(DOY) %>%
        summarize_all(list(median=~median(., na.rm=TRUE),
            quant25=~quantile(., na.rm=TRUE)[2],
            quant75=~quantile(., na.rm=TRUE)[4]))

    if(standalone){
        pdf(file=outfile, width=10, height=10)
        par(mfrow=c(2, 1), oma=c(1, 1, 0, 0))
    }
    par(mar=c(0, 3, 1, 1), lend=2)

    # gpplim = c(0, max(smry$GPP_C_filled_quant75, na.rm=TRUE))
    # erlim = c(min(smry$ER_C_filled_quant25, na.rm=TRUE), 0)
    # gpplim = c(0, max(gpplim[2], abs(erlim[1])))
    # erlim = c(min(abs(gpplim[2]) * -1, erlim[1]), 0)
    gpplim=c(0, 5)
    erlim=c(-5, 0)

    plot(smry$DOY, smry$GPP_C_filled_median, ylab='', yaxs='i', type='l',
        bty='n', lwd=4, xlab='', ylim=gpplim, xaxs='i', xaxt='n', yaxt='n',
        col='gray30')
    polygon(x=c(smry$DOY, rev(smry$DOY)),
        y=c(smry$GPP_C_filled_quant25, rev(smry$GPP_C_filled_quant75)),
        border=NA, col=alpha('forestgreen', alpha=0.6))
    axis(2, las=2, line=0, xpd=NA, tck=-.02, labels=FALSE,
        at=round(seq(0, gpplim[2], length.out=5), 1))
    axis(2, las=2, line=-0.5, xpd=NA, tcl=0, col='transparent',
        at=round(seq(0, gpplim[2], length.out=5), 1))
    abline(h=0, lty=2, lwd=2, col='gray60')
    medsums = round(colSums(select(smry, contains('median'))), 1)

    if(standalone){
        mtext(expression(paste(bold("gC") ~ bold("m") ^ bold("-2") ~
                bold(" d") ^ bold("-1"))), side=2,
            line=-0.5, outer=TRUE)
        # mtext(expression(paste("GPP (gC"~"m"^"-2"~" d"^"-1"*')')), side=2,
        #     line=2.5)
        # legend('topleft', legend=c('Median', '', '25-75%', '', 'NEP Median'),
        #     col=c('forestgreen', 'sienna4', alpha('forestgreen', alpha=0.6),
        #         alpha('sienna', alpha=0.6), 'black'),
        #     bty='n', lty=1, lwd=c(4, 4, 10, 10, 4))
        if(filter_label){
            legend('topright', title='Filters', bty='n', title.col='gray30',
                lty=1, seg.len=0.2, lwd=2, legend=c(..., var_quant_filt))
        }
    }

    legend('right', title='Cumulative\nMedian Sums', bty='n',
        legend=c(paste('GPP:', medsums[1]), paste('ER:', medsums[2]),
            paste('NEP:', medsums[3])), title.col='gray30')
    legend('left', paste('Sites included:', nsites_included), bty='n')

    par(mar=c(3, 3, 0, 1))

    plot(smry$DOY, smry$ER_C_filled_median, ylab='', yaxs='i', type='l',
        bty='n', lwd=4, xlab='', ylim=erlim, xaxs='i', xaxt='n', yaxt='n')
    polygon(x=c(smry$DOY, rev(smry$DOY)),
        y=c(smry$ER_C_filled_quant25, rev(smry$ER_C_filled_quant75)),
        border=NA, col=alpha('sienna', alpha=0.6))
    axis(2, las=2, line=0, xpd=NA, tck=-.02, labels=FALSE,
        at=round(seq(0, erlim[1], length.out=5), 1))
    axis(2, las=2, line=-0.5, xpd=NA, tcl=0, col='transparent',
        at=round(seq(0, erlim[1], length.out=5), 1))
    axis(1, line=0, tck=-.02, labels=FALSE, at=seq(0, max(smry$DOY), 30))
    axis(1, line=-0.5, tcl=0, col='transparent', at=seq(0, max(smry$DOY), 30))
    lines(smry$DOY, smry$NEP_C_filled_median, col='black', lwd=4, xpd=NA, lend=1)

    if(standalone){
        # mtext(expression(paste("ER (gC"~"m"^"-2"~" d"^"-1"*')')), side=2, line=2.5)
        mtext('DOY', side=1, line=2, font=2)
        dev.off()
    }
}


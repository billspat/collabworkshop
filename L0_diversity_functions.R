## TITLE:         NEON Organismal Data: functions for computing diversity measures
## AUTHOR:        Phoebe Zarnetske, Quentin Read, Jonathan Knott
## COLLABORATORS: Sydne Record (Bryn Mawr), Ben Baiser (UFL), Angela Strecker (PSU), 
##                John M. Grady (MSU/Bryn Mawr), Jonathan Belmaker (Tel Aviv U), Mao-Ning Tuanmu (Academia Sinica),
##                Lydia Beaudrot (Rice U), Kate Thibault 
## DATA:          NEON organismal data: all species, all years, all sites
## PROJECT:       "NEON's continental-scale biodiversity"
## DATE:          initiated: July 7, 2020; last updated: July 9, 2020


### This is a copy from https://github.com/NEON-biodiversity/neonbiodiversity/blob/master/code/L0toL1/neon/L0_diversity_functions.R
### copied for demonstrationin a workrshop for the MSB project code collaboration

### This script contains functions for calculating various diversity measures 
### for the NEON taxonomic groups. These are consistently used across all taxa,
### and this script is sourced at the top of any scripts using these functions.


library(lme4)
library(dplyr)
library(purrr)
library(reshape2)
library(lubridate)
library(iNEXT)

# Get rid of anything not listed as target taxa (some taxa do not have this distinction). 
# For example for mammals that would be bycatch that are not nocturnal rodents
# For now just get rid of anything not identified to species
# Use the following rule: if there is no species-level identification within a given genus, keep those individuals because they can all be treated as a single species
# But if there are any individuals in a genus that are identified to species, we have to get rid of all the un-ID'd ones because they could be part of that species
# If anything is identified to a level even coarser than genus, get rid of it

keep_taxa <- function(dat, column = 'taxonID') {
    not_to_sp <- grepl('sp.', dat$scientificName, fixed = TRUE) | grepl('spp.', dat$scientificName, fixed = TRUE) | dat$scientificName == ''
    if (with(dat, exists('specificEpithet'))) {
        not_to_sp[is.na(dat$specificEpithet) | dat$specificEpithet == ''] <- TRUE
    }
    genera_with_species <- unique(dat$genus[!not_to_sp])
    na.omit(unique(dat[!(((dat$genus %in% genera_with_species) & not_to_sp) | dat$taxonRank %in% c('family','order','phylum','class','kingdom','superorder')), column]))
}

# Added by QDR, 19 June 2018
# Modified by QDR, 20 June 2018: Create functions for cumulative richness by year and for reshaping data to plot.

# Chao1 richness estimator
# Takes as input a vector x of species IDs
# Requires observed abundances.
estimator_chao1 <- function(x) {
    xcomm <- table(x)
    S_obs <- length(xcomm) # Number of species observed
    f1 <- sum(xcomm == 1) # Number of singletons
    f2 <- sum(xcomm == 2) # Number of doubletons
    chao1 <- S_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1)) # Calculate chao1 estimator
    var_chao1 <- f2 * ( ((f1/f2)/4)^4 + (f1/f2)^3 + ((f1/f2)/2)^2 ) # Variance of estimator
    if (!is.finite(var_chao1)) var_chao1 <- 0 # If no doubletons, variance is zero
    return(data.frame(chao1 = chao1, 
                      chao1_var = var_chao1,
                      chao1_CImin = max(S_obs, chao1 - 1.96 * sqrt(var_chao1)),
                      chao1_CImax = chao1 + 1.96 * sqrt(var_chao1)))
}

# Asymptotic richness estimator, using iNEXT package.
# Takes as input a vector x of species IDs
# Requires observed abundances.

# For now, just get the asymptotic estimator for richness, and the bounds of its 95% conf int
# Later we can extract even more output from the iNEXT output object.
estimator_asymp <- function(x) {
    require(iNEXT)
    xcomm <- table(x)
    inext_out <- iNEXT(x = list(as.numeric(xcomm)), q = 0, datatype='abundance', nboot = 99) # run iNEXT on the community
    richness_est <- subset(inext_out$AsyEst, Diversity == 'Species richness') # Extract only richness info from output
    return(with(richness_est, data.frame(asymp_est = Estimator, 
                                         asymp_est_stderr = s.e., 
                                         asymp_est_CImin = LCL, 
                                         asymp_est_CImax = UCL)))
}

# Function to get cumulative richness estimators by site and year
# Sequentially add years and see what happens to the cumulative observed richness and richness estimators.
# Each year's result represents all data up to that point in time.

# Modified 26 June: add option to do by site or by plot and to keep either the final observed values or all years
# (group argument and by_year = TRUE or FALSE)
# Modified 29 June: add option for different plot name for the aquatic taxa (uses namedLocation)
richness_cumulative <- function(dat, column = 'taxonID', group = 'site', plot_name = 'plotID', by_year = TRUE) {
    dat <- dat %>%
        rename(sp = !!column, plotID = !!plot_name) %>%
        dplyr::select(siteID, plotID, year, sp)
    if (group == 'site') dat <- dat %>% group_by(siteID)
    if (group == 'plot') dat <- dat %>% group_by(siteID, plotID)
    out <- dat %>%
        do(cbind(year = min(.$year):max(.$year),
                 richness = map_int(min(.$year):max(.$year), function(yr) length(unique(.$sp[.$year <= yr]))),
                 map_dfr(min(.$year):max(.$year), function(yr) estimator_chao1(.$sp[.$year <= yr])),
                 map_dfr(min(.$year):max(.$year), function(yr) estimator_asymp(.$sp[.$year <= yr]))
        ))
    if (!by_year) {
        out %>% filter(year == max(year)) %>% select(-year)
    } else {
        out
    }
}

# Simplified function for phylogenetic and functional diversity, with null models (added 29 June)
# Modified from nasa code, and required for function below
# Only does it based on incidence (no abundance weighting)
pd_fd <- function(sp_list, pddist, fddist, nnull = 99, phylo_spp = NULL, func_spp = NULL) {
    
    require(picante)
    
    # convert species list to a table
    m <- t(as.matrix(table(sp_list)))
    
    # Get rid of species that aren't in phylogenetic and functional diversity.
    if (!is.null(phylo_spp)) mphy <- m[, dimnames(m)[[2]] %in% phylo_spp, drop = FALSE] else mphy <- m
    if (!is.null(func_spp)) mfunc <- m[, dimnames(m)[[2]] %in% func_spp, drop = FALSE] else mfunc <- m
    
    if (dim(mphy)[2] > 1) {
        # Calculate PD and do null models
        MPD <- ses.mpd(mphy, pddist, null.model = 'taxa.labels', abundance.weighted = FALSE, runs = nnull)
        MNTD <- ses.mntd(mphy, pddist, null.model = 'taxa.labels', abundance.weighted = FALSE, runs = nnull)
        # Get all z-scores
        MPD_z <- (MPD$mpd.obs[1] - mean(MPD$mpd.rand.mean, na.rm=TRUE))/sd(MPD$mpd.rand.mean, na.rm=TRUE)
        MNTD_z <- (MNTD$mntd.obs[1] - mean(MNTD$mntd.rand.mean, na.rm=TRUE))/sd(MNTD$mntd.rand.mean, na.rm=TRUE)
        MPD_obs <- MPD$mpd.obs[1]
        MNTD_obs <- MNTD$mntd.obs[1]
    } else {
        MPD_z <- MNTD_z <- MPD_obs <- MNTD_obs <- NA
    }
    if (dim(mfunc)[2] > 1) {
        # Calculate FD and do null models
        MPDfunc <- ses.mpd(mfunc, fddist, null.model = 'taxa.labels', abundance.weighted = FALSE, runs = nnull)
        MNTDfunc <- ses.mntd(mfunc, fddist, null.model = 'taxa.labels', abundance.weighted = FALSE, runs = nnull)
        # Get all z-scores
        MPDfunc_z <- (MPDfunc$mpd.obs[1] - mean(MPDfunc$mpd.rand.mean, na.rm=TRUE))/sd(MPDfunc$mpd.rand.mean, na.rm=TRUE)
        MNTDfunc_z <- (MNTDfunc$mntd.obs[1] - mean(MNTDfunc$mntd.rand.mean, na.rm=TRUE))/sd(MNTDfunc$mntd.rand.mean, na.rm=TRUE)
        MPDfunc_obs <- MPDfunc$mpd.obs[1]
        MNTDfunc_obs <- MNTDfunc$mntd.obs[1]
    } else {
        MPDfunc_z <- MNTDfunc_z <- MPDfunc_obs <- MNTDfunc_obs <- NA
    }
    # Concatenate observed values and z-scores into a vector and return.
    return(data.frame(MPD_phy = MPD_obs, MPD_phy_z = MPD_z, 
                      MNTD_phy = MNTD_obs, MNTD_phy_z = MNTD_z,
                      MPD_func = MPDfunc_obs, MPD_func_z = MPDfunc_z, 
                      MNTD_func = MNTDfunc_obs, MNTD_func_z = MNTDfunc_z))
    
}

# Function to do functional and phylogenetic cumulatively, instead of taxonomic (added 29 June)
# Requires p_dist and f_dist which are matrices with row and column names corresponding to the taxonID, and entries being the pairwise distances
# Also set number of null iterations
PDFD_cumulative <- function(dat, p_dist, f_dist, column = 'taxonID', group = 'site', plot_name = 'plotID', by_year = TRUE, n_iter = 99) {
    
    physpp <- dimnames(p_dist)[[1]]
    funspp <- dimnames(f_dist)[[1]]
    
    dat <- dat %>%
        rename(sp = !!column, plotID = !!plot_name) %>%
        select(siteID, plotID, year, sp)
    if (group == 'site') dat <- dat %>% group_by(siteID)
    if (group == 'plot') dat <- dat %>% group_by(siteID, plotID)
    out <- dat %>%
        do(cbind(year = min(.$year):max(.$year),
                 map_dfr(min(.$year):max(.$year), function(yr) pd_fd(.$sp[.$year <= yr], p_dist, f_dist, n_iter, physpp, funspp))
        ))
    if (!by_year) {
        out %>% filter(year == max(year)) %>% select(-year)
    } else {
        out
    }
}

# 


# Abundance weighted taxonomic diversity metrics (accepts as input a vector of basal areas)
taxonomic_diversity <- function(x) {
    diversity_q1 <- d(x, q = 1)
    diversity_q2 <- d(x, q = 2)
    data.frame(diversity_q1 = diversity_q1,
               diversity_q2 = diversity_q2,
               diversity_shannon = log(diversity_q1),
               diversity_simpson = 1 - 1/diversity_q2)
}


pd_fd_trees <- function(sp_list, abundances = NULL, pddist, fddist, nnull = 99, phylo_spp = NULL, func_spp = NULL) {
    
    require(picante)
    
    ab_weighted <- !is.null(abundances) # Only do abundance-weighted if supplied.
    
    # if no abundances, convert species list to a table
    if (!ab_weighted) {
        m <- t(as.matrix(table(sp_list))) 
    } else {
        m <- t(as.matrix(abundances))
        dimnames(m) <- list(NULL, sp_list)
    }
    
    # Get rid of species that aren't in phylogenetic and functional diversity.
    if (!is.null(phylo_spp)) mphy <- m[, dimnames(m)[[2]] %in% phylo_spp, drop = FALSE] else mphy <- m
    if (!is.null(func_spp)) mfunc <- m[, dimnames(m)[[2]] %in% func_spp, drop = FALSE] else mfunc <- m
    
    if (dim(mphy)[2] > 1) {
        # Calculate PD and do null models
        MPD <- ses.mpd(mphy, pddist, null.model = 'taxa.labels', abundance.weighted = ab_weighted, runs = nnull)
        MNTD <- ses.mntd(mphy, pddist, null.model = 'taxa.labels', abundance.weighted = ab_weighted, runs = nnull)
        # Get all z-scores
        MPD_z <- (MPD$mpd.obs[1] - mean(MPD$mpd.rand.mean, na.rm=TRUE))/sd(MPD$mpd.rand.mean, na.rm=TRUE)
        MNTD_z <- (MNTD$mntd.obs[1] - mean(MNTD$mntd.rand.mean, na.rm=TRUE))/sd(MNTD$mntd.rand.mean, na.rm=TRUE)
        MPD_obs <- MPD$mpd.obs[1]
        MNTD_obs <- MNTD$mntd.obs[1]
    } else {
        MPD_z <- MNTD_z <- MPD_obs <- MNTD_obs <- NA
    }
    if (dim(mfunc)[2] > 1) {
        # Calculate FD and do null models
        MPDfunc <- ses.mpd(mfunc, fddist, null.model = 'taxa.labels', abundance.weighted = ab_weighted, runs = nnull)
        MNTDfunc <- ses.mntd(mfunc, fddist, null.model = 'taxa.labels', abundance.weighted = ab_weighted, runs = nnull)
        # Get all z-scores
        MPDfunc_z <- (MPDfunc$mpd.obs[1] - mean(MPDfunc$mpd.rand.mean, na.rm=TRUE))/sd(MPDfunc$mpd.rand.mean, na.rm=TRUE)
        MNTDfunc_z <- (MNTDfunc$mntd.obs[1] - mean(MNTDfunc$mntd.rand.mean, na.rm=TRUE))/sd(MNTDfunc$mntd.rand.mean, na.rm=TRUE)
        MPDfunc_obs <- MPDfunc$mpd.obs[1]
        MNTDfunc_obs <- MNTDfunc$mntd.obs[1]
    } else {
        MPDfunc_z <- MNTDfunc_z <- MPDfunc_obs <- MNTDfunc_obs <- NA
    }
    # Concatenate observed values and z-scores into a vector and return.
    res <- data.frame(MPD_phy = MPD_obs, MPD_phy_z = MPD_z, 
                      MNTD_phy = MNTD_obs, MNTD_phy_z = MNTD_z,
                      MPD_func = MPDfunc_obs, MPD_func_z = MPDfunc_z, 
                      MNTD_func = MNTDfunc_obs, MNTD_func_z = MNTDfunc_z)
    
    # Specify if abundance weighted or not in the names of the result
    if (!ab_weighted) {
        names(res) <- paste(names(res), 'incidence', sep = '_')
    } else {
        names(res) <- paste(names(res), 'abundance', sep = '_')
    }
    
    return(res)
    
}




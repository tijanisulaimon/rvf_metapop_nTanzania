# program to simulate stochastic SEIR model on simulated movement network 
# aim: simulate RVFV transmission for 1+ years, test interventions
# Tijani Sulaimon
# 27 August 2022

#########
# TO DO #
#########
#### Improvements to the movement and disease model (future work?) ###
# use cattle+sheep+goats, not just cattle
# validate movement data entry
# use data on slaughterhouses
# include main roads in model. e.g. indicator for pairs joined by roads, and/or presents of roads (GIS?)
# include international border as ward characteristic

# clear memory and switch off graphics devices
rm(list = ls())
graphics.off()


#setwd() if needed

# load packages
library(tidyverse)
library(dplyr)
library(igraph)
library(SimInf)
set_num_threads(1)
library(parallel)
library(scales)
library(RColorBrewer)
library(dirmult)
library(corrplot)
library(sf)

# functions

# load sim.hurdle.move function
source('sim.functions.R')


# the date the epidemic starts
start.date <- as.Date("2015-01-01")

# runs through which months?
n.months <- 12
n.days <- n.months/12 * 365
all.months <- formatC(1:n.months, width = 2, flag = "0")

# which day of the year (1:365) is the first of each month?
all.months.day1 <-
  sapply(all.months, function(mm) {
    print(mm)
    yyyy <- 2014 + ceiling(as.numeric(mm)/12)
    mm <- ((as.numeric(mm) - 1) %% 12) + 1
    as.Date(paste(yyyy, mm, "01", sep = "-")) - start.date + 1
  })

# select species
sp <- c("animals", "cattle", "caprine")[2]


# import wards data - required to identify the nearest primary market
wards <- 
  read.csv("WardSpatialData.csv", 
           stringsAsFactors = FALSE, quote = "")

# drop tanga
wards <- droplevels(wards[wards$region %in% c("Arusha", "Manyara", "Kilimanjaro"), ])
all.ward.names <- wards$fullname


# make risk score for RVF spread (as opposed to emergence, which is defined as proportion of
# ward with NDVI between 0.15 and 0.4 and exists in wards DF as wards$rvf.risk).
# there is a negative correlation 
# Will de G find a strong and negative relationship with continuous NDVI
# for seroprevalence risk" 
wards$rvf.risk2 <- 1 - wards$ndvi

# wards to include in the analysis: all wards enough animals 
### check sensitivity to choice of minimum number (500, 1000, etc)
min.ward.pop <- 1000
wards[wards[, paste0("pop.", sp)] < min.ward.pop, paste0("pop.", sp)] <- min.ward.pop
ward.names <- wards$fullname
rownames(wards) <- ward.names

all(all.ward.names == ward.names)

# which wards contain secondary markets?
is.mkt <- ward.names %in% c("arusha/monduli/meserani", "arusha/arusha/bwawani", 
                            "kilimanjaro/hai/machamekusini")
mkts <- ward.names[is.mkt]
not.mkts <- ward.names[!is.mkt]


# import list of contiguous wards, to allow spread via local diffusion as well as 
# via the movement network
contig.wards.all <- 
  read.csv("contiguous.wards.csv", 
           stringsAsFactors = FALSE, quote = "")
rownames(contig.wards.all) <- paste(contig.wards.all$from, contig.wards.all$to, sep = ".")

# import movement parameters

# load predicted probability of movement
mu.z.list <-
  lapply(all.months, function(month) {
    m <- formatC(1 + (as.numeric(month) - 1) %% 12, width = 2, flag = "0")
    mu.z <- 
      as.matrix(
        read.csv(paste0(sp,"/", sp, ".month", m, ".mvt.matrix.prob.2018-06-07-00.csv"), 
                 row.names = 1))
    # make rownames and colnames consistent                  
    colnames(mu.z) <- rownames(mu.z)
    # keep movements to and from included wards
    mu.z[ward.names, ward.names]
  })


# load predicted number moved given movement
mu.c.list <-
  lapply(all.months, function(month) {
    m <- formatC(1 + (as.numeric(month) - 1) %% 12, width = 2, flag = "0")
    mu.c <- 
      as.matrix(
        read.csv(paste0(sp,"/", sp, ".month", m, ".mvt.matrix.rate.2018-06-07-00.csv"), 
                 row.names = 1))
    # make rownames and colnames consistent                  
    colnames(mu.c) <- rownames(mu.c)
    # keep movements to and from included wards
    mu.c[ward.names, ward.names]
  })

# load dispersion (theta) parameter for the ZTNB distribution
theta.list <- 
  lapply(all.months, function(month) {
    m <- formatC(1 + (as.numeric(month) - 1) %% 12, width = 2, flag = "0")
    read.csv(paste0(sp,"/", sp, ".month", m, ".mvt.matrix.rate.theta.2018-06-07-00.csv"))$theta
  })

names(mu.z.list) <- names(mu.c.list) <- names(theta.list) <- all.months


# calculate expected livestock flow for the whole time period (all.months)
mu.yr <- Reduce("+", lapply(all.months, function(m) mu.z.list[[m]] * mu.c.list[[m]]))

# set weights less than 0.2 to 0
#mu.mth <- mu.yr/n.months
 #mu.yr[mu.yr < 0.2] <- 0
 #mu.mth[mu.mth < 0.2] <- 0
# define a market as a ward with expected outflow > 0 for the year
mkts.all <- unique(c(mkts, rownames(mu.yr)[rowSums(mu.yr) != 0]))

# remove secondary market wards from the contiguous wards network
# because these wards are assumed to have no standing population of animals
contig.wards <- contig.wards.all[contig.wards.all$from %in% not.mkts & contig.wards.all$to %in% not.mkts, ]

### create spatial matrix from contiguous ward dataframe
spatial.mat <- as.matrix(as_adjacency_matrix( graph_from_data_frame(contig.wards)))
spatial.mat <- spatial.mat[match(rownames(mu.yr), rownames(spatial.mat)), ]
spatial.mat <- spatial.mat[ ,match(colnames(mu.yr), colnames(spatial.mat))]
spatial.mat[is.na(spatial.mat)] <- 0

# Now add spatial layer to market matrix to create a multiplex matrix

mp.mu.yr.mat <- mu.yr + spatial.mat*12
# mp.mu.mth.df <- bind_rows(
#   as.data.frame(as.table(mu.mth)),
#   as.data.frame(as.table(spatial.mat)),
#   .id = "layer"
# ) %>%
#   filter(Freq != 0) %>%
#   select(
#     from = Var1,
#     to = Var2,
#     weight = Freq,
#     layer
#   )
# 
# mp.mu.mth.df$layer[mp.mu.mth.df$layer == 1] <- "market"
# mp.mu.mth.df$layer[mp.mu.mth.df$layer == 2] <- "spatial"
# ##
# mp.mu.net <- graph_from_data_frame(
#   mp.mu.mth.df, directed = TRUE
# )

bal.mu.net <- graph_from_adjacency_matrix(
  mp.mu.yr.mat, mode = "directed", weighted = TRUE) %>% 
  simplify()

bal.mth.mu.net <- graph_from_adjacency_matrix(
  mu.mth+spatial.mat, mode = "directed", weighted = TRUE) %>% 
  simplify()

yr.mu.net <- graph_from_adjacency_matrix(
  mu.yr, mode = "directed", weighted = TRUE) %>% 
  simplify()
hist(degree(yr.mu.net)); mean(degree(yr.mu.net))

mth.mu.net <- graph_from_adjacency_matrix(
  mu.mth, mode = "directed", weighted = TRUE) %>% 
  simplify()
hist(degree(yr.mu.net)); mean(degree(yr.mu.net))

spatial.net <- graph_from_adjacency_matrix(
  spatial.mat, mode = "directed", weighted = TRUE)

#hist(degree(bal.mu.net)); mean(degree(bal.mu.net))
deg.dist <- degree.distribution(bal.mu.net)
deg.seq <- degree(bal.mu.net)
dens <- edge_density(bal.mu.net)
n.nodes <- vcount(bal.mu.net)
n.edges <- ecount(bal.mu.net)


###################
###################
###################
# GLOBAL SETTINGS #
###################
###################
###################

# how many subnodes to split the wards into, to simulate heterogeneity within wards 
# subnodes form a square grid
n.subnodes <- 8^2 # 1 means don't split # 3^2 means 3x3 grid, etc

# strength of coupling between subnodes (used if n.subnodes > 1)
coupling <- 0.02
#coupling <- 0.0001

# SEIR model parameters
inf.pars <- 
  list(
    beta = 4/7,    # S -> E rate
    epsilon = 1/7,   # E -> I rate
    gamma = 1/7,     # I -> R rate
    i0 = 20,         # number of infectious animals to seed into ward(s) at the start of the epidemic
    var.beta = TRUE)         

# R0
(R0 <- inf.pars$beta / inf.pars$gamma)
inf.pars <- 
  c(inf.pars, 
    cov = 1 - 1/R0) #max(0.5, min(0.1 + 1 - 1/R0, 0.9)))

# exclude seed ward from intervention
exclude.seed <- FALSE     # move to global settings section

# reporting threshold
rep.thresh <- 0.005

# number of animals moving between contiguous wards
contig.move.num <- 
  0.25/n.subnodes * coupling * mean(wards[, paste0("pop.", sp)]) * 365/12 * inf.pars$gamma



#hist(wards[ward.names, "ndvi"])

# allow beta to vary with RVF risk
if(inf.pars$var.beta > 0) {
  wards$beta <- 
    inf.pars$beta * 
    (1 - wards[ward.names, "ndvi"]) / mean(1 - wards[ward.names, "ndvi"])
} else wards$beta <- inf.pars$beta

# R0
wards$R0 <- wards$beta / inf.pars$gamma
hist(wards$R0, xlab = "R0", main = "R0 distribution across 398 wards")
mean(wards$R0)


wards$cov <- 1 - 1/wards$R0
wards$cov[wards$cov < 0] <- 0

# how many vaccines are available? Inf means unlimited
n.vax.dose <- c(Inf, round(sum(wards[, paste0("pop.", sp)]) / 10))[2]



# initial node states
# these inital states are common to all simulations (outside the loop)
n <- length(ward.names)

# how many months does the movement ban apply to? (assume it is reactive)
mban.months <- all.months[-1]

# how many replicates of each scenario to run?
nrep <- 400 # delete this line, it has been defined above?
# Generate theoretical networks
n.net <- nrep

net.list <- list(
  er = replicate(n.net, {network.fn(n.nodes = n.nodes, type = "er")}, simplify = FALSE),
  pa = replicate(n.net, {network.fn(n.nodes = n.nodes, type = "pa")}, simplify = FALSE),
  sw = replicate(n.net, {network.fn(n.nodes = n.nodes, type = "sw")}, simplify = FALSE)
  )

save(net.list,
     file = paste0("sim.output/net.list.", n.net, "_", Sys.Date(), ".Rdata"))
#load(paste0("sim.output/net.list.", 300, "_2022-07-25", ".Rdata"))
####################################################
### Calculate and plot correlation between ranks ###
##### Averaged over all synthetic networks #########

ave.meas.list <-
  lapply(c("er", "pa", "sw"),
         function(y){
           data.frame(
             deg = rank(
               rowMeans(
                 do.call("cbind",
                         lapply(lapply(net.list[[y]],
                                       FUN = igraph::degree),
                                         as.vector)),
                 na.rm = TRUE )
             ),
             btw = rank(
               rowMeans(
                 do.call("cbind",
                         lapply(lapply(net.list[[y]],
                                       FUN = igraph::betweenness),
                                as.vector)),
                 na.rm = TRUE )
             ),
             pag = rank(
               rowMeans(
                 do.call("cbind",
                         lapply(lapply(net.list[[y]], function(x){
                   page_rank(x)$vector}),
                   as.vector)),
                 na.rm = TRUE)
             ),
             risk = rank(wards$rvf.risk),
             row.names = V(net.list[[y]][[1]])$name
           )
         })

names(ave.meas.list) <- c("er", "pa", "sw")

# add tz network measure df to list
ave.meas.list$tz <-
  data.frame(
    deg = rank(degree(delete.edges(yr.mu.net, E(yr.mu.net)[E(yr.mu.net)$weight < 1]))),
    btw = rank(betweenness(delete.edges(yr.mu.net, E(yr.mu.net)[E(yr.mu.net)$weight < 1]))),
    pag = rank(page.rank(delete.edges(yr.mu.net, E(yr.mu.net)[E(yr.mu.net)$weight < 1]))$vector),
    risk = rank(wards$rvf.risk),
    row.names = V(delete.edges(yr.mu.net, E(yr.mu.net)[E(yr.mu.net)$weight < 1]))$name
  )

# calc correlation btw measures for each network and visualise

ave.cor.list <-
  lapply(ave.meas.list,
         cor, method = "spearman")

col1 <- colorRampPalette(c('#7F0000', 'red', '#FF7F00', 'yellow', 'white',
                                  'cyan', '#007FFF', 'blue','#00007F'))

#this is my only addition
par(mfrow=c(2,2))

corrplot(ave.cor.list$er, method="color", cl.pos = "n", type = "lower",
         col=col1(10), cl.length=21,order = "original", addCoef.col="grey")
corrplot(ave.cor.list$pa, method="color", cl.pos = "n",type = "lower",
         col=col1(20), cl.length=21,order = "original", addCoef.col="grey")
corrplot(ave.cor.list$sw, method="color", cl.pos = "n",type = "lower",
         col=col1(20), cl.length=21,order = "original", addCoef.col="grey")
corrplot(ave.cor.list$tz, method="color", cl.pos = "n",type = "lower",
         col=col1(20), cl.length=21,order = "original", addCoef.col="grey")
par(mfrow=c(1,1))
###
# 2. Get a table of network measures

net.summary <-  data.frame(
  network = c("er", "pa", "sw", "tz"),
  N = c(vcount(net.list[["er"]][[1]]),
        vcount(net.list[["pa"]][[1]]),
        vcount(net.list[["sw"]][[1]]),
        vcount(delete.edges(bal.mth.mu.net, E(bal.mth.mu.net)[E(bal.mth.mu.net)$weight < 1]))),

  edge.dens = c(mean(unlist(lapply(net.list[["er"]], edge_density))),
                mean(unlist(lapply(net.list[["pa"]], edge_density))),
                mean(unlist(lapply(net.list[["sw"]], edge_density))),
                edge_density(delete.edges(bal.mth.mu.net, E(bal.mth.mu.net)[E(bal.mth.mu.net)$weight < 1]))),

   mean.deg =  c(mean(unlist(lapply(lapply(net.list[["er"]], degree), mean))),
                 mean(unlist(lapply(lapply(net.list[["pa"]], degree), mean))),
                 mean(unlist(lapply(lapply(net.list[["sw"]], degree), mean))),
                 mean(degree(delete.edges(bal.mth.mu.net, E(bal.mth.mu.net)[E(bal.mth.mu.net)$weight < 1])))),

   sd.deg = c(mean(unlist(lapply(lapply(net.list[["er"]], degree), sd))),
              mean(unlist(lapply(lapply(net.list[["pa"]], degree), sd))),
              mean(unlist(lapply(lapply(net.list[["sw"]], degree), sd))),
              sd(degree(delete.edges(bal.mth.mu.net, E(bal.mth.mu.net)[E(bal.mth.mu.net)$weight < 1])))),

  mean.distance = c(mean(unlist(lapply(net.list[["er"]], mean_distance))),
                mean(unlist(lapply(net.list[["pa"]], mean_distance))),
                mean(unlist(lapply(net.list[["sw"]], mean_distance))),
                mean_distance(bal.mth.mu.net)),

  Transitivity = c(mean(unlist(lapply(net.list[["er"]], transitivity))),
                   mean(unlist(lapply(net.list[["pa"]], transitivity))),
                   mean(unlist(lapply(net.list[["sw"]], transitivity))),
                   transitivity(bal.mth.mu.net)),

  GWCC = c(round(mean(unlist(
    lapply(net.list[["er"]], function(x){max(components(x, "weak")$csize)})))),
    round(mean(unlist(
      lapply(net.list[["pa"]], function(x){max(components(x, "weak")$csize)})))),
    round(mean(unlist(
      lapply(net.list[["sw"]], function(x){max(components(x, "weak")$csize)})))),
           max(components(bal.mth.mu.net, "weak")$csize) ),

  GSCC = c(round(mean(unlist(
    lapply(net.list[["er"]], function(x){max(components(x, "strong")$csize)})))),
    round(mean(unlist(
      lapply(net.list[["pa"]], function(x){max(components(x, "strong")$csize)})))),
    round(mean(unlist(
      lapply(net.list[["sw"]], function(x){max(components(x, "strong")$csize)})))),
           max(components(bal.mth.mu.net, "strong")$csize) )
)

net.summary

####################################################

intervention <- 
  expand.grid(
    network = c("tz", "er", "pa", "sw"),
    vax = c("none", "deg", "btw", "pag", "risk", "rand"),
    mban = "none", 
    prop.vax = seq(0, 0.5, 0.1),
    prop.mban = 0,
    cor = seq(0, 1, 0.1),
    beta = inf.pars$beta,
    gamma = inf.pars$gamma,
    rep = 1:nrep, 
    stringsAsFactors = FALSE)

# intervention <- rbind(intervention.vax, intervention.mban)
dim(intervention)

intervention$prop.vax[intervention$vax %in% "none"] <- 0
intervention$vax[intervention$prop.vax == 0] <- "none"
intervention$prop.mban[intervention$mban %in% "none"] <- 0
intervention$mban[intervention$prop.mban == 0] <- "none"
intervention$cor[intervention$vax %in% "none" & 
                   intervention$mban %in% "none" ] <- 1
intervention$cor[intervention$vax %in% "rand"] <- 1
intervention$cor[intervention$mban %in% "rand" ] <- 1
# intervention$cor[intervention$prop.vax != 0.2] <- 1

intervention <- intervention[!duplicated(intervention), ]
dim(intervention)

if(all(intervention$vax == "none")) 
  intervention$plan <- intervention$mban else
    if(all(intervention$mban == "none")) intervention$plan <- intervention$vax else
      intervention$plan <- paste0("VX: ", intervention$vax, "; MB: ", intervention$mban)
dim(intervention)

intervention$k <- 1:nrow(intervention)

intervention$seed.ward <- "balanced"


########## Function to create event lists for simulation ######################
create.event.tab  <- function(k, intervention.df){

  print(paste("Simulating events for scenario", k, "of", nrow(intervention.df)))
  print(intervention.df[k, ])
  
  # this code doesn't work if intervention.df$seed.ward is e.g. "balanced" or "rand"
  ward.names.ns <-
    if(exclude.seed) ward.names[!ward.names %in% intervention.df$seed.ward[k]] else ward.names
  
  if(intervention.df$network[k] == "tz"){
    #net <- bal.mu.net
    net <- yr.mu.net
    sec.markets <- c("arusha/monduli/meserani", "arusha/arusha/bwawani", 
                  "kilimanjaro/hai/machamekusini")
    is.market <- ward.names %in% sec.markets
    markets <- ward.names[is.market]
    not.markets <- ward.names[!is.market]
    markets.all <- unique(c(markets, rownames(mu.yr)[rowSums(mu.yr) != 0]))
    
    intervention.df$contig.move.n <- contig.move.num
    
    par.list <-
      lapply(all.months, function(m) {
        print(m)
        
        # predicted probability of movement
        mu.z <- mu.z.list[[m]]
        
        # predicted number moved given movement
        mu.c <- mu.c.list[[m]]
        
        # dispersion (theta) parameter for the ZTNB distribution
        theta <- theta.list[[m]]
        
        # problem:
        # because the pink slips record only onward journeys from markets (n = 111), 
        # wards with no market (n = 287) appear to produce no animals. 
        # a simple solution is to assume that they send all their animals to
        # the nearest market. this could be improved by taking account of Gemma's
        # survey data, and also the trade-off between distance and size.
        # go with the simple solution for now.
        # identify wards that sent out no animals
        zero.out <- data.frame(ward = rownames(mu.z)[rowSums(mu.z) == 0], stringsAsFactors = FALSE)
        rowSums(mu.z[zero.out$ward, ]) # zero probability of outflow
        rowSums(mu.c[zero.out$ward, ]) # but nonzero rate given any outflow
        # for those wards that send out no animals 
        # assume that this is due to missing data and assume that all
        # animals from this ward went to the nearest active primary wards (active = outflow > 0) 
        zero.out$nearest.active.primary <- 
          sapply(zero.out$ward, function(w) {
            nap.index <- 
              which.nearest(
                unlist(wards[w, c("lat", "long")]),
                wards[!wards$fullname %in% zero.out$ward & wards$market == "Primary", 
                      c("lat", "long")])
            names(nap.index)
          })
        # how many primary markets have outflow now, but no inflow?
        rownames(mu.z)[rowSums(mu.z) > 0 & !rownames(mu.z) %in% zero.out$nearest.active.primary]
        wards[rownames(mu.z)[rowSums(mu.z) > 0 & !rownames(mu.z) %in% zero.out$nearest.active.primary], "prodsys"]
        # assume that all of these wards are supplied internally.
        # now distribute the outflow from each "nearest active primary" over 
        # the upstream wards with no recorded outflow
        for(nap in unique(zero.out$nearest.active.primary)) {
          mu.z[zero.out$ward[zero.out$nearest.active.primary == nap], nap] <- 1/sum(zero.out$nearest.active.primary == nap)
        }; rm(nap)
        sort(rowSums(mu.z))
        # diagonals should all be zero
        all(diag(mu.z) == 0)
        
        # apply movement ban via the movement probability matrix mu.z
        if( intervention.df$mban[k] != "none" & m %in% mban.months) {
         # keep.wards <- mban.all$ward[mban.all$ward %in% ward.names.ns]
          keep.wards <- 
            select.nodes.measure(net = net, measure = intervention.df$mban[k],
                                 prop = intervention.df$prop.mban[k], 
                                 cor.coeff = intervention.df$cor[k], 
                                 markets = markets.all)
          mu.z[keep.wards, ] <- 0
          mu.z[, keep.wards] <- 0
        }
        
        # multiply matrices to give expected flow
        mu <- mu.z * mu.c
        # check all wards now have outputs
        print(table(rowSums(mu) == 0))
        # output
        list(mu.z = mu.z, mu.c = mu.c, mu = mu, theta = theta)
      })
    names(par.list) <- all.months
    
    # simulate balanced moves
    sim.moves.list <- 
      lapply(par.list, 
             function(mpar) 
               sim.hurdle.move(mu.z = mpar$mu.z, mu.c = mpar$mu.c, alpha = mpar$theta, is.mkt = is.market))
  } else{
    sec.markets <- NULL
    is.market <- ward.names %in% sec.markets
    markets <- ward.names[is.market]
    not.markets <- ward.names[!is.market]
    markets.all <- NULL
    
    intervention.df$contig.move.n <- 0
    
    # generate network for simulating monthly movements and apply movement ban
    net <- net.list[[intervention$network[k]]][[intervention$rep[k]]]
    n.cattle.moved <- 100
    sim.move <- as.matrix(as_adjacency_matrix( net )) * n.cattle.moved
    sim.moves.list <- replicate(n.months, sim.move, simplify = FALSE)
    names(sim.moves.list) <- all.months
    
    # apply movement ban via the theoretical adjacency matrix sim.moves.list
    lapply(all.months, function(m){
      if( intervention.df$mban[k] != "none" & m %in% mban.months) {
  
        keep.wards <- 
          select.nodes.measure(net = net, measure = intervention.df$mban[k],
                               prop = intervention.df$prop.mban[k], 
                               cor.coeff = intervention.df$cor[k])
        sim.moves.list[[m]][keep.wards, ] <- 0
        sim.moves.list[[m]][, keep.wards] <- 0
      }
    })
    }
  
  # remove secondary market wards from the contiguous wards network
  # because these wards are assumed to have no standing population of animals
  contig.wards <- contig.wards.all[contig.wards.all$from %in% not.markets & contig.wards.all$to %in% not.markets, ]

  # create events
  # start with event type "extTrans": moves between nodes
  
  # turn the move matrix list into a table of pairwise moves
  contig.move.n <- intervention.df$contig.move.n[k]
  evts.list <-
    lapply(names(sim.moves.list), 
           function(mth) {
             print(mth)
             #             date1 <- as.Date(paste("2015", mth, "01", sep = "-"))
             date1 <- all.months.day1[mth]
             move <- sim.moves.list[[mth]] 
             
             # make a data frame of movements, one row per movement
             move.pair <- 
               expand.grid(fr.ward = rownames(move), 
                           to.ward = colnames(move), 
                           time = NA, n.moved = NA,
                           stringsAsFactors = FALSE)
             rownames(move.pair) <- paste(move.pair$fr.ward, move.pair$to.ward, sep = ".")
             
             
             # remove rows where no animals were moved
             # plus contiguous ward pairs, if any neighbour-neighbour movement 
             dim(move.pair)
             keep.pairs <- 
               apply(which(move > 0, arr.ind = TRUE), 1, 
                     function(i) paste(rownames(move)[i[1]], colnames(move)[i[2]], sep = "."))
             if(contig.move.n > 0) keep.pairs <- unique(c(keep.pairs, rownames(contig.wards)))
             move.pair <- move.pair[keep.pairs, ]
             dim(move.pair)
             
             # separate into the four different types of movement
             move.order <- c("wm", "mm", "mw", "ww")
             move.pair$move.type <- 
               factor(tolower(paste0(c("w", "m")[(move.pair$fr.ward %in% sec.markets) + 1], 
                                     c("w", "m")[(move.pair$to.ward %in% sec.markets) + 1])), 
                      move.order)
             
             for(m in move.order) {
               move.pair$time[move.pair$move.type == m] <- match(m, move.order) + date1
               move.pair$n.moved[move.pair$move.type == m] <- 
                 diag(
                   as.matrix(
                     move[move.pair$fr.ward[move.pair$move.type == m], 
                          move.pair$to.ward[move.pair$move.type == m]]))
               move.pair
             }; rm(m)
             move.pair <- move.pair[order(move.pair$time), ]
             
             # check that the same total number of animals is being moved
             sum(move.pair$n.moved) == sum(move)
             
             # simulate movements between contiguous wards
             # by moving an extra contig.move.n animals between all contiguous wards each month
             if(contig.move.n > 0) {
               move.pair[rownames(contig.wards), "n.moved"] <- 
                 move.pair[rownames(contig.wards), "n.moved"] + 
                 rpois(nrow(contig.wards), lambda = contig.move.n)
             }
             
             print(paste(round(nrow(contig.wards) * contig.move.n), sp, 
                         "moved between neighbours in month", mth))
             
             # output data frame of movement events
             if(nrow(move.pair) > 0.5) {
               move.df <- 
                 data.frame(
                   event = "extTrans", # options: exit (death), enter (birth), intTrans (e.g. S -> I), extTrans (movement)
                   # Events scheduled at same time processed inorder: exit, enter, internal trans, external trans
                   time = move.pair$time,
                   node = match(move.pair$fr.ward, ward.names),
                   dest = match(move.pair$to.ward, ward.names),
                   n = move.pair$n.moved,
                   proportion = 0,
                   select = 2, # chooses the column (compartment) of E the event operates on
                   shift = 0,
                   stringsAsFactors = FALSE)
             } else {
               move.df <- 
                 data.frame(event = character(0), time = integer(0), node = integer(0),
                            dest = integer(0), n = integer(0), proportion = integer(0),
                            select = integer(0), shift = integer(0)) 
             }
             
             # what proportion of movements are of 1 animal 
             print(paste0(round(100 * mean(move.df$n == 1)), "% of market movements have batch size = 1 in month ", mth))
             # measure changes in ward populations for each node
             active.nodes <- sort(unique(c(move.df$node, move.df$dest)))
             active.nodes <- active.nodes[!active.nodes %in% which(is.market)]
             inout.df.list <- 
               lapply(active.nodes, function(j) {
                 # n out
                 n.out <- sum(move.df$n[move.df$node == j & move.df$event == "extTrans"])
                 # n in
                 n.in <- sum(move.df$n[move.df$dest == j & move.df$event == "extTrans"])
                 imbalance <- n.in - n.out
                 # create births or deaths to correct imbalance in in/out-flow
                 if(n.out != n.in) {
                   demog.df <- 
                     data.frame(
                       event = ifelse(imbalance > 0, "exit", "enter"), # options: exit (death), enter (birth), intTrans (e.g. S -> I), extTrans (movement)
                       time = abs(max(sign(imbalance) * move.pair$time)) + sign(imbalance),
                       node = j,
                       dest = 0,
                       n = abs(imbalance),
                       proportion = 0,
                       select = ifelse(imbalance > 0, 2, 1), # chooses the column (compartment) of E the event operates on
                       shift = 0,
                       stringsAsFactors = FALSE)
                   return(demog.df)
                 } else return(NULL)
               })
             output.df <- rbind(move.df, do.call("rbind", inout.df.list))
             # check balance between inflow, outflow, births and deaths 
             imbalances <- 
               sapply(1:n, function(i) {
                 sum(output.df$n[output.df$event == "extTrans" & output.df$dest == i]) -
                   sum(output.df$n[output.df$event == "extTrans" & output.df$node == i]) +
                   sum(output.df$n[output.df$event == "enter" & output.df$node == i]) - 
                   sum(output.df$n[output.df$event == "exit" & output.df$node == i])
               })
             stopifnot(all(imbalances == 0))
             if(nrow(output.df) > 0.5) output.df$month <- mth
             output.df
           })
  sapply(evts.list, dim)
  
  # rbind the movements into a data frame
  evts <- do.call("rbind", evts.list)
  # table(evts$time)
  evts <- evts[order(evts$time), ]

  # output events table
  return(evts) 
}

###########  ###############

## Pre-split the data into m/n chunks
# split k into 100, 
nrow(intervention)/100
nsplit <- 100
k.split <- split(intervention$k,
                           rep(1:nsplit, each = nrow(intervention)/nsplit))
intervention$k.split <- rep(k.split[[1]], nsplit)
#intervention$k.split <- intervention$k
length(k.split[[1]])
start.time <- Sys.time()

for (i in names(k.split)){
  print(paste0("Creating event table for split ", i))
  assign(paste0("input.list", i),
    mclapply(k.split[[i]],
             FUN = create.event.tab, intervention.df = intervention,
             mc.cores = detectCores())
    )
  output <- paste0("sim.output/input.list", i, ".Rdata")
  save(list = paste0("input.list", i), file = output)
  remove(list = paste0("input.list", i))
}

End.time <- Sys.time()
End.time - start.time


# run SEIR model

sim.disease.trans <- 
  function(k, intervention.df, event.list){
    print(paste("Simulating disease transmission for scenario", k, "of", nrow(intervention.df)))
    print(intervention.df[k, ])
    
    if(intervention.df$network[k] == "tz"){
      sec.markets <- c("arusha/monduli/meserani", "arusha/arusha/bwawani", 
                    "kilimanjaro/hai/machamekusini")
      is.market <- ward.names %in% sec.markets
      markets <- ward.names[is.market]
      not.markets <- ward.names[!is.market]
      cattle.pop <- wards[ward.names[!is.market], paste0("pop.", sp)]
      
      net <- yr.mu.net
    } else{
      sec.markets <- NULL
      is.market <- ward.names %in% sec.markets
      markets <- ward.names[is.market]
      not.markets <- ward.names[!is.market]
      cattle.pop <- rep(10000, length(ward.names))
      net <- net.list[[intervention.df$network[k]]][[intervention.df$rep[k]]]
    }
    
    n <- length(ward.names)
    
    u0.outer <- 
      data.frame(
        S = rep(0, n),
        E = rep(0, n),
        I = rep(0, n),
        R = rep(0, n))
    rownames(u0.outer) <- ward.names
    
    # add susceptible animals
    u0.outer$S[!is.market] <- cattle.pop
    
    # add cumulative count of infecteds and vax column
    u0.outer$Ic <- u0.outer$I
    
    seed.ward <- sample(not.markets, 10,
                        prob = wards[not.markets, "rvf.risk"])
    
   # seed.ward <- sample(not.markets, 10)
    
    # add infected animals to u0
    u0 <- u0.outer
    u0[seed.ward, "I"] <- u0[seed.ward, "Ic"] <- inf.pars$i0 
    u0$S <- u0$S - u0$I
    
    # get vaccination plan for scenario k 
    if(intervention.df$vax[k] == "none"){
      vax.df <- data.frame(
        ward = "none",
        day = Inf,
        cov = 0,
        stringsAsFactors = FALSE)
    } else{
      vax.df <- 
        data.frame(
          ward = select.nodes.measure(
            net = net,
            measure = intervention.df$vax[k],
            prop = intervention.df$prop.vax[k],
            cor.coeff = intervention.df$cor[k]), 
          day = 0, 
          cov = inf.pars$cov, 
          stringsAsFactors = FALSE,
          row.names = NULL)
    }
    vax.per.ward <- 0
    
    if(all(vax.df$day == 0) & !(intervention.df$vax[k] %in% c("none"))) {
      #ward.names.ns <- ward.names[!ward.names %in% seed.ward]
      ward.names.ns <- ward.names
      
      vax.per.ward <- round(u0[vax.df$ward, "S"] * vax.df$cov)
      u0[vax.df$ward, "R"] <- u0[vax.df$ward, "R"] + vax.per.ward
      u0[vax.df$ward, "S"] <- u0[vax.df$ward, "S"] - vax.per.ward
    }
    n.ward.vax <- nrow(vax.df[vax.df$ward != "none", ])
    vax.used <- sum(vax.per.ward)
    if(is.na(vax.used)) vax.used <- 0
    print(intervention.df$vax[k])
    print(paste(n.ward.vax, "wards vaccinated.", vax.used, "doses used."))
    
    ## set up SEIR model ##
    
    # split nodes into subnodes 
    mod.input <- split.nodes.sir(n.subnodes = n.subnodes, u0 = u0, seir = TRUE)
    
    # create N matrix, which determines how internal and external transfer events
    # shift indiviudals between compartments. here the only event that uses shift 
    # is vaccination
    N.mat <- cbind(mod.input$E[, 1])
    # susceptibles move to recovered/removed
    N.mat[substr(rownames(N.mat), 1, 1) == "S", ] <- 
      which(substr(rownames(N.mat), 1, 1) == "R") - which(substr(rownames(N.mat), 1, 1) == "S")
    # exposed stay exposed (vax too late)
    N.mat[substr(rownames(N.mat), 1, 1) == "E", ] <- 0
    # infectious stay infectious (vax too late), and cumulative incidence should never move
    N.mat[substr(rownames(N.mat), 1, 1) == "I", ] <- 0
    # recovereds shouldn't move, as they're already effectively vaxed
    N.mat[substr(rownames(N.mat), 1, 1) == "R", ] <- 0
    
    # set up model 
    model <- mparse(transitions = mod.input$transitions,
                    compartments = mod.input$compartments,
                    gdata = 
                      c(beta = intervention.df$beta[k],  # transmission rate from susceptible to exposed
                        epsilon = inf.pars$epsilon, # rate from exposed to infected 
                        gamma = intervention.df$gamma[k], # recovery rate from infected to recovered
                        coupling = coupling),  # strength of coupling between neighbouring subnodes
                    u0 = mod.input$u0,
                    tspan = 1:n.days,
                    events = event.list[[intervention.df$k.split[k]]],
                    #ldata = matrix(rep(intervention.df$beta[k], nrow(wards)), nrow = 1, dimnames = list("beta", NULL)),
                    E = mod.input$E,
                    N = N.mat)
    
    # add node-specific betas
    #model@ldata <- matrix(wards$beta, nrow = 1)
    
    # run the simulation
    result <- run(model, threads = 1)
    
    u0[u0$Ic > 0, ]
    # plot(result, index =which(u0$I > 0))
    
    event.tab <- table(model@events@event)
    event.report <- 
      paste("There have been", event.tab[1], "deaths,", event.tab[2], "births, and",
            event.tab[3], "movements")
    print(event.report)
    
    # extract the total cumulative incidence across all wards (except markets) 
    # and the total number of animals
    use.wards <- match(not.markets, rownames(u0))
    sir.comp <- mod.input$compartments[mod.input$compartments != "Ic"]
    i.comp <- sir.comp[substring(sir.comp, 1, 1) %in% c("E", "I")]
    cumInc.full <- sapply(use.wards, function(w) trajectory(result, index =w)$Ic)
    cumInc.full.prop <- 
      sapply(use.wards, function(w) {
        trj <- trajectory(result, index =w)
        out <- trj$Ic / rowSums(trj[, sir.comp])
        out[is.na(out) | is.infinite(out)] <- 0
        out
        })
    Inc.full <- sapply(use.wards, function(w) rowSums(cbind(trajectory(result, index =w)[, i.comp])))
    Inc.full.prop <- 
      sapply(use.wards, function(w) {
        trj <- trajectory(result, index =w)
        out <- rowSums(cbind(trj[, i.comp])) / rowSums(trj[, sir.comp])
        out[is.na(out) | is.infinite(out)] <- 0
        out
      })
    Inc <- rowSums(Inc.full)
    cumInc <- rowSums(cumInc.full)
    N <- 
      rowSums(sapply(use.wards, function(w)
        rowSums(trajectory(result, index =w)[, sir.comp])))
    out <-
      cbind(
        intervention.df[k, ],
        day = 1:length(cumInc),
        cumInc = cumInc, 
        Inc = Inc, 
        N = N, 
        cumInc.ward = rowSums(cumInc.full.prop > rep.thresh),
        Inc.ward = rowSums(Inc.full.prop > rep.thresh),
        N.ward = n,
        cumInc.district = 
          sapply(1:n.days, 
                 function(d) length(unique(wards$district[cumInc.full[d, ] > 0]))),
        N.district = length(unique(wards$district)),
        row.names = NULL)
    attr(out, "cumInc.full.prop") <- cumInc.full.prop
    attr(out, "Inc.full.prop") <- Inc.full.prop
    attr(out, "n.ward.vax") <- n.ward.vax
    attr(out, "vax.used") <- vax.used
    attr(out, "event.report") <- event.report
    return(out)
    }

start.time <- Sys.time()

for (i in names(k.split)[1:nsplit]) {
  print(paste0("simulating split ", i))
  load(paste0("sim.output/input.list", i, ".Rdata"))
  # check all input events objects are data frames
  stopifnot(sapply(get(paste0("input.list", i)), class) == "data.frame")
  assign(paste0("sim.res.list", i),
         mclapply(k.split[[i]],
                FUN = sim.disease.trans, intervention.df = intervention,
                event.list = get(paste0("input.list", i)),
                mc.cores = detectCores())
  )
  output <- paste0("sim.output/sim.res.list", i, ".Rdata")
  save(list = paste0("sim.res.list", i), file = output)
  remove(list = paste0("input.list", i))
  remove(list = paste0("sim.res.list", i))
}


End.time <- Sys.time()
End.time - start.time


sim.res.df <- data.frame()
for (i in names(k.split)[1:nsplit]){
  print(paste0("fetching result of split ", i))
  load(paste0("sim.output/sim.res.list", i, ".Rdata"))
  # identify and select successful runs, exclude failed runs
  sim.res.list <- get(paste0("sim.res.list", i))
  sim.res.list <- 
    sim.res.list[which(sapply(sim.res.list, class) == "data.frame")]
  sim.res.df <- 
    rbind.data.frame(sim.res.df, do.call("rbind", sim.res.list),
          stringsAsFactors = FALSE)
  remove(sim.res.list)
  remove(list = paste0("sim.res.list", i))
}


sim.res <- sim.res.df

save(sim.res,
file = paste0("sim.output/sim.res", nrep, "runs.", Sys.Date(), ".Rdata"))

# load("sim.output/sim.res400runs.2022-07-29.Rdata")

############# TJ: summarise results and use ggplot

## average cum inc ward, with confidence interval
cumInc.ward.res <- sim.res %>%
# select(network, vax, prop.vax, cor, day, cumInc.ward, N.ward) %>%
  group_by(network, vax, prop.vax, cor, day) %>%
  summarise(mean.cumInc.ward = mean(cumInc.ward, na.rm = TRUE),
            median.cumInc.ward = median(cumInc.ward, na.rm = TRUE),
            sd.cumInc.ward = sd(cumInc.ward, na.rm = TRUE),
            N.ward = mean(N.ward, na.rm = TRUE),
            n.rep = n(),
            min.cumInc.ward = min(cumInc.ward, na.rm = TRUE),
            max.cumInc.ward = max(cumInc.ward, na.rm = TRUE)) %>%
  mutate(se.cumInc.ward= sd.cumInc.ward / sqrt(n.rep),
         lci.cumInc.ward =
           mean.cumInc.ward - 1.96 * se.cumInc.ward,
         uci.cumInc.ward =
           mean.cumInc.ward + 1.96 * se.cumInc.ward) %>%
  ungroup()



net.types <- c(
  er = "Erdos-Renyi random",
  pa = "Scale-free",
  sw = "Small-world",
  tz = "Data-driven"
)

library("ggthemes")

pdf("sim.output/meanCumInc.pdf")
# Creating a plot
cumInc.ward.res %>%
  filter( cor == 1, prop.vax == 0 ) %>%
  ggplot( mapping = aes(x = day, y = mean.cumInc.ward/N.ward,
                        group = network, color = network) )+
  geom_line()+
  geom_ribbon(aes(ymin = lci.cumInc.ward/N.ward,
                  ymax = uci.cumInc.ward/N.ward, group = network, fill = network),
              alpha=0.5)+
  # facet_wrap(vars(network),
  #            labeller = labeller(
  #              network = net.types))+
  labs(x = "Time (days)",
       y = "Mean cumulative incidence")+
  scale_fill_discrete(name = "Network", labels = net.types)+
  scale_color_discrete(name = "Network", labels = net.types)+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.margin=margin(-5,-5,-5,-5),
        legend.position = "top",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(face = "italic", size = 12),
        strip.text = element_text(face = "bold", size = 12),
        strip.background = element_rect(fill = "white"))

# Closing the graphical device
dev.off()

# get outbreak size at the end of day 365
cumInc.ward.lastday <- cumInc.ward.res %>%
  filter(day == n.days)

cumInc.novax <- cumInc.ward.lastday %>%
  filter(vax == "none")
# cumInc.novax$prop.cumIncWard <- cumInc.novax$mean.cumInc.ward/cumInc.novax$N.ward

cumInc.novax.exp <- expand_grid(
  cumInc.novax,
  vax0 = c("deg", "btw", "pag", "rand", "risk")
) %>% mutate(vax = vax0) %>% select(-vax0)
# check to see all combinations are present
# unique(cumInc.novax.exp$vax)

cumInc.ward.lastday <- rbind(
  cumInc.ward.lastday, cumInc.novax.exp
)

# create new columns of relative difference in cumulative incidence
cumInc.ward.lastday <- cumInc.ward.lastday %>%
  left_join(cumInc.novax %>%
              select(network,
                     mean.novax = mean.cumInc.ward,
                     lci.novax = lci.cumInc.ward,
                     uci.novax = uci.cumInc.ward),
            by =  "network")

cumInc.ward.lastday$RC.cumInc.ward <-
  (cumInc.ward.lastday$mean.novax - cumInc.ward.lastday$mean.cumInc.ward)/
  cumInc.ward.lastday$mean.novax

col1 <- c( "blue", "red", "darkorchid1", "black", "darkgoldenrod3")
col2 <- c( "blue", "red", "darkorchid1", "darkgoldenrod3")

pdf("sim.output/CumIncByPropVac.pdf")          # Paper size

cumInc.ward.lastday %>% filter(cor == 1, vax != "none") %>%
  ggplot(aes(x = prop.vax,
             y = RC.cumInc.ward*100,
             group = vax, color = vax))+
  geom_line(size = 1)+
  geom_point()+
  # geom_errorbar(aes(x=prop.vax, ymin=lci.cumInc.ward, ymax= uci.cumInc.ward,
  #                    group = vax, color = vax))+
  facet_wrap(vars(network),
             labeller = labeller(
               network = net.types))+
  scale_color_manual(name="Vax",
                     labels=c("Betweenness", "Degree",
                              "PageRank", "Random", "Risk"),
                     values=col1)+
  labs(x = "Proportion of vaccinated nodes",
       y = "% reduction in cumulative incidence Vs no-vaccination")+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.margin=margin(-5,-5,-5,-5),
        legend.position = "top",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(face = "italic", size = 12),
        strip.text = element_text(face = "bold", size = 12),
        strip.background = element_rect(fill = "white"))

# Closing the graphical device
dev.off()


# compare other strategies with random vaccination
# FInd better ways of doing this to avoid repetitions?
outsize.rand10 <- cumInc.ward.lastday %>%
  filter(vax == "rand", prop.vax == 0.1)

outsize.rand20 <- cumInc.ward.lastday %>%
  filter(vax == "rand", prop.vax == 0.2)

outsize.rand30 <- cumInc.ward.lastday %>%
  filter(vax == "rand", prop.vax == "0.3")

outsize.rand40 <- cumInc.ward.lastday %>%
  filter(vax == "rand", prop.vax == 0.4)

outsize.rand50 <- cumInc.ward.lastday %>%
  filter(vax == "rand", prop.vax == 0.5)

cumInc.ward.lastday$rand10.cumInc.ward <- 0

cumInc.ward.lastday$rand10.cumInc.ward[cumInc.ward.lastday$network == "er"] <-
  outsize.rand10$mean.cumInc.ward[outsize.rand10$network == "er" &
                                    outsize.rand10$cor == 1]

cumInc.ward.lastday$rand10.cumInc.ward[cumInc.ward.lastday$network == "pa"] <-
  outsize.rand10$mean.cumInc.ward[outsize.rand10$network == "pa"&
                                    outsize.rand10$cor == 1]

cumInc.ward.lastday$rand10.cumInc.ward[cumInc.ward.lastday$network == "sw"] <-
  outsize.rand10$mean.cumInc.ward[outsize.rand10$network == "sw"&
                                    outsize.rand10$cor == 1]

cumInc.ward.lastday$rand10.cumInc.ward[cumInc.ward.lastday$network == "tz"] <-
  outsize.rand10$mean.cumInc.ward[outsize.rand10$network == "tz"&
                                    outsize.rand10$cor == 1]

cumInc.ward.lastday$RCrand10.cumInc.ward <-
  (cumInc.ward.lastday$rand10.cumInc.ward - cumInc.ward.lastday$mean.cumInc.ward)/
  cumInc.ward.lastday$rand10.cumInc.ward

cumInc.ward.lastday$rand20.cumInc.ward <- 0

cumInc.ward.lastday$rand20.cumInc.ward[cumInc.ward.lastday$network == "er"] <-
  outsize.rand20$mean.cumInc.ward[outsize.rand20$network == "er"&
                                    outsize.rand20$cor == 1]

cumInc.ward.lastday$rand20.cumInc.ward[cumInc.ward.lastday$network == "pa"] <-
  outsize.rand20$mean.cumInc.ward[outsize.rand20$network == "pa"&
                                    outsize.rand20$cor == 1]

cumInc.ward.lastday$rand20.cumInc.ward[cumInc.ward.lastday$network == "sw"] <-
  outsize.rand20$mean.cumInc.ward[outsize.rand20$network == "sw"&
                                    outsize.rand20$cor == 1]

cumInc.ward.lastday$rand20.cumInc.ward[cumInc.ward.lastday$network == "tz"] <-
  outsize.rand20$mean.cumInc.ward[outsize.rand20$network == "tz"&
                                    outsize.rand20$cor == 1]

cumInc.ward.lastday$RCrand20.cumInc.ward <-
  (cumInc.ward.lastday$rand20.cumInc.ward - cumInc.ward.lastday$mean.cumInc.ward)/
  cumInc.ward.lastday$rand20.cumInc.ward

cumInc.ward.lastday$rand30.cumInc.ward <- 0

cumInc.ward.lastday$rand30.cumInc.ward[cumInc.ward.lastday$network == "er"] <-
  outsize.rand30$mean.cumInc.ward[outsize.rand30$network == "er"&
                                    outsize.rand30$cor == 1]

cumInc.ward.lastday$rand30.cumInc.ward[cumInc.ward.lastday$network == "pa"] <-
  outsize.rand30$mean.cumInc.ward[outsize.rand30$network == "pa"&
                                    outsize.rand30$cor == 1]

cumInc.ward.lastday$rand30.cumInc.ward[cumInc.ward.lastday$network == "sw"] <-
  outsize.rand30$mean.cumInc.ward[outsize.rand30$network == "sw"&
                                    outsize.rand30$cor == 1]

cumInc.ward.lastday$rand30.cumInc.ward[cumInc.ward.lastday$network == "tz"] <-
  outsize.rand30$mean.cumInc.ward[outsize.rand30$network == "tz"&
                                    outsize.rand30$cor == 1]

cumInc.ward.lastday$RCrand30.cumInc.ward <-
  (cumInc.ward.lastday$rand30.cumInc.ward - cumInc.ward.lastday$mean.cumInc.ward)/
  cumInc.ward.lastday$rand30.cumInc.ward


cumInc.ward.lastday$rand40.cumInc.ward <- 0

cumInc.ward.lastday$rand40.cumInc.ward[cumInc.ward.lastday$network == "er"] <-
  outsize.rand40$mean.cumInc.ward[outsize.rand40$network == "er"&
                                    outsize.rand40$cor == 1]

cumInc.ward.lastday$rand40.cumInc.ward[cumInc.ward.lastday$network == "pa"] <-
  outsize.rand40$mean.cumInc.ward[outsize.rand40$network == "pa"&
                                    outsize.rand40$cor == 1]

cumInc.ward.lastday$rand40.cumInc.ward[cumInc.ward.lastday$network == "sw"] <-
  outsize.rand40$mean.cumInc.ward[outsize.rand40$network == "sw"&
                                    outsize.rand40$cor == 1]

cumInc.ward.lastday$rand40.cumInc.ward[cumInc.ward.lastday$network == "tz"] <-
  outsize.rand40$mean.cumInc.ward[outsize.rand40$network == "tz"&
                                    outsize.rand40$cor == 1]

cumInc.ward.lastday$RCrand40.cumInc.ward <-
  (cumInc.ward.lastday$rand40.cumInc.ward - cumInc.ward.lastday$mean.cumInc.ward)/
  cumInc.ward.lastday$rand40.cumInc.ward


cumInc.ward.lastday$rand50.cumInc.ward <- 0

cumInc.ward.lastday$rand50.cumInc.ward[cumInc.ward.lastday$network == "er"] <-
  outsize.rand50$mean.cumInc.ward[outsize.rand50$network == "er"&
                                    outsize.rand50$cor == 1]

cumInc.ward.lastday$rand50.cumInc.ward[cumInc.ward.lastday$network == "pa"] <-
  outsize.rand50$mean.cumInc.ward[outsize.rand50$network == "pa"&
                                    outsize.rand50$cor == 1]

cumInc.ward.lastday$rand50.cumInc.ward[cumInc.ward.lastday$network == "sw"] <-
  outsize.rand50$mean.cumInc.ward[outsize.rand50$network == "sw"&
                                    outsize.rand50$cor == 1]

cumInc.ward.lastday$rand50.cumInc.ward[cumInc.ward.lastday$network == "tz"] <-
  outsize.rand50$mean.cumInc.ward[outsize.rand50$network == "tz"&
                                    outsize.rand50$cor == 1]

cumInc.ward.lastday$RCrand50.cumInc.ward <-
  (cumInc.ward.lastday$rand50.cumInc.ward - cumInc.ward.lastday$mean.cumInc.ward)/
  cumInc.ward.lastday$rand50.cumInc.ward


pdf("sim.output/cumIncByCorCoef10.pdf")          # Paper size

cumInc.ward.lastday %>% filter(prop.vax == 0.1, vax != "rand") %>%
  ggplot(aes(x = cor,
             y = RCrand10.cumInc.ward*100,
             group = vax, color = vax))+
  geom_line(size = 1)+
  geom_point()+
  facet_wrap(vars(network),
             scales = "free_y",
             labeller = labeller(
               network = net.types))+
  scale_color_manual(name="Strategy",
                     labels=c("Betweenness", "Degree",
                              "PageRank", "Risk"),
                     values=col2)+
  labs(x = "Correlation coefficient",
       y = "% reduction in cummulative incidence Vs random vaccination")+
  scale_x_reverse()+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.margin=margin(-5,-5,-5,-5),
        legend.position = "top",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(face = "italic", size = 12),
        strip.text = element_text(face = "bold", size = 12),
        strip.background = element_rect(fill = "white"))

dev.off()

pdf("sim.output/cumIncByCorCoef20.pdf")

cumInc.ward.lastday %>% filter(prop.vax == 0.2, vax != "rand") %>%
  ggplot(aes(x = cor,
             y = RCrand20.cumInc.ward*100,
             group = vax, color = vax))+
  geom_line(size = 1)+
  geom_point()+
  facet_wrap(vars(network),
             scales = "free_y",
             labeller = labeller(
               network = net.types))+
  scale_color_manual(name="Strategy",
                     labels=c("Betweenness", "Degree",
                              "PageRank", "Risk"),
                     values=col2)+
  labs(x = "Correlation coefficient",
       y = "% reduction in cummulative incidence Vs random vaccination")+
  scale_x_reverse()+
  theme_bw()+
  theme(
        legend.margin=margin(-5,-5,-5,-5),
        legend.position = "top",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(face = "italic", size = 10),
        strip.text = element_text(face = "bold", size = 12),
        strip.background = element_rect(fill = "white"))

dev.off()

pdf("sim.output/cumIncByCorCoef30.pdf")

cumInc.ward.lastday %>% filter(prop.vax == "0.3", vax != "rand") %>%
  ggplot(aes(x = cor,
             y = RCrand20.cumInc.ward*100,
             group = vax, color = vax))+
  geom_line(size = 1)+
  geom_point()+
  facet_wrap(vars(network),
             scales = "free_y",
             labeller = labeller(
               network = net.types))+
  scale_color_manual(name="Strategy",
                     labels=c("Betweenness", "Degree",
                              "PageRank", "Risk"),
                     values=col2)+
  labs(x = "Correlation coefficient",
       y = "% reduction in cummulative incidence Vs random vaccination")+
  scale_x_reverse()+
  theme_bw()+
  theme(
    legend.margin=margin(-5,-5,-5,-5),
    legend.position = "top",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(face = "italic", size = 10),
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "white"))

dev.off()

pdf("sim.output/cumIncByCorCoef40.pdf")

cumInc.ward.lastday %>% filter(prop.vax == 0.4, vax != "rand") %>%
  ggplot(aes(x = cor,
             y = RCrand20.cumInc.ward*100,
             group = vax, color = vax))+
  geom_line(size = 1)+
  geom_point()+
  facet_wrap(vars(network),
             scales = "free_y",
             labeller = labeller(
               network = net.types))+
  scale_color_manual(name="Strategy",
                     labels=c("Betweenness", "Degree",
                              "PageRank", "Risk"),
                     values=col2)+
  labs(x = "Correlation coefficient",
       y = "% reduction in cummulative incidence Vs random vaccination")+
  scale_x_reverse()+
  theme_bw()+
  theme(
    legend.margin=margin(-5,-5,-5,-5),
    legend.position = "top",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(face = "italic", size = 10),
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "white"))

dev.off()

pdf("sim.output/cumIncByCorCoef50.pdf")

cumInc.ward.lastday %>% filter(prop.vax == 0.5, vax != "rand") %>%
  ggplot(aes(x = cor,
             y = RCrand20.cumInc.ward*100,
             group = vax, color = vax))+
  geom_line(size = 1)+
  geom_point()+
  facet_wrap(vars(network),
             scales = "free_y",
             labeller = labeller(
               network = net.types))+
  scale_color_manual(name="Strategy",
                     labels=c("Betweenness", "Degree",
                              "PageRank", "Risk"),
                     values=col2)+
  labs(x = "Correlation coefficient",
       y = "% reduction in cummulative incidence Vs random vaccination")+
  scale_x_reverse()+
  theme_bw()+
  theme(
    legend.margin=margin(-5,-5,-5,-5),
    legend.position = "top",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(face = "italic", size = 10),
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "white"))

dev.off()
# 
# 
# ############ Results summary table #####################
# summary.tab <- cumInc.ward.res %>% filter(day == n.days) %>% 
#   select(network, vax, prop.vax, cor, 
#          mean.cumInc.ward, lci.cumInc.ward, uci.cumInc.ward) %>% 
#   mutate(mean.cumInc.ward = round(mean.cumInc.ward),
#          lci.cumInc.ward = round(lci.cumInc.ward),
#          uci.cumInc.ward = round(uci.cumInc.ward))
# summary.tab$mean.ci = paste0(summary.tab$mean.cumInc.ward, " (",  
#                              summary.tab$lci.cumInc.ward, 
#                              ",", " ", summary.tab$uci.cumInc.ward,")")
# summary.tab <- summary.tab %>% 
#   select(-mean.cumInc.ward, -lci.cumInc.ward, -uci.cumInc.ward)
# prop.vaxtab <- pivot_wider(summary.tab %>% filter(cor==1), 
#                            names_from = prop.vax, values_from = mean.ci)
# xtable::print.xtable(xtable::xtable(prop.vaxtab %>% select(-cor)),
#                      include.rownames = FALSE,
#                      booktabs = TRUE,
#                      scalebox = 0.8)
# ##################################################
# 
# #### Plot RVF risks on map ############
# 
ward.poly <- sf::st_read("2012 Wards Shapefiles/TZwards.shp")
district.poly <- sf::st_read("2012 Wards Shapefiles/districts.geojson")
# 
tz.map.dat <- wards %>% select(
  name, district, region, long, lat, market,
  area.km2, pop2012, status, prodsys,
  rvf.risk, pop.cattle, popdens.cattle,
  ndvi
)

dat.sf <- sf::st_as_sf(x=tz.map.dat,
                   coords = c("long", "lat"))
sf::st_crs(dat.sf) <- sf::st_crs(ward.poly)

### combine both sf data types
join.dat <- dat.sf %>% as.data.frame() %>%
  left_join(ward.poly %>% as.data.frame(),
             by = c( "name" = "Ward_Name" ))

join.dat.sf <- sf::st_as_sf(join.dat)

# randomly assign them
join.dat.sf$rvf.risk.rand1 <- sample(join.dat.sf$rvf.risk)
join.dat.sf$rvf.risk.rand2 <- sample(join.dat.sf$rvf.risk)
join.dat.sf$rvf.risk.rand3 <- sample(join.dat.sf$rvf.risk)

join.dat.sf <- join.dat.sf %>%
  sf::st_sf(sf_column_name = 'geometry.y')

region.cord <-
  data.frame(
    Region = c("Arusha", "Kilimanjaro", "Manyara"),
    long = c(36, 37.63, 36.83),
    lat = c(-3, -3.75, -4.5))

# Creating a plot
plot(join.dat.sf["rvf.risk"],
     key.pos = 4, axes = TRUE,
     extent = st_bbox(c(xmin = 34.8, xmax = 38.5,
                        ymax = -1.9, ymin = -5.5)),
     main =
     "Inter-ward variation in RVF risk")
# dev.off()

# extract latitude
join.dat.sf <- join.dat.sf %>%
  mutate(lat = sf::st_coordinates(geometry.x)[,1],
              lon = sf::st_coordinates(geometry.x)[,2])

tz.map.dat <- rownames_to_column(tz.map.dat, var = "fullname")
ward.coords <- tz.map.dat %>%
  select(fullname, long, lat)

wards <- wards %>%
  mutate(ward_type = case_when(
    fullname %in% mkts.all ~ "Origin",
    TRUE ~ "Destination"
  ),
  deg = degree(bal.mu.net),
  indeg = degree(bal.mu.net, mode = "in"),
  outdeg = degree(bal.mu.net, mode = "out"),
  btw = betweenness(bal.mu.net),
  pag = page_rank(bal.mu.net)$vector
  )

wards$mkts <- ifelse(
  wards$fullname %in% mkts.all, "Market", "Non-market"
)


# Multiplex network consisting of both market movement and spatial layer
bal.mu.net.df <- as_data_frame(bal.mu.net)
bal.mu.net.df <- bal.mu.net.df %>%
  left_join(ward.coords, by = c("from" = "fullname"))

bal.mu.net.df <- bal.mu.net.df %>%
  left_join(ward.coords, by = c("to" = "fullname"))

#pdf("sim.output/bal.map.pdf")

ggplot(data = join.dat.sf) +
  geom_sf(aes(fill = market), color = "grey") +
  geom_curve(data = bal.mu.net.df %>% filter(weight >= 1),
             aes(x = long.x, xend = long.y,
                 y = lat.x, yend = lat.y),
             size = 0.25,
             alpha = 0.3, curvature = 0.5)+
  geom_point(data = wards,
             aes(long, lat,
                 color = ward_type,
                 shape = ward_type), alpha = 0.8)+
  geom_point(data = wards,
             aes(long, lat,
                 color = ward_type,
                 shape = ward_type, size = deg), show.legend = FALSE)+
  scale_color_manual(values = c("cyan", "red"))+
 # scale_color_ipsum()+
  scale_shape_manual(values = c(19, 15))+
  scale_fill_pander()+
  geom_text(data = region.cord,
             aes(long, lat, label = Region),
             size = 5, color = "white")+
  # scale_size_continuous(limits = c(1, max(ward.coords$deg)))+
  coord_sf(xlim = c(34.7, 38.45), ylim = c(-6, -1.65), expand = FALSE)+
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())+
  theme_void()

#dev.off()

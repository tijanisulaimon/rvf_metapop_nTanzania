# This is an attempt to clean/organise Paul's code so we have functions
# in one file and source it into another file where we do explorations

# function to simulate movements (hurdle model) ------------------------------------------
# inputs:
#   - inter-node movement probabilities matrix, mu.z
#   - inter-node movement rate matrix (assuming zero-truncated 
#     negative binomial), mu.c
#   - alpha, the standard dispersion parameter (aka "size") for the
#     negative binomial distribution
#   - a logical vector indicating which nodes (rows/columns) are markets 

# the major change in v3 (2017-12-01) is to allow the mu matrices not to
# be square, i.e. not all wards can generate movements, but all wards 
# can receive movements, so nrow(mu) < ncol(mu)
library(dirmult)
sim.hurdle.move <-
  function(mu.z, mu.c, alpha, is.mkt) {
    
    # if the probability of movement is zero for all pairs of nodes, return 
    # a matrix of zero moves
    if(all(mu.z == 0)) return(mu.z)
    
    # simulate numbers being moved backwards because movements
    # out of markets are more likely to be recorded than movements into markets
    
    # create empty matrix to store one round of moves
    move.c <- move.z <- rep(NA, length = prod(dim(mu.z)))
    attributes(move.c) <- attributes(move.z) <- attributes(mu.z)
    
    # expected number moved
    (mu <- mu.z * mu.c)
    
    # which pairs had any movement at all
    move.z[, ] <- rbinom(prod(dim(mu.z)), 1, c(mu.z))
    # mu.z
    # move.z
    
    # move 4: ward-ward movements
    print("simulating ward-ward movements")
    # first, which wards have any movement
    # move.z[!is.mkt, !is.mkt] 
    
    # and how many movements are there given that movement has happened
    move.c[!is.mkt, !is.mkt] <- 
      rztnbinom.mu(n = sum(!is.mkt)^2, 
                   mu = c(mu.c[!is.mkt, !is.mkt]), 
                   size = alpha) #* c(mu.z[!is.mkt, !is.mkt])
    
    # move 3: movements from markets to wards
    print("simulating market-ward movements")
    # rates:
    # mu.z[is.mkt, !is.mkt]
    # mu.c[is.mkt, !is.mkt]
    # move.z[is.mkt, !is.mkt]
    move.c[is.mkt, !is.mkt] <- 
      rztnbinom.mu(n = sum(is.mkt) * sum(!is.mkt), mu = c(mu.c[is.mkt, !is.mkt]), 
                   size = alpha)
    (move <- move.c * move.z)
    
    # the number that moved out of each market
    (n.move3 <- sum(move[is.mkt, !is.mkt]))
    # rowSums(move[is.mkt, !is.mkt])
    
    # move 2: moves among markets
    print("simulating market-market movements")
    # assume that the between market records are as good as the out of market records,
    # so a proportion of the moves out of markets (move 3), had moved between markets
    # e.g. ...  
    # sum(move["M1", !is.mkt]) # ...animals moved out of market 1. 
    # on average, this number...
    # mu[is.mkt, "M1"]
    # came into M1 from other markets
    
    # the expected number moved among markets:
    # sum(mu[is.mkt, is.mkt])
    # expected proportion of movements from markets to wards that had moved between markets
    # sum(mu[is.mkt, is.mkt]) / sum(mu[is.mkt, !is.mkt])  #### constraint! #### must be < 1 
    
    n.move2 <- 
      if(n.move3 < 0.5) 0 else
        rbinom(1, n.move3, min(1, sum(mu[is.mkt, is.mkt]) / sum(mu[is.mkt, !is.mkt])))
    
    if(n.move2 == 0) {
      move[is.mkt, is.mkt] <- 0
      n.move1.tot <-
        rowSums(move[is.mkt, !is.mkt]) + rowSums(move[is.mkt, is.mkt]) - colSums(move[is.mkt, is.mkt])
    } else {
      repeat {
        move[is.mkt, is.mkt] <-
          rztmultinom.dir(n = 1, size = n.move2, prob = prop.table(c(mu[is.mkt, is.mkt])), 
                          alpha = alpha * length(c(mu[is.mkt, is.mkt])))

        # to recap, this number of animals
        # n.move3
        # have moved from market to wards.
        # this number of them
        # n.move2
        # moved among markets.
        # a total of 
        # n.move3
        # animals must move from wards to markets in move 1.
        # the number moved from each market is 
        # rowSums(move[is.mkt, !is.mkt])
        # however inter-market transfers mean that the numbers moving into
        # the markets must be adjusted
        
        # inflow to each market in move 2
        # colSums(move[is.mkt, is.mkt])
        # outflow from each market in move 2
        # rowSums(move[is.mkt, is.mkt])
        
        print("simulating ward-market movements")
        # to have a net flow of zero through the markets, the number put in from wards must be
        n.move1.tot <-
          rowSums(move[is.mkt, !is.mkt]) + rowSums(move[is.mkt, is.mkt]) - colSums(move[is.mkt, is.mkt])
        
        # due to the patchy nature of the data, it's possible for some of these input totals to be negative
        # which is not possible. this might not happen with real data, but if it does, re-draw until all >= 0
        if(all(n.move1.tot >= 0)) break
      }
    }
    
    # these can be drawn from a multinomial distribution
    if(sum(mu.z[!is.mkt, is.mkt]) > 0) {
      move[!is.mkt, is.mkt] <-
        sapply(rownames(mu.z)[is.mkt], 
               function(mkt) {
                 print(mkt)
                 lambda <- c(mu[!is.mkt, mkt])
                 rztmultinom.dir(n = 1, size = n.move1.tot[mkt], prob = lambda / sum(lambda), alpha = alpha)
               })
    } else {
      move[!is.mkt, is.mkt] <- 0
    }
    
    stopifnot(sum(move * move.z) == sum(move * move.z))
    move
  }


# function to simulate from zero-truncated negative binomial --------------
rztnbinom.mu <-
  function(n, size = Inf, prob, mu) {
    require(actuar)
    if(missing(prob)) prob <- size / (size + mu)
    rztnbinom(n = n, size = size, prob = prob)
  }



# function to approximate a zero-truncated negative binomial --------
# This function approximates a zero-truncated negative binomial using multinomial-dirichlet
# so that the total can be fixed
rztmultinom.dir <- 
  function(n, size, prob, alpha) {
    require(dirmult)
    # if the total is zero, the output elements must all be zero
    if(size == 0) return(rep(0, length(prob)))
    which.zero <- prob == 0
    prob.nz <- prob[!which.zero]
    # because this is a zero-truncated count distribution, all the elements with non-zero mu must generate 1s 
    # so total must be at least the number of non-zero elements
    # if not, output 1s for the highest probability elements, up to a total of "size"
    if(size <= length(prob.nz)) {
      out.nz <- rep(0, length(prob.nz))
      out.nz[(1 + length(prob.nz) - rank(prob.nz)) <= size] <- 1
    } else {
      # only need to simulate for the nonzero elements, so select these
      p <- rdirichlet(n = n, alpha = alpha * prob.nz / sum(prob.nz))
      print(paste("size", size))
      out.nz <- apply(p, 1, function(pp) 1 + rmultinom(n = 1, size = size - length(pp), prob = pp))[, 1]
    }
    out <- rep(NA, length(prob))
    out[which.zero] <- 0
    # replace nonzero elements with out.nz. first check they are the same length
    stopifnot(length(out[!which.zero]) == length(out.nz))
    out[!which.zero] <- out.nz
    out
  }



# Function to calculate the geodesic distance between two points ---------
# specified by radian latitude/longitude using the
# Spherical Law of Cosines (slc)
# adapted from http://www.r-bloggers.com/great-circle-distance-calculations-in-r/
gcd.slc <- function(long, lat) {
  R <- 6371 # Earth mean radius [km]
  latr <- lat*pi/180
  longr <- long*pi/180
  x <- prod(sin(latr)) + prod(cos(latr)) * cos(diff(longr))
  x <- min(x, 1)
  x <- max(x, -1)
  d <- acos(x) * R
  return(d) # Distance in km
}


# function to find nearest location ---------------------------------------
which.nearest <-
  function(latlon, latlon.mat, show.dist = FALSE) {
    d <- 
      apply(latlon.mat, 1, 
            function(x) 
              gcd.slc(long = c(x[2], latlon[2]), lat = c(x[1], latlon[1])))
    i <- which.min(d)
    if(show.dist) print(paste(round(d[i], 2), "km"))
    attributes(i)$min.dist.km <- d[i]
    i
  }


# function to split nodes into subnodes to create within-node heterogeneity --------
# the subnodes are arranged on a square grid, so sqrt(n.subnodes) must be an integer
split.nodes.sir<- 
  function(n.subnodes, u0, adj = NULL, seir = TRUE) {
    
    # n.subnodes <- 4
    # u0 <- data.frame(S = 1000, E = 10, I = 0, R = 0, Ic = 10)
    # adj <- NULL
    # seir <- TRUE
    
    # check that n.subnodes subnodes can be arranged in a square grid
    stopifnot(sqrt(n.subnodes) == round(sqrt(n.subnodes)))
    # define compartments
    compartments <- c("S", "E", "I", "R")
    if(!seir) compartments <- compartments[compartments != "E"]
    # split the compartments into subnodes
    sub.comp <- 
      paste0(rep(compartments, each = n.subnodes), rep(1:n.subnodes, length(compartments)))
    n.sub.comp <- length(sub.comp)
    # create adjacency matrix, if none is supplied
    if(is.null(adj)) {
      d <- sqrt(n.subnodes)
      adj <- 
        structure(rep(0, n.subnodes^2), dim = c(n.subnodes, n.subnodes), 
                  dimnames = list(1:n.subnodes, 1:n.subnodes))
      for(i in 1:nrow(adj)) {
        for(j in 1:ncol(adj)) {
          # detect if subnodes i and j are adjacent on a d X d square matrix. if so, indicate on adjacency matrix with 1,
          # otherwise leave at 0. the spatial arrangement of subnodes is matrix(1:n.subnodes, nrow = d)
          # and adjacency is defined as being one sideways move away, so diagonal neighbours are not adjacent
          if(((i - j) == 1 && j%%d != 0) | ((i - j) == -1 && i%%d != 0) | abs(i - j) == d) adj[i, j] <- 1
        }; rm(j)
      }; rm(i)
    }
    stopifnot(all(dim(adj) == n.subnodes))
    # create the list of transitions between compartments, including contamination from coupled (adjacent)
    # subnodes.
    # for an SIR model the transitions for a single node are
    # "S1 -> (S1+I1+R1) > 0 ? beta*S1*(I1)/(S1+I1+R1) : 0 -> I1 + Ic" 
    # "I1 -> gamma*I1 -> R1" 
    # for an SEIR model the transitions for a single node are
    # "S1 -> (S1+I1+R1) > 0 ? beta*S1*(I1)/(S1+I1+R1) : 0 -> E1 + Ec" 
    # "E1 -> epsilon*E1 -> I1" 
    # "I1 -> gamma*I1 -> R1" 
    # these compartments are split across subnodes
    transitions.tab <-
      sapply(1:n.subnodes, function(sn) {
        nbrs <- rownames(adj)[adj[sn, ] == 1]
        env.contam <- 
          if(length(nbrs) > 0) paste0(" + (", paste0("I", nbrs, collapse = "+"), ")*coupling") else NULL
        N <- paste0("(", paste0(compartments, sn, collapse = "+"), ")")
        s.trans <-
          paste0("S", sn, 
                 " -> ", N, " > 0 ? beta*S", sn, "*(I", sn, # "?" is ifelse in C. E.g. "if a > 0 ? b : c" means ifelse(a > 0, b, c)  
                 env.contam,
                 ")/", N,
                 " : 0 -> ", compartments[2], sn, " + Ic")
        e.trans <-
          if(seir) paste0("E", sn, " -> epsilon*E", sn, " -> I", sn) else NULL
        i.trans <-
          paste0("I", sn, " -> gamma*I", sn, " -> R", sn)
        c(s.trans, e.trans, i.trans)
      })
    transitions <- c(t(transitions.tab))
    # split the initial numbers in each compartment for each node across the subnodes.
    # numbers from each node are randomly (mulitnomially) allocated to subnodes
    u0.list <- 
      lapply(compartments, 
             function(cmp) {
               out <- t(sapply(1:nrow(u0), function(i) rmultinom(1, u0[i, cmp], rep(1, n.subnodes))))
               if(n.subnodes == 1) out <- t(out)
               colnames(out) <- paste0(cmp, 1:n.subnodes)
               out
             })
    u0.out <- data.frame(do.call("cbind", u0.list))
    stopifnot(all(names(u0.out) %in% sub.comp) & all(sub.comp %in% names(u0.out)))
    u0.out <- u0.out[, sub.comp]
    u0.out$Ic <- u0$Ic
    # make the events table. this allows different types of event (e.g. movement, vaccination, etc)
    # to happen to different sets of animals. for example, movements might be taken from all subnodes 
    # and all compartments, while a movement across a border might just happen to edge subnodes.
    E <- matrix(rep(0, n.sub.comp * 2),
                nrow = n.sub.comp,
                dimnames = list(sub.comp))
    E[substr(rownames(E), 1, 1) == "S", 1] <- 1
    E[, 2] <- 1
    E <- 
      rbind(
        cbind(E, do.call("rbind", rep(list(diag(n.subnodes)), length(compartments)))), 
        Ic = 0)
    stopifnot(all(rownames(E) %in% names(u0.out)) & all(names(u0.out) %in% rownames(E)))
    # output all elements as a list
    list(compartments = rownames(E), transitions = transitions, E = E, u0 = u0.out, adj = adj)
  }


# Create a PA network using the Albert-Barabasi model. Notice the GSCC
# of graphs produced by this model is 1. Aim is to increase it by 
# simply allowing wards have outgoing edges
sample_pa_withHighGSCC <- function(n, power, m){
  pa <- sample_pa(n = n, power = power, m = m, directed = TRUE)
  n.df <- igraph::as_data_frame(pa) # convert to df
  # Get all edgelist with incoming edges
  inc.edge.present <- n.df[n.df$from %in% which(degree(pa, mode = "in") != 0), ]
  # Get those without incoming edges
  no.inc.edge <- n.df[n.df$from %in% which(degree(pa, mode = "in") ==0), ]
  # Swap the direction of one of their edges, selected at random
  to.swap <- no.inc.edge %>% group_by(from) %>% 
    sample_n(1)
  to.keep <- no.inc.edge[!(paste(no.inc.edge$from, no.inc.edge$to) %in% 
                             paste(to.swap$from, to.swap$to)  
  ), ]
  swaped <- tibble(
    from = to.swap$to,
    to = to.swap$from
  )
  new.df <- rbind(
    inc.edge.present, to.keep, swaped 
  )
  return(graph_from_data_frame(new.df, directed = TRUE))
}


# function to generate another vector for a given r and a vector x --------
correlated <- function(x, r) { # generate another vector for a given r and vector x
  if(r == 0){
    return( sample(x) )
  } else{
    return(x * sign(r) + rnorm(length(x), sd=sqrt(var(x)/r^2 - var(x))) )
  }
}

# generate networks
network.fn <- function(n.nodes, type){
  stopifnot(type %in% c("er", "pa", "sw"))
  if(type == "er") {
    net <- sample_gnm(n = n.nodes, m = 2*n.nodes, directed = TRUE)
  } else if (type == "pa"){
    net <- sample_pa_withHighGSCC(n = n.nodes, power = 1, m = 2)
    } else if (type == "sw"){
    net1 <- make_lattice(dimvector = n.nodes, nei = 2, directed = TRUE, 
                         circular = TRUE)
    net <- rewire(net1, with = each_edge(prob = 0.1, multiple = FALSE))
  }
  vertex_attr(net, "name", index = V(net)) <- V(bal.mu.net)$name
  return(net)
}


select.nodes.measure <- function(net, measure, prop, cor.coeff, markets = NULL){
  stopifnot(measure %in% c("deg", "btw", "eig", "pag", "rand", "risk", 
                           "indeg", "outdeg", "pag.alt"))
  
  n.nodes <- length( V(net)$name ) 
  n.keep <- round(n.nodes * prop)
  measure.df <- data.frame( node = V(net)$name )

  if(measure == "deg") measure.df$cent <- degree(net)
  if(measure == "indeg") measure.df$cent <- degree(net, mode = "in")
  if(measure == "outdeg") measure.df$cent <- degree(net, mode = "out")
  if(measure == "btw") measure.df$cent <- betweenness(net)
  if(measure == "eig") measure.df$cent <- eigen_centrality(net)$vector
  if(measure == "pag") measure.df$cent <- page_rank(net)$vector
  if(measure == "risk") measure.df$cent <- wards$rvf.risk
  
  # for movement ban, list of all markets have to be selected, important for TZ network
  if( length(markets) != 0 ){
    measure.df <- subset(measure.df, node %in% markets)
    n.nodes <- nrow(measure.df)
    n.keep <- round(n.nodes * prop)
  }
  
  if(measure == "rand") return( as.vector(sample(measure.df$node)[1:n.keep]) )
  
  measure.df$tr.cent.rk <- jitter(rank(measure.df$cent), factor = 0.001)
  measure.df$tr.cent.rk.norm <- (measure.df$tr.cent.rk-1)/(n.nodes-1)
  measure.df$noise.cent.rk <- correlated(x = measure.df$tr.cent.rk.norm, r = cor.coeff)
  measure.df$noise.cent.rk.norm <- 
    (measure.df$noise.cent.rk - min(measure.df$noise.cent.rk))/
    (max(measure.df$noise.cent.rk) - min(measure.df$noise.cent.rk)) 
  measure.df <- measure.df[order(-measure.df$noise.cent.rk.norm)[1:n.keep], ]
  return( as.vector(measure.df$node) )
}




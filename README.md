# Supporting data and code for 'Modelling the effectiveness of targeting vaccination for controlling Rift Valley fever virus using imperfect network information'

This repository contains partial data and code, which were used to produced the manuscript titled 'Modelling the effectiveness of targeting vaccination for controlling Rift Valley fever virus using imperfect network information'


In the manuscript, we used a mathematical model to describe the spread of RVF within and between 398 wards in Northern Tanzania, connected by cattle movements, on which we evaluated the impact of targeting vaccination using imperfect movement data. We show that pre-emptive vaccination guided by only market movement permit data is sufficient to prevent large outbreaks. Targeted control (either by the risk of RVF introduction or onward transmission) at any level of imperfect movement information is preferred over random vaccination, and any improvement in information reliability is advantageous to their effectiveness. 

## Data

The study exploits the multiplex network of cattle movements in Northern Tanzania generated by Chaters et al. (2019). The multiplex contains two layers: a movement network of cattle through markets, and a network connecting adjacent wards. In the first layer, used market movement permit data to generate a static network of cattle movements between 398 wards within three regions (Arusha, Kilimanjaro, and Manyara) in Northern Tanzania, where the number of cattle moved is based on the estimates of the number moved in a month. The study defined a ward as an administrative unit of a mean area $243 \ \text{km}^2$ containing a mean human population of $12\ 000$ and a mean cattle population of $9\ 000$ across all $398$ wards. 

The market permit data were collected as part of the SEEDZ (Social, Economic and Environmental Drivers of Zoonoses in Tanzania) project. This dataset belongs to the Ministry of Livestock and Fisheries of Tanzania, and we do not have permission to share it publicly. This also applies to the data on the number of cattle in each ward. For information on how to access this data, please contact the corresponding author.

## Code
ALl codes, including network generation, the metapopulation model simulation, visualisations were performed in R. The script sim.funtions.R contains all major functions, which were sourced for in the sim.exploration.R script.

#### The functions written in the sim.functions.R script include:

- `sim.hurdle.move`: this function simulate movements based on the hurdle model described in Chaters et al. (2019), and its inputs are:
  - inter-node movement probabilities matrix, mu.z
  - inter-node movement rate matrix (assuming zero-truncated negative binomial), mu.c
  - alpha, the standard dispersion parameter (aka "size") for the negative binomial distribution
  - a logical vector indicating which nodes (rows/columns) are markets 

- `rztnbinom.mu`: function to simulate from zero-truncated negative binomial.

- `rztmultinom.dir`: function to approximate a zero-truncated negative binomial. This function approximates a zero-truncated negative binomial using multinomial-dirichlet so that the total movements can be fixed.

- `gcd.slc`: function to calculate the geodesic distance between two points, specified by radian latitude/longitude using the Spherical Law of Cosines (slc) adapted from http://www.r-bloggers.com/great-circle-distance-calculations-in-r/

- `which.nearest`: function to find nearest location

- `split.nodes.sir`: function to split nodes into subnodes to create within-node heterogeneity the subnodes are arranged on a square grid, so \sqrt{n.subnodes} must be an integer.

- `sample_pa_withHighGSCC': function to create a PA network using the Albert-Barabasi model. Notice the GSCC of graphs produced by this model is 1. Aim is to increase it by simply allowing wards have outgoing edges

- `correlated`: function to generate another vector y for a given r (correlation coefficient) and a vector x 

- `network.fn`: function to generate theoretical networks (ER random, scale-free and small-world)

- `select.nodes.measure`: function to select a proportion of highly ranked nodes in a network, based on a measure of interest, and a correlation coefficient value

#### The sim.explorations.R script 
This script simulate stochastic SEIR model of RVFV transmission on simulated movement networks for $>1$ years, and test interventions. Visulations are produced from this script


 




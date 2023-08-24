#' Create a Transmission Cloud
#'
#' @description Generate a transmission cloud, a collection of all
#' plaussible transmission pairs.
#'
#' @details There are many options to model the DNA evolution;
#' see \code{ape::dist.dna}.
#'
#' @param outbreak object of class "simOutbreak"
#' @param downsampling_prob probability of sampling infected cases
#' @param cutoff genomic distance cutoff
#' @param dna_model DNA evolution model; see Details
#' @param seed seed number
#'
#' @import epicontacts
#' @import igraph
#' @import dplyr
#'
#' @return The function returns the following outputs:\cr
#' \itemize{
#'     \item \code{tt} a linelist of the true transmission pairs.
#'     \item \code{tc} a data frame of all plausible transmission pairs.
#' }
#'
#' The \code{tc} data frame consists of these following columns:\cr
#' \itemize{
#'     \item \code{inf.source} IDs of the potential infectors
#'     \item \code{inf.ID} IDs of the infectees
#'     \item \code{si} serial interval data
#'     \item \code{type} transmission type: coprimary, non-coprimary, neither
#'     \item \code{M} number of unsampled intermediaries in the non-coprimary transmission
#'     \item \code{d} genomic distance
#' }
#'
#' @export
#'
#' @examples
#' p <- .7 # sampling proportion
#' epsilon <- 15 # genomic distance cutoff
#' outbreak <- simOutbreak(genint_params = c(4.5,2))
#' mytc <- createTC(outbreak, p, epsilon)
#' plot(mytc$tt,
#'      thin = F,
#'      x_axis = "inf.times",
#'      node_color = "group",
#'      node_value = 1,
#'      col_pal = c(unsampled = "black", sampled = "gold"),
#'      arrow_size = .5,
#'      node_size = 5,
#'      edge_width = .5,
#'      node_width = .5,
#'      height = 800,
#'      width = 600,
#'      label = F)

createTC.simOutbreak <- function(outbreak,
                     downsampling_prob,
                     cutoff,
                     dna_model = "N",
                     seed = 1){

  set.seed(seed)
  # downsampling the true transmission tree
  unsampled <- outbreak$epidata %>%
    slice_sample(prop = 1-downsampling_prob) %>%
    pull(inf.ID) #get unsampled cases' name
  sampled <- outbreak$epidata %>%
    filter(!(inf.ID %in% unsampled)) %>%
    pull(inf.ID)


  # find SNP distance
  pdist <- ape::dist.dna(outbreak$aligndata[sampled], dna_model, as.matrix = T)


  # infectee whose infector is sampled
  newepi <- outbreak$epidata %>%
    filter(!(inf.ID %in% unsampled), !(inf.source %in% unsampled), !is.na(si)) %>%
    select(inf.source, inf.ID, si) %>%
    mutate(type = "non-coprimary", M = 0)
  if(nrow(newepi) > 1) {
    newepi$d <- diag(pdist[newepi$inf.source, newepi$inf.ID])
  } else newepi$d <- pdist[newepi$inf.source, newepi$inf.ID]


  # find infectee that has unsampled infector
  findInf <- sampled[which(!(sampled %in% newepi$inf.ID))]


  # create a transmission cloud
  ## find all potential infectors for each intectee whose infector is unsampled
  newtc <- lapply(findInf, function(x){
    infector <- names(which(pdist[,x] <= cutoff))
    infector <- infector[infector != x]
    d <- pdist[infector, x]
    si <- outbreak$epidata$inf.times[outbreak$epidata$inf.ID == x] - outbreak$epidata$inf.times[outbreak$epidata$inf.ID %in% infector]

    data.frame(inf.ID = rep(x, length(infector)), inf.source = infector, si = si, d = d) %>%
      filter(si >= 0)
  })
  newtc <- bind_rows(newtc)
  row.names(newtc) <- NULL

  ## get the transmission types and exclude the coprimary with intermediaries
  ### define graph to find the transmission path
  gg <- igraph::graph_from_data_frame(outbreak$epidata %>%
                                        select(inf.source, inf.ID, si) %>%
                                        filter(!is.na(si)),
                                      directed = F,
                                      outbreak$epidata %>%
                                        select(inf.ID))
  ## get the transmission type
  getTType <- function(case1, case2, graph, outbreak){

    # find the shortest path
    path <- igraph::shortest_paths(graph, case1, case2)
    path <- labels(path$vpath[[1]])



    # find the mrca
    # mrca must be the one having the earliest infectiousness time in the path
    if(length(path) > 0){
      mrca <- outbreak$epidata %>%
        filter(inf.ID %in% path) %>%
        filter(inf.times == inf.times[which.min(inf.times)]) %>%
        pull(inf.ID)
    } else{
      mrca <- outbreak$epidata %>%
        filter(inf.ID %in% c(case1, case2)) %>%
        filter(inf.times == inf.times[which.min(inf.times)]) %>%
        pull(inf.ID)
    }


    # get the transmission type
    k <- length(path)-2
    if(!(mrca %in% c(case1, case2))) {
      if(k == 1){
        type <- "coprimary"
        M <- NA
      } else {
        type <- "neither"
        M <- NA
      }
    } else {
      type <- "non-coprimary"
      M <- k
    }


    # output
    return(data.frame(from = case1, to = case2, type = type, M = M))
  }

  cl <- makeCluster(detectCores()-1, "SOCK")
  registerDoSNOW(cl)

  tmp <- foreach(i=1:nrow(newtc),
                 .packages = c("dplyr"),
                 #.export = c("getTType"),
                 .combine = bind_rows,
                 .errorhandling = "pass") %dopar% {
                   .GlobalEnv$getTType <- getTType

                   getTType(newtc$inf.source[i], newtc$inf.ID[i], gg, outbreak)
                 }
  stopCluster(cl)

  newtc <- newtc %>%
    left_join(tmp, by = c("inf.source"="from", "inf.ID"="to")) #%>%
  #filter(type != "neither") # exclude coprimary with intermediate

  ## the transmission cloud
  newepi <- newepi %>% bind_rows(newtc)
  newepi$type <- factor(newepi$type, c("coprimary", "non-coprimary", "neither"))



  # plot true transmission tree
  outbreak$epidata$group <- ifelse(!(outbreak$epidata$inf.ID %in% unsampled), "sampled", "unsampled")
  epicont <- epicontacts::make_epicontacts(outbreak$epidata,
                                           contacts = outbreak$epidata %>%
                                             select(inf.source, inf.ID, nmut) %>%
                                             filter(!is.na(inf.source)),
                                           id = "inf.ID",
                                           from = "inf.source",
                                           to = "inf.ID",
                                           directed = T)



  # output
  list(tt = epicont, tc = newepi)
}




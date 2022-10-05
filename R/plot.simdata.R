#' Plot Transmission Network
#'
#' A function to visualize the generated transmission network by \code{simulateData()}.
#'
#' @param x object of class \code{simdata}.
#'
#' @return A \code{visNetwork} visualization plot of the generated transmission network. \cr\cr
#'
#' @details The node size is scaled proportionally to the amount of direct/indirect infectees
#' connected to the case, and the edge length is proportional to the period of symptom onset
#' times between two connected cases.
#'
#' @export
#'
#' @import visNetwork
#' @import dplyr
#'
#' @examples
#' data <- simulateData(n = 100, mu = 4, sigma = 2, R0 = .7)
#' plot(data)
plot.simdata <- function(x){
  tc <- x$transnet


  # create nodes and edge data frames
  epi.nodes <- data.frame(id = unique(c(tc$case_i, tc$case_j))) %>%
    left_join((tc %>% group_by(case_i) %>% summarise(value = n())), by = c("id"="case_i")) %>%
    mutate(value = ifelse(is.na(value), 1, value+1), font.size = 0, color = "coral")
  epi.edges <- data.frame(from = tc$case_i,
                          to = tc$case_j,
                          length = tc$onset_diff+5,
                          arrows = "to",
                          color = "black")
  #



  # draw visnetwork
  myplot <- visNetwork(epi.nodes, epi.edges, width = "100%", height = "500px") %>%
    visOptions(highlightNearest = list(enabled = TRUE, algorithm = "hierarchical")) %>%
    visNodes(scaling = list(min = 10, max = 50)) %>%
    visEdges(scaling = list(min = 1, max = 10)) %>%
    visLegend() %>%
    visLayout(randomSeed = 1234567)

  return(myplot)
}

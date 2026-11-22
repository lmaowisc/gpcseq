#' Four patient histories at two calendar looks
#'
#' @description Return a small, fully specified data cut in the format accepted by \code{\link{gpc_fit}} and \code{\link{gpc_update}}. It provides independently calculable answers for all three methods.
#' @param look Integer 1 or 2. In this example the calendar time of each look equals its index. All patients enter at calendar time 0.
#' @details Treated patient T1 remains alive with no admission through time 2; T2 dies at 1.5. Control patient C1 dies at time 1; C2 is hospitalized at 0.5 and remains alive. For the state analysis, \code{"hospitalized"} denotes a history of admission, rather than current hospital occupancy; the best-to-worst order is \code{"well"}, \code{"hospitalized"}, \code{"dead"}. Death is absorbing.
#'
#' At look 1 the treatment member wins in every pair using survival first and hospitalization count second. All four scores are identical and no monitoring test is performed. At look 2 there are three wins and one loss, giving a mean signed score of 0.5 and an estimated variance of 0.25. Through horizon 2 the win odds and win ratio both equal 3, and the restricted net win time equals 0.875 time units. These are calculation examples, not evidence about trial operating characteristics.
#' @returns A list of three data frames: \code{patients} (four enrolled records), \code{events} (clinical events for win--loss comparison), and \code{transitions} (observed state changes for Aalen--Johansen), restricted to the selected look. All times use the same arbitrary unit.
#' @seealso \code{\link{gpc_fit}}, \code{\link{gpc_update}}.
#' @examples
#' cut <- gpc_example(2)
#' cut$patients
#' plan <- gpc_design(2, 0.05, "ipcw", horizons = 2)
#' fit <- gpc_fit(plan, 1, cut$patients, cut$events)
#' gpc_influence(fit)
#' @export
gpc_example <- function(look=1L) {
  if (!is.numeric(look) || length(look)!=1L || is.na(look) || !look %in% 1:2)
    stop("For the example, look is 1 or 2.",call.=FALSE)
  time <- c(Inf,1.5,1,Inf)
  patients <- data.frame(id=c("T1","T2","C1","C2"),arm=c(1L,1L,0L,0L),
    entry=0, followup=pmin(look,time),death=time<=look,state0="well")
  events <- data.frame(id="C2",time=0.5,type="hospitalization",value=NA_real_)
  transitions <- data.frame(id=c("C2","C1","T2"),time=c(.5,1,1.5),from="well",
                            to=c("hospitalized","dead","dead"))
  list(patients=patients,events=events,transitions=transitions[transitions$time<=look,,drop=FALSE])
}

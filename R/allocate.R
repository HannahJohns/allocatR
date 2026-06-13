#' Randomise a participant
#'
#' @description Randomise a participant according to the MSB, CS-MSB or CSSize-MSB algorithms.
#'
#' @usage allocate( method = c("simple","csmsb","msb"),
#'                  arm,
#'                  data = NULL,
#'                  newData = NULL,
#'                  settings = NULL
#'                  allocationSequence = NULL,
#'                  implementation = c("cpp", "R")
#'                  )
#'
#' @param method A string indicating the randomisation method to be used
#' @param arm An integer vector containing the current treatment assignments
#' @param data A data frame with \code{length(arm)} rows containing relevant data (e.g. covariates) for currently assigned participants
#' @param newData A one-row data frame with the same columns as \code{data} containning relevant data (e.g. covariates) for the new participant to be assigned
#' @param settings A list containing settings to be used by a specific randomisation method
#' @param allocationSequence An object of type \code{Alloc_seq}
#' @param implementation A string indicating which implementation of a given method should be used
#'
#' @returns A list with class \code{Allocation} containing the following:
#' \describe{
#' \item assignment
#' \item details a list containing internal steps taken to reach this assignment result
#' }
#'
#' @details
#' Additional details...
#'
#'
allocate <- function(  method = c("simple","csmsb","msb"),
                       arm,
                       data = NULL,
                       newData = NULL,
                       settings = NULL,
                       randNumber = NULL,
                       implementation = c("cpp", "R")
){




}






allocate_simple <- function()


allocate_msb <- function()


allocate_msb <- function()



#' data, arm, newRecord,
#'                   method=c('CSMSB', 'MSB'),
#'                   voteMethod=NULL,
#'                   biasMethod=c("static", "linear"),
#'                   weights = NULL,
#'                   alpha,
#'                   biasParameters,
#'                   groupImbalanceThreshold=Inf,
#'                   imbalancePriority=c("abstain", "covariate", "group")
#'                   trng = NULL,
#'                   jump_trng = 0
#'                   )
#' @param method
#' @param data A data frame of covariates for current records. 0 rows if first randomisation.
#' @param arm  A vector of arm assignments for current records. Length 0 if first randomisation.
#' @param newRecord A 1 by \code{ncol(data)} data frame of covariates for the new record.
#' @param method A character string indicating the type of method to use.  Defaults to CSMSB.
#' @param voteMethod A character string indicating how votes should be combined. Defaults to "majority" for MSB and "pooled" for CSMSB.
#' @param biasMethod A character string indicating how a bias coin probability should be generated. Defaults to "static"
#' @param weights A vector of length \code{ncol(data)} giving the weights for each column in order. Defaults to even weights if not supplied.
#' @param alpha A number indicating the level of significance for when MSB/CSMSB intervenes.
#' @param biasParameters A numeric vector providing parameters to the indicated bias method
#' @param groupImbalanceThreshold Group allocation imbalance threshold required to trigger CSSize-MSB allocation rules. By default will never trigger.
#' @param imbalancePriority Priority rule for CSSize-MSB. Defaults to "abstain"
#' @param trng an optional parallel random number generator from the rTRNG package. If an integer is supplied instead, a yarn3 engine will be generated using this integer as the seed.
#' @param jump_trng a nonnegative integer specifying how many steps the trng engine should be jumped forward before it is used.
#'
#' @details
#'
#' Variable type is read from the class of each column in \code{data}.
#' Factors are assumed categorical, ordered factors are assumed ordinal
#' and numeric are assumed continuous.
#'
#' When \code{method=="MSB"}, the Minimal Sufficient Balance algorithm is used.
#' As implemented here, it uses the Chi-square test to evaluate imbalance on
#' categorical variables, and the two group z-test (assuming equal variance) to
#' evaluate imbalance on continuous variables. This is done even when the
#' chi-square test may be inaccurate due to low cell counts or where the sample variance is large.
#' Ordinal variables are not supported by this algorithm.
#'
#' If the MSB algorithm is used, \code{voteMethod} must be \code{"majority"} and
#' \code{biasMethod} must be \code{"static"}
#'
#' When \code{method=="CSMSB"}, the Common Scale Minimal Sufficient Balance
#' algorithm is used. It always uses the Wilcoxon-Mann-Whitney odds ratio to
#' evaluate imbalance. For categorical variables, it uses a binary representation
#' of the variable, which is 1 if it matches the value in \code{newRecord} or 0 otherwise.
#'
#' \code{biasParameters} should be either a single number or a vector of numbers
#' depending on \code{biasMethod}
#'
#' \enumerate{
#' \item static - a single number between 0.5 and 1 giving the strength of the biased coin
#' \item linear - a single number between 1 and infinity that scales the imbalance measure to calculate a biased coin strength according to the formula 0.5+imbalance*biasParameters
#' }
#'
#' When using majority or weighted voting, there are a few quirks to this implementation.
#' \enumerate{
#'     \item When checking each variable, we assume that a vote is always valid and log the vote as if it is.
#'     \item We remove invalid votes based on \code{alpha}, immediately before we tally votes.
#' }
#' This is mainly done so we can log small imbalances for understanding trial behaviour over time.
#'
#' The \code{details} list returned by this function varies in contents depending on the method used.
#'
#' If \code{method=="MSB"} then it will contain the following:
#'
#' \itemize{
#'  \item testStatistic - the test statistics used for each
#'  \item vote - the votes placed by each variable
#' }
#'
#' If \code{method=="CSMSB"} then it will contain the following:
#'
#' \itemize{
#'  \item testStatistic - The test statistic(s) used by the CSMSB algorithm
#'  \item vote - The votes placed by each variable
#'  \item imbalance - The current imbalance on each variable
#'  \item seImbalance - The standard error of the current imbalance on each variable
#'  \item prospectiveImbalance0 - The imalance on each variable if the new participant is assigned to control
#'  \item prospectiveImbalance1 - The imalance on each variable if the new participant is assigned to intervention
#'  \item BD - The pooled measure of imbalance
#'  \item seBD - The standard error of the pooled measure of imbalance
#'  }
#'
#' @return  A list containing the following:
#' \itemize{
#'     \item decision - A randomisation decision (0 or 1)
#'     \item prob - The probability of assignment to arm 0
#'     \item biasDirection - The direction of the assignment bias (0 or 1)
#'     \item details - A list containing further information about the randomisation decision. The specific contents of this list depend on the method. See Details above.
#'}
allocate <- function(  method = c("simple","csmsb","msb"),
                       arm,
                       data = NULL,
                       newData = NULL,
                       settings = NULL
                       randNumber = NULL,
                       implementation = c("cpp", "R")
){

}

#' The Params virtual class
#'
#' Virtual S4 class that all other Params classes inherit from.
#'
#' @section Parameters:
#'
#' The Params class defines the following parameters:
#'
#' \describe{
#'     \item{\code{nGenes}}{The number of genes to simulate.}
#'     \item{\code{nCells}}{The number of cells to simulate.}
#'     \item{\code{[seed]}}{Seed to use for generating random numbers.}
#' }
#'
#' The parameters not shown in brackets can be estimated from real data.
#'
#' @name Params
#' @rdname Params
#' @aliases Params-class
setClass("Params",
         contains = "VIRTUAL",
         slots = c(nGenes = "numeric",
                   nCells = "numeric",
                   seed = "numeric"),
         prototype = prototype(nGenes = 10000, nCells = 100,
                               seed = sample(seq_len(1e6), 1)))


#' The SimpleParams class
#'
#' S4 class that holds parameters for the simple simulation.
#'
#' @section Parameters:
#'
#' The simple simulation uses the following parameters:
#'
#' \describe{
#'     \item{\code{nGenes}}{The number of genes to simulate.}
#'     \item{\code{nCells}}{The number of cells to simulate.}
#'     \item{\code{[seed]}}{Seed to use for generating random numbers.}
#'     \item{\code{mean.shape}}{The shape parameter for the mean gamma
#'     distribution.}
#'     \item{\code{mean.rate}}{The rate parameter for the mean gamma
#'     distribution.}
#'     \item{\code{[count.disp]}}{The dispersion parameter for the counts
#'     negative binomial distribution.}
#' }
#'
#' The parameters not shown in brackets can be estimated from real data using
#' \code{\link{simpleEstimate}}. For details of the simple simulation
#' see \code{\link{simpleSimulate}}.
#'
#' @name SimpleParams
#' @rdname SimpleParams
#' @aliases SimpleParams-class
#' @exportClass SimpleParams
setClass("SimpleParams",
         contains = "Params",
         slots = c(mean.shape = "numeric",
                   mean.rate = "numeric",
                   count.disp = "numeric"),
         prototype = prototype(mean.shape = 0.4, mean.rate = 0.3,
                               count.disp = 0.1))

#' The SplatParams class
#'
#' S4 class that holds parameters for the Splat simulation.
#'
#' @section Parameters:
#'
#' The Splat simulation requires the following parameters:
#'
#' \describe{
#'     \item{\code{nGenes}}{The number of genes to simulate.}
#'     \item{\code{nCells}}{The number of cells to simulate.}
#'     \item{\code{[seed]}}{Seed to use for generating random numbers.}
#'     \item{\emph{Batch parameters}}{
#'         \describe{
#'             \item{\code{[nBatches]}}{The number of batches to simulate.}
#'             \item{\code{[batchCells]}}{Vector giving the number of cells in
#'             each batch.}
#'             \item{\code{[batch.facLoc]}}{Location (meanlog) parameter for the
#'             batch effect factor log-normal distribution. Can be a vector.}
#'             \item{\code{[batch.facScale]}}{Scale (sdlog) parameter for the
#'             batch effect factor log-normal distribution. Can be a vector.}
#'             \item{\code{[batch.rmEffect]}}{Logical, removes the batch effect
#'             and continues with the simulation when TRUE. This allows the 
#'             user to test batch removal algorithms without having to calculate
#'             the new expected cell means with batch removed.}
#'         }
#'     }
#'     \item{\emph{Mean parameters}}{
#'         \describe{
#'             \item{\code{mean.shape}}{Shape parameter for the mean gamma
#'             distribution.}
#'             \item{\code{mean.rate}}{Rate parameter for the mean gamma
#'             distribution.}
#'         }
#'     }
#'     \item{\emph{Library size parameters}}{
#'         \describe{
#'             \item{\code{lib.loc}}{Location (meanlog) parameter for the
#'             library size log-normal distribution, or mean parameter if a
#'             normal distribution is used.}
#'             \item{\code{lib.scale}}{Scale (sdlog) parameter for the library
#'             size log-normal distribution, or sd parameter if a normal
#'             distribution is used.}
#'             \item{\code{lib.norm}}{Logical. Whether to use a normal
#'             distribution for library sizes instead of a log-normal.}
#'         }
#'     }
#'     \item{\emph{Expression outlier parameters}}{
#'         \describe{
#'             \item{\code{out.prob}}{Probability that a gene is an expression
#'             outlier.}
#'             \item{\code{out.facLoc}}{Location (meanlog) parameter for the
#'             expression outlier factor log-normal distribution.}
#'             \item{\code{out.facScale}}{Scale (sdlog) parameter for the
#'             expression outlier factor log-normal distribution.}
#'         }
#'     }
#'     \item{\emph{Group parameters}}{
#'         \describe{
#'             \item{\code{[nGroups]}}{The number of groups or paths to
#'             simulate.}
#'             \item{\code{[group.prob]}}{Probability that a cell comes from a
#'             group.}
#'         }
#'     }
#'     \item{\emph{Differential expression parameters}}{
#'         \describe{
#'             \item{\code{[de.prob]}}{Probability that a gene is differentially
#'             expressed in a group. Can be a vector.}
#'             \item{\code{[de.downProb]}}{Probability that a differentially
#'             expressed gene is down-regulated. Can be a vector.}
#'             \item{\code{[de.facLoc]}}{Location (meanlog) parameter for the
#'             differential expression factor log-normal distribution. Can be a
#'             vector.}
#'             \item{\code{[de.facScale]}}{Scale (sdlog) parameter for the
#'             differential expression factor log-normal distribution. Can be a
#'             vector.}
#'         }
#'     }
#'     \item{\emph{Biological Coefficient of Variation parameters}}{
#'         \describe{
#'             \item{\code{bcv.common}}{Underlying common dispersion across all
#'             genes.}
#'             \item{\code{bcv.df}}{Degrees of Freedom for the BCV inverse
#'             chi-squared distribution.}
#'         }
#'     }
#'     \item{\emph{Dropout parameters}}{
#'         \describe{
#'             \item{\code{dropout.type}}{The type of dropout to simulate.
#'             "none" indicates no dropout, "experiment" is global dropout using
#'             the same parameters for every cell, "batch" uses the same
#'             parameters for every cell in each batch, "group" uses the same
#'             parameters for every cell in each groups and "cell" uses a
#'             different set of parameters for each cell.}
#'             \item{\code{dropout.mid}}{Midpoint parameter for the dropout
#'             logistic function.}
#'             \item{\code{dropout.shape}}{Shape parameter for the dropout
#'             logistic function.}
#'         }
#'     }
#'     \item{\emph{Differentiation path parameters}}{
#'         \describe{
#'             \item{\code{[path.from]}}{Vector giving the originating point of
#'             each path. This allows path structure such as a cell type which
#'             differentiates into an intermediate cell type that then
#'             differentiates into two mature cell types. A path structure of
#'             this form would have a "from" parameter of c(0, 1, 1) (where 0 is
#'             the origin). If no vector is given all paths will start at the
#'             origin.}
#'             \item{\code{[path.nSteps]}}{Vector giving the number of steps to
#'             simulate along each path. If a single value is given it will be
#'             applied to all paths. This parameter was previously called
#'             \code{path.length}.}
#'             \item{\code{[path.skew]}}{Vector giving the skew of each path.
#'             Values closer to 1 will give more cells towards the starting
#'             population, values closer to 0 will give more cells towards the
#'             final population. If a single value is given it will be applied
#'             to all paths.}
#'             \item{\code{[path.nonlinearProb]}}{Probability that a gene
#'             follows a non-linear path along the differentiation path. This
#'             allows more complex gene patterns such as a gene being equally
#'             expressed at the beginning an end of a path but lowly expressed
#'             in the middle.}
#'             \item{\code{[path.sigmaFac]}}{Sigma factor for non-linear gene
#'             paths. A higher value will result in more extreme non-linear
#'             variations along a path.}
#'     }
#'   }
#' }
#'
#' The parameters not shown in brackets can be estimated from real data using
#' \code{\link{splatEstimate}}. For details of the Splat simulation
#' see \code{\link{splatSimulate}}.
#'
#' @name SplatParams
#' @rdname SplatParams
#' @aliases SplatParams-class
#' @exportClass SplatParams
setClass("SplatParams",
         contains = "Params",
         slots = c(nBatches = "numeric",
                   batchCells = "numeric",
                   batch.facLoc = "numeric",
                   batch.facScale = "numeric",
                   batch.rmEffect = "logical",
                   mean.shape = "numeric",
                   mean.rate = "numeric",
                   lib.loc = "numeric",
                   lib.scale = "numeric",
                   lib.norm = "logical",
                   out.prob = "numeric",
                   out.facLoc = "numeric",
                   out.facScale = "numeric",
                   nGroups = "numeric",
                   group.prob = "numeric",
                   de.prob = "numeric",
                   de.downProb = "numeric",
                   de.facLoc = "numeric",
                   de.facScale = "numeric",
                   bcv.common = "numeric",
                   bcv.df = "numeric",
                   dropout.type = "character",
                   dropout.mid = "numeric",
                   dropout.shape = "numeric",
                   path.from = "numeric",
                   path.nSteps = "numeric",
                   path.skew = "numeric",
                   path.nonlinearProb = "numeric",
                   path.sigmaFac = "numeric"),
         prototype = prototype(nBatches = 1,
                               batchCells = 100,
                               batch.facLoc = 0.1,
                               batch.facScale = 0.1,
                               batch.rmEffect = FALSE,
                               mean.rate = 0.3,
                               mean.shape = 0.6,
                               lib.loc = 11,
                               lib.scale = 0.2,
                               lib.norm = FALSE,
                               out.prob = 0.05,
                               out.facLoc = 4,
                               out.facScale = 0.5,
                               nGroups = 1,
                               group.prob = 1,
                               de.prob = 0.1,
                               de.downProb = 0.5,
                               de.facLoc = 0.1,
                               de.facScale = 0.4,
                               bcv.common = 0.1,
                               bcv.df = 60,
                               dropout.type = "none",
                               dropout.mid = 0,
                               dropout.shape = -1,
                               path.from = 0,
                               path.nSteps = 100,
                               path.skew = 0.5,
                               path.nonlinearProb = 0.1,
                               path.sigmaFac = 0.8))

#' Distance sampling
#'
#' `Distance` is a simple way to fit detection functions and estimate
#' abundance using distance sampling methodology.
#'
#' Underlying `Distance` is the package `mrds`, for more advanced
#' analyses (such as those involving double observer surveys) one may find it
#' necessary to use `mrds`.
#'
#' Examples of distance sampling analyses are available at
#' <https://distancesampling.org/resources/vignettes.html>.
#'
#' For help with distance sampling and this package, there is a Google Group
#' <https://groups.google.com/forum/#!forum/distance-sampling>.
#' 
#' Bugs can be reported at <https://github.com/DistanceDevelopment/Distance/issues>.
#'
#' @name Distance-package
#' @import mrds
#' @aliases Distance-package Distance
#' @author David L. Miller <dave@@ninepointeightone.net>
#' @references
#' "_PACKAGE"
#'
#' Key References:
#'
#' Miller D.L., E. Rexstad, L. Thomas, L. Marshall and J.L. Laake. 2019.
#'   Distance Sampling in R. Journal of Statistical Software, 89(1), 1-28.
#'   \doi{10.18637/jss.v089.i01}
#'
#' Background References:
#'
#' Laake, J.L. and D.L. Borchers. 2004. Methods for incomplete
#'   detection at distance zero. In: Advanced Distance Sampling, eds. S.T.
#'   Buckland, D.R.Anderson, K.P. Burnham, J.L. Laake, D.L. Borchers, and L.
#'   Thomas. Oxford University Press.
#'
#' Marques, F.F.C. and S.T. Buckland. 2004. Covariate models for the detection
#'   function. In: Advanced Distance Sampling, eds. S.T. Buckland,
#'   D.R.Anderson, K.P. Burnham, J.L. Laake, D.L. Borchers, and L. Thomas.
#'   Oxford University Press.
#'
#' @keywords Statistical Models
#'
NULL

#' The flatfile data format
#'
#' `Distance` allows loading data as a "flat file" and analyse data (and
#' obtain abundance estimates) straight away, provided that the format of the
#' flat file is correct. One can provide the file as, for example, an Excel
#' spreadsheet using [`readxl::read_xls`][readxl::read_xls] in or CSV using
#' [`read.csv`][utils::read.csv].
#'
#' Each row of the data table corresponds to either: (1) an observation or (2)
#' a sample (transect) without observations. In either case the following
#' columns must be present:
#'   * `distance` observed distance to object
#'   * `object` a unique identifier for each observation (only required when
#'   using [`dht2`][dht2])
#'   * `Sample.Label` identifier for the sample (transect id)
#'   * `Effort` effort for this transect (e.g. line transect length or number
#'   of times point transect was visited)
#'   * `Region.Label` label for a given stratum (see below)
#'   * `Area` area of the strata`
#' When the row represents a transect without observations, `distance` and any
#' other observation-specific covariates (including `size` and detection
#' function covariates) take the value `NA`.
#'
#' Note that in the simplest case (one area surveyed only once) there is only
#' one `Region.Label` and a single corresponding `Area` duplicated for each
#' observation.
#'
#' The example given below was provided by Eric Rexstad. Additional examples
#' can be found at  <https://distancesampling.org/resources/vignettes.html>.
#'
#' @name flatfile
#' @docType methods
#' @examples
#' \dontrun{
#' library(Distance)
#' # Need to have the readxl package installed from CRAN
#' require(readxl)
#'
#' # Need to get the file path first
#' minke.filepath <- system.file("minke.xlsx", package="Distance")
#'
#' # Load the Excel file, note that col_names=FALSE and we add column names after
#' minke <- read_xlsx(minke.filepath, col_names=FALSE)
#' names(minke) <- c("Region.Label", "Area", "Sample.Label", "Effort",
#'                   "distance")
#' # One may want to call edit(minke) or head(minke) at this point
#' # to examine the data format
#'
#' ## perform an analysis using the exact distances
#' pooled.exact <- ds(minke, truncation=1.5, key="hr", order=0)
#' summary(pooled.exact)
#'
#'
#' ## Try a binned analysis
#' # first define the bins
#' dist.bins <- c(0,.214, .428,.643,.857,1.071,1.286,1.5)
#' pooled.binned <- ds(minke, truncation=1.5, cutpoints=dist.bins, key="hr",
#'                     order=0)
#'
#' # binned with stratum as a covariate
#' minke$stratum <- ifelse(minke$Region.Label=="North", "N", "S")
#' strat.covar.binned <- ds(minke, truncation=1.5, key="hr",
#'                          formula=~as.factor(stratum), cutpoints=dist.bins)
#'
#' # Stratified by North/South
#' full.strat.binned.North <- ds(minke[minke$Region.Label=="North",],
#'                   truncation=1.5, key="hr", order=0, cutpoints=dist.bins)
#' full.strat.binned.South <- ds(minke[minke$Region.Label=="South",],
#'                      truncation=1.5, key="hr", order=0, cutpoints=dist.bins)
#'
#' ## model summaries
#' model.sel.bin <- data.frame(name=c("Pooled f(0)", "Stratum covariate",
#'                                    "Full stratification"),
#'                             aic=c(pooled.binned$ddf$criterion,
#'                                   strat.covar.binned$ddf$criterion,
#'                                   full.strat.binned.North$ddf$criterion+
#'                                   full.strat.binned.South$ddf$criterion))
#'
#' # Note model with stratum as covariate is most parsimonious
#' print(model.sel.bin)
#' }
NULL


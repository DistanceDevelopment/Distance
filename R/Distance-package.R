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
#' <http://examples.distancesampling.org/>.
#'
#' For help with distance sampling and this package, there is a Google Group
#' <https://groups.google.com/forum/#!forum/distance-sampling>.
#'
#' @name Distance-package
#' @import mrds
#' @aliases Distance-package Distance
#' @docType package
#' @author David L. Miller <dave@@ninepointeightone.net>
#' @references
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


#' Simulated minke whale data
#'
#' Data simulated from models fitted to 1992/1993 Southern Hemisphere minke
#' whale data collected by the International Whaling Commission. See Branch and
#' Butterworth (2001) for survey details (survey design is shown in figure
#' 1(e)). Data simulated by David Borchers.
#'
#' Data are included here as both R data and as an Excel spreadsheet to
#' illustrate the "flat file" input method. See [`flatfile`][flatfile] for how
#' to load this data and an example analysis.
#'
#' @references Branch, T.A. and D.S. Butterworth (2001) Southern Hemisphere
#' minke whales: standardised abundance estimates from the 1978/79 to 1997/98
#' IDCR-SOWER surveys. Journal of Cetacean Research and Management 3(2):
#' 143-174
#'
#' Hedley, S.L., and S.T. Buckland. Spatial Models for Line Transect Sampling.
#' Journal of Agricultural, Biological, and Environmental Statistics 9, no. 2
#' (2004): 181-199. \doi{10.1198/1085711043578}.
#'
#' @name minke
#' @keywords datasets
#' @source Shipped with the Distance for Windows.
#' @docType data
#' @format `data.frame` with 99 observations of 5 variables:
#'   * `Region.Label` stratum label (`"North"` or `"South"`)
#'   * `Area` stratum area
#'   * `Sample.Label` transect identifier
#'   * `Effort` transect length
#'   * `distance` observed distance
#' @examples
#' data(minke)
#' head(minke)
NULL

#' The flatfile data format
#'
#' `Distance` allows loading data as a "flat file" and analyse data (and
#' obtain abundance estimates) straight away, provided that the format of the
#' flat file is correct. One can provide the file as, for example, an Excel
#' spreadsheet using [`readxl::read_xls`][readxl::read_xls] in or CSV using
#' [`read.csv`][utils::read.csv].
#'
#' Each row of the data table corresponds to one observation and must have a
#' the following columns:
#'   * `distance` observed distance to object
#'   * `Sample.Label` Identifier for the sample (transect id)
#'   * `Effort` effort for this transect (e.g. line transect length or number
#'   of times point transect was visited)
#'   * `Region.Label` label for a given stratum (see below)
#'   * `Area` area of the strata`
#'
#' Note that in the simplest case (one area surveyed only once) there is only
#' one `Region.Label` and a single corresponding `Area` duplicated for each
#' observation.
#'
#' The example given below was provided by Eric Rexstad.
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


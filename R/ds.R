#' Fit detection functions and calculate abundance from line or point transect
#' data
#'
#' This function fits detection functions to line or point transect data and
#' then (provided that survey information is supplied) calculates abundance and
#' density estimates. The examples below illustrate some basic types of
#' analysis using `ds()`.
#'
#' @param data a `data.frame` containing at least a column called `distance` or
#' a numeric vector containing the distances.  NOTE!  If there is a column
#' called `size` in the data then it will be interpreted as group/cluster size,
#' see the section "Clusters/groups", below. One can supply data as a "flat
#' file" and not supply `region_table`, `sample_table` and `obs_table`, see
#' "Data format", below and [`flatfile`][flatfile].
#' @param truncation either truncation distance (numeric, e.g. 5) or percentage
#' (as a string, e.g. "15%"). Can be supplied as a `list` with elements `left`
#' and `right` if left truncation is required (e.g.  `list(left=1,right=20)` or
#' `list(left="1%",right="15%")` or even `list(left="1",right="15%")`).  By
#' default for exact distances the maximum observed distance is used as the
#' right truncation. When the data is binned, the right truncation is the
#' largest bin end point. Default left truncation is set to zero.
#' @param transect indicates transect type "line" (default) or "point".
#' @param formula formula for the scale parameter. For a CDS analysis leave
#' this as its default `~1`.
#' @param key key function to use; `"hn"` gives half-normal (default), `"hr"`
#' gives hazard-rate and `"unif"` gives uniform. Note that if uniform key is
#' used, covariates cannot be included in the model.
#' @param adjustment adjustment terms to use; `"cos"` gives cosine (default),
#' `"herm"` gives Hermite polynomial and `"poly"` gives simple polynomial. A
#' value of `NULL` indicates that no adjustments are to be fitted.
#' @param nadj the number of adjustment terms to fit. The default value
#' (`NULL`) will select via AIC (using a sequential forward selection
#' algorithm) up to `max.adjustment` adjustments (unless `order` is specified).
#' A non-negative integer value will cause the specified number of adjustments
#' to be fitted. The order of adjustment terms used will depend on the `key`
#' and `adjustment`. For `key="unif"`, adjustments of order 1, 2, 3, ... are
#' fitted when `adjustment = "cos"` and order 2, 4, 6, ... otherwise. For
#' `key="hn"` or `"hr"` adjustments of order 2, 3, 4, ... are fitted when
#' `adjustment = "cos"` and order 4, 6, 8, ... otherwise. See Buckland et al.
#' (2001, p. 47) for details.
#' @param order order of adjustment terms to fit. The default value (`NULL`)
#' results in `ds` choosing the orders to use - see `nadj`. Otherwise a scalar
#' positive integer value can be used to fit a single adjustment term of the
#' specified order, and a vector of positive integers to fit multiple
#' adjustment terms of the specified orders. For simple and Hermite polynomial
#' adjustments, only even orders are allowed. The number of adjustment terms
#' specified here must match `nadj` (or `nadj` can be the default `NULL` value).
#' @param scale the scale by which the distances in the adjustment terms are
#' divided. Defaults to `"width"`, scaling by the truncation distance. If the
#' key is uniform only `"width"` will be used. The other option is `"scale"`:
#' the scale parameter of the detection
#' @param cutpoints if the data are binned, this vector gives the cutpoints of
#' the bins. Ensure that the first element is 0 (or the left truncation
#' distance) and the last is the distance to the end of the furthest bin.
#' (Default `NULL`, no binning.) Note that if `data` has columns `distbegin`
#' and `distend` then these will be used as bins if `cutpoints` is not
#' specified. If both are specified, `cutpoints` has precedence.
#' @param monotonicity should the detection function be constrained for
#' monotonicity weakly (`"weak"`), strictly (`"strict"`) or not at all
#' (`"none"` or `FALSE`). See Monotonicity, below. (Default `"strict"`). By
#' default it is on for models without covariates in the detection function,
#' off when covariates are present.
#' @param dht_group should density abundance estimates consider all groups to
#' be size 1 (abundance of groups) `dht_group=TRUE` or should the abundance of
#' individuals (group size is taken into account), `dht_group=FALSE`. Default
#' is `FALSE` (abundance of individuals is calculated).
#' @param region_table `data_frame` with two columns:
#'   * `Region.Label` label for the region
#'   * `Area` area of the region
#'   * `region_table` has one row for each stratum. If there is no
#'   stratification then `region_table` has one entry with `Area` corresponding
#'   to the total survey area. If `Area` is omitted density estimates only are
#'   produced.
#' @param sample_table `data.frame` mapping the regions to the samples
#' (i.e. transects). There are three columns:
#'   * `Sample.Label` label for the sample
#'   * `Region.Label` label for the region that the sample belongs to.
#'   * `Effort` the effort expended in that sample (e.g. transect length).
#' @param obs_table `data.frame` mapping the individual observations
#' (objects) to regions and samples. There should be three columns:
#'   * `object` unique numeric identifier for the observation
#'   * `Region.Label` label for the region that the sample belongs to
#'   * `Sample.Label` label for the sample
#' @param convert_units conversion between units for abundance estimation, see
#' "Units", below. (Defaults to 1, implying all of the units are "correct"
#' already.)
#' @param er_var encounter rate variance estimator to use when abundance
#' estimates are required. Defaults to "R2" for line transects and "P3" for
#' point transects. See [`dht2`][dht2] for more information and if more
#' complex options are required.
#' @param method optimization method to use (any method usable by
#' [`optim`][stats::optim] or [`optimx`][optimx::optimx]). Defaults to
#' `"nlminb"`.
#' @param debug_level print debugging output. `0`=none, `1-3` increasing levels
#' of debugging output.
#' @param quiet suppress non-essential messages (useful for bootstraps etc).
#' Default value `FALSE`.
#' @param initial_values a `list` of named starting values, see
#' [`mrds-opt`][mrds::mrds-opt]. Only allowed when AIC term selection is not
#' used.
#' @param max_adjustments maximum number of adjustments to try (default 5) only
#' used when `order=NULL`.
#' @param er_method encounter rate variance calculation: default = 2 gives the method of Innes et al, using expected counts in the encounter rate. Setting to 1 gives observed counts (which matches Distance for Windows) and 0 uses binomial variance (only useful in the rare situation where study area = surveyed area). See [`dht.se`][mrds::dht.se] for more details.
#' @param dht_se should uncertainty be calculated when using `dht`? Safe to leave as `TRUE`, used in `bootdht`.
#' @param dht.group deprecated, see same argument with underscore, above.
#' @param region.table deprecated, see same argument with underscore, above.
#' @param sample.table deprecated, see same argument with underscore, above.
#' @param obs.table deprecated, see same argument with underscore, above.
#' @param convert.units deprecated, see same argument with underscore, above.
#' @param er.var deprecated, see same argument with underscore, above.
#' @param debug.level deprecated, see same argument with underscore, above.
#' @param initial.values deprecated, see same argument with underscore, above.
#' @param max.adjustments deprecated, see same argument with underscore, above.
#' @return a list with elements:
#'   * `ddf` a detection function model object.
#'   * `dht` abundance/density information (if survey region data was supplied,
#'   else `NULL`)
#'
#' @section Details:
#' If abundance estimates are required then the `data.frame`s `region_table`
#' and `sample_table` must be supplied. If `data` does not contain the columns
#' `Region.Label` and `Sample.Label` then the `data.frame` `obs_table` must
#' also be supplied. Note that stratification only applies to abundance
#' estimates and not at the detection function level. Density and abundance
#' estimates, and corresponding estimates of variance and confidence intervals,
#' are calculated using the methods described in Buckland et al. (2001)
#' sections 3.6.1 and 3.7.1 (further details can be found in the documentation
#' for [`dht`][mrds::dht]).
#'
#' For more advanced abundance/density estimation please see the
#' [`dht`][mrds::dht] and [`dht2`][dht2] functions.
#'
#' Examples of distance sampling analyses are available at
#' <http://examples.distancesampling.org/>.
#'
#' Hints and tips on fitting (particularly optimisation issues) are on the
#' [`mrds-opt`][mrds::mrds-opt] manual page.
#'
#' @section Clusters/groups:
#'  Note that if the data contains a column named `size`, cluster size will be
#'  estimated and density/abundance will be based on a clustered analysis of
#'  the data. Setting this column to be `NULL` will perform a non-clustered
#'  analysis (for example if "`size`" means something else in your dataset).
#'
#' @section Truncation:
#' The right truncation point is by default set to be largest observed distance
#' or bin end point. This is a default will not be appropriate for all data and
#' can often be the cause of model convergence failures. It is recommended that
#' one plots a histogram of the observed distances prior to model fitting so as
#' to get a feel for an appropriate truncation distance. (Similar arguments go
#' for left truncation, if appropriate). Buckland et al (2001) provide
#' guidelines on truncation.
#'
#' When specified as a percentage, the largest `right` and smallest `left`
#' percent distances are discarded. Percentages cannot be supplied when using
#' binned data.
#'
#' For left truncation, there are two options: (1) fit a detection function to
#' the truncated data as is (this is what happens when you set `left`).  This
#' does not assume that g(x)=1 at the truncation point. (2) manually remove
#' data with distances less than the left truncation distance -- effectively
#' move the centre line out to be the truncation distance (this needs to be
#' done before calling `ds`). This then assumes that detection is certain at
#' the left truncation distance. The former strategy has a weaker assumption,
#' but will give higher variance as the detection function close to the line
#' has no data to tell it where to fit -- it will be relying on the data from
#' after the left truncation point and the assumed shape of the detection
#' function. The latter is most appropriate in the case of aerial surveys,
#' where some area under the plane is not visible to the observers, but their
#' probability of detection is certain at the smallest distance.
#'
#' @section Binning:
#' Note that binning is performed such that bin 1 is all distances greater or
#' equal to cutpoint 1 (>=0 or left truncation distance) and less than cutpoint
#' 2. Bin 2 is then distances greater or equal to cutpoint 2 and less than
#' cutpoint 3 and so on.
#'
#' @section Monotonicity:
#' When adjustment terms are used, it is possible for the detection function to
#' not always decrease with increasing distance. This is unrealistic and can
#' lead to bias. To avoid this, the detection function can be constrained for
#' monotonicity (and is by default for detection functions without covariates).
#'
#' Monotonicity constraints are supported in a similar way to that described
#' in Buckland et al (2001). 20 equally spaced points over the range of the
#' detection function (left to right truncation) are evaluated at each round
#' of the optimisation and the function is constrained to be either always
#' less than it's value at zero (`"weak"`) or such that each value is
#' less than or equal to the previous point (monotonically decreasing;
#' `"strict"`). See also [`check.mono`][mrds::check.mono].
#'
#' Even with no monotonicity constraints, checks are still made that the
#' detection function is monotonic, see [`check.mono`][mrds::check.mono].
#'
# THIS IS STOLEN FROM mrds, sorry Jeff!
#' @section Units:
#'  In extrapolating to the entire survey region it is important that the unit
#'  measurements be consistent or converted for consistency. A conversion
#'  factor can be specified with the `convert_units` argument. The values of
#'  `Area` in `region_table`, must be made consistent with the units for
#'  `Effort` in `sample_table` and the units of `distance` in the `data.frame`
#'  that was analyzed. It is easiest if the units of `Area` are the square of
#'  the units of `Effort` and then it is only necessary to convert the units of
#'  `distance` to the units of `Effort`. For example, if `Effort` was entered
#'  in kilometres and `Area` in square kilometres and `distance` in metres then
#'  using `convert_units=0.001` would convert metres to kilometres, density
#'  would be expressed in square kilometres which would then be consistent with
#'  units for `Area`. However, they can all be in different units as long as
#'  the appropriate composite value for `convert_units` is chosen. Abundance
#'  for a survey region can be expressed as: `A*N/a` where `A` is `Area` for
#'  the survey region, `N` is the abundance in the covered (sampled) region,
#'  and `a` is the area of the sampled region and is in units of `Effort *
#'  distance`. The sampled region `a` is multiplied by `convert_units`, so it
#'  should be chosen such that the result is in the same units as `Area`.  For
#'  example, if `Effort` was entered in kilometres, `Area` in hectares (100m x
#'  100m) and `distance` in metres, then using `convert_units=10` will convert
#'  `a` to units of hectares (100 to convert metres to 100 metres for distance
#'  and .1 to convert km to 100m units).
#'
#' @section Data format:
#' One can supply `data` only to simply fit a detection function. However, if
#' abundance/density estimates are necessary further information is required.
#' Either the `region_table`, `sample_table` and `obs_table` `data.frame`s can
#' be supplied or all data can be supplied as a "flat file" in the `data`
#' argument. In this format each row in data has additional information that
#' would ordinarily be in the other tables. This usually means that there are
#' additional columns named: `Sample.Label`, `Region.Label`, `Effort` and
#' `Area` for each observation. See [`flatfile`][flatfile] for an example.
#'
#' @section Density estimation:
#' If column `Area` is omitted, a density estimate is generated but note that
#' the degrees of freedom/standard errors/confidence intervals will not match
#' density estimates made with the `Area` column present.
#'
#' @author David L. Miller
#' @seealso [`flatfile`][flatfile], [`AIC.ds`][AIC.ds], [`ds.gof`][ds.gof],
#' [`p_dist_table`][p_dist_table], [`plot.ds`][plot.ds],
#' [`add_df_covar_line`][add_df_covar_line]
#' @export
#'
#' @importFrom stats quantile as.formula
#' @importFrom methods is
#' @references
#' Buckland, S.T., Anderson, D.R., Burnham, K.P., Laake, J.L., Borchers, D.L.,
#' and Thomas, L. (2001). Distance Sampling. Oxford University Press. Oxford,
#' UK.
#'
#' Buckland, S.T., Anderson, D.R., Burnham, K.P., Laake, J.L., Borchers, D.L.,
#' and Thomas, L. (2004). Advanced Distance Sampling. Oxford University Press.
#' Oxford, UK.
#'
#' @examples
#'
#' # An example from mrds, the golf tee data.
#' library(Distance)
#' data(book.tee.data)
#' tee.data <- subset(book.tee.data$book.tee.dataframe, observer==1)
#' ds.model <- ds(tee.data, 4)
#' summary(ds.model)
#' plot(ds.model)
#'
#' \dontrun{
#' # same model, but calculating abundance
#' # need to supply the region, sample and observation tables
#' region <- book.tee.data$book.tee.region
#' samples <- book.tee.data$book.tee.samples
#' obs <- book.tee.data$book.tee.obs
#'
#' ds.dht.model <- ds(tee.data, 4, region_table=region,
#'                    sample_table=samples, obs_table=obs)
#' summary(ds.dht.model)
#'
#' # specify order 2 cosine adjustments
#' ds.model.cos2 <- ds(tee.data, 4, adjustment="cos", order=2)
#' summary(ds.model.cos2)
#'
#' # specify order 2 and 3 cosine adjustments, turning monotonicity
#' # constraints off
#' ds.model.cos23 <- ds(tee.data, 4, adjustment="cos", order=c(2, 3),
#'                    monotonicity=FALSE)
#' # check for non-monotonicity -- actually no problems
#' check.mono(ds.model.cos23$ddf, plot=TRUE, n.pts=100)
#'
#' # include both a covariate and adjustment terms in the model
#' ds.model.cos2.sex <- ds(tee.data, 4, adjustment="cos", order=2,
#'                         monotonicity=FALSE, formula=~as.factor(sex))
#' # check for non-monotonicity -- actually no problems
#' check.mono(ds.model.cos2.sex$ddf, plot=TRUE, n.pts=100)
#'
#' # truncate the largest 10% of the data and fit only a hazard-rate
#' # detection function
#' ds.model.hr.trunc <- ds(tee.data, truncation="10%", key="hr",
#'                         adjustment=NULL)
#' summary(ds.model.hr.trunc)
#'
#' # compare AICs between these models:
#' AIC(ds.model)
#' AIC(ds.model.cos2)
#' AIC(ds.model.cos23)
#'}
ds <- function(data, truncation=ifelse(is.null(cutpoints),
                                     ifelse(is.null(data$distend),
                                            max(data$distance),
                                            max(data$distend)),
                                     max(cutpoints)),
             transect="line",
             formula=~1, key=c("hn", "hr", "unif"),
             adjustment=c("cos", "herm", "poly"),
             nadj=NULL,
             order=NULL, scale=c("width", "scale"),
             cutpoints=NULL, dht_group=FALSE,
             monotonicity=ifelse(formula==~1, "strict", "none"),
             region_table=NULL, sample_table=NULL, obs_table=NULL,
             convert_units=1, er_var=ifelse(transect=="line", "R2", "P3"),
             method="nlminb", quiet=FALSE, debug_level=0,
             initial_values=NULL, max_adjustments=5, er_method=2, dht_se=TRUE,
             # deprecated below here:
             dht.group,
             region.table,
             sample.table,
             obs.table,
             convert.units,
             er.var,
             debug.level,
             initial.values,
             max.adjustments){

  # capture the call
  this_call <- match.call(expand.dots = FALSE)

  # deprecation work around for dsims
  if("er.var" %in% names(this_call)){
    this_call$er_var <- this_call$er.var
    er_var <- er.var
    this_call$er.var <- NULL
    warning(paste0("Argument: er.var is deprecated, check documentation."))
  }
  if("max.adjustments" %in% names(this_call)){
    this_call$max_adjustments <- this_call$max.adjustments
    max_adjustments <- max.adjustments
    this_call$max.adjustments <- NULL
    warning(paste0("Argument: max.adjustments is deprecated, check documentation."))
  }
  # check for deprecated arguments
  .deprecated_args(c("dht.group", "region.table", "sample.table", "obs.table",
                     "convert.units", "er.var", "debug.level", "initial.values",
                     "max.adjustments"), this_call)

  # this routine just creates a call to mrds, it's not very exciting
  # or fancy, it does do a lot of error checking though

  # check the data, format into the correct tables if we have a flat file
  data <- checkdata(data, region_table, sample_table, obs_table, formula)
  region_table <- data$region.table
  sample_table <- data$sample.table
  obs_table    <- data$obs.table
  data         <- data$data

  # setup left and right truncation (width)
  truncation <- get_truncation(truncation, cutpoints, data)
  left <- truncation$left
  width <- truncation$width

  ### binning
  if(is.null(cutpoints)){
    if(any(names(data)=="distend") & any(names(data)=="distbegin")){
      message("Columns \"distbegin\" and \"distend\" in data: performing a binned analysis...")
      binned <- TRUE
      breaks <- sort(unique(c(data$distend, data$distbegin)))
      data$distance <- (data$distend + data$distbegin)/2
    }else{
      binned <- FALSE
      breaks <- NULL
    }
  }else{
    # make sure that the first bin starts 0 or left
    if(!is.null(left)){
      if(cutpoints[1]!=left){
        stop("The first cutpoint must be 0 or the left truncation distance!")
      }
    }else if(cutpoints[1]!=0){
      stop("The first cutpoint must be 0 or the left truncation distance!")
    }

    # remove distbegin and distend if they already exist
    if(any(names(data)=="distend") & any(names(data)=="distbegin")){
      message("data already has distend and distbegin columns, removing them and appling binning as specified by cutpoints.")
      data$distend <- NULL
      data$distbegin <- NULL
    }
    # send off to create_bins to make the correct columns in data
    data <- create_bins(data, cutpoints)
    binned <- TRUE
    breaks <- cutpoints
  }

  # transect type
  point <- switch(transect,
                  "line"=FALSE,
                  "point"=TRUE,
                  stop("Only \"point\" or \"line\" transects may be supplied.")
                 )

  # key and adjustments
  key <- match.arg(key)
  # keep the name for the key function
  key.name <- switch(key,
                     hn   = "half-normal",
                     hr   = "hazard-rate",
                     unif = "uniform"
                    )

  # no covariates with uniform
  if((as.formula(formula)!=~1) & key=="unif"){
    stop("Can't use uniform key with covariates.")
  }
  # uniform key must use width scaling
  scale <- match.arg(scale)
  if(key=="unif"){
    scale <- "width"
  }

  # check we have an allowed adjustment
  if(!is.null(adjustment)){
    adjustment <- match.arg(adjustment)
  }

  # if the user supplied order=0, that's equivalent to adjustment=NULL
  if((!is.null(order) & all(order==0)) | max_adjustments==0){
    adjustment <- NULL
  }

  if(!is.null(adjustment)){
    # by default don't do AIC selection, will turn on later
    aic.search <- FALSE
    if(!is.null(order)){

      # if both nadj and order are supplied
      if(!is.null(nadj)){
        if(length(order) != nadj) stop("The number of adjustment orders specified in 'order' must match 'nadj'")
      }

      if(any(order != ceiling(order))){
          stop("Adjustment orders must be integers.")
      }

      order <- sort(order)

    }else if(!is.null(nadj)){
      order <- get_adj_orders(nadj, key, adjustment)
    }else{
      # if there are covariates then don't do the AIC search
      if(formula != ~1){
        aic.search <- FALSE
        message("Model contains covariate term(s): no adjustment terms will be included.")
      }else{
      # otherwise go ahead and set up the candidate adjustment orders
        aic.search <- TRUE
        order <- get_adj_orders(max_adjustments, key, adjustment)
      }
    }

    # keep the name for the adjustments
    adj.name <- switch(adjustment,
                       cos  = "cosine",
                       herm = "Hermite",
                       poly = "simple polynomial"
                      )
  }else{
    aic.search <- FALSE
  }

  # monotonicity
  if(is.logical(monotonicity)){
    if(!monotonicity){
      mono <- FALSE
      mono.strict <- FALSE
    }
  }else if(monotonicity=="none"){
    mono <- FALSE
    mono.strict <- FALSE
  }else if(monotonicity=="weak"){
    mono <- TRUE
    mono.strict <- FALSE
  }else if(monotonicity=="strict"){
    mono <- TRUE
    mono.strict <- TRUE
  }else{
    stop("monotonicity must be one of \"none\", FALSE, \"weak\" or \"strict\".")
  }

  # can't do monotonicity and covariates, fail!
  if(mono & formula!=as.formula("~1")){
    stop("Monotonicity cannot be enforced with covariates.")
  }

  # set up the control options
  control <- list(optimx.method=method, showit=debug_level)

  # if initial values were supplied, pass them on
  if(!is.null(initial_values) & !aic.search){
    control$initial <- initial_values
  }else if(!is.null(initial_values) & aic.search){
    stop("Cannot supply initial values when using AIC term selection")
  }

  ### Actually fit some models here

  # construct the meta data object...
  meta.data <- list(width = width,point = point,binned = binned,
                    mono=mono, mono.strict=mono.strict)
  if(!is.null(left)){
    meta.data$left <- left
  }
  if(binned){
    meta.data$breaks <- breaks
  }

  # if we are doing an AIC-based search then, create the indices for the
  # for loop to work along, else just give the length of the order object
  if(aic.search){
    for.ind <- c(0, seq(along=order))
    message("Starting AIC adjustment term selection.")
  }else if(!is.null(adjustment)){
    for.ind <- length(order)
  }else{
    for.ind <- 1
  }

  # dummy last model
  last.model <- list(criterion=Inf)

  # loop over the orders of adjustments
  for(i in for.ind){
    # construct model formulae
    # CDS model
    if(formula==as.formula("~1")){
      model.formula <- paste("~cds(key =\"", key,"\", formula = ~1",sep="")
    # MCDS model
    }else{
      model.formula <- paste("~mcds(key = \"",key,"\",",
                                 "formula =~", as.character(formula)[2], sep="")
    }

    # build a message to let the user know what is being fitted
    this.message <- paste("Fitting ", key.name, " key function", sep="")

    # adjustments?
    # this handles the case when we have adjustments but are doing AIC search
    # so want to fit a key function alone to begin with.
    if(!is.null(adjustment) & i!=0){
      if(length(order[1:i])==1){
        order.str <- order[1:i]
      }else{
        order.str <- paste("c(", paste(order[1:i], collapse=","), ")", sep="")
      }

      model.formula <- paste(model.formula, ",",
                           "adj.series=\"", adjustment,
                           "\", adj.order=", order.str, ",",
                           "adj.scale=\"", scale, "\"", sep="")

      this.message <- paste(this.message,
                            " with ", adj.name,"(",
                            paste(order[1:i],collapse=","),
                            ") adjustments", sep="")

      # use the last parameter values as starting values if
      # we are doing AIC search and we're at step 2 and onwards
      if(aic.search && length(order[1:i]) > 1){
        lastpar <- last.model$par
        control$initial <- list()
        if(key == "hr"){
          control$initial$shape <- lastpar[1]
          lastpar <- lastpar[-1]
        }
        if(key != "unif"){
          control$initial$scale <- lastpar[1]
          lastpar <- lastpar[-1]
        }
        if(length(lastpar)>0){
          control$initial$adjustment <- lastpar
        }

        # add a space for a new parameter
        # rep() here to handle the case where the previous iteration
        # failed and we are using initial values from more than 1 model ago
        control$initial$adjustment <- c(control$initial$adjustment,
                                        rep(0,
                                            i -
                                            length(control$initial$adjustment)))
      }
    }

    model.formula <- paste(model.formula,")",sep="")

    # tell the user what is being fitted
    message(this.message)

    # turn-off monotonicity if we have a key only model
    if(i==0){
      mono.save <- meta.data$mono
      mono.strict.save <- meta.data$mono.strict
      meta.data$mono <- FALSE
      meta.data$mono.strict <- FALSE
    }

    # actually fit a model
    # wrap everything around this so we don't print out a lot of useless
    # stuff...
    model <- suppressPackageStartupMessages(
               try(
                                  ddf(dsmodel = as.formula(model.formula),
                                      data = data, method = "ds",
                                      control=control,
                                      meta.data = meta.data), silent=quiet))

    # turn-off monotonicity if we have a key only model
    if(i==0){
      meta.data$mono <- mono.save
      meta.data$mono.strict <- mono.strict.save
    }

    # if that worked
    if(!inherits(model, "try-error")){
      if(model$ds$converge==0){

        model$name.message <- sub("^Fitting ","",this.message)

        # need this to get plotting to work!
        model$call$dsmodel <- as.formula(model.formula)

        message(paste("AIC=",round(model$criterion,3)))

        if(aic.search){
          # if this models AIC is worse (bigger) than the last
          # return the last model and stop looking.
          if(model$criterion > last.model$criterion){
            model <- last.model
            # capitalise!
            model_name <- model$name.message
            model_name <- paste0(toupper(substring(model_name, 1, 1)),
                                substring(model_name, 2))
            message(paste0("\n", model_name, " selected."))
            break
          }else{
            # otherwise keep this, best model
            last.model <- model
          }
        } # end AIC selection fiddling
      }else{
        # if the model didn't converge warn the user
        message("  Model failed to converge.")
        # return the last model
        if(aic.search & !is.infinite(last.model$criterion)){
          model <- last.model
          # capitalise!
          model_name <- model$name.message
          model_name <- paste0(toupper(substring(model_name, 1, 1)),
                              substring(model_name, 2))
          message(paste0("\n", model_name, " selected."))
          break
        }else{
          # if that was the only model just need to return NULL,
          # no last.model to fall back on
          model <- NULL
          break
        }
      } # end model convergence check
    }else{
      if(last.model$criterion == Inf & length(last.model)==1){
        message("\n\nAll models failed to fit!\n")
        model <- NULL
        break
      }else{
        message(paste0("\n\nError in model fitting, returning: ",
                       sub("^Fitting ","",last.model$name.message)))
        message(paste0("\n  Error: ",model[1],"\n"))
        model <- NULL
        model <- last.model
      }
      break
    } # end try error check
  } # end for() over adjustments

  if(is.null(model)){
    stop("No models could be fitted.")
  }

  # check that hazard models have a reasonable scale parameter
  if(key=="hr" && model$par[1] < sqrt(.Machine$double.eps)){
    warning("Estimated hazard-rate scale parameter close to 0 (on log scale). Possible problem in data (e.g., spike near zero distance).")
  }

  # check to see if resulting function is monotonic
  mono.chk <- mrds::check.mono(model, n.pts=20)

  ## Now calculate abundance/density using dht()
  if(!is.null(region_table) & !is.null(sample_table)){
    # if obs_table is not supplied, then data must have the Region.Label and
    # Sample.Label columns

    # setup dht options
    dht_options <- list(group         = dht_group,
                        ervar         = er_var,
                        varflag       = er_method,
                        convert.units = convert_units)

    # if no obs_table
    if(is.null(obs_table)){
      if(all(c("Region.Label", "Sample.Label") %in% names(data))){

        if(any(is.na(model$hessian))){
          message("Some variance-covariance matrix elements were NA, possible numerical problems; only estimating detection function.\n")
          dht.res <- NULL
        }else{
          dht.res <- dht(model, region_table, sample_table,
                         options=dht_options, se=dht_se)
        }
      }else{
        message("No obs_table supplied nor does data have Region.Label and Sample.Label columns, only estimating detection function.\n")
        dht.res <- NULL
      }
    }else{
      # from ?dht:
      # For animals observed in tight clusters, that estimator gives the
      # abundance of groups (group=TRUE in options) and the abundance of
      # individuals is estimated as s_1/p_1 + s_2/p_2 + ... + s_n/p_n, where
      # s_i is the size (e.g., number of animals in the group) of each
      # observation(group=FALSE in options).

      if(any(is.na(model$hessian))){
        message("Some variance-covariance matrix elements were NA, possible numerical problems; only estimating detection function.\n")
        dht.res <- NULL
      }else{
        dht.res <- dht(model, region_table, sample_table, obs_table,
                       options=dht_options, se=dht_se)
      }
    }
  }else{
    # if no information on the survey area was supplied just return
    # the detection function stuff
    dht.res <- NULL

    if(!quiet){
      message("No survey area information supplied, only estimating detection function.\n")
    }
  }

  # construct return object
  ret.obj <- list(ddf  = model,
                  dht  = dht.res,
                  call = this_call)

  # give it some class
  class(ret.obj) <- "dsmodel"

  return(ret.obj)
}

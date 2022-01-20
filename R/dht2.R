#' Abundance estimation for distance sampling models
#'
#' Once a detection function is fitted to data, this function can be used to
#' compute abundance estimates over required areas. The function also allows
#' for stratification and variance estimation via various schemes (see below).
#'
#' @param ddf model fitted by [`ds`][Distance::ds] or [`ddf`][mrds::ddf]
#' @param strat_formula a formula giving the stratification structure (see
#' "Stratification" below). Currently only one level of stratification is
#' supported.
#' @param observations `data.frame` to link detection function data (indexed by
#' `object` column IDs) to the transects (indexed by `Sample.Label` column
#' IDs). See "Data" below.
#' @param transects `data.frame` with information about samples (points or
#' line transects). See "Data" below.
#' @param geo_strat `data.frame` with information about any geographical
#' stratification. See "Data" below.
#' @param flatfile data in the flatfile format, see [`flatfile`][flatfile].
#' @param convert_units conversion factor between units for the distances,
#' effort and area. See "Units" below.
#' @param er_est encounter rate variance estimator to be used. See "Variance"
#' below and [`varn`][mrds::varn].
#' @param multipliers `list` of `data.frame`s. See "Multipliers" below.
#' @param sample_fraction proportion of the transect covered (e.g., 0.5 for
#' one-sided line transects). May be specified as either a single number or a
#' `data.frame` with 2 columns `Sample.Label` and `fraction` (if fractions are
#' different for each transect).
#' @param stratification what do strata represent, see "Stratification" below.
#' @param ci_width for use with confidence interval calculation (defined as
#' 1-alpha, so the default 95 will give a 95% confidence interval).
#' @param innes logical flag for computing encounter rate variance using either
#' the method of Innes et al (2002) where estimated abundance per transect
#' divided by effort is used as the encounter rate, vs. (when `innes=FALSE`)
#' using the number of observations divided by the effort (as in Buckland et
#' al., 2001)
#' @param total_area for options `stratification="effort_sum"` and
#' `stratification="replicate"` the area to use as the total for combined,
#' weighted final estimates.
#' @param binomial_var if we wish to estimate abundance for the covered area
#' only (i.e., study area = surveyed area) then this must be set to be
#' `TRUE` and use the binomial variance estimator of Borchers et al.
#' (1998). This is only valid when objects are not clustered. (This situation
#' is rare.)
#' @return a `data.frame` (of class `dht_result` for pretty printing) with
#' estimates and attributes containing additional information, see "Outputs"
#' for information on column names.
#'
#' @export
#' @importFrom rlang .data
#' @importFrom stats qt na.omit predict terms var qnorm
#' @importFrom dplyr group_by group_by_at mutate ungroup select distinct
#' mutate_if if_else summarize_all "%>%" filter_at inner_join anti_join
#' bind_rows left_join arrange vars
#' @importFrom mrds DeltaMethod
#' @section Data:
#' The data format allows for complex stratification schemes to be set-up. Three
#' objects are always required:
#'   * `ddf` the detection function (see [`ds`][Distance::ds] or
#'   [`ddf`][mrds::ddf] for information on the format of their inputs).
#'   * `observations` has one row per observation and links the observations to
#'   the transects. Required columns:
#'     * `object` (unique ID for the observation, which must match with the
#'     data in the detection function)
#'     * `Sample.Label` (unique ID for the transect).
#'     * Additional columns for strata which are not included in the detection
#'     function are required (stratification covariates that are included in
#'     the detection function do not need to be included here). The important
#'     case here is group size, which must have column name `size` (but does
#'     not need to be in the detection function).
#'   * `transects` has one row per sample (point or line transect). At least
#'   one row is required. Required columns: `Sample.Label` (unique ID for the
#'   transect), `Effort` (line length for line transects, number of visits for
#'   point transects), if there is more than one geographical stratum.
#'
#' With only these three arguments, abundance can only be calculated for the
#' covered area. Including additional information on the area we wish to
#' extrapolate to (i.e., the study area), we can obtain abundance estimates:
#'   * `geo_strat` has one row for each stratum that we wish to estimate
#'   abundance for. For abundance in the study area, at least one row is
#'   required. Required columns: `Area` (the area of that stratum). If there
#'   is >1 row, then additional columns, named in `strat_formula`.`
#'
#' Note that if the `Area` column is set to all 0, then only density estimates
#' will be returned.
#'
#' @section Multipliers:
#' It is often the case that we cannot measure distances to individuals or
#' groups directly, but instead need to estimate distances to something they
#' produce (e.g., for whales, their blows; for elephants their dung) -- this is
#' referred to as indirect sampling. We may need to use estimates of production
#' rate and decay rate for these estimates (in the case of dung or nests) or
#' just production rates (in the case of songbird calls or whale blows). We
#' refer to these conversions between "number of cues" and "number of animals"
#' as "multipliers".
#'
#' The `multipliers` argument is a `list`, with 2 possible elements (`creation`
#' and `decay`). Each element of which is a `data.frame` and must have at least
#' a column named `rate`, which abundance estimates will be divided by (the
#' term "multiplier" is a misnomer, but kept for compatibility with Distance
#' for Windows). Additional columns can be added to give the standard error and
#' degrees of freedom for the rate if known as `SE` and `df`, respectively. You
#' can use a multirow `data.frame` to have different rates for different
#' geographical areas (for example). In this case the rows need to have a
#' column (or columns) to `merge` with the data (for example `Region.Label`).
#'
#' @section Stratification:
#' The `strat_formula` argument is used to specify a column to use to stratify
#' the results, using the form `~column.name` where `column.name` is the column
#' name you wish to use.
#'
#' The `stratification` argument is used to specify which of four types of
#' stratification are intended:
#'  * `"geographical"` if each stratum represents a different geographical
#'  areas and you want the total over all the areas
#'  * `"effort_sum"` if your strata are in fact from replicate
#'  surveys (perhaps using different designs) but you don't have many
#'  replicates and/or want an estimate of "average variance"
#'  * `"replicate"` if you have replicate surveys but have many of them, this
#'  calculates the average abundance and the variance between those many
#'  surveys (think of a population of surveys)
#'  * `"object"` if the stratification is really about the type of object
#'  observed, for example sex, species or life stage and what you want is the
#'  total number of individuals across all the classes of objects. For example,
#'  if you have stratified by sex and have males and females, but also want a
#'  total number of animals, you should use this option.
#'
#' A simple example of using `stratification="geographical"` is given below.
#' Further examples can be found at <http://examples.distancesampling.org/>
#' (see, e.g., the deer pellet survey).
#'
#' @section Variance:
#' Variance in the estimated abundance comes from multiple sources. Depending
#' on the data used to fit the model and estimate abundance, different
#' components will be included in the estimated variances. In the simplest
#' case, the detection function and encounter rate variance need to be
#' combined. If group size varies, then this too must be included. Finally, if
#' multipliers are used and have corresponding standard errors given, this are
#' also included. Variances are combined by assuming independence between the
#' measures and adding variances. A brief summary of how each component is
#' calculated is given here, though see references for more details.
#'   * *detection function*: variance from the detection function parameters is
#'   transformed to variance about the abundance via a sandwich estimator (see
#'   e.g., Appendix C of Borchers et al (2002)).
#'   * *encounter rate*: for strata with >1 transect in them, the encounter
#'   rate estimators given in Fewster et al (2009) can be specified via the
#'   `er_est` argument. If the argument `innes=TRUE` then calculations use the
#'   estimated number of individuals in the transect (rather than the
#'   observed), which was give by Innes et al (2002) as a superior estimator.
#'   When there is only one transect in a stratum, Poisson variance is assumed.
#'   Information on the Fewster encounter rate variance estimators are given in
#'   [`varn`][mrds::varn]
#'   * *group size*: if objects occur in groups (sometimes "clusters"), then
#'   the empirical variance of the group sizes is added to the total variance.
#'   * *multipliers*: if multipliers with standard errors are given, their
#'   corresponding variances are added. If no standard errors are supplied,
#'   then their contribution to variance is assumed to be 0.
#'
#' @section Units:
#' It is often the case that distances are recorded in one convenient set of
#' units, whereas the study area and effort are recorded in some other units.
#' To ensure that the results from this function are in the expected units, we
#' use the `convert_units` argument to supply a single number to convert the
#' units of the covered area to those of the study/stratification area (results
#' are always returned in the units of the study area). For line transects, the
#' covered area is calculated as `2 * width * length` where `width` is the
#' effective (half)width of the transect (often referred to as w in the
#' literature) and `length` is the line length (referred to as L). If `width`
#' and `length` are measured in kilometres and the study area in square
#' kilometres, then all is fine and `convert_units` is 1 (and can be ignored).
#' If, for example, line length and distances were measured in metres, we
#' instead need to convert this to be kilometres, by dividing by 1000 for each
#' of distance and length, hence `convert_units=1e-6`. For point transects,
#' this is slightly easier as we only have the radius and study area to
#' consider, so the conversion is just such that the units of the truncation
#' radius are the square root of the study area units.
#'
#' @section Output:
#' On printing the output from call to `dht2`, three tables are produced. Below is a guide to the output columns names, per table.
#'
#' - Summary statistics table
#'   - `Region.Label` Stratum name (this first column name depends on the `formula` supplied)
#'   - `Area` Size of stratum
#'   - `CoveredArea` Surveyed area in stratum (2 x w x L)
#'   - `Effort` Transect length or number of point visits per stratum
#'   - `n` Number of detections
#'   - `k` Number of replicate transects
#'   - `ER` Encounter rate
#'   - `se.ER` Standard error of encounter rate
#'   - `cv.ER` Coefficient of variation of encounter rate
#' - Abundance or density estimates table:
#'   - `Region.Label` As above
#'   - `Estimate` Point estimate of abundance or density
#'   - `se` Standard error
#'   - `cv` Coefficient of variation
#'   - `LCI` Lower confidence bound
#'   - `UCI` Upper confidence bound
#'   - `df` Degrees of freedom used for confidence interval computation
#' - Components percentage of variance:
#'   - `Region.Label` As above
#'   - `Detection` Percent of variance in abundance/density associated with
#'   detection function uncertainty
#'   - `ER` Percent of variance in abundance/density associated with
#'   variability in encounter rate
#'   - `Multipliers` Percent of variance in abundance/density associated with
#'   uncertainty in multipliers
#'
#' @references
#'
#' Borchers, D.L., S.T. Buckland, P.W. Goedhart, E.D. Clarke, and S.L. Hedley.
#' 1998. Horvitz-Thompson estimators for double-platform line transect surveys.
#' *Biometrics* 54: 1221-1237.
#'
#' Borchers, D.L., S.T. Buckland, and W. Zucchini. 2002 *Estimating Animal
#' Abundance: Closed Populations*. Statistics for Biology and Health. Springer
#' London.
#'
#' Buckland, S.T., E.A. Rexstad, T.A. Marques, and C.S. Oedekoven. 2015
#' *Distance Sampling: Methods and Applications*. Methods in Statistical
#' Ecology. Springer International Publishing.
#'
#' Buckland, S.T., D.R. Anderson, K. Burnham, J.L. Laake, D.L. Borchers, and L.
#' Thomas. 2001 *Introduction to Distance Sampling: Estimating Abundance of
#' Biological Populations*. Oxford University Press.
#'
#' Innes, S., M. P. Heide-Jorgensen, J.L. Laake, K.L. Laidre, H.J. Cleator, P.
#' Richard, and R.E.A. Stewart. 2002 Surveys of belugas and narwhals in the
#' Canadian high arctic in 1996. *NAMMCO Scientific Publications* 4, 169-190.
#'
#' @name dht2
#' @examples
#' \dontrun{
#' # example of simple geographical stratification
#' # minke whale data, with 2 strata: North and South
#' data(minke)
#' # first fitting the detection function
#' minke_df <- ds(minke, truncation=1.5, adjustment=NULL)
#' # now estimate abundance using dht2
#' # stratum labels are in the Region.Label column
#' minke_dht2 <- dht2(minke_df, flatfile=minke, stratification="geographical",
#'                    strat_formula=~Region.Label)
#' # could compare this to minke_df$dht and see the same results
#' minke_dht2
#' # can alternatively report density
#' print(minke_dht2, report="density")
#'}
dht2 <- function(ddf, observations=NULL, transects=NULL, geo_strat=NULL,
                 flatfile=NULL,
                 strat_formula, convert_units=1, er_est=c("R2", "P2"),
                 multipliers=NULL, sample_fraction=1,
                 ci_width=0.95, innes=FALSE,
                 stratification="geographical",
                 total_area=NULL, binomial_var=FALSE){

  # just get the ds model if we have Distance::ds output
  if(inherits(ddf, "dsmodel")){
    ddf <- ddf$ddf
  }

  # get default variance estimation
  if(all(er_est == c("R2", "P2"))){
    if(ddf$ds$aux$point){
      er_est <- "P2"
    }else{
      er_est <- "R2"
    }
  }

  # check we have a valid stratification option
  if(!(stratification %in% c("geographical", "object",
                             "effort_sum", "replicate"))){
    stop("'stratification' must be one of: \"geographical\", \"object\", \"effort_sum\" or \"replicate\"")
  }

  # check ci width
  if(ci_width > 1 | ci_width < 0){
    stop("ci_width must be between 0 and 1!")
  }else{
    ci_width <- (1-ci_width)/2
  }

  # can't do both innes and binomial var
  if(innes & binomial_var){
    stop("Only 'innes' or 'binomial_var' may be TRUE")
  }

  # grouped estimation
  # time to recurse
  if(!is.null(ddf$data$size) & !all(ddf$data$size==1) ){
    mc <- match.call(expand.dots = FALSE)
    dddf <- ddf
    dddf$data$size <- 1
    mc$ddf <- dddf
    if(!is.null(flatfile)){
      ff <- flatfile
      ff$size <- 1
      mc$flatfile <- ff
    }
    grouped <- eval.parent(mc)
    rm(dddf, mc)
  }else{
    grouped <- NULL
  }
  # ^^ we'll save this output later on

  # holder of transect type
  transect_type <- if(ddf$ds$aux$point) "point" else "line"

  # what are the stratum labels specicied in strat_formula?
  stratum_labels <- attr(terms(strat_formula), "term.labels")

  # TODO: currently break if >1 stratum is defined
  # https://github.com/DistanceDevelopment/Distance/issues/46
  if(length(stratum_labels) > 1){
    stop("Only one level of stratification is currently supported")
  }

  if(!is.null(observations) & !is.null(transects)){
    if(!is.null(geo_strat)){
      # what if there were as.factor()s in the formula?
      geo_strat <- safe_factorize(strat_formula, geo_strat)
      # which occur at the geo level?
      geo_stratum_labels <- stratum_labels[stratum_labels %in%
                                           colnames(geo_strat)]
    }else{
      geo_stratum_labels <- NULL
    }


    # what if there were as.factor()s in the formula?
    transects <- safe_factorize(strat_formula, transects)
    observations <- safe_factorize(strat_formula, observations)

    # TODO: do some data checking at this point
    # - check duplicate column names (e.g., obs$sex and df$data$sex)
    # - check column names don't conflict with columns created below
    #    - list of protected column names that can't be in the data
    #         protected <- c("p", ".Label", )

    # check the data
    dht2_checkdata(ddf, observations, transects, geo_strat, strat_formula,
                   stratum_labels, geo_stratum_labels)

    # drop unused levels of factors
    ddf$data <- droplevels(ddf$data)
    observations <- droplevels(observations)
    transects <- droplevels(transects)
    geo_strat <- droplevels(geo_strat)

    # prepare data
    obj_keep <- ddf$data$object[ddf$data$distance <= ddf$ds$aux$width &
                                ddf$data$distance >= ddf$ds$aux$left]
    bigdat <- ddf$data[ddf$data$object %in% obj_keep, ]
    observations <- observations[observations$object %in% obj_keep, ]


    # get probabilities of detection
    bigdat$p <- predict(ddf)$fitted

    bigdat <- merge(bigdat, observations, all.x=TRUE, by="object",
                    suffixes=c("DUPLICATE", ""))

    # remove column duplicates
    if(any(grepl("DUPLICATE", names(bigdat)))){
      bigdat[, grepl("DUPLICATE", names(bigdat))] <- NULL
    }

    # merge onto transects
    join_labs <- intersect(names(bigdat), names(transects))
    join_labs <- join_labs[c("Sample.Label", geo_stratum_labels) %in% join_labs]
    bigdat <- merge(bigdat, transects, all.x=TRUE, all.y=TRUE,
                    by=join_labs,
                    suffixes=c("DUPLICATE", ""))
    # remove Sample.Label dupes
    if(any(grepl("DUPLICATE", names(bigdat)))){
      bigdat[, grepl("DUPLICATE", names(bigdat))] <- NULL
    }

    # merge on the geographical strata
    if(!is.null(geo_strat)){
      # if we have ~1 we ignore stratification
      if(strat_formula==~1){
        geo_strat$Area <- sum(geo_strat$Area)
        geo_strat$.Label <- "Total"
        geo_strat <- unique(geo_strat[, c(".Label", "Area")])

        bigdat$.Label <- "Total"
        stratum_labels <- ".Label"
        geo_stratum_labels <- ".Label"
      }
      # TODO: do something here with strat_formula
      bigdat <- merge(bigdat, geo_strat, all.x=TRUE, by=geo_stratum_labels)
    }else{
      bigdat$Area <- NA
      bigdat$Label <- "Total"
      stratum_labels <- c("Label", stratum_labels)
    }
  }else if(!is.null(flatfile)){
    # if we have a flatfile
    # TODO: check flatfile format here
    # should this just run Distance:::checkdata(flatfile)?

    # make a dummy distance column for use later on
    # this overwrites the column that's there BUT that's okay
    # since we need to make sure it's consistent with the bins
    if(is.null(flatfile$distance)){
      if(!(c("distend", "distbegin") %in% names(flatfile))){
        stop("flatfile must include columns named either 'distance' or 'distbegin' and 'distend'")
      }
      flatfile$distance <- (flatfile$distend+flatfile$distbegin)/2
    }

    # check regular columns exist
    flatfile_labels <- c("distance", "Sample.Label", "Effort", "Area")

    if(!all(flatfile_labels %in% names(flatfile))){
      stop(paste("Column(s):",
                 paste(flatfile_labels[!(flatfile_labels %in% names(flatfile))],
                       collapse=", "),
                 "not in `flatfile`"))
    }
    # safely truncate the data, respecting the data structure
    flatfile <- safetruncate(flatfile, ddf$meta.data$width, ddf$meta.data$left)
    bigdat <- flatfile

    # what if there were as.factor()s in the formula?
    bigdat <- safe_factorize(strat_formula, bigdat)

    # check strat columns are in the data
    if(!all(stratum_labels %in% names(bigdat))){
      stop(paste("Column(s):",
                 paste(stratum_labels[!(stratum_labels %in% names(bigdat))],
                       collapse=", "),
                 "not in `flatfile`"))
    }

    # make object column if not present
    if(is.null(bigdat$object)){
      # ensure that there isn't a size in the data if this is a
      # placeholder row for a sample unit
      bigdat$object <- NA
      bigdat$object[!is.na(bigdat$distance)] <- 1:sum(!is.na(bigdat$distance))
    }else if(!all(is.na(bigdat$distance) == is.na(bigdat$object))){
      stop("NAs in distance column do not match those in the object column, check data")
    }
    # sort by object ID
    bigdat <- bigdat[order(bigdat$object), ]
    # remove the rows where there were no observations
    bigdat_nona <- bigdat[!is.na(bigdat$object), ]
    # get probabilities of detection
    pp <- predict(ddf, newdata=bigdat_nona, compute=TRUE)$fitted
    bigdat$p <- NA
    bigdat$p[!is.na(bigdat$object)] <- pp

    if(strat_formula==~1){
      bigdat$Area <- sum(unique(bigdat[, c("Area", "Region.Label")])$Area)
      bigdat$Region.Label <- NULL
      bigdat$.Label <- "Total"

      stratum_labels <- ".Label"
      geo_stratum_labels <- ".Label"
    }else{
      geo_stratum_labels <- c()
    }

    # TODO: check when implementing multiple stratification
    # https://github.com/DistanceDevelopment/Distance/issues/46
    if(stratification=="object"){
      # what are all the possible combinations of obs level stratum
      # levels and sample labels?
      ex <- expand.grid(lapply(bigdat[, c("Sample.Label", stratum_labels)],
                               function(x) unique(na.omit(x))),
                        stringsAsFactors=FALSE)
      # which are not represented in the data?
      aj <- anti_join(ex, bigdat, by=c("Sample.Label", stratum_labels))
      # join the unrepresented sample combinations to the extra cols
      # (i.e., add Area, Effort data to aj)
      aj <- left_join(aj, unique(bigdat[,c("Sample.Label","Effort","Area")]),
                      by="Sample.Label")

      # remove the transects with no stratum data
      bigdat2 <- filter_at(bigdat, stratum_labels, function(x) !is.na(x))

      # rbind that onto the original data
      bigdat <- bind_rows(bigdat2, aj)
    }

    # TODO: this needs to be checked for the multi-strata case
    # https://github.com/DistanceDevelopment/Distance/issues/46
    # check that Area-stratum combinations are coherent
    ucomb <- unique(bigdat[, c("Area", stratum_labels)])
    if(length(na.omit(ucomb[,stratum_labels])) >
       length(na.omit(unique(bigdat[,stratum_labels])))){
      stop(">1 Area defined for a single stratum label, check data")
    }

  }else{
    stop("Need to supply either observations, transects and geo_strat OR flatfile")
  }

  # stop if any of the transects has zero or negative effort
  if(any(is.na(bigdat[["Effort"]])) || any(bigdat[["Effort"]] <= 0)){
    stop("Some transects have <=0 or NA Effort")
  }

  # handle multipliers
  bigdat <- dht2_multipliers(multipliers, bigdat)
  mult <- TRUE
  if(attr(bigdat, "multipliers")) mult <- TRUE else mult <- FALSE

  # make group size 1 if not in the data
  if(is.null(bigdat$size)){
    # ensure that there isn't a size in the data if this is a
    # placeholder row for a sample unit
    bigdat$size <- NA
    bigdat$size[!is.na(bigdat$distance)] <- 1
  }
  # make object column if not present
  if(is.null(bigdat$object)){
    # ensure that there isn't a size in the data if this is a
    # placeholder row for a sample unit
    bigdat$object <- NA
    bigdat$object[!is.na(bigdat$distance)] <- 1:sum(!is.na(bigdat$distance))
  }else if(!all(is.na(bigdat$distance) == is.na(bigdat$object))){
    stop("NAs in distance column do not match those in the object column, check data")
  }

  # now do some calculations
  bigdat$Nhat <- bigdat$size/bigdat$p

  df_width <- ddf$ds$aux$width*convert_units
  df_left <- ddf$ds$aux$left*convert_units

  # handle sample fractions
  bigdat <- dht2_sample_fraction(sample_fraction, bigdat)

  # TODO: clean-up the data, removing stratum labels with zero observations
  #       or NA labels (and warning)

  # dplyr cheatsheet:
  # - group_by : subsequent commands operate per group
  # - mutate   : adds a new column
  # - distinct : select only the unique row combinations
  # - select   : retain only these columns
  # - .data$   : get the column from the current data, don't go to global env
  # note that this all turns out to be non-standard dplyr code because
  # you can't have unquoted variable names (NSE) in CRAN-submitted packages
  # see https://dplyr.tidyverse.org/articles/programming.html for "why"s
  # so we use rlang::.data to get around this

  # first do transect level calculations
  res <- bigdat %>%
    group_by_at(vars("Sample.Label"))

  # make sure summaries are made at the sample-stratum level here
  # (since sample labels may be duplicated between strata)
  res <- res %>%
    group_by_at(.vars=stratum_labels, .add=TRUE)

  res <- res %>%
      # *observations* per transect
      mutate(transect_n = sum(.data$size, na.rm=TRUE),
             transect_n_observations = length(na.omit(unique(.data$object))),
             # abundance estimate per transect in covered area
             transect_Nc = sum(.data$Nhat, na.rm=TRUE)) %>%
    ungroup()
  # undo second grouping
  if(stratification=="object"){
    res <- res %>%
      ungroup()
  }

  # save the sample-level stats
  res_sample <- res
  res_sample$distance <- res_sample$object <- NULL
  res_sample <- unique(res_sample)

  ## now at at the stratum level, calculate the abundance in the covered area,
  ## number of observations
  res <- res %>%
    group_by_at(.vars=stratum_labels) %>%
      # calculate various summary stats
      mutate(
             # individuals and observations per stratum
             n_individuals  = sum(.data$size, na.rm=TRUE),
             n_observations = length(na.omit(unique(.data$object))),
             # abundance estimate per stratum in covered area
             Nc             = sum(.data$Nhat, na.rm=TRUE),
             # covered area per transect
             Covered_area   = area_calc(df_width, .data$Effort,
                                         transect_type, .data$sample_fraction),
             # get group size stats
             group_var      = if_else(.data$n_observations>1,
                                      var(.data$size, na.rm=TRUE)/
                                       sum(!is.na(.data$size)),
                                      0),
             group_mean     = mean(.data$size, na.rm=TRUE)) %>%
      # report n as n_observations
      mutate(n = .data$n_observations)
  # if we didn't have any areas, then set to 1 and estimate density
  est_density <- FALSE
  if(all(res$Area == 0) | all(is.na(res$Area))){
    res$Area <- 1
    est_density <- TRUE
  }

# TODO: make this more elegant
# TODO: include cue count summary info?
if(mult){
  res <- res %>%
    mutate(Nc_cuecorrected = .data$Nc/.data$rate)
}else{
  res <- res %>%
    mutate(Nc_cuecorrected = NA)
}

  # detection function uncertainty
  # do this sans-pipe, as we need cross-terms and to get the matrix
  df_unc <- varNhat(res, ddf)
  df_Nhat_unc <- df_unc$Nhat

  # extract groupings
  vardat <- attr(df_unc, "vardat_str")
  vardat$.rows <- NULL
  vardat$df_var <- diag(df_Nhat_unc$variance)

  # detection function p uncertainty
  ddf_summary <- summary(ddf)
  vardat$p_var <- ddf_summary$average.p.se[1,1]^2
  vardat$p_average <- ddf_summary$average.p

  # we interrupt your regularly-scheduled grouping to bring you...
  # detection function uncertainty
  res <- res %>% ungroup()

  # merge variances back into res
  res <- merge(res, vardat, all.x=TRUE)

  # normal programming resumed
  res <- res %>%
    group_by_at(.vars=stratum_labels) %>%
      # now we're done with observation level stuff, so remove the columns
      # that we don't need
      select(!!stratum_labels, "Sample.Label", "Area", "n", "Nc", "transect_n",
             "Effort", "Covered_area", "df_var", "transect_Nc", "group_var",
             "group_mean", "Nc_cuecorrected", "rate_var",  "rate", "rate_df",
             "rate_CV", "p_var", "p_average", "transect_n_observations") %>%
      # keep only unique rows
      distinct()

# TODO: fix
if(mult){
  res <- res %>%
    mutate(Nc = .data$Nc_cuecorrected) %>%
    mutate(transect_Nc = .data$transect_Nc/.data$rate)
}
  # save ungrouped version for later calculations
  resT <- ungroup(res)

  # calculate ER variance
  res <- ER_var_f(res, innes=innes, er_est=er_est, binomial_var=binomial_var)

  # calculate final summaries
  res <- res %>%
    # number of transects, total effort and covered area per stratum
    mutate(k            = length(.data$Sample.Label),
           Effort       = sum(.data$Effort),
           Covered_area = sum(.data$Covered_area)) %>%
      mutate(group_var_Nhat = (.data$Area/.data$Covered_area *
                               .data$Nc)^2*
                              .data$group_var/.data$group_mean^2) %>%
      mutate(rate_var_Nhat = (.data$Area/.data$Covered_area *
                               .data$Nc)^2*
                              .data$rate_var/.data$rate^2) %>%
    ## keep only these columns
    select(!!stratum_labels, "Area", "Nc", "n", "ER_var", "Effort", "k",
           "Covered_area", "df_var", "group_var", "group_mean",
           "group_var_Nhat", "ER_var_Nhat", "rate_var", "rate_var_Nhat", "rate",
           "rate_df", "rate_CV", "p_var", "p_average") %>%
    ## now just get the distinct cases
    distinct()

  # we calculated n differently above, so reconcile this in the
  # encounter rate calculation
  if(stratification=="object"){
    res <- res %>%
      mutate(ER = sum(.data$n)/.data$Effort)
  }else{
    res <- res %>%
      mutate(ER = .data$n/.data$Effort)
  }
  res <- res %>%
    mutate(ER_CV = sqrt(.data$ER_var)/.data$ER,
           ER_df = compute_df(.data$k, type=er_est)) %>%
    # calculate stratum abundance estimate
    mutate(Abundance = (.data$Area/.data$Covered_area) * .data$Nc) %>%
    mutate(df_CV = sqrt(.data$df_var)/.data$Abundance) %>%
    mutate(group_CV = if_else(.data$group_var==0, 0,
                              sqrt(.data$group_var)/.data$group_mean))

  # se and CV
  res <- res %>%
    # first add the ER+group size and detfct variance components
    mutate(Abundance_CV = sqrt(sum(c(.data$ER_var_Nhat,
                                     .data$df_var),
                               na.rm=TRUE))/
                           .data$Abundance) %>%
    # now add in the multiplier rate CV
    mutate(Abundance_CV = sqrt(sum(c(.data$Abundance_CV^2,
                                     .data$group_CV^2,
                                     .data$rate_CV^2),
                                   na.rm=TRUE))) %>%
    mutate(Abundance_se = .data$Abundance_CV * .data$Abundance) %>%
    distinct()

  # total degrees of freedom and CI calculation
  if(binomial_var){
    # normal approximation for binomial_var
    res <- res %>%
      mutate(bigC = exp((abs(qnorm(ci_width)) *
                     sqrt(log(1 + .data$Abundance_CV^2))))) %>%
      mutate(df = 0)
  }else{
    res <- res %>%
      mutate(df = .data$Abundance_CV^4/
                    sum(c(if_else(.data$k==1, 0, .data$ER_CV^4/.data$ER_df),
                          .data$df_CV^4/(length(ddf$fitted) - length(ddf$par)),
                          .data$group_CV^4/(.data$n-1),
                          .data$rate_CV^4/.data$rate_df),
                     na.rm=TRUE)) %>%
      # adjust if df is too small
      mutate(df = if_else(.data$Abundance_CV==0, 1, .data$df)) %>%
      mutate(df = if_else(.data$df<1, 1, .data$df)) %>%
      # big C for Satterthwaite
      mutate(bigC = exp((abs(qt(ci_width, .data$df)) *
                     sqrt(log(1 + .data$Abundance_CV^2)))))
  }

  # actually calculate the CIs
  res <- res %>%
    mutate(LCI = if_else(.data$Abundance_CV==0,
                         .data$Abundance, .data$Abundance / .data$bigC),
           UCI = if_else(.data$Abundance_CV==0,
                         .data$Abundance, .data$Abundance * .data$bigC)) %>%
    # done!
    ungroup() %>%
    distinct()


  # make a summary
  res <- as.data.frame(res)
  # TODO: this is untested for multiple strata
  # https://github.com/DistanceDevelopment/Distance/issues/46
  # here we loop over the different stratification variables in the below
  # terminology, "row" indicates a summary row (so for a given stratification
  # variable, we have mutiple values ("strata": North, South etc) and then
  # summarize to one row total at the end. This is for generalizability for to
  # multiple stratification variables later)
  for(this_stratum in stratum_labels){
    dat_row <- res

    # remove NA labels
    iind <- which(is.na(dat_row[, this_stratum]))
    if(length(iind)>0) dat_row <- dat_row[-iind, ]

    # save row names
    stra_row <- dat_row[, stratum_labels, drop=FALSE]
    # remove labels
    dat_row[, stratum_labels] <- NULL

    # don't do anything unless this stratum has rows attached to it
    if(length(stra_row[[this_stratum]]) > 1){
      this_stra_row <- stra_row[1, , drop=FALSE]
      this_stra_row[] <- NA
      this_stra_row[[this_stratum]] <- "Total"

      dat_row <- dat_row %>%
        mutate(n  = sum(.data$n))

      #### stratification options
      ## in Distance for Windows these are in terms of density
      if(stratification=="geographical"){
        # "weighting by area"
        #  which is adding up the abundances
        dat_row <- dat_row %>%
          # for the density case weight by covered area
          mutate(weight       = if_else(rep(est_density, nrow(dat_row)),
                                        .data$Covered_area/
                                         sum(.data$Covered_area), 1),
                 Area         = sum(.data$Area),
                 Covered_area = sum(.data$Covered_area),
                 Effort       = sum(.data$Effort),
                 k            = sum(.data$k)) %>%
          # now summarize ER variance and degrees of freedom
          mutate(ER_var       = sum(.data$weight^2 * .data$ER_var,
                                    na.rm=TRUE)) %>%
          mutate(ER_var_Nhat  = sum(.data$weight^2 * .data$ER_var_Nhat,
                                    na.rm=TRUE)) %>%
          mutate(ER_df = .data$ER_var_Nhat^2/
                         sum((res$ER_var_Nhat^2/.data$ER_df)))
      }else if(stratification %in% c("effort_sum", "replicate")){
        # check that all areas are the same value
        if(length(unique(dat_row$Area))>1 &
           is.null(total_area)){
          stop(paste0("More than 1 Area value in data, need a single Area for stratification=\"",
                      stratification, "\", fix or supply \"total_area\""))
        }
        # if the user didn't supply total_area, but the areas are the same
        # use that as the area
        if(is.null(total_area)){
          total_area <- dat_row$Area[1]
          xt <- 1
        }else{
          xt <- total_area/dat_row$Area
        }
        if(stratification == "replicate"){
          xt <- 1
        }

        dat_row <- dat_row %>%
          mutate(weight       = xt * .data$Effort/sum(.data$Effort)) %>%
          mutate(ER_var       = sum((.data$Effort/sum(.data$Effort))^2*.data$ER_var,
                                    na.rm=TRUE)) %>%
          mutate(ER_var_Nhat  = sum(.data$weight^2*.data$ER_var_Nhat,
                                    na.rm=TRUE)) %>%
          mutate(ER_df        = .data$ER_var_Nhat^2/sum(
                                  ((.data$weight^2 * res$ER_var_Nhat)^2/
                                   .data$ER_df))) %>%
          mutate(Area         = total_area,
                 Covered_area = sum(.data$Covered_area),
                 Effort       = sum(.data$Effort),
                 k            = sum(.data$k))
      }else if(stratification=="object"){
        # things you want to add up like object type
        dat_row <- dat_row %>%
          mutate(ER_var       = sum(.data$ER_var, na.rm=TRUE)) %>%
          mutate(ER_var_Nhat  = sum(.data$ER_var_Nhat, na.rm=TRUE)) %>%
          mutate(weight       = 1,
                 Covered_area = .data$Covered_area[1],
                 Area         = .data$Area[1],
                 Effort       = .data$Effort[1],
                 k            = .data$k[1]) %>%
          mutate(ER_df        = .data$ER_var_Nhat^2/sum((res$ER_var_Nhat^2/
                                                         .data$ER_df)))
      }

      # calculate other summaries
      dat_row <- dat_row %>%
        mutate(ER         = .data$n/.data$Effort) %>%
        mutate(ER_CV      = if_else(.data$ER==0,
                                    0,
                                    sqrt(.data$ER_var)/.data$ER)) %>%
        mutate(group_mean     = mean(.data$group_mean),
               group_var      = sum(.data$group_var),
               group_var_Nhat = sum(.data$group_var_Nhat)) %>%
        mutate(group_CV   = if_else(.data$group_var==0, 0,
                                    sqrt(.data$group_var)/.data$group_mean))

      # calculate mean abundance, unless we are doing replicate, where we need
      # the per-stratum abundances for variance later
      if(stratification != "replicate"){
        dat_row <- dat_row %>%
          mutate(Abundance  = sum(.data$weight*.data$Abundance))
      }

      # calculate total variance for detection function
      vcov <- df_Nhat_unc$variance
      df_tvar <- matrix(dat_row$weight, nrow=1) %*%
                  vcov %*%
                  matrix(dat_row$weight, ncol=1)
      dat_row <- dat_row %>%
        mutate(df_CV  = sqrt(df_tvar[1, 1])/dat_row$Abundance[1])


      # calculate total variance
      if(stratification=="replicate"){
        # Buckland 2001, 3.84-3.87
        nrep <- nrow(res)
        # get "between" variance (empirical variance of strata)
        rvar <- sum(dat_row$weight*(dat_row$Abundance -
                     sum(dat_row$weight*dat_row$Abundance))^2)/
                    (sum(dat_row$weight) * (nrep-1))
        # now calculate abundance
        dat_row <- dat_row %>%
          mutate(Abundance  = sum(.data$weight*.data$Abundance),
                 ER_df      = nrep-1)
        # add the detection function variance
        tvar <- rvar + df_tvar
      }else if(stratification=="effort_sum"){
        # add the pre-weighted CVs
        tvar <- dat_row$Abundance[1]^2 *
                 sum(c(dat_row$ER_CV[1]^2,
                       dat_row$df_CV[1]^2,
                       dat_row$rate_CV[1]^2,
                       dat_row$group_CV[1]^2),
                     na.rm=TRUE)
      }else{
        # add sources of variance
        tvar <- dat_row$ER_var_Nhat[1] +
                df_tvar +
                dat_row$Abundance[1]^2*dat_row$rate_CV[1]^2
        # add-in group size component if not doing Innes et al
        if(!innes){
          tvar <- dat_row$group_var_Nhat[1] +
                  tvar
        }
      }

      dat_row <- dat_row %>%
        mutate(Nc     = sum(.data$weight*.data$Nc),
               df_var = df_tvar[1,1]) %>%
        mutate(Abundance_se = sqrt(tvar)) %>%
        mutate(Abundance_CV = .data$Abundance_se/.data$Abundance) %>%
        mutate(df   = NA,
               bigC = NA,
               LCI  = NA,
               UCI  = NA)
      # drop weights as we don't need them any more
      dat_row$weight <- NULL
      # drop unwanted rows
      dat_row <- dat_row %>%
        distinct()

      # compute degrees of freedom
      if(stratification == "replicate"){
        dat_row <- dat_row %>%
          mutate(wtcv = sum(c((sqrt(rvar)/.data$Abundance[1])^4/
                               .data$ER_df[1],
                              (df_tvar/.data$Abundance[1]^2)^2/
                                (length(ddf$fitted) - length(ddf$par)),
                              if_else(.data$df==0, 0 ,
                                      (.data$rate_var_Nhat[1]/
                                       .data$Abundance[1]^2)^2/
                                        .data$rate_df[1]),
                              (.data$group_var_Nhat[1]/
                               .data$Abundance[1]^2)^2/
                               (length(ddf$fitted)-1)
                             ),
                            na.rm=TRUE)) %>%
          # calculate Satterthwaite df
          mutate(df = sum(c((sqrt(rvar)/.data$Abundance[1])^2,
                            (df_tvar/.data$Abundance[1]^2),
                            if_else(.data$df==0, 0 ,
                                    (.data$rate_var_Nhat[1]/
                                     .data$Abundance[1]^2)),
                            (.data$group_var_Nhat[1]/.data$Abundance[1]^2)
                           ),
                          na.rm=TRUE)^2) %>%
          mutate(df = .data$df/.data$wtcv) %>%
          #mutate(df = .data$ER_df + (nrep - length(ddf$par))) %>%
          mutate(bigC = exp((abs(qt(ci_width, .data$df)) *
                           sqrt(log(1 + .data$Abundance_CV^2)))))
          # drop weight column
          dat_row$wtcv <- NULL
      }else{
        if(binomial_var){
          # normal approximation for binomial_var
          dat_row <- dat_row %>%
            mutate(bigC = exp((abs(qnorm(ci_width)) *
                           sqrt(log(1 + .data$Abundance_CV^2))))) %>%
            mutate(df = 0)
        }else{
          df_tvar <- df_tvar[1, 1]
          dat_row <- dat_row %>%
            # average multiplier
            mutate(rate          = mean(.data$rate),
                   rate_df       = sum(.data$rate_df),
                   rate_var      = sum(.data$rate_var),
                   rate_var_Nhat = .data$Abundance^2 * .data$rate_CV^2,
                   rate_CV       = sqrt(sum(.data$rate_var))/mean(.data$rate))
          dat_row <- dat_row %>%
            # CV weights for Satterthwaite df
            mutate(wtcv = sum(c((sqrt(.data$ER_var_Nhat[1])/
                                      .data$Abundance[1])^4/
                                 .data$ER_df[1],
                                (df_tvar/.data$Abundance[1]^2)^2/
                                  (length(ddf$fitted) - length(ddf$par)),
                                if_else(.data$df==0, 0 ,
                                        (.data$rate_var_Nhat[1]/
                                         .data$Abundance[1]^2)^2/
                                          .data$rate_df[1]),
                                (.data$group_var_Nhat[1]/
                                 .data$Abundance[1]^2)^2/
                                 (length(ddf$fitted)-1)
                               ),
                              na.rm=TRUE)) %>%
            # calculate Satterthwaite df
            mutate(df = sum(c((sqrt(.data$ER_var_Nhat[1])/.data$Abundance[1])^2,
                              (df_tvar/.data$Abundance[1]^2),
                              if_else(.data$df==0, 0 ,
                                      (.data$rate_var_Nhat[1]/
                                       .data$Abundance[1]^2)),
                              (.data$group_var_Nhat[1]/.data$Abundance[1]^2)
                             ),
                            na.rm=TRUE)^2) %>%
            mutate(df = .data$df/.data$wtcv) %>%
            mutate(bigC = exp((abs(qt(ci_width, .data$df)) *
                           sqrt(log(1 + .data$Abundance_CV^2)))))
          # drop weight column
          dat_row$wtcv <- NULL
        }
      }
      # actually calculate the CIs
      dat_row <- dat_row %>%
        mutate(LCI = .data$Abundance / .data$bigC,
               UCI = .data$Abundance * .data$bigC)

      this_stra_row <- cbind.data.frame(this_stra_row, dat_row, row.names=NULL)

      res <- rbind(res, this_stra_row)
    }
  }


  ## last bit of formatting
  res <- as.data.frame(res)
  res <- unique(res)

  # warn if we only had one transect in one or more strata
  if(any(res$k == 1)){
    warning("One or more strata have only one transect, cannot calculate empirical encounter rate variance")
  }

  # fix area == covered area for compatibility with mrds::dht
  if(est_density){
    #res$Area <- res$Covered_area
    res <- res %>%
      mutate(Area = .data$Covered_area)
  }else{
    # generate density results too!
    dens_res <- res %>%
      mutate(Density = .data$Abundance/.data$Area,
             df_var = .data$df_var/.data$Area^2) %>%
      mutate(Density_se = sqrt(.data$Abundance_se^2/.data$Area^2)) %>%
      mutate(Density_CV = .data$Density_se/.data$Density) %>%
      mutate(bigC = exp((abs(qt(ci_width, .data$df)) *
                     sqrt(log(1 + .data$Density_CV^2))))) %>%
      mutate(LCI   = .data$Density / .data$bigC,
             UCI   = .data$Density * .data$bigC,
             Area  = .data$Covered_area) %>%
      select(!!stratum_labels, "Area", "n", "ER_var", "Effort", "k",
             "Density", "Density_CV", "Density_se", "UCI", "LCI", "df",
             "Covered_area", "df_var", "group_var", "group_mean",
             "rate_var", "rate", "rate_df", "rate_CV", "p_var", "p_average")

    # store this
    attr(res, "density") <- dens_res
  }

  # save stratum labels
  attr(res, "stratum_labels") <- stratum_labels
  # were we really estimating density?
  attr(res, "density_only") <- est_density
  # save the sample level estimates
  attr(res, "sample_res") <- res_sample
  # detection function variance data
  attr(res, "df_var") <- df_Nhat_unc
  # save the variance proportions
  attr(res, "prop_var") <- variance_contributions(res)
  # save grouped analysis (might be NULL)
  attr(res, "grouped") <- grouped
  # save enounter rate variance information
  attr(res, "ER_var") <- c(er_est, innes, binomial_var)
  # save stratification type
  attr(res, "stratification") <- stratification
  # save multiplier info
  attr(res, "multipliers") <- names(multipliers)
  # save sample_fraction
  attr(res, "sample_fraction") <- sample_fraction

  class(res) <- c("dht_result", "data.frame")
  return(res)
}

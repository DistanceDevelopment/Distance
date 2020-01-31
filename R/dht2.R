#' Abundance estimation for distance sampling models
#'
#' Once a detection function is fitted to data, this function can be used to compute abundance estimates over required areas. The function also allows for stratification and variance estimation via various schemes (see below).
#'
#' @param ddf model fitted by \code{\link[Distance]{ds}} or \code{\link[mrds]{ddf}}
#' @param strat_formula a formula giving the stratification structure (see "Stratification" below)
#' @param observations \code{data.frame} to link detection function data (indexed by \code{object} column IDs) to the transects (indexed by \code{Sample.Label} column IDs). See "Data" below.
#' @param transects \code{data.frame} with information about samples (points or line transects). See "Data" below.
#' @param geo_strat \code{data.frame} with information about any geographical stratification. See "Data" below.
#' @param flatfile data in the flatfile format, see \code{\link[Distance]{flatfile}}
#' @param convert_units conversion factor between units for the distances, effort and area. See "Units" below.
#' @param er_est encounter rate variance estimator to be used. See "Variance" below and \code{\link{varn}}.
#' @param multipliers \code{list} of \code{data.frame}s. See "Multipliers" below.
#' @param sample_fraction what proportion of the transects was covered (e.g., 0.5 for one-sided line transects)
#' @param stratification what do strata represent, see "Stratification" below.
#' @param ci_width for use with confidence interval calculation (defined as 1-alpha, so the default 95 will give a 95\% confidence interval).
#' @param innes logical flag for computing encounter rate variance using either the method of Innes et al (2002) where estimated abundance per transect divided by effort is used as the encounter rate, vs. (when \code{innes=FALSE}) using the number of observations divided by the effort (as in Buckland et al., 2001)
#' @param total_area for options \code{stratification="effort_sum"} and \code{stratification="replicate"} the area to use as the total for combined, weighted final estimates.
#' @return a \code{data.frame} with estimates and attributes containing additional information
#'
#' @export
#' @importFrom rlang .data
#' @importFrom stats qt na.omit predict terms var
#' @importFrom dplyr group_by group_by_at mutate ungroup select distinct mutate_if if_else summarize_all "%>%" filter_at inner_join anti_join bind_rows left_join arrange vars
#' @importFrom mrds DeltaMethod
#' @section Data:
#' The data format allows for complex stratification schemes to be set-up. Before going into this detail, three objects are always required:
#' \describe{
#' \item{ddf}{the detection function (see \code{Distance::ds} or \code{mrds::ddf} for information on the format of their inputs).}
#' \item{\code{observations}}{has one row per observation and links the observations to the transects. Required columns: \code{object} (unique ID for the observation, which must match with the data in the detection function) and \code{Sample.Label} (unique ID for the transect). Additional columns for strata which are not included in the detection function are required (stratification covariates that are included in the detection function do not need to be included here). The important case here is group size, which must have column name \code{size} (but does not need to be in the detection function).}
#' \item{\code{transects}}{has one row per sample (point or line transect). At least one row is required. Required columns: \code{Sample.Label} (unique ID for the transect), \code{Effort} (line length for line transects, number of visits for point transects), if there is more than one geographical stratum.}
#' }
#' With only these three arguments, abundance can only be calculated for the covered area. Including additional information on the area we wish to extrapolate to (i.e., the study area), we can obtain abundance estimates:
#' \describe{
#' \item{\code{geo_strat}}{has one row for each stratum that we wish to estimate abundance for. For abundance in the study area, at least one row is required. Required columns: \code{Area} (the area of that stratum). If there is >1 row, then additional columns, named in \code{strat_formula}.}
#' }
#' @section Multipliers:
#' It is often the case that we cannot measure distances to individuals or groups directly, but instead need to estimate distances to something they produce (e.g., for whales, their blows; for elephants their dung) -- this is referred to as indirect sampling. We may need to use estimates of production rate and decay rate for these estimates (in the case of dung or nests) or just production rates (in the case of songbird calls or whale blows). We refer to these conversions between "number of cues" and "number of animals" as "multipliers".
#' The \code{multipliers} argument is a \code{list}, with 2 possible elements (\code{creation} and \code{decay} Each element of which is a \code{data.frame} and must have at least a column named \code{rate}, which abundance estimates will be divided by (the term "multiplier" is a misnomer, but kept for compatibility with Distance for Windows). Additional columns can be added to give the standard error and degrees of freedom for the rate if known as \code{SE} and \code{df}, respectively.
#' @section Stratification:
#' There are four stratification options:
#' \describe{
#'  \item{geographical}{if each stratum represents a different geographical areas and you want the total over all the areas}
#'  \item{effort_sum}{if your strata are in fact from replicate surveys (perhaps using different designs) but you don't have many replicates and/or want an estimate of "average variance".}
#'  \item{replicate}{if you have replicate surveys but have many of them, this calculates the average abundance and the variance between those many surveys (think of a population of surveys)}
#'  \item{object}{if the stratification is really about the type of object observed, for example sex, species or life stage and what you want is the total number of individuals across all the classes of objects. For example, if you have stratified by sex and have males and females, but also want a total number of animals, you should use this option.}
#' }
#' @section Variance:
#' Variance in the estimated abundance comes from multiple sources. Depending on the data used to fit the model and estimate abundance, different components will be included in the estimated variances. In the simplest case, the detection function and encounter rate variance need to be combined. If group size varies, then this too must be included. Finally, if multipliers are used and have corresponding standard errors given, this are also included. Variances are combined by assuming independence between the measures and adding variances. A brief summary of how each component is calculated is given here, though see references for more details.
#' \describe{
#' \item{detection function}{variance from the detection function parameters is transformed to variance about the abundance via a sandwich estimator (see e.g., Appendix C of Borchers et al (2002)).}
#' \item{encounter rate}{for strata with >1 transect in them, the encounter rate estimators given in Fewster et al (2009) can be specified via the \code{er_est} argument. If the argument \code{innes=TRUE} then calculations use the estimated number of individuals in the transect (rather than the observed), which was give by Innes et al (2002) as a superior estimator. When there is only one transect in a stratum, Poisson variance is assumed. Information on the Fewster encounter rate variance estimators are given in \code{\link{varn}}}
#' \item{group size}{if objects occur in groups (sometimes "clusters"), then the empirical variance of the group sizes is added to the total variance.}
#' \item{multipliers}{if multipliers with standard errors are given, their corresponding variances are added. If no standard errors are supplied, then their contribution to variance is assumed to be 0.}
#' }
#' @section Units:
#' It is often the case that distances are recorded in one convenient set of units, whereas the study area and effort are recorded in some other units. To ensure that the results from this function are in the expected units, we use the \code{convert_units} argument to supply a single number to convert the units of the covered area to those of the study/stratification area (results are always returned in the units of the study area). For line transects, the covered area is calculated as \code{2 * width * length} where \code{width} is the effective (half)width of the transect (often referred to as w in the literature) and \code{length} is the line length (referred to as L). If \code{width} and \code{length} are measured in kilometres and the study area in square kilometres, then all is fine and \code{convert_units} is 1 (and can be ignored). If, for example, line length and distances were measured in metres, we instead need to convert this to be kilometres, by dividing by 1000 for each of distance and length, hence \code{convert_units=1e-6}. For point transects, this is slightly easier as we only have the radius and study area to consider, so the conversion is just such that the units of the truncation radius are the square root of the study area units.
#'
#' @references
#' Borchers, D.L., S.T. Buckland, and W. Zucchini. 2002 \emph{Estimating Animal Abundance: Closed Populations}. Statistics for Biology and Health. Springer London.
#'
#' Buckland, S.T., E.A. Rexstad, T.A. Marques, and C.S. Oedekoven. 2015 \emph{Distance Sampling: Methods and Applications}. Methods in Statistical Ecology. Springer International Publishing.
#'
#' Buckland, S.T., D.R. Anderson, K. Burnham, J.L. Laake, D.L. Borchers, and L. Thomas. 2001 \emph{Introduction to Distance Sampling: Estimating Abundance of Biological Populations}. Oxford University Press.
#'
#' Innes, S., M. P. Heide-Jorgensen, J.L. Laake, K.L. Laidre, H.J. Cleator, P. Richard, and R.E.A. Stewart. 2002 Surveys of belugas and narwhals in the Canadian high arctic in 1996. \emph{NAMMCO Scientific Publications} \bold{4}, 169–-190.
#' @name dht2
dht2 <- function(ddf, observations=NULL, transects=NULL, geo_strat=NULL,
                 flatfile=NULL,
                 strat_formula, convert_units=1, er_est=c("R2", "P2"),
                 multipliers=NULL, sample_fraction=1,
                 ci_width=0.95, innes=FALSE,
                 stratification="geographical",
                 total_area=NULL){

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

  # check ci width
  if(ci_width > 1 | ci_width < 0){
    stop("ci_width must be between 0 and 1!")
  }else{
    ci_width <- (1-ci_width)/2
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
    grouped <- eval(mc, parent.frame())
    rm(dddf, mc)
  }else{
    grouped <- NULL
  }
  # ^^ we'll save this output later on

  # holder of transect type
  transect_type <- if(ddf$ds$aux$point) "point" else "line"

  # what are the stratum labels specicied in strat_formula?
  stratum_labels <- attr(terms(strat_formula), "term.labels")

  if(!is.null(observations) & !is.null(transects)){
    if(!is.null(geo_strat)){
      # which occur at the geo level?
      geo_stratum_labels <- stratum_labels[stratum_labels %in%
                                           colnames(geo_strat)]
    }else{
      geo_stratum_labels <- NULL
    }

# list of protected column names that can't be in the data
# protected <- c("p", ".Label", )


    # TODO: do some data checking at this point
    # - check duplicate column names (e.g., obs$sex and df$data$sex)
    # - check column names don't conflict with columns created below
    # check the data
    dht2_checkdata(ddf, observations, transects, geo_strat, strat_formula,
                   stratum_labels, geo_stratum_labels)

    # prepare data
    obj_keep <- ddf$data$object[ddf$data$distance <= ddf$ds$aux$width &
                                ddf$data$distance >= ddf$ds$aux$left]
    bigdat <- ddf$data[ddf$data$object %in% obj_keep, ]
    observations <- observations[observations$object %in% obj_keep, ]


    # get probabilities of detection
    bigdat$p <- predict(ddf)$fitted

    bigdat <- merge(bigdat, observations, all.x=TRUE, by="object")

    # remove Sample.Label dupes
    if(!is.null(bigdat[["Sample.Label.x"]])){
      bigdat[["Sample.Label"]] <- bigdat[["Sample.Label.x"]]
      bigdat[["Sample.Label.x"]] <- NULL
      bigdat[["Sample.Label.y"]] <- NULL
    }

    # merge onto transects
    bigdat <- merge(bigdat, transects, all.x=TRUE, all.y=TRUE,
                    by="Sample.Label")

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
    # TODO: check that unqiue(Area, stratum_labels) make sense

    flatfile <- safetruncate(flatfile, ddf$meta.data$width, ddf$meta.data$left)
    bigdat <- flatfile

    # get probabilities of detection
    pp <- predict(ddf, bigdat, compute=TRUE)$fitted
    bigdat$p <- pp

    if(strat_formula==~1){
      bigdat$Area <- sum(unique(bigdat[, c("Area", "Region.Label")])$Area)
      bigdat$Region.Label <- NULL
      bigdat$.Label <- "Total"

      stratum_labels <- ".Label"
      geo_stratum_labels <- ".Label"
    }else{
      geo_stratum_labels <- c()
    }

    # TODO: unchecked for multiple stratification
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

    # TODO: this needs to be fixed for the multi-strata case
    # check that Area-stratum combinations are coherent
    ucomb <- unique(bigdat[, c("Area", stratum_labels)])
    if(length(na.omit(ucomb[,stratum_labels])) >
       length(na.omit(unique(bigdat[,stratum_labels])))){
      stop(">1 Area defined for a single stratum label, check data")
    }

  }else{
    stop("Need to supply either observations, transects and geo_strat OR flatfile")
  }

  # merge-in multipliers
  if(!is.null(multipliers)){
    if(!is.list(multipliers)){
      stop("multipliers must be a list")
    }
    if(!all(names(multipliers) %in% c("creation", "decay")) |
       is.null(names(multipliers))){
      stop("Multipliers must be named \"creation\" and \"decay\"")
    }
    if(length(multipliers)>2){
      stop("Only one creation and one decay rate may be provided")
    }

    # base multipliers
    bigmult <- data.frame(rate = 1,
                          # need to set df to zero here as we will
                          # add later on...
                          rate_df = 0,
                          rate_SE = 0,
                          rate_CV = 0)

    for(ii in names(multipliers)){
      if(!is.data.frame(multipliers[[ii]])){
        stop(paste0("multipliers[[", ii, "]] must be a data.frame"))
      }
      # check multipliers has at least a rate column
      if(!("rate" %in% names(multipliers[[ii]]))){
        stop(paste("You need at least a column named \"rate\" in",
                   names(multipliers)[ii], "multiplier"))
      }

      if(is.null(multipliers[[ii]]$df)){
        # this deals with the no df case
        multipliers[[ii]]$df <- Inf
      }
      if(is.null(multipliers[[ii]]$SE)){
        multipliers[[ii]]$SE <- 0
      }

      bigmult$rate_CV <- bigmult$rate_CV^2 +
                          (multipliers[[ii]]$SE/
                           multipliers[[ii]]$rate)^2
      # since we are dividing, use the sandwich estimator,
      # var(1/x) = 1/x^2 * var(x) * 1/x^2 => se(1/x) = se(x)/x^2
      if(ii=="decay" && multipliers[[ii]]$SE!=0){
        multipliers[[ii]]$SE <- multipliers[[ii]]$SE/multipliers[[ii]]$rate^2
      }

      bigmult$rate <- bigmult$rate*multipliers[[ii]]$rate
      bigmult$rate_SE <- sqrt(bigmult$rate_SE^2 + multipliers[[ii]]$SE^2)
      bigmult$rate_df <- bigmult$rate_df + multipliers[[ii]]$df

    }
    bigmult$rate_CV <- sqrt(bigmult$rate_CV)
    bigdat <- merge(bigdat, bigmult, all.x=TRUE)
    mult <- TRUE
  }else{
    # setup "fake" data for when we don't have multipliers
    # this makes the calculations cleaner below
    mult <- FALSE
    bigdat <- bigdat %>%
      mutate(Nc_cuecorrected = NA,
             rate = 1,
             rate_df = 1,
             rate_SE = 0,
             rate_CV = 0)
  }

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

  df_width <- (ddf$ds$aux$width - ddf$ds$aux$left)*convert_units

  # TODO: check that sample_fraction is positive and a single number
  area_calc <- function(width, effort, transect_type, sample_fraction){
    if(transect_type=="point"){
      return(effort*pi*width^2*sample_fraction)
    }else{
      return(effort*2*width*sample_fraction)
    }
  }

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

  # first do transect level calculations
  res <- bigdat %>%
    group_by_at(vars("Sample.Label"))

  # if we are stratifying by object-level covariates, need to
  # make sure summaries are made at the sample-stratum level here
  if(stratification=="object"){
    res <- res %>%
      group_by_at(.vars=stratum_labels, .add=TRUE)
  }

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
             n_individuals = sum(.data$size, na.rm=TRUE),
             n_observations = length(na.omit(unique(.data$object))),
             # abundance estimate per stratum in covered area
             Nc = sum(.data$Nhat, na.rm=TRUE),
             # covered area per transect
             Covered_area = area_calc(df_width, .data$Effort,
                                      transect_type, sample_fraction),
             # get group size stats
             group_var  = if_else(.data$n_observations>1,
                                  var(.data$size, na.rm=TRUE)/
                                   sum(!is.na(.data$size)),
                                  0),
             group_mean = mean(.data$size, na.rm=TRUE)) %>%
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
  vardat$p_var <- summary(ddf)$average.p.se[1,1]^2
  vardat$p_average <- summary(ddf)$average.p

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
      select(!!stratum_labels, "Sample.Label", "Area", "n", "Nc", "transect_n", "Effort",
             "Covered_area", "df_var", "transect_Nc", "group_var", "group_mean",
             "Nc_cuecorrected", "rate_SE", "rate", "rate_df", "rate_CV", "p_var", "p_average",
             "transect_n_observations") %>%
      # keep only unique rows
      distinct()

# TODO: fix
if(mult){
  res <- res %>%
    mutate(Nc = .data$Nc_cuecorrected) %>%
    mutate(transect_Nc = .data$transect_Nc/.data$rate)
}
  # save ungrouped version for summary calculation later
  resT <- ungroup(res)

  # calculate ER variance
  res <- ER_var_f(res, innes=innes, er_est=er_est, est_density)

  # calculate final summaries
  res <- res %>%
    # number of transects, total effort and covered area per stratum
    mutate(k = length(.data$Sample.Label),
           Effort = sum(.data$Effort),
           Covered_area = sum(.data$Covered_area)) %>%
    ## keep only these columns
    select(!!stratum_labels, "Area", "Nc", "n", "ER_var", "Effort", "k",
           "Covered_area", "df_var", "group_var", "group_mean", "ER_var_Nhat",
           "rate_SE", "rate", "rate_df", "rate_CV", "p_var", "p_average") %>%
    ## now just get the distinct cases
    distinct()
    # calculate stratum encounter rate

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

#if(!is.null(grouped)){
#  res <- res %>%
#    mutate(Abundance_CV = sqrt(sum(c(ER_CV^2,
#                                     df_CV^2))))
#  if(nrow(grouped)>1){
#    gr <- grouped[-nrow(grouped),,drop=FALSE]
#  }
#  gr <- gr[, c(stratum_labels, "Abundance_CV", "Abundance", "Effort")]
#  names(gr) <- c(stratum_labels, "group_Abundance_CV",
#                 "group_Abundance", "group_Effort")
#  res <- merge(res, gr, by=stratum_labels)
#
#  # get var info
#  gr_var <- attr(grouped, "df_var")
#  vcov <- solvecov(ddf$hessian)$inv
#
#  res <- res %>%
#    mutate(scale     = .data$Area/.data$Covered_area) %>%
#    mutate(group_cov = covn(.data$group_Effort/(scale*sum(.data$Effort)),
#                            .data$group_Abundance/scale,
#                            .data$Abundance/scale,
#                            er_est)) %>%
#    # cov.Nc.Ncs <- as.vector(by(obs$size*(1 - obs$pdot)/obs$pdot^2,
#    #                            obs$Region.Label, sum))
#    mutate(group_cov  = (.data$group_cov +
#                         diag(t(gr_var$partial)%*%
#                          vcov%*%
#                          df_Nhat_unc$partial))^2) %>%
#    mutate(group_mean = .data$group_Abundance_CV^2 +
#                         .data$Abundance_CV^2 -
#                         2*.data$group_cov/
#                        (.data$Abundance * .data$group_Abundance)) %>%
#    mutate(group_CV = if_else(.data$group_var==0, 0,
#                              sqrt(.data$group_var)/.data$group_mean))
#  res$scale <- NULL
#  res$group_Abundance_CV <- NULL
#  res$group_Abundance <- NULL
#}

  # se and CV
  res <- res %>%
    mutate(Abundance_CV = sqrt(sum(c(.data$ER_CV^2,
                                     .data$df_CV^2,
                                     .data$rate_CV^2,
                                     .data$group_CV^2),
                                   na.rm=TRUE)))%>%
    mutate(Abundance_se = .data$Abundance_CV*.data$Abundance) %>%
    distinct()

  # total degrees of freedom and CI calculation
  res <- res %>%
    mutate(df = .data$Abundance_CV^4/
                  sum(c(if_else(.data$k==1, 0, .data$ER_CV^4/.data$ER_df),
                        .data$df_CV^4/(length(ddf$fitted) - length(ddf$par)),
                        .data$group_CV^4/(.data$n-1),
                        .data$rate_CV^4/.data$rate_df),
                   na.rm=TRUE)) %>%
    # adjust if df is too small
    mutate(df = if_else(.data$Abundance_CV==0, 1, .data$df)) %>%
    # big C for Satterthwaite
    mutate(bigC = exp((abs(qt(ci_width, .data$df)) *
                   sqrt(log(1 + .data$Abundance_CV^2))))) %>%
    # actually calculate the cis
    mutate(LCI = if_else(.data$Abundance_CV==0, .data$Abundance, .data$Abundance / .data$bigC),
           UCI = if_else(.data$Abundance_CV==0, .data$Abundance, .data$Abundance * .data$bigC)) %>%
    # done!
    ungroup() %>%
    distinct()

## TODO: summary stuff, this is BAD code
  # make a summary
  res <- as.data.frame(res)
# TODO: this is untested for multiple labels
  for(this_stratum in stratum_labels){
    dat_row <- res

    # remove NA labels
    iind <- which(is.na(dat_row[, this_stratum]))
    if(length(iind)>0) dat_row <- dat_row[-iind, ]

    # save row names
    stra_row <- dat_row[, stratum_labels, drop=FALSE]
    # remove labels
    dat_row[, stratum_labels] <- NULL

    if(length(stra_row[[this_stratum]]) > 1){
      this_stra_row <- stra_row[1, , drop=FALSE]
      this_stra_row[] <- NA
      this_stra_row[[this_stratum]] <- "Total"

# from mrds:
    # df for total estimate assuming sum of indep region estimates; uses
    # variances instead of cv's because it is a sum of means for encounter
    # rate portion of variance (df.total)

      dat_row <- dat_row %>%
        mutate(n  = sum(.data$n))
      #### stratification options
      ## in Distance for Windows these are in density
      if(stratification=="geographical"){
        # "weighting by area"
        #  which is adding up the abundances
        dat_row <- dat_row %>%
          mutate(ER_var       = sum(.data$ER_var, na.rm=TRUE)) %>%
          mutate(ER_var_Nhat  = sum(.data$ER_var_Nhat, na.rm=TRUE)) %>%
          mutate(weight       = 1,
                 Area         = sum(.data$Area),
                 Covered_area = sum(.data$Covered_area),
                 Effort       = sum(.data$Effort),
                 k            = sum(.data$k)) %>%
          mutate(ER_df = .data$ER_var_Nhat^2/sum((res$ER_var_Nhat^2/.data$ER_df)))
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
        # replicates, 2 different ways
        dat_row <- dat_row %>%
          mutate(weight       = xt * .data$Effort/sum(.data$Effort)) %>%
          mutate(ER_var       = sum((.data$Effort/sum(.data$Effort))^2*.data$ER_var,
                                    na.rm=TRUE)) %>%
          mutate(ER_var_Nhat  = sum(.data$weight^2*.data$ER_var_Nhat,
                                    na.rm=TRUE)) %>%
          mutate(ER_df        = .data$ER_var_Nhat^2/sum(
                                  ((.data$weight^2 * res$ER_var_Nhat)^2/.data$ER_df)))%>%
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
          mutate(ER_df        = .data$ER_var_Nhat^2/sum((res$ER_var_Nhat^2/.data$ER_df)))
      }else{
        # TODO: move to top
        stop("Invalid weighting option")
      }

      # calculate the weighted abundance
      dat_row <- dat_row %>%
        mutate(ER         = .data$n/.data$Effort) %>%
        mutate(ER_CV      = if_else(.data$ER==0, 0, sqrt(.data$ER_var)/.data$ER)) %>%
        mutate(Abundance  = sum(.data$weight*.data$Abundance)) %>%
        mutate(group_mean = mean(.data$group_mean),
               group_var  = sum(.data$group_var)) %>%
        mutate(group_CV   = if_else(.data$group_var==0, 0,
                                    sqrt(.data$group_var)/.data$group_mean))


      # calculate total variance for detection function
      vcov <- df_Nhat_unc$variance
      # vvv could do this if we wanted to ignore covariance
      #vcov[] <- 0
      #diag(vcov) <- diag(df_Nhat_unc$variance)
      df_tvar <- matrix(dat_row$weight, nrow=1) %*%
                  vcov %*%
                  matrix(dat_row$weight, ncol=1)
      dat_row <- dat_row %>%
        mutate(df_CV  = sqrt(df_tvar[1,1])/dat_row$Abundance[1])

      if(stratification=="replicate"){
        # get "between" variance (empirical variance of strata)
        tvar <- sum((dat_row$Abundance[-nrow(dat_row)] -
                     dat_row$Abundance[nrow(dat_row)])^2)/(nrow(dat_row)-2)
      }else if(stratification=="effort_sum"){
        # add the pre-weighted CVs
        tvar <- dat_row$Abundance[1]^2 *
                 sum(c(dat_row$ER_CV[1]^2,
                       dat_row$df_CV[1]^2,
                       dat_row$rate_CV[1]^2,
                       dat_row$group_CV[1]^2),
                     na.rm=TRUE)
      }else{
        # add all sources of variance (weighted as above)
        # here we subtract the detection function variance and add-in the total
        # variance, which includes a covar term
        tvar <- sum(dat_row$weight^2*(
                      (dat_row$Abundance_se^2-
                       res$Abundance^2*res$rate_CV^2 -
                        dat_row$df_var)),
                    na.rm=TRUE) +
                    df_tvar +
                    dat_row$Abundance[1]^2*dat_row$rate_CV[1]^2
      }


      dat_row <- dat_row %>%
        mutate(Nc = sum(.data$weight*.data$Nc),
               df_var = df_tvar[1,1]) %>%
        mutate(Abundance_se = sqrt(tvar)) %>%
        mutate(Abundance_CV = .data$Abundance_se/.data$Abundance) %>%
        mutate(df=NA,
               bigC = NA,
               LCI = NA,
               UCI = NA)
      # drop weights as we don't need them any more
      dat_row$weight <- NULL
      # drop unwanted rows
      dat_row <- dat_row %>%
        distinct()

      dat_row <- dat_row %>%
        # compute degrees of freedom
        # CV weights for Satterthwaite df
        mutate(wtcv = sum(c((sqrt(.data$ER_var_Nhat)/.data$Abundance)^4/.data$ER_df,
                            (df_tvar/.data$Abundance^2)^2/(length(ddf$fitted) - length(ddf$par)),
                            (.data$rate_SE/.data$Abundance)^4/.data$rate_df,
                            (.data$group_var/.data$Abundance^2)^2/(.data$n-1)),
                          na.rm=TRUE)) %>%
        # calculate Satterthwaite df
        mutate(df = sum(c((sqrt(.data$ER_var_Nhat)/.data$Abundance)^2,
                          (sqrt(df_tvar)/.data$Abundance)^2,
                          (.data$rate_SE/.data$Abundance)^2,
                          (sqrt(.data$group_var)/.data$Abundance)^2),
                        na.rm=TRUE)^2) %>%
        mutate(df = .data$df/.data$wtcv)
      # drop weight column
      dat_row$wtcv <- NULL

      # actually calculate the CIs
      dat_row <- dat_row %>%
        mutate(bigC = exp((abs(qt(ci_width, .data$df)) *
                       sqrt(log(1 + .data$Abundance_CV^2))))) %>%
        mutate(LCI = .data$Abundance / .data$bigC,
               UCI = .data$Abundance * .data$bigC)

      this_stra_row <- as.data.frame(c(this_stra_row, dat_row))

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
    res$Area <- res$Covered_area
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
             "rate_SE", "rate", "rate_df", "rate_CV", "p_var", "p_average")

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

  class(res) <- c("dht_result", "data.frame")
  return(res)
}

# calculate encounter rate variance

# this includes group size for all options, given the Innes et al estimator
#  will include this component

# note that as in mrds::dht we assume independence between strata
#  so the vcov matric for ER is diagonal (excluding the "total" row/col)
# erdat is assumed to be grouped by stratum labels at this point
#' @importFrom dplyr group_vars across all_of
ER_var_f <- function(erdat, innes, binomial_var=FALSE){
  if(binomial_var){
    # "varflag=0"
    # do the binomial var if A=a
    erdat <- erdat %>%
      mutate(pdot = .data$transect_n/.data$transect_Nc) %>%
      mutate(ER_var = sum(((1-.data$pdot)/.data$pdot^2)*.data$transect_n,
                          na.rm=TRUE)) %>%
      mutate(ER_var = if_else(is.infinite(.data$ER_var), 0, .data$ER_var)) %>%
      mutate(ER_var_Nhat = .data$ER_var/sum(.data$Covered_area)^2) %>%
      # if any stratum only had one transect:
      mutate(ER_var_Nhat = ifelse(length(unique(.data$Sample.Label))>1,
                                  .data$ER_var_Nhat, 0))
      erdat$pdot <- NULL
  }else{

    # sort the data if we use O2/O3 estimators
    if(any(erdat$er_est %in% c("O1", "O2", "O3", "S1", "S2"))){
      warning("Using one of O1, O2, O3, S1 or S2 encounter rate variance estimators, assuming that sorting on Sample.Label is meaningful.", immediate. = TRUE, call. = FALSE)
      if(!is.numeric(erdat$Sample.Label)){
        warning("Additionally, Sample.Label is not numeric, this may cause additional issues", immediate. = TRUE, call. = FALSE)
      }
      # Add in a column indexing the original order
      erdat$.orig.order <- 1:nrow(erdat)
      erdat <- erdat %>%
        # Sort by sample label
        arrange(.data$Sample.Label)
    }

#<<<<<<< HEAD
#    erdat <- erdat %>%
#      mutate(ER_var = varn(.data$Effort, .data$transect_n_observations,
#                           type=er_est)) %>%
#      mutate(ER_var = ifelse(length(unique(.data$Sample.Label))>1,
#                             .data$ER_var, 0))
#
#    # should the estimator of Innes et al be used?
#    if(innes){
#      # this is the "varflag=2"
#      erdat <- erdat %>%
#        # note that the ER reported in summary is using n/L but we need to use
#        # the N/L estimate for the total variance calculation
#        mutate(ER_var_Nc = varn(.data$Effort, .data$transect_Nc, type=er_est)) %>%
#        mutate(ER_var_Nc = ifelse(length(unique(.data$Sample.Label))>1,
#                               .data$ER_var_Nc,
#                               0)) %>%
#=======
    # save the current variables that the data is grouped by
    groupings_saved <- group_vars(erdat)

    # Innes et al estimator
    if(innes){
      # this is the "varflag=2"
      erdat <- erdat %>%
        # note that the ER reported in summary is using n/L but we need to use
        # the N/L estimate for the total variance calculation
        mutate(ER_var = varn(.data$Effort, .data$transect_Nc, type=unique(.data$er_est))) %>%
        mutate(ER_var = ifelse(length(unique(.data$Sample.Label))>1,
                               .data$ER_var,
                               0))# %>%
      erdat <- erdat %>%
        group_by(across("er_est"), .add=TRUE) %>%
#>>>>>>> multiddf-dht2
        # put ER var on the Nhat scale
        mutate(ER_var_Nhat = varn(.data$Effort/(sum(.data$Effort)*
                                                unique(.data$Area)/
                                                sum(.data$Covered_area)),
                                  .data$transect_Nc,
                                  type=unique(.data$er_est))) %>%
        ungroup() %>%
        # this is completely deranged hieroglyphics, replaces
        # group_by_at from dplyr 1.0.0
        group_by(across(all_of(groupings_saved)))

      erdat <- erdat %>%
        # if any strata only had one transect:
        mutate(ER_var_Nhat = ifelse(length(unique(.data$Sample.Label))>1,
                                    .data$ER_var_Nhat,
                                    0))
    }else{
    # "classic" ER estimator, see e.g. Fewster et al
      # this is the "varflag=1"
      erdat <- erdat %>%
#<<<<<<< HEAD
#=======
        group_by(across("er_est"), .add=TRUE) %>%
        mutate(ER_var = varn(.data$Effort, .data$transect_n_observations,
                             type=unique(.data$er_est))) %>%
        ungroup() %>%
        # this is completely deranged hieroglyphics, replaces
        # group_by_at from dplyr 1.0.0
        group_by(across(all_of(groupings_saved)))

      erdat <- erdat %>%
        mutate(ER_var = ifelse(length(unique(.data$Sample.Label))>1,
                               .data$ER_var, 0)) %>%
#>>>>>>> multiddf-dht2
        # put ER var on the Nhat scale
        mutate(ER_var_Nhat = ((.data$Area/sum(.data$Covered_area))*
                              .data$Nc*sum(.data$Effort))^2 *
                              .data$ER_var/
                              sum(.data$transect_n_observations)^2) %>%
        mutate(ER_var_Nhat = ifelse(length(unique(.data$Sample.Label))>1,
                                    .data$ER_var_Nhat,
                                    0))
    }
  }

  # let the Nhat estimate be 0 if the ER_var was 0
  erdat <- erdat %>%
    mutate(ER_var_Nhat = ifelse(is.na(.data$ER_var_Nhat) |
                                is.nan(.data$ER_var_Nhat),
                                0, .data$ER_var_Nhat))

# TODO: fix this for multi
# do we need to sort????
#  # put the data back in the order it entered the function,
#  # if we re-ordered it
#  if(er_est %in% c("O2", "O3")){
#    erdat <- erdat %>%
#      arrange(.data$.originalorder)
#    erdat$.originalorder <- NULL
#  }

  return(erdat)
}

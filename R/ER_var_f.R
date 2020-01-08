# calculate encounter rate variance
# note that as in mrds::dht we assume independence between strata
#  so the vcov matric for ER is (block?) diagonal
# encounter rate variance estimation function
ER_var_f <- function(erdat, innes, er_est, est_density){
  if(est_density){
    # "varflag=0"
    # do the binomial var if A=a
    erdat <- erdat %>%
      mutate(pdot = .data$n/.data$Nc) %>%
      mutate(ER_var = sum(.data$size^2*(1-.data$pdot)/.data$pdot^2) +
                           .data$Nc^2 * .data$group_var/.data$group_mean^2)
      erdat$pdot <- NULL

# TODO: this might not be right
    erdat <- erdat %>%
      mutate(ER_var_Nhat = (.data$Nc*sum(.data$Effort))^2 *
                             .data$ER_var/sum(.data$transect_n)^2) %>%
      # if any stratum only had one transect:
      mutate(ER_var_Nhat = ifelse(length(unique(.data$Sample.Label))>1,
                                  .data$ER_var_Nhat, 0))
  }else{

    # sort the data if we use O2/O3 estimators
    if(er_est %in% c("O2", "O3")){
      warning(paste("Using the", er_est, "encounter rate variance estimator, assuming that sorting on Sample.Label is meaningful"))
      if(!is.numeric(erdat$Sample.Label)){
        warning("Additionally, Sample.Label is not numeric, this may cause additional issues")
      }
      erdat <- erdat %>%
        mutate(.originalorder = 1:nrow(erdat)) %>%
        arrange(.data$Sample.Label)
    }

    if(innes){
      # this is the "varflag=2"
      erdat <- erdat %>%
        mutate(ER_var = varn(.data$Effort, .data$transect_Nc, type=er_est)) %>%
        mutate(ER_var = ifelse(length(unique(.data$Sample.Label))>1,
                               .data$ER_var, 0)) %>%
        # put ER var on the Nhat scale
        mutate(ER_var_Nhat = (.data$Area/sum(.data$Covered_area))^2 *
                             (.data$Nc*sum(.data$Effort))^2 *
                               .data$ER_var/sum(.data$transect_n)^2) %>%
        # if any strata only had one transect:
        mutate(ER_var_Nhat = ifelse(length(unique(.data$Sample.Label))>1,
                                    .data$ER_var_Nhat, 0))
#                                    ER_var_Nhat,
#                                    ((Area/sum(Covered_area))*Nc^2)*
#                                     (1/((Area/sum(Covered_area))*Nc) +
#                                      group_var/group_mean^2)))
    }else{
      # this is the "varflag=1"
      erdat <- erdat %>%
        mutate(ER_var = varn(.data$Effort, .data$transect_n_observations, type=er_est)) %>%
        mutate(ER_var = ifelse(length(unique(.data$Sample.Label))>1,
                               .data$ER_var, 0)) %>%
        # put ER var on the Nhat scale
        mutate(ER_var_Nhat = ((.data$Area/sum(.data$Covered_area))*.data$Nc*sum(.data$Effort))^2 *
                                .data$ER_var/
                                  sum(.data$transect_n_observations)^2 +
                               ((.data$Area/sum(.data$Covered_area))*.data$Nc)^2 * .data$group_var/.data$group_mean^2) %>%
        mutate(ER_var_Nhat = ifelse(length(unique(.data$Sample.Label))>1,
                                    .data$ER_var_Nhat, 0))
#                                    ((Area/sum(Covered_area))*Nc)^2 *
#                                     (1/transect_n_observations +
#                                      group_var/group_mean^2)))
#
    }
  }

  # let the Nhat estimate be 0 if the ER_var was 0
  erdat <- erdat %>%
    mutate(ER_var_Nhat = ifelse(is.na(.data$ER_var_Nhat) |
                                is.nan(.data$ER_var_Nhat),
                                0, .data$ER_var_Nhat))

  if(er_est %in% c("O2", "O3")){
    erdat <- erdat %>%
      arrange(.data$.originalorder)
    erdat$.originalorder <- NULL
  }

  return(erdat)
}

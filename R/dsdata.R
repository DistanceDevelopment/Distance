#' Simulated minke whale data with cluster size
#'
#' Data simulated from models fitted to 1992/1993 Southern Hemisphere minke whale data collected by the International Whaling Commission. See Branch and Butterworth (2001) for survey details (survey design is shown in figure 1(e)). Data simulated by David Borchers.
#'
#' @references Branch, T.A. and D.S. Butterworth. (2001) Southern Hemisphere minke whales: standardised abundance estimates from the 1978/79 to 1997/98 IDCR-SOWER surveys. Journal of Cetacean Research and Management 3(2): 143-174
#'
#' Hedley, S.L., and S.T. Buckland. (2004) Spatial models for line transect sampling. Journal of Agricultural, Biological, and Environmental Statistics 9: 181-199. doi:10.1198/1085711043578.
#'
#' @name ClusterExercise
#' @aliases ClusterExercise_units
#' @keywords datasets
#' @docType data
#' @format \code{data.frame} with 99 observations of 9 variables:
#' \tabular{ll}{
#' \code{Region.Label} \tab stratum label (\code{"North"} or \code{"South"}) \cr
#' \code{Area} \tab stratum area  (square nautical mile)\cr
#' \code{Sample.Label} \tab transect identifier\cr
#' \code{Effort} \tab transect length  (nautical mile)\cr
#' \code{distance} \tab observed distance  (nautical mile)\cr
#' \code{Cluster.strat} \tab strata based on cluster size: 1, 2 and 3+\cr 
#' \code{size} \tab {cluster size}\cr
#' \code{Study.Area} \tab {name of study area}}
NULL

#' Cue counts of whale blows
#'
#' Cues are treated as an indirect count, requiring the use of multipliers.
#'
#' Because whale blows disappear instantaneously, there is no need to measure a decay rate.
#' However a cue production rate (blows per individual per unit time) is required,
#' as is a measure of variability of that rate.
#'
#' @name CueCountingExample
#' @aliases CueCountingExample_units
#' @docType data
#' @format A data frame with 109 rows and 15 variables.
#'   \tabular{ll}{ 
#'     \code{Region.Label} \tab stratum labels \cr
#'     \code{Area} \tab size (km^2) of each stratum \cr
#'     \code{Sample.Label} \tab transect labels \cr
#'     \code{Cue.rate} \tab rate of blows per animal per hour \cr
#'     \code{Cue.rate.SE} \tab variability in cue rate \cr
#'     \code{Cue.rate.df} \tab degrees of freedom (number of animals sampled for cues) \cr
#'     \code{object} \tab object ID \cr
#'     \code{distance} \tab perpendicular distance (km) \cr
#'     \code{Sample.Fraction} \tab proportion of full circle scanned (radians) \cr
#'     \code{Sample.Fraction.SE} \tab variability in sampling fraction (0) \cr
#'     \code{Search.time} \tab Duration of scanning effort (hr) \cr
#'     \code{bss} \tab Beaufort sea state \cr
#'     \code{sp} \tab Species detected (all observations W in these data) \cr
#'     \code{size} \tab Number of animals in group (all 1 in these data) \cr
#'     \code{Study.Area} \tab study area name \cr
#'   }
#' @note  There are two other nuances in this survey.  Even though the survey is taking place
#' on a moving ship, effort is measured as amount of time scanning for blows.  In some instances,
#' it is not possible for the observer to scan the sea all around them as view may be restricted
#' by the ship's superstructure.  Here a \code{sampling fraction} multiplier is employed 
#' to deal with restricted vision.  Units of measure of \code{cue.rate} and \code{Search.time}
#' must be equal.
#' @keywords datasets
NULL
#' Eastern Tropical Pacific spotted dolphin survey
#'
#' Observers aboard tuna vessels detecting dolphin schools along with a number
#' of possibly useful covariates for modelling the detection function.
#'
#'
#' @name ETP_Dolphin
#' @aliases ETP_Dolphin_units
#' @docType data
#' @format A data frame with 1090 rows and 13 variables.
#'   \tabular{ll}{  
#'     \code{Region.Label} \tab stratum labels (only one)\cr
#'     \code{Area} \tab size (nmi) of each stratum \cr
#'     \code{Sample.Label} \tab transect labels \cr
#'     \code{Effort} \tab transect length (nmi) \cr
#'     \code{object} \tab object ID \cr
#'     \code{distance} \tab perpendicular distance (nmi) \cr
#'     \code{LnCluster} \tab natural log of cluster size \cr
#'     \code{Month} \tab month of detection \cr
#'     \code{Beauf.class} \tab Beaufort sea state \cr
#'     \code{Cue.type} \tab initial cue triggering detection \cr
#'     \code{Search.method} \tab observer method making the detection \cr
#'     \code{size} \tab cluster size \cr
#'     \code{Study.Area} \tab study area name \cr
#'   }
#' @keywords datasets
#' @source Inter-American Tropical Tuna Commission
#' @note Several different search methods included in these data
#'   \tabular{ll}{
#'   \code{0} \tab binoculars from crows nest \cr
#'   \code{2} \tab binoculars from elsewhere on ship \cr
#'   \code{3} \tab helicopter searching ahead of ship \cr
#'   \code{5} \tab radar detects of seabirds above dolphin schools \cr
#'   }
#' @seealso Several cue types were also recorded by observers.
#'   \tabular{ll}{
#'   \code{1} \tab seabirds above the school \cr
#'   \code{2} \tab water splashes \cr
#'   \code{3} \tab unspecified \cr
#'   \code{4} \tab floating objects such as logs \cr
#'   }
NULL
#' Simulated line transect survey data
#'
#' Simulated line transect survey. Twelve transects, detection function is half-normal. 
#' True object density is 79.8 animals per km^2.
#'
#' @name LTExercise
#' @aliases LTExercise_units
#' @docType data
#' @format A data frame with 106 rows and 7 variables
#'   \tabular{ll}{
#'     \code{Region.Label} \tab strata names (single stratum) \cr
#'     \code{Area} \tab size of study area (1 in this case, making abundance and density equal) \cr
#'     \code{Sample.Label} \tab transect ID \cr
#'     \code{Effort} \tab length of transects (km) \cr
#'     \code{object} \tab object ID  \cr
#'     \code{distance} \tab perpendicular distance (m) \cr
#'     \code{Study.Area} \tab name of study area \cr
#'   }
#' @note  There is no unit object associated with this dataset 
#' @source Simulated data, from the distance sampling introductory course, Centre for
#'  Research into Ecological & Environmental Modelling, University of St Andrews.
#' @keywords datasets
NULL
#' Simulated point transect survey data
#'
#' Simulated point transect survey. Thirty point transects, detection function is half-normal. 
#' True object density is 79.6 animals per hectare.
#'
#' @name PTExercise
#' @aliases PTExercise_units
#' @docType data
#' @format A data frame with 144 rows and 7 variables
#'   \tabular{ll}{
#'     \code{Region.Label} \tab strata names (single stratum) \cr
#'     \code{Area} \tab size of study area (0 in this case) \cr
#'     \code{Sample.Label} \tab transect ID \cr
#'     \code{Effort} \tab number of visits to point \cr
#'     \code{object} \tab object ID  \cr
#'     \code{distance} \tab radial distance (m) \cr
#'     \code{Study.Area} \tab name of study area \cr
#'   }
#' @source Simulated data, from the distance sampling introductory course, Centre for
#'  Research into Ecological & Environmental Modelling, University of St Andrews.
#' @keywords datasets
NULL
#' Savanna sparrow point transects 
#'
#' Point transect data collected in Colorado 1980/81 to examine effect of
#' agricultural practices upon avian community.  
#' 
#' Design consisted of point transects placed in multiple pastures 
#' (3 in 1980, 4 in 1981).  While many species were observed, only data for 
#' Savannah sparrows (Passerculus sandwichensis) are included here.
#'
#' @references Knopf, F.L., J.A. Sedgwick, and R.W. Cannon. (1988) Guild 
#' structure of a riparian avifauna relative to seasonal cattle grazing. 
#' The Journal of Wildlife Management 52 (2): 280–290. https://doi.org/10.2307/3801235.
#'
#' @name Savannah_sparrow_1980
#' @aliases Savannah_sparrow_1981 Savannah_sparrow_1980_units Savannah_sparrow_1981_units
#' @keywords datasets
#' @docType data
#' @format \code{data.frame} with 468 observations of 7 variables:
#' \tabular{ll}{
#' \code{Region.Label} \tab stratum label (pasture ID) \cr
#' \code{Area} \tab stratum area (set to 1 so density is reported) \cr
#' \code{Sample.Label} \tab transect identifier \cr
#' \code{Effort} \tab number of visits \cr
#' \code{object} \tab object ID \cr
#' \code{distance} \tab radial distance (m) \cr
#' \code{Study.Area} \tab name of study area
#' }
#' 
#' @note Data structure for 1981 data set is identical, but there are 
#' 448 observations across 4 pastures.
NULL

#' Simulated minke whale data 
#'
#' Data simulated from models fitted to 1992/1993 Southern Hemisphere minke whale data collected by the International Whaling Commission. See Branch and Butterworth (2001) for survey details (survey design is shown in figure 1(e)). Data simulated by David Borchers.
#'
#' @references Branch, T.A. and D.S. Butterworth. (2001) Southern Hemisphere minke whales: standardised abundance estimates from the 1978/79 to 1997/98 IDCR-SOWER surveys. Journal of Cetacean Research and Management 3(2): 143-174
#'
#' Hedley, S.L., and S.T. Buckland. (2004) Spatial models for line transect sampling. Journal of Agricultural, Biological, and Environmental Statistics 9: 181-199. doi:10.1198/1085711043578.
#'
#' @name Stratify_example
#' @aliases Stratify_example_units
#' @keywords datasets
#' @docType data
#' @format \code{data.frame} with 99 observations of 7 variables:
#' \tabular{ll}{
#' \code{Region.Label} \tab stratum label (\code{"North"} or \code{"South"}) \cr
#' \code{Area} \tab stratum area  (square nautical mile)\cr
#' \code{Sample.Label} \tab transect identifier\cr
#' \code{Effort} \tab transect length  (nautical mile)\cr
#' \code{object} \tab object ID \cr
#' \code{distance} \tab observed distance (nautical mile)\cr
#' \code{Study.Area} \tab {name of study area}}
NULL

#' Hawaiian amakihi point transect data
#'
#' Also known as the Common 'Amakihi, a type of Hawaiian honeycreeper 
#'
#' @name amakihi
#' @aliases amakihi_units
#' @docType data
#' @format A data frame with 1487 rows and 12 variables
#'   \tabular{ll}{
#'     \code{Region.Label} \tab strata names (seven strata) \cr
#'     \code{Area} \tab size of study area (set to 0) \cr
#'     \code{Sample.Label} \tab transect ID \cr
#'     \code{Effort} \tab number of visits to point \cr
#'     \code{object} \tab object ID  \cr
#'     \code{distance} \tab radial distance (m) \cr
#'     \code{Month} \tab month survey conducted (not used) \cr
#'     \code{OBs} \tab observer ID (note capitalisation of variable name) \cr
#'     \code{Sp} \tab species code (COAM) for all detections \cr
#'     \code{MAS} \tab Time after sunrise (min) \cr
#'     \code{HAS} \tab Time after sunrise (hours) \cr
#'     \code{Study.Area} \tab name of study area \cr
#'   }
#' @note Example for investigating covariates in the detection function.  Note high 
#' colinearity between two measures of time since sunrise.  Convergence problems can
#' result from models with several factor covariates.
#' @references Marques, T.A., L. Thomas, S.G. Fancy and S.T. Buckland. (2007)
#'  Improving estimates of bird density using multiple-covariate distance sampling. 
#'  The Auk 124 (4): 1229–1243. 
#'  https://doi.org/10.1642/0004-8038(2007)124[1229:IEOBDU]2.0.CO;2. 
#' @keywords datasets
NULL
#' Capercaillie in Monaughty Forest
#'
#' Data from a line transect survey of capercaillie in Monaughty Forest, Moray, 
#' Scotland.
#'
#' @name capercaillie
#' @docType data
#' @keywords datasets
#' @aliases capercaillie_units
#' 
#' @format   A data frame with 112 observations on the following 9 variables.
#' \tabular{ll}{
#'   \code{Sample.Label} \tab name of single transect \cr
#'   \code{Effort} \tab transect length (km) \cr
#'   \code{distance} \tab perpendicular distance (m) \cr
#'   \code{object} \tab object ID \cr
#'   \code{size} \tab only individual birds detected \cr
#'   \code{detected} \tab whether detected \cr
#'   \code{observer} \tab single observer data \cr
#'   \code{Region.Label} \tab stratum name \cr
#'   \code{Area} \tab size of Monaughty Forest (ha) \cr
#' }
NULL
#' Ducknest line transect survey data
#'
#' Simulated line transect survey of duck nests, designed to reproduce the data of Figure 2 in Anderson and Pospahala (1970).
#'
#' The Monte Vista National Wildlife Refuge is in southern Colorado in the USA at an altitude
#' of roughly 2400m.
#'
#' @name ducknest
#' @aliases ducknest_units
#' @docType data
#' @format A data frame with 534 rows and 7 variables
#'   \tabular{ll}{
#'     \code{Region.Label} \tab strata names (single stratum in this instance) \cr
#'     \code{Area} \tab size of refuge (0 in this case, actual size 60km^2) \cr
#'     \code{Sample.Label} \tab transect ID \cr
#'     \code{Effort} \tab length of transects (km) \cr
#'     \code{object} \tab nest ID  \cr
#'     \code{distance} \tab perpendicular distance (m) \cr
#'     \code{Study.Area} \tab name of wildlife refuge \cr
#'   }
#' @references Anderson, D. R., and R. S. Pospahala. 1970. Correction of bias in belt transect studies of immotile objects. The Journal of Wildlife Management 34 (1): 141–146. https://doi.org/10.2307/3799501
#' @source Simulated data, from the distance sampling introductory course, Centre for Research into Ecological & Environmental Modelling, University of St Andrews.
#' @keywords datasets
#' @aliases ducknests_units
NULL
#' Golf tee data
#'
#' The data are from independent surveys by eight observers of a population of 250
#' groups (760 individuals) of golf tees.  The tees, of two colours, were placed 
#' in groups of between 1 and 8 in a survey region of 1680 m^2^, either exposed 
#' above the surrounding grass, or at least partially hidden by it.  They were 
#' surveyed by the 1999 statistics honours class at the Univ of St Andrews.
#' 
#' We treat each group of golf tees as a single animal with size equal to the 
#' number of tees in the group; yellow tees are male, green are female; tees 
#' exposed above the surrounding grass are classified as exposed, others as 
#' unexposed.  We are grateful to Miguel Bernal for making these data available; 
#' they were collected by him as part of a masters project.
#'
#' @name golftees
#' @aliases golftees_units
#' @docType data
#' @format The format is:
#' \describe{ List of 4 
#' $ book.tee.dataframe:'data.frame': 
#'   \item{$ object}{object ID}
#'   \item{$ observer}{observer ID}
#'   \item{$ detected}{detected or not detected}
#'   \item{$ distance}{perpendicular distance}
#'   \item{$ size}{group size}
#'   \item{$ sex}{number of tees in group}
#'   \item{$ exposure}{tee height above ground}  
#'   $ book.tee.region   :'data.frame': 2 obs. of 2 variables: ..
#'   \item{$ Region.Label}{stratum name}
#'   \item{$ Area}{stratum size}
#'   $ book.tee.samples   :'data.frame': 11 obs. of 3 variables: ..
#'   \item{$ Sample.Label}{transect label}
#'   \item{$ Region.Label}{stratum name}
#'   \item{$ Effort}{transect length}
#'   $ book.tee.obs :'data.frame': 162 obs. of 3 variables:
#'   \item{$ object}{object ID}
#'   \item{$ Region.Label}{stratum in which it was detected}
#'   \item{$ Sample.Label}{transect on which it was detected}
#'   }
#' @keywords datasets
#' @references 
#' Borchers, D. L., S.T. Buckland, and W. Zucchini. 2002. Estimating Animal Abundance: Closed Populations. Statistics for Biology and Health. London: Springer-Verlag. https://www.springer.com/gp/book/9781852335601.
#' 
#' Buckland, S.T., D.R. Anderson, K.P. Burnham, J.L. Laake, D.L. Borchers, and L. Thomas. Advanced Distance Sampling: Estimating Abundance of Biological Populations. OUP Oxford, 2004.
NULL
#' Sika deer pellet data from southern Scotland
#'
#' Because sika deer spend most of their time in woodland areas, abundance
#' estimates are based on pellet group counts.  Line transect methods were applied
#' to estimate deer pellet group density by geographic block.
#'
#' Data presented here are from the Peebleshire portion of the study described by 
#' Marques et al. (2001).
#'
#' @name sikadeer
#' @aliases sikadeer_units
#' @docType data
#' @format A data frame with 1923 rows and 11 variables.
#'   \tabular{ll}{  
#'     \code{Region.Label} \tab stratum labels \cr
#'     \code{Area} \tab size (ha) of each stratum \cr
#'     \code{Sample.Label} \tab transect labels \cr
#'     \code{Defecation.rate} \tab rate of dung production per individual per day \cr
#'     \code{Defecation.rate.SE} \tab variability in defecation rate \cr
#'     \code{Decay.rate} \tab time (days) for dung to become undetectable \cr
#'     \code{Decay.rate.SE} \tab variability in decay rate \cr
#'     \code{Effort} \tab transect length (km) \cr
#'     \code{object} \tab object ID \cr
#'     \code{distance} \tab perpendicular distance (cm) \cr
#'     \code{Study.Area} \tab study area name \cr
#'   }
#' @references Marques, F.F.C., S.T. Buckland, D. Goffin, C.E. Dixon, D.L. Borchers, B.A. Mayle, and A.J. Peace. (2001).
#'  Estimating deer abundance from line transect surveys of dung: 
#'  sika deer in southern Scotland. Journal of Applied Ecology 38 (2): 349–363.
#'   https://doi.org/10.1046/j.1365-2664.2001.00584.x
#' @keywords datasets
NULL
#' Simulation of encounter rate variance 
#'
#' Simulated line transect data with large differences in transect length.
#' In \code{systematic_var_2} that transect length gradient is coupled with a
#' strong animal gradient; exaggerating encounter rate variance between transects.  
#' 
#' True population size is 1000 objects in the study area of size 0.5 km^2; such that 
#' true density is 2000 objects per km.
#' 
#' @references Fewster, R.M., S.T. Buckland, K.P. Burnham, 
#' D.L. Borchers, P.E. Jupp, J.L. Laake and L. Thomas. (2009) Estimating the 
#' encounter rate variance in distance sampling. Biometrics 65 (1): 225–236.
#'  https://doi.org/10.1111/j.1541-0420.2008.01018.x.
#'
#' @name Systematic_variance_1
#' @aliases Systematic_variance_2 Systematic_variance_1_units Systematic_variance_2_units
#' @keywords datasets
#' @docType data
#' @format \code{data.frame} with 253 observations of 7 variables:
#' \tabular{ll}{
#' \code{Region.Label} \tab stratum label (default) \cr
#' \code{Area} \tab stratum area (0.5 km^2) \cr
#' \code{Sample.Label} \tab transect identifier \cr
#' \code{Effort} \tab transect length (km) \cr
#' \code{object} \tab object ID \cr
#' \code{distance} \tab perpendicular distance (m) \cr
#' \code{Study.Area} \tab name of study area
#' }
#' 
#' @note Data structure for \code{systematic_var_2} is identical, but there are 
#' 256 observations and a strong animal gradient.
NULL

#' Simulated line transect survey data with covariates
#'
#' Simulated line transect survey. Only eight line transects, detection function is 
#' half-normal. 
#'
#' @name unimak
#' @aliases unimak_units
#' @docType data
#' @format A data frame with 60 rows and 9 variables
#'   \tabular{ll}{
#'     \code{Region.Label} \tab strata names (single stratum) \cr
#'     \code{Area} \tab size of study area (mi^2) \cr
#'     \code{Sample.Label} \tab transect ID \cr
#'     \code{Effort} \tab transect length (mi) \cr
#'     \code{object} \tab object ID  \cr
#'     \code{distance} \tab perpendicular distance (km) \cr
#'     \code{MSTDO} \tab time since medication taken by observer (min) \cr
#'     \code{Hour} \tab time of day of sighting (hour) \cr
#'     \code{Study.Area} \tab name of study area \cr
#'   }
#' @note \code{Hour} is covariate that has no effect on detection function, while 
#' \code{MSTDO} does affect the detection function.  Examine the ability of model 
#' selection to choose the correct model.
#' @source Simulated data, from the distance sampling introductory course, Centre for
#'  Research into Ecological & Environmental Modelling, University of St Andrews.
#' @keywords datasets
NULL
#' Steve Buckland's winter wren surveys
#'
#' Observations of winter wren (Troglodytes troglodytes L.) collected by Steve Buckland in woodland/parkland at Montrave Estate near Leven, Fife, Scotland.
#'
#' Four different surveys were carried out:
#' \describe{
#'   \item{\code{wren_5min}}{5-minute point count}
#'   \item{\code{wren_snapshot}}{snapshot method}
#'   \item{\code{wren_cuecount}}{cue count}
#'   \item{\code{wren_lt}}{line transect survey}
#' }
#'
#' @name wren
#' @docType data
#' @references Buckland, S. T. (2006) Point-transect surveys for songbirds: robust methodologies. The Auk 123 (2): 345–357.
#' @source Steve Buckland
#' @keywords datasets
#' @aliases wren_5min wren_snapshot wren_cuecount wren_lt wren_5min_units wren_snapshot_units wren_cuecount_units wren_lt_units
#' @note wren_5min is  data frame with 134 observations of 8 variables
#'   \tabular{ll}{
#'     \code{Region.Label} \tab stratum name (single stratum) \cr
#'     \code{Area} \tab size (ha) of Montrave study area \cr
#'     \code{Sample.Label} \tab point label \cr
#'     \code{Effort} \tab Number of visits to point \cr
#'     \code{object} \tab Object ID \cr
#'     \code{distance} \tab radial distance (m) \cr
#'     \code{direction} \tab direction of detection from point \cr
#'     \code{Study.Area} \tab Montrave Estate \cr
#'   }
#' @note wren_snapshot is  data frame with 119 observations of 7 variables
#'   \tabular{ll}{
#'     \code{Region.Label} \tab stratum name (single stratum) \cr
#'     \code{Area} \tab size (ha) of Montrave study area \cr
#'     \code{Sample.Label} \tab point label \cr
#'     \code{Effort} \tab Number of visits to point \cr
#'     \code{object} \tab Object ID \cr
#'     \code{distance} \tab radial distance (m) \cr
#'     \code{Study.Area} \tab Montrave Estate \cr
#'   }
#' @note wren_cuecount is  data frame with 774 observations of 9 variables
#'   \tabular{ll}{
#'     \code{Region.Label} \tab stratum name (single stratum) \cr
#'     \code{Area} \tab size (ha) of Montrave study area \cr
#'     \code{Sample.Label} \tab point label \cr
#'     \code{Cue.rate} \tab Production rate (per min) of cues \cr
#'     \code{Cue.rate.SE} \tab  SE of cue production rate \cr
#'     \code{object} \tab Object ID \cr
#'     \code{distance} \tab radial distance (m) \cr
#'     \code{Search.time} \tab Time (min) listening for cues \cr
#'     \code{Study.Area} \tab Montrave Estate \cr
#'   }
#' @note wren_lt is  data frame with 156 observations of 8 variables
#'   \tabular{ll}{
#'     \code{Region.Label} \tab stratum name (single stratum) \cr
#'     \code{Area} \tab size (ha) of Montrave study area \cr
#'     \code{Sample.Label} \tab transect label \cr
#'     \code{Effort} \tab  transect length (km) \cr
#'     \code{object} \tab Object ID \cr
#'     \code{distance} \tab perpendicular distance (m) \cr
#'     \code{Study.Area} \tab Montrave Estate \cr
#'   }
NULL

#' Duiker camera trap survey
#'
#' Study took place in Tai National Park Cote d'Ivoire in 2014.  Filmed Maxwell's duikers (Philantomba maxwellii) were assigned to distance intervals; recorded distances are the midpoints of the intervals. This data includes only observations recorded at times of peak activity.
#'
#' @name DuikerCameraTraps
#' @aliases DuikerCameraTraps_units
#' @docType data
#' @format A data frame with 6277 rows and 6 variables
#'   \tabular{ll}{
#'     \code{Region.Label} \tab strata names (single stratum) \cr
#'     \code{Area} \tab size of study area (40.37 km^2) \cr
#'     \code{multiplier} \tab spatial effort, as the proportion of a circle covered by the angle of view of the camera (42 degrees for these cameras) \cr
#'     \code{Sample.Label} \tab camera station identifier (21 functioning cameras in this data set) \cr
#'     \code{Effort} \tab temporal effort, i.e. the number of 2-second time-steps over which the camera operated\cr
#'     \code{distance} \tab radial distance (m) to interval midpoint\cr
#'   }
#' @source Howe, E.J., Buckland, S.T., Després-Einspenner, M.-L. and Kühl, H.S. (2017), Distance sampling with camera traps. Methods Ecol Evol, 8: 1558-1565. doi:10.1111/2041-210X.12790
#'
#' Howe, Eric J. et al. (2018), Data from: Distance sampling with camera traps, Dryad, Dataset, https://doi.org/10.5061/dryad.b4c70
#' @keywords datasets
NULL


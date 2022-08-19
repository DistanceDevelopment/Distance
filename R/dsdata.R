#' Simulated minke whale data with cluster size
#'
#' Data simulated from models fitted to 1992/1993 Southern Hemisphere minke
#' whale data collected by the International Whaling Commission. See Branch and
#' Butterworth (2001) for survey details (survey design is shown in figure
#' 1(e)). Data simulated by David Borchers.
#'
#' @references Branch, T.A. and D.S. Butterworth. (2001) Southern Hemisphere
#' minke whales: standardised abundance estimates from the 1978/79 to 1997/98
#' IDCR-SOWER surveys. Journal of Cetacean Research and Management 3(2):
#' 143-174
#'
#' Hedley, S.L., and S.T. Buckland. (2004) Spatial models for line transect
#' sampling. Journal of Agricultural, Biological, and Environmental Statistics
#' 9: 181-199. \doi{10.1198/1085711043578}.
#'
#' @name ClusterExercise
#' @aliases ClusterExercise_units
#' @keywords datasets
#' @docType data
#' @format `data.frame` with 99 observations of 9 variables:
#'   * `Region.Label` stratum label (`"North"` or `"South"`)
#'   * `Area` stratum area  (square nautical mile)
#'   * `Sample.Label` transect identifier
#'   * `Effort` transect length  (nautical mile)
#'   * `distance` observed distance  (nautical mile)
#'   * `Cluster.strat` strata based on cluster size: 1, 2 and 3+
#'   * `size` cluster size
#'   * `Study.Area` name of study area
NULL

#' Cue counts of whale blows
#'
#' Cues are treated as an indirect count, requiring the use of multipliers.
#'
#' Because whale blows disappear instantaneously, there is no need to measure a
#' decay rate. However a cue production rate (blows per individual per unit
#' time) is required, as is a measure of variability of that rate.
#'
#' @name CueCountingExample
#' @aliases CueCountingExample_units
#' @docType data
#' @format A `data.frame` with 109 rows and 15 variables.
#'   * `Region.Label stratum labels
#'   * `Area` size (km^2) of each stratum
#'   * `Sample.Label` transect labels
#'   * `Cue.rate` rate of blows per animal per hour
#'   * `Cue.rate.SE` variability in cue rate
#'   * `Cue.rate.df` degrees of freedom (number of animals sampled for cues)
#'   * `object` object ID
#'   * `distance` perpendicular distance (km)
#'   * `Sample.Fraction` proportion of full circle scanned (radians)
#'   * `Sample.Fraction.SE` variability in sampling fraction (0)
#'   * `Search.time` Duration of scanning effort (hr)
#'   * `bss` Beaufort sea state
#'   * `sp` Species detected (all observations W in these data)
#'   * `size` Number of animals in group (all 1 in these data)
#'   * `Study.Area` study area name
#' @note  There are two other nuances in this survey.  Even though the survey
#' is taking place on a moving ship, effort is measured as amount of time
#' scanning for blows.  In some instances, it is not possible for the observer
#' to scan the sea all around them as view may be restricted by the ship's
#' superstructure.  Here a `sampling fraction` multiplier is employed to deal
#' with restricted vision.  Units of measure of `cue.rate` and `Search.time`
#' must be equal.
#' @keywords datasets
NULL

#' Eastern Tropical Pacific spotted dolphin survey
#'
#' Observers aboard tuna vessels detecting dolphin schools along with a number
#' of possibly useful covariates for modelling the detection function.
#'
#' Several different search methods included in these data
#'   * `0` binoculars from crows nest
#'   * `2` binoculars from elsewhere on ship
#'   * `3` helicopter searching ahead of ship
#'   * `5` radar detects of seabirds above dolphin schools
#'
#' Several cue types were also recorded by observers.
#'   * `1` seabirds above the school
#'   * `2` water splashes
#'   * `3` unspecified
#'   * `4` floating objects such as logs
#'
#' @name ETP_Dolphin
#' @aliases ETP_Dolphin_units
#' @docType data
#' @format A `data.frame` with 1090 rows and 13 variables:
#'   * `Region.Label` stratum labels (only one)
#'   * `Area` size (nmi) of each stratum
#'   * `Sample.Label` transect labels
#'   * `Effort` transect length (nmi)
#'   * `object` object ID
#'   * `distance` perpendicular distance (nmi)
#'   * `LnCluster` natural log of cluster size
#'   * `Month` month of detection
#'   * `Beauf.class` Beaufort sea state
#'   * `Cue.type` initial cue triggering detection
#'   * `Search.method` observer method making the detection
#'   * `size` cluster size
#'   * `Study.Area` study area name
#' @keywords datasets
#' @source Inter-American Tropical Tuna Commission
NULL

#' Simulated line transect survey data
#'
#' Simulated line transect survey. Twelve transects, detection function is
#' half-normal.  True object density is 79.8 animals per km^2.
#'
#' @name LTExercise
#' @aliases LTExercise_units
#' @docType data
#' @format A `data.frame` with 106 rows and 7 variables
#'   * `Region.Label` strata names (single stratum)
#'   * `Area` size of study area (1 in this case, making abundance and density
#'   equal)
#'   * `Sample.Label` transect ID
#'   * `Effort` length of transects (km)
#'   * `object` object ID 
#'   * `distance` perpendicular distance (m)
#'   * `Study.Area` name of study area
#' @note  There is no unit object associated with this dataset 
#' @source Simulated data, from the distance sampling introductory course,
#' Centre for Research into Ecological & Environmental Modelling, University of
#' St Andrews.
#' @keywords datasets
NULL

#' Simulated point transect survey data
#'
#' Simulated point transect survey. Thirty point transects, detection function
#' is half-normal.  True object density is 79.6 animals per hectare.
#'
#' @name PTExercise
#' @aliases PTExercise_units
#' @docType data
#' @format A `data.frame` with 144 rows and 7 variables
#'    * `Region.Label` strata names (single stratum)
#'    * `Area` size of study area (0 in this case)
#'    * `Sample.Label` transect ID
#'    * `Effort` number of visits to point
#'    * `object` object ID
#'    * `distance` radial distance (m)
#'    * `Study.Area` name of study area
#' @source Simulated data, from the distance sampling introductory course,
#' Centre for Research into Ecological & Environmental Modelling, University of
#' St Andrews.
#' @keywords datasets
NULL

#' Savanna sparrow point transects
#'
#' Point transect data collected in Colorado 1980/81 to examine effect of
#' agricultural practices upon avian community.
#'
#' Design consisted of point transects placed in multiple pastures (3 in 1980
#' and 4 in 1981). While many species were observed, only data for Savannah
#' sparrows (*Passerculus sandwichensis*) are included here.
#'
#' Data given here are different from the Distance for Windows example project.
#' Here each individual sighting is treated as an independent observation. This
#' corresponds to the analysis in Buckland et al. (2001) Section 8.7.  In the
#' Distance for Windows project objects are clusters of individuals. This
#' should not affect the results too greatly as most clusters were of size 1,
#' and so the results obtained should not be too far out.
#'
#' @references
#'
#' Knopf, F.L., J.A. Sedgwick, and R.W. Cannon. (1988) Guild structure of a
#' riparian avifauna relative to seasonal cattle grazing.  The Journal of
#' Wildlife Management 52 (2): 280–290.  \doi{10.2307/3801235}
#'
#' @name Savannah_sparrow_1980
#' @aliases Savannah_sparrow_1981 Savannah_sparrow_1980_units
#' Savannah_sparrow_1981_units
#' @keywords datasets
#' @docType data
#' @format `data.frame` with 468 observations (1980) and 448 observations
#' (1981) of 7 variables:
#'   * `Region.Label` stratum label (pasture ID)
#'   * `Area` stratum area (set to 1 so density is reported)
#'   * `Sample.Label` transect identifier
#'   * `Effort` number of visits
#'   * `object` object ID
#'   * `distance` radial distance (m)
#'   * `Study.Area` name of study area
NULL

#' Simulated minke whale data
#'
#' Data simulated from models fitted to 1992/1993 Southern Hemisphere minke
#' whale data collected by the International Whaling Commission. See Branch and
#' Butterworth (2001) for survey details (survey design is shown in figure
#' 1(e)). Data simulated by David Borchers.
#'
#' @references Branch, T.A. and D.S. Butterworth. (2001) Southern Hemisphere
#' minke whales: standardised abundance estimates from the 1978/79 to 1997/98
#' IDCR-SOWER surveys. Journal of Cetacean Research and Management 3(2):
#' 143-174
#'
#' Hedley, S.L., and S.T. Buckland. (2004) Spatial models for line transect
#' sampling. Journal of Agricultural, Biological, and Environmental Statistics
#' 9: 181-199. \doi{10.1198/1085711043578}.
#'
#' @name Stratify_example
#' @aliases Stratify_example_units
#' @keywords datasets
#' @docType data
#' @format `data.frame` with 99 observations of 7 variables:
#' `Region.Label` stratum label (`"North"` or `"South"`)
#' `Area` stratum area  (square nautical mile)
#' `Sample.Label` transect identifier
#' `Effort` transect length  (nautical mile)
#' `object` object ID
#' `distance` observed distance (nautical mile)
#' `Study.Area` name of study area
NULL

#' Hawaiian amakihi point transect data
#'
#' Also known as the Common 'Amakihi, a type of Hawaiian honeycreeper
#'
#' @name amakihi
#' @aliases amakihi_units
#' @docType data
#' @format A `data.frame` with 1487 rows and 12 variables
#'    * `Region.Label` strata names (seven strata)
#'    * `Area` size of study area (set to 0)
#'    * `Sample.Label` transect ID
#'    * `Effort` number of visits to point
#'    * `object` object ID
#'    * `distance` radial distance (m)
#'    * `Month` month survey conducted (not used)
#'    * `OBs` observer ID (note capitalisation of variable name)
#'    * `Sp` species code (COAM) for all detections
#'    * `MAS` Time after sunrise (min)
#'    * `HAS` Time after sunrise (hours)
#'    * `Study.Area` name of study area
#' @note Example for investigating covariates in the detection function.  Note
#' high colinearity between two measures of time since sunrise.  Convergence
#' problems can result from models with several factor covariates.
#' @references Marques, T.A., L. Thomas, S.G. Fancy and S.T. Buckland. (2007)
#' Improving estimates of bird density using multiple-covariate distance
#' sampling.  The Auk 124 (4): 1229–1243.
#' \doi{10.1642/0004-8038(2007)124[1229:IEOBDU]2.0.CO;2}
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
#' @format A `data.frame` with 112 observations on the following 9 variables.
#'    * `Sample.Label` name of single transect
#'    * `Effort` transect length (km)
#'    * `distance` perpendicular distance (m)
#'    * `object` object ID
#'    * `size` only individual birds detected
#'    * `detected` whether detected
#'    * `observer` single observer data
#'    * `Region.Label` stratum name
#'    * `Area` size of Monaughty Forest (ha)
NULL

#' Ducknest line transect survey data
#'
#' Simulated line transect survey of duck nests, designed to reproduce the data
#' of Figure 2 in Anderson and Pospahala (1970).
#'
#' The Monte Vista National Wildlife Refuge is in southern Colorado in the USA
#' at an altitude of roughly 2400m.
#'
#' @name ducknest
#' @aliases ducknest_units
#' @docType data
#' @format A `data.frame` with 534 rows and 7 variables
#'    * `Region.Label` strata names (single stratum in this instance)
#'    * `Area` size of refuge (0 in this case, actual size 60km^2)
#'    * `Sample.Label` transect ID
#'    * `Effort` length of transects (km)
#'    * `object` nest ID
#'    * `distance` perpendicular distance (m)
#'    * `Study.Area` name of wildlife refuge
#' @references Anderson, D. R., and R. S. Pospahala. 1970. Correction of bias
#' in belt transect studies of immotile objects. The Journal of Wildlife
#' Management 34 (1): 141–146. \doi{10.2307/3799501}
#' @source Simulated data, from the distance sampling introductory course,
#' Centre for Research into Ecological & Environmental Modelling, University of
#' St Andrews.
#' @keywords datasets
#' @aliases ducknests_units
NULL

#' Golf tee data
#'
#' The data are from independent surveys by eight observers of a population of
#' 250 groups (760 individuals) of golf tees.  The tees, of two colours, were
#' placed in groups of between 1 and 8 in a survey region of 1680 m^2, either
#' exposed above the surrounding grass, or at least partially hidden by it.
#' They were surveyed by the 1999 statistics honours class at the University of
#' St Andrews.
#'
#' We treat each group of golf tees as a single animal with size equal to the
#' number of tees in the group; yellow tees are male, green are female; tees
#' exposed above the surrounding grass are classified as exposed, others as
#' unexposed.  We are grateful to Miguel Bernal for making these data
#' available; they were collected by him as part of a masters project.
#'
#' @name golftees
#' @aliases golftees_units
#' @docType data
#' @format Data is a `list` with 4 elements each of which is a `data.frame`:
#'   * `book.tee.dataframe`
#'     * `object` object ID
#'     * `observer` observer ID
#'     * `detected` detected or not detected
#'     * `distance` perpendicular distance
#'     * `size` group size
#'     * `sex` number of tees in group
#'     * `exposure` tee height above ground
#'   * `book.tee.region`
#'     * `Region.Label` stratum name
#'     * `Area` stratum size
#'   * `book.tee.samples`
#'     * `Sample.Label` transect label
#'     * `Region.Label` stratum name
#'     * `Effort` transect length
#'   * `book.tee.obs`
#'     * `object` object ID
#'     * `Region.Label` stratum in which it was detected
#'     * `Sample.Label` transect on which it was detected
#' @keywords datasets
#' @references
#' Borchers, D. L., S.T. Buckland, and W. Zucchini. 2002. Estimating Animal
#' Abundance: Closed Populations. Statistics for Biology and Health. London:
#' Springer-Verlag. <https://link.springer.com/book/10.1007/978-1-4471-3708-5>
#'
#' Buckland, S.T., D.R. Anderson, K.P. Burnham, J.L. Laake, D.L. Borchers, and
#' L. Thomas. Advanced Distance Sampling: Estimating Abundance of Biological
#' Populations. Oxford University Press. Oxford, 2004.
NULL

#' Sika deer pellet data from southern Scotland
#'
#' Because sika deer spend most of their time in woodland areas, abundance
#' estimates are based on pellet group counts.  Line transect methods were
#' applied to estimate deer pellet group density by geographic block.
#'
#' Data presented here are from the Peebleshire portion of the study described
#' by Marques et al. (2001).
#'
#' @name sikadeer
#' @aliases sikadeer_units
#' @docType data
#' @format A `data.frame` with 1923 rows and 11 variables.
#'    * `Region.Label` stratum labels
#'    * `Area` size (ha) of each stratum
#'    * `Sample.Label` transect labels
#'    * `Defecation.rate` rate of dung production per individual per day
#'    * `Defecation.rate.SE` variability in defecation rate
#'    * `Decay.rate` time (days) for dung to become undetectable
#'    * `Decay.rate.SE` variability in decay rate
#'    * `Effort` transect length (km)
#'    * `object` object ID
#'    * `distance` perpendicular distance (cm)
#'    * `Study.Area` study area name
#' @references
#'
#' Marques, F.F.C., S.T. Buckland, D. Goffin, C.E. Dixon, D.L.  Borchers, B.A.
#' Mayle, and A.J. Peace. (2001). Estimating deer abundance from line transect
#' surveys of dung: sika deer in southern Scotland. Journal of Applied Ecology
#' 38 (2): 349–363.  \doi{10.1046/j.1365-2664.2001.00584.x}
#' @keywords datasets
NULL

#' Simulation of encounter rate variance
#'
#' `systematic_var_1` consists of simulated line transect data with large
#' differences in transect length. In `systematic_var_2` that transect length
#' gradient is coupled with a strong animal gradient; exaggerating encounter
#' rate variance between transects.
#'
#' True population size is 1000 objects in the study area of size 0.5 km^2;
#' such that true density is 2000 objects per km.
#'
#' @references Fewster, R.M., S.T. Buckland, K.P. Burnham, D.L. Borchers, P.E.
#' Jupp, J.L. Laake and L. Thomas. (2009) Estimating the encounter rate
#' variance in distance sampling. Biometrics 65 (1): 225–236.
#' \doi{10.1111/j.1541-0420.2008.01018.x}
#'
#' @name Systematic_variance_1
#' @aliases Systematic_variance_2 Systematic_variance_1_units
#' Systematic_variance_2_units
#' @keywords datasets
#' @docType data
#' @format `data.frame` with 253 observations (`systematic_var_1`) or 256
#' observations (`systematic_var_2`) of 7 variables:
#' `Region.Label` stratum label (default)
#' `Area` stratum area (0.5 km^2)
#' `Sample.Label` transect identifier
#' `Effort` transect length (km)
#' `object` object ID
#' `distance` perpendicular distance (m)
#' `Study.Area` name of study area
NULL

#' Simulated line transect survey data with covariates
#'
#' Simulated line transect survey. Only eight line transects, detection
#' function is half-normal.
#'
#' @name unimak
#' @aliases unimak_units
#' @docType data
#' @format A `data.frame` with 60 rows and 9 variables
#'    * `Region.Label` strata names (single stratum)
#'    * `Area` size of study area (mi^2)
#'    * `Sample.Label` transect ID
#'    * `Effort` transect length (mi)
#'    * `object` object ID
#'    * `distance` perpendicular distance (km)
#'    * `MSTDO` time since medication taken by observer (min)
#'    * `Hour` time of day of sighting (hour)
#'    * `Study.Area` name of study area
#' @note `Hour` is covariate that has no effect on detection function,
#' while `MSTDO` does affect the detection function.  Examine the ability
#' of model selection to choose the correct model.
#' @source Simulated data, from the distance sampling introductory course,
#' Centre for Research into Ecological & Environmental Modelling, University of
#' St Andrews.
#' @keywords datasets
NULL

#' Steve Buckland's winter wren surveys
#'
#' Observations of winter wren (*Troglodytes troglodytes L.*) collected by Steve
#' Buckland in woodland/parkland at Montrave Estate near Leven, Fife, Scotland.
#'
#' Four different surveys were carried out:
#'   * `wren_5min` 5-minute point count
#'   * `wren_snapshot` snapshot method
#'   * `wren_cuecount` cue count
#'   * `wren_lt` line transect survey
#'
#' @name wren
#' @docType data
#' @references Buckland, S. T. (2006) Point-transect surveys for songbirds:
#' robust methodologies. The Auk 123 (2): 345–357.
#' @source Steve Buckland
#' @keywords datasets
#' @aliases wren_5min wren_snapshot wren_cuecount wren_lt wren_5min_units
#' wren_snapshot_units wren_cuecount_units wren_lt_units
#' @note `wren_5min`: 134 observations of 8 variables
#'    * `Region.Label`  stratum name (single stratum)
#'    * `Area`  size (ha) of Montrave study area
#'    * `Sample.Label`  point label
#'    * `Effort`  Number of visits to point
#'    * `object`  Object ID
#'    * `distance`  radial distance (m)
#'    * `direction`  direction of detection from point
#'    * `Study.Area`  Montrave Estate
#' @note `wren_snapshot`: 119 observations of 7 variables
#'    * `Region.Label`  stratum name (single stratum)
#'    * `Area`  size (ha) of Montrave study area
#'    * `Sample.Label`  point label
#'    * `Effort`  Number of visits to point
#'    * `object`  Object ID
#'    * `distance`  radial distance (m)
#'    * `Study.Area`  Montrave Estate
#' @note `wren_cuecount`: 774 observations of 9 variables
#'    * `Region.Label`  stratum name (single stratum)
#'    * `Area`  size (ha) of Montrave study area
#'    * `Sample.Label`  point label
#'    * `Cue.rate`  Production rate (per min) of cues
#'    * `Cue.rate.SE`   SE of cue production rate
#'    * `object`  Object ID
#'    * `distance`  radial distance (m)
#'    * `Search.time`  Time (min) listening for cues
#'    * `Study.Area`  Montrave Estate
#' @note `wren_lt`: 156 observations of 8 variables
#'    * `Region.Label`  stratum name (single stratum)
#'    * `Area`  size (ha) of Montrave study area
#'    * `Sample.Label`  transect label
#'    * `Effort`   transect length (km)
#'    * `object`  Object ID
#'    * `distance`  perpendicular distance (m)
#'    * `Study.Area`  Montrave Estate
NULL

#' Duiker camera trap survey
#'
#' Study took place in Tai National Park Cote d'Ivoire in 2014.  Filmed
#' Maxwell's duikers (Philantomba maxwellii) were assigned to distance
#' intervals; recorded distances are the midpoints of the intervals. This data
#' includes only observations recorded at times of peak activity.
#'
#' @name DuikerCameraTraps
#' @aliases DuikerCameraTraps_units
#' @docType data
#' @format A `data.frame` with 6277 rows and 6 variables
#'   * `Region.Label`  strata names (single stratum)
#'   * `Area`  size of study area (40.37 km^2)
#'   * `multiplier`  spatial effort, as the proportion of a circle covered by
#'   the angle of view of the camera (42 degrees for these cameras)
#'   * `Sample.Label`  camera station identifier (21 functioning cameras in
#'   this data set)
#'   * `Effort`  temporal effort, i.e. the number of 2-second time-steps over
#'   which the camera operated
#'   * `distance`  radial distance (m) to interval midpoint
#' @source Howe, E.J., Buckland, S.T., Després-Einspenner, M.-L. and Kühl, H.S.
#' (2017), Distance sampling with camera traps. Methods Ecol Evol, 8:
#' 1558-1565. \doi{10.1111/2041-210X.12790}
#'
#' Howe, Eric J. et al. (2018), Data from: Distance sampling with camera traps,
#' Dryad, Dataset, \doi{10.5061/dryad.b4c70}
#' @keywords datasets
NULL

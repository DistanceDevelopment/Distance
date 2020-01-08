#' Generate table of unit conversions
#'
#' Returns a table of conversions between the units used in Distance for Windows. This is extracted from the \code{DistIni.mdb} default database.
#' @export
#'
#' @author David L Miller
#' @importFrom utils read.delim
# Taken from readdst, wheere it is in turn taken from Distance for Windows
units_table <- function(){

  # here is a lookup table
  # this trick thanks to Noam Ross
  unit_tab <- read.delim(textConnection(
  '                 Unit     |     Conversion
# ------------------------ | -----------------------------
                     Degree  | 1.745329251994329547437e-02
                        Gon  | 1.570796326794896696777e-02
                       Grad  | 1.570796326794896696777e-02
                Microradian  | 9.999999999999999547481e-07
                   Mil 6400  | 9.817477042468100018047e-04
                     Minute  | 2.908882086657215803975e-04
          Minute Centesimal  | 1.570796326794900068646e-04
                     Radian  | 1.000000000000000000000e+00
                     Second  | 4.848136811095359842698e-06
          Second Centesimal  | 1.570796326794900051705e-06
                       Acre  | 4.046873000000000047294e+03
                    Hectare  | 1.000000000000000000000e+04
                    Section  | 2.589998000000000000000e+06
          Square centimeter  | 1.000000000000000047922e-04
          Square centimetre  | 1.000000000000000047922e-04
                Square foot  | 9.290304000000000617110e-02
                Square inch  | 6.451599999999999826214e-04
           Square kilometer  | 1.000000000000000000000e+06
           Square kilometre  | 1.000000000000000000000e+06
               Square meter  | 1.000000000000000000000e+00
               Square metre  | 1.000000000000000000000e+00
                Square mile  | 2.589988000000000000000e+06
          Square millimeter  | 9.999999999999999547481e-07
          Square millimetre  | 9.999999999999999547481e-07
       Square nautical mile  | 3.429904000000000000000e+06
                Square yard  | 8.361269999999999535945e-01
                 Centimeter  | 1.000000000000000020817e-02
                 Centimetre  | 1.000000000000000020817e-02
           Chain (Benoit A)  | 2.011678240000000172927e+01
             Chain (Benoit)  | 2.011678249437587240323e+01
             Chain (Clarke)  | 2.011661949000000149113e+01
              Chain (Sears)  | 2.011676512155263196746e+01
          Chain (US Survey)  | 2.011684023368049878400e+01
                     Fathom  | 1.828799999999999981171e+00
                       Foot  | 3.048000000000000153655e-01
                Foot (1865)  | 3.048008333333332986470e-01
            Foot (Benoit A)  | 3.047997333333329894600e-01
            Foot (Benoit B)  | 3.047997347632709908005e-01
              Foot (Clarke)  | 3.047972651150995804237e-01
          Foot (Gold Coast)  | 3.047997101815089759924e-01
         Foot (Indian 1937)  | 3.047984100000000196040e-01
         Foot (Indian 1962)  | 3.047996000000000038632e-01
         Foot (Indian 1975)  | 3.047995000000000009877e-01
              Foot (Indian)  | 3.047995179900422901831e-01
   Foot (Modified American)  | 3.048122529845059824893e-01
               Foot (Sears)  | 3.047994715386762032416e-01
           Foot (US Survey)  | 3.048006096012192411848e-01
                    Furlong  | 2.011680000000000063665e+02
                       Inch  | 2.539999999999999896749e-02
                  Kilometer  | 1.000000000000000000000e+03
                  Kilometre  | 1.000000000000000000000e+03
                       Link  | 2.011661949759657175285e-01
            Link (Benoit A)  | 2.011678239999999950882e-01
              Link (Benoit)  | 2.011678249437587051585e-01
               Link (Sears)  | 2.011676512155262941395e-01
           Link (US Survey)  | 2.011684023368049967218e-01
                      Meter  | 1.000000000000000000000e+00
                      Metre  | 1.000000000000000000000e+00
             Metre (German)  | 1.000001359699999925468e+00
                       Mile  | 1.609344000000000050932e+03
           Mile (US Survey)  | 1.609347218694440016407e+03
                 Millimeter  | 1.000000000000000020817e-03
                 Millimetre  | 1.000000000000000020817e-03
              Nautical Mile  | 1.852000000000000000000e+03
                        Rod  | 5.029200000000000336797e+00
                       Yard  | 9.143999999999999905853e-01
            Yard (Benoit A)  | 9.143991999999999675808e-01
            Yard (Benoit B)  | 9.143992042898120287120e-01
              Yard (Clarke)  | 9.143917950000000072208e-01
         Yard (Indian 1937)  | 9.143952300000000033009e-01
         Yard (Indian 1962)  | 9.143987999999999560785e-01
         Yard (Indian 1975)  | 9.143985000000000029630e-01
              Yard (Indian)  | 9.143985539701268150381e-01
               Yard (Sears)  | 9.143984146160286652361e-01'),

  sep='|', comment.char="#", strip.white=TRUE)
  names(unit_tab) <- c("Unit", "Conversion")
  unit_tab$Conversion <- as.numeric(unit_tab$Conversion)

  return(unit_tab)
}

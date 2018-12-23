\name{nhanes2007}
\alias{nhanes2007}
\docType{data}
\title{
    National Health and Nutrition Examination Survey data, 2007-2008
}
\description{
Demographic and dietary intake variables from a U.S. national household survey
}

\usage{data(nhanes2007)}
\format{
  A data frame with 4,329 person-level observations on the following 26 variables measuring 24-hour dietary recall. See \url{https://wwwn.cdc.gov/nchs/nhanes/2013-2014/DR2IFF_H.htm} for more details about the variables.
  \describe{
    \item{\code{SEQN}}{Identification variable}
    \item{\code{SDMVSTRA}}{Stratum}
    \item{\code{SDMVPSU}}{Primary sampling unit, numbered within each stratum (1,2)}
    \item{\code{WTDRD1}}{Dietary day 1 sample weight}
    \item{\code{GENDER}}{Gender (0 = female; 1 = male)}
    \item{\code{RIDAGEYR}}{Age in years at the time of the screening interview; reported for survey participants between the ages of 1 and 79 years of age. All responses of participants aged 80 years and older are coded as 80.}
    \item{\code{RIDRETH1}}{Race/Hispanic origin
        (1 = Mexican American;
        2 = Other Hispanic;
        3 = Non-Hispanic White;
        4 = Non-Hispanic Black;
        5 = Other Race including multiracial)}
    \item{\code{BMXWT}}{Body weight (kg)}
    \item{\code{BMXBMI}}{Body mass Index ((weight in kg) / (height in meters)**2)}
    \item{\code{DIET}}{On any diet (0 = No; 1 = Yes)}
    \item{\code{CALDIET}}{On a low-calorie diet (0 = No; 1 = Yes)}
    \item{\code{FATDIET}}{On a low-fat diet (0 = No; 1 = Yes)}
    \item{\code{CARBDIET}}{On a low-carbohydrate diet (0 = No; 1 = Yes)}
    \item{\code{DR1DRSTZ}}{Dietary recall status that indicates quality and completeness of survey participant's response to dietary recall section.
        (1 = Reliable and met the minimum criteria;
        2 = Not reliable or not met the minimum criteria;
        4 = Reported consuming breast-milk (infants and children only))}
    \item{\code{DR1TKCAL}}{Energy (kcal)}
    \item{\code{DR1TPROT}}{Protein (gm)}
    \item{\code{DR1TCARB}}{Carbohydrate (gm)}
    \item{\code{DR1TSUGR}}{Total sugars (gm)}
    \item{\code{DR1TFIBE}}{Dietary fiber (gm)}
    \item{\code{DR1TTFAT}}{Total fat (gm)}
    \item{\code{DR1TSFAT}}{Total saturated fatty acids (gm)}
    \item{\code{DR1TMFAT}}{Total monounsaturated fatty acids (gm)}
    \item{\code{DR1TPFAT}}{Total polyunsaturated fatty acids (gm)}
    \item{\code{DR1TCAFF}}{Caffeine (mg)}
    \item{\code{DR1TALCO}}{Alcohol (gm)}
    \item{\code{DR1_320Z}}{Total plain water drank yesterday (gm)}
  }
}
\details{
The National Health and Nutrition Examination Survey (NHANES) is a program of studies designed to assess the health and nutritional status of adults and children in the United States. The survey is unique in that it combines interviews and physical examinations. The \code{nhis2007} data set contains observations for 4,329 persons collected in 2007-2008.
}
\source{
National Health and Nutrition Examination Survey of 2007-2008 conducted by the U.S. National Center for Health Statistics.
\url{http://www.cdc.gov/nchs/nhanes.htm}
}

\examples{
data(nhanes2007)
str(nhanes2007)
summary(nhanes2007)
}
\keyword{datasets}

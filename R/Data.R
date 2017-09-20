#' Oxygen Consumption of Leukocytes
#'
#' A dataset containing measurements on the oxygen consumption of leukocytes in the presence and absence of inactivated staphylococci.
#'
#' @format A data frame with 144 rows and 5 variables:
#' \describe{
#'   \item{O2}{oxygen consumption of leukocytes in \eqn{\mu}l}
#'   \item{Staphylococci}{whether or not inactivated staphylococci were added, 1 denotes yes, 0 no}
#'   \item{Time}{the measurements were taken after 6, 12 and 18 minutes}
#'   \item{Group}{the treatment group, either P for Placebo or V for Verum}
#'   \item{Subject}{the subject id}
#' }
#' 
#' @usage data(o2cons)
#' 
#' @source Friedrich, S., Brunner, E. and Pauly, M. (2016). Permuting longitudinal data despite all the dependencies. arXiv preprint arXiv:1509.05570v2
"o2cons"

#' EEG Measurements in Patients with Alzheimer's Disease (long format)
#'
#' At the Department of Neurology, University Clinic of Salzburg, 160 patients were diagnosed
#' with either AD, MCI, or SCC, based on neuropsychological diagnostics. This data set contains z-scores for brain rate and Hjorth complexity,
#' each measured at frontal, temporal and central electrode positions and averaged across hemispheres. In addition to standardization, complexity
#' values were multiplied by -1 in order to make them more easily comparable to brain rate
#' values: For brain rate we know that the values decrease with age and pathology, while
#' Hjorth complexity values are known to increase with age and pathology.
#' The three between-subjects factors considered were sex (men vs. women), diagnosis (AD
#' vs. MCI vs. SCC), and age (\eqn{< 70} vs. \eqn{>= 70} years). Additionally, the within-subjects factors region (frontal, temporal, central) and 
#' feature (brain rate, complexity) structure the response vector.
#'
#' @format A data frame with 960 rows and 7 variables:
#' \describe{
#'   \item{resp}{EEG measurements}
#'   \item{sex}{sex of the patient}
#'   \item{age}{age of the patient, coded as 0 for less than 70 years and 1 for \eqn{>= 70} years}
#'   \item{diagnosis}{neuropsychological diagnosis, AD for Alzheimer's Disease, MCI for mild cognitive impairment or SCC for subjective cognitive complaints without clinically significant deficits}
#'   \item{region}{brain region of the EEG measurements, one of "temporal", "frontal" and "central"}
#'   \item{feature}{feature of the EEG measurements, either "brainrate" or "complexity"}
#'   \item{id}{Subject id}
#' }
#' 
#' @usage data(EEG)
#' 
#' @source  Bathke, A., Friedrich, S., Konietschke, F., Pauly, M., Staffen, W., Strobl, N. and Hoeller, Y. (2016). Using EEG, SPECT, and Multivariate Resampling Methods
#' to Differentiate Between Alzheimer's and other Cognitive Impairments. arXiv preprint arXiv:1606.09004.
"EEG"


#' EEG Measurements in Patients with Alzheimer's Disease (wide format)
#'
#' At the Department of Neurology, University Clinic of Salzburg, 160 patients were diagnosed
#' with either AD, MCI, or SCC, based on neuropsychological diagnostics. This data set contains z-scores for brain rate and Hjorth complexity,
#' each measured at frontal, temporal and central electrode positions and averaged across hemispheres. In addition to standardization, complexity
#' values were multiplied by -1 in order to make them more easily comparable to brain rate
#' values: For brain rate we know that the values decrease with age and pathology, while
#' Hjorth complexity values are known to increase with age and pathology.
#' The three between-subjects factors considered were sex (men vs. women), diagnosis (AD
#' vs. MCI vs. SCC), and age (\eqn{< 70} vs. \eqn{>= 70} years). Additionally, the within-subjects factors region (frontal, temporal, central) and 
#' feature (brain rate, complexity) structure the response vector.
#'
#' @format A data frame with 160 rows and 9 variables:
#' \describe{
#'   \item{brainrate_temporal}{EEG measurements for brainrate in temporal regions}
#'   \item{brainrate_frontal}{EEG measurements for brainrate in frontal regions}
#'   \item{brainrate_central}{EEG measurements for brainrate in central regions}
#'   \item{complexity_temporal}{EEG measurements for complexity in temporal regions}
#'   \item{complexity_frontal}{EEG measurements for complexity in frontal regions}
#'   \item{complexity_central}{EEG measurements for complexity in central regions}
#'   \item{sex}{sex of the patient}
#'   \item{age}{age of the patient}
#'   \item{diagnosis}{neuropsychological diagnosis, AD for Alzheimer's Disease, MCI for mild cognitive impairment or SCC for subjective cognitive complaints without clinically significant deficits}
#'   \item{AgeGroup}{categorized age, coded as 0 for less than 70 years and 1 for \eqn{>= 70} years}
#' }
#' 
#' @usage data(EEGwide)
#' 
#' @source  Bathke, A., Friedrich, S., Konietschke, F., Pauly, M., Staffen, W., Strobl, N. and Hoeller, Y. (2016). Using EEG, SPECT, and Multivariate Resampling Methods
#' to Differentiate Between Alzheimer's and other Cognitive Impairments. arXiv preprint arXiv:1606.09004.
"EEGwide"
#' Moisture data
#' 
#' Moisture level at two stages in a food manufacturing process, in-process and final.
#' These data are given in Example 14.1 
#' 
#' 
#' @name moisture.df
#' @docType data
#' @format A data frame with 25 observations on the following 6 variables.
#' \itemize{ \item{batch. the batch number of the measurement}
#' \item{proc.level. the in-process moisture level} \item{final.level. natural
#' the final moisture level of the batch} \item{ls.fit the least squares fitted value
#' of final.level given proc.level} \item{residiual. the least squares residual}
#' \item{residiual.sq. the squaredvleast squares residual}}
#' @keywords datasets
#' @examples
#' 
#' data(moisture.df)
#' plot(final.level~proc.level, data = moisture.df)
#' 
#' 
"moisture.df"

#' Slug data
#' 
#' Lengths and weights of 100 slugs from the species Limax maximus collected
#' around Hamilton, New Zealand.
#' 
#' 
#' @name slug
#' @docType data
#' @format A data frame with 100 observations on the following 4 variables.
#' \itemize{ \item{length. length (mm) of the slug}
#' \item{weight. weight (g) of the slug} \item{log.len. natural
#' logarithm of the \code{length}} \item{log.wt. natural logarithm of
#' the \code{weight}} }
#' @references Barker, G. and McGhie, R. (1984). The Biology of Introduced
#' Slugs (Pulmonata) in New Zealand: Introduction and Notes on Limax Maximus,
#' NZ Entomologist 8, 106--111.
#' @keywords datasets
#' @examples
#' 
#' data(slug)
#' plot(weight~length, data = slug)
#' plot(log.wt~log.len, data = slug)
#' 
#' 
"slug"

#' Data for simple random sampling, stratified sampling, and clusting sampling
#' experiments
#' 
#' A simulated population made up of 100 individuals. The individuals come from
#' three ethnic groups with population proportions of 40\%, 40\%, and 20\%,
#' respectively. There are twenty neighborhoods, and five individuals live in
#' each one.  Now, the income distribution may be different for the three
#' ethnic groups.  Also, individuals in the same neighborhood tend to be more
#' similar than individuals in different neighborhoods.
#' 
#' 
#' @name sscsample.data
#' @docType data
#' @format A data frame with 100 observations on the following 3 variables.
#' \itemize{ \item{income. Simulated income in $10,000}
#' \item{ethnicity. A numerical vector indicating the ethnic group of
#' the observation} \item{neighborhood. A numeric vector indicating the
#' neighborhood of the observation} }
#' @keywords datasets
#' @examples
#' 
#' data(sscsample.data)
#' plot(income~ethnicity, data = sscsample.data)
#' 
"sscsample.data"

#' bears
#' 
#' Body measurements for 143 wild bears.
#' 
#' Wild bears were anesthetized, and their bodies were measured and weighed.
#' One goal of the study was to make a table (or perhaps a set of tables) for
#' hunters, so they could estimate the weight of a bear based on other
#' measurements. This would be used because in the forest it is easier to
#' measure the length of a bear, for example, than it is to weigh it.
#' 
#' @name bears
#' @docType data
#' @format A data frame with 143 observations on the following 12 variables.
#' \itemize{ \item{ID. Indentification number}
#' \item{Age. Bear's age, in months. Note, wild bears are always born
#' in January, so an expert can estimate the bear's age without directly asking
#' it how old it is. } \item{Month. Month when the measurement was
#' made. 1 = Jan., 12 = Dec. Since bears hibernate in the winter, their body
#' shape probably depends on the season. } \item{Sex. 1 = male 2 =
#' female } \item{Head.L. Length of the head, in inches }
#' \item{Head.W. Width of the head, in inches }
#' \item{Neck.G. Girth (distance around) the neck, in inches }
#' \item{Length. Body length, in inches } \item{Chest.G. Girth
#' (distance around) the chest, in inches } \item{Weight. Weight of the
#' bear, in pounds } \item{Obs.No. Observation number for this bear.
#' For example, the bear with ID = 41 (Bertha) was measured on four occasions,
#' in the months coded 7, 8, 11, and 5. The value of Obs.No goes from 1 to 4
#' for these observations. } \item{Name. The names of the bears given
#' to them by the researchers} }
#' @references This data set was supplied by Gary Alt. Entertaining references
#' are in Reader's Digest April, 1979, and Sports Afield September, 1981.
#' @source This data is in the example data set Bears.MTW distributed with
#' Minitab
#' @keywords datasets
#' @examples
#' 
#' data(bears)
#' boxplot(Weight~Sex, data = bears)
#' 
NULL





#' Bolstad Functions
#' 
#' A set of R functions and data sets for the book Introduction to Bayesian
#' Statistics, Bolstad, W.M. (2007), John Wiley & Sons ISBN 0-471-27020-2. Most
#' of the package functions replicate the Minitab macros that are provided with
#' the book. Some additional functions are provided to simplfy inference about
#' the posterior distribution of the parameters of interest.
#' 
#' \tabular{ll}{ Package: \tab Bolstad\cr Type: \tab Package\cr Version: \tab
#' 0.2-26\cr Date: \tab 2015-05-01\cr License: \tab GPL 2\cr }
#' 
#' @name Bolstad-package
#' @aliases Bolstad-package Bolstad
#' @docType package
#' @author James Curran Maintainer: James Curran <j.curran@@auckland.ac.nz> ~~
#' The author and/or maintainer of the package ~~
#' @references Bolstad, W.M. (2007), Introduction to Bayesian Statistics, John
#' Wiley & Sons.
#' @keywords package
"_PACKAGE"

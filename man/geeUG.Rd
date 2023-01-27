\name{geeUG}
\alias{geeUG}
\title{
Fit Unit Gamma Generalized Estimating Equations (GEEUG)
}
\description{
The geeUG function fits generalized estimating equations using the reparametrized unit gamma distribution (Mousa et al., 2016).
}
\usage{
geeUG(formula.mu, formula.phi, data, id, tol = 0.001,
                 maxiter = 25, corstr = "independence", linkmu = "logit",
                 linkphi = "log", silence = FALSE)
}
\arguments{
\item{formula.mu}{
an object of class "formula" (or one that can be coerced to that class) related to the mean predictor: a expression of the form "response ~ predictors". See the documentation of \code{\link{glm}} and \code{\link{formula}} for details.
}
\item{formula.phi}{
an object of class "formula" (or one that can be coerced to that class) related to the precision predictor: a expression of the form "~ predictors".
}
  \item{data}{
an optional data frame, list or environment containing the variables in the formula.
}
  \item{id}{
a vector which identifies the cluster of observations. The length of id should be the same as the number of observations.
}
  \item{tol}{
the tolerance used in the fitting algorithm.
}
  \item{maxiter}{
the maximum number of iterations.
}
  \item{corstr}{
a character string specifying the correlation structure for the mean model. The following are permitted: "exchangeable", "AR-1", "unstructured".
}
  \item{linkmu}{
The link function for mean model. The following are permitted: "logit", "cloglog", "probit" and "identity".
}
  \item{linkphi}{
The link function for precision model. The following are permitted: "log", "invsquare", "sqrt" and "identity".
}
  \item{silence}{
a logical argument: if FALSE the coefficients in each iteration will be printed.
}
}
\details{
It is important that the input data be sorted by the argument called "id" so that the repeated measures structure is correctly specified.
}
\value{
An object of class "geeUG" representing the fit.
}
\references{
MOUSA, A. M., EL-SHEIKH, A. A. and ABDEL-FATTAH, M. A. (2016). A gamma regression for bounded continuous variables. Advances and Applications in Statistics, 49, 305–326.

FREITAS, J.V., NOBRE, J., ESPINHEIRA, P., and RÊGO, L. (2023). Unit gamma regression models to correlated bounded data.
}
\author{
João de Freitas, jvbfreitas@ime.unicamp.br.
}

\examples{
dad = read.table("http://www.leg.ufpr.br/lib/exe/fetch.php/publications:papercompanions:dataset.txt")
dad$TRIM = factor(dad$TRIM)
dad$LOCAL = factor(dad$LOCAL)
dad$LOCAL = factor(dad$LOCAL, levels = levels(dad$LOCAL)[c(2,3,1)])
idna = as.numeric(rownames(dad[is.na(dad$y),]))
dados = dad[-idna,]
dados$LOCALL = dados$LOCAL
levels(dados$LOCALL) = c("updown","Reservoir","updown")
aux = geeUG(y~TRIM+LOCAL,~TRIM+LOCALL,dados,as.numeric(factor(dados$id)),
            corstr = "unstructured",
            tol = 0.00001^2,maxiter = 300,linkphi = "log",
            linkmu = "cloglog")
aux
summary(aux)
plot(aux)
dadaux = data.frame(aux$resq)
library(car)
library(qqplotr)
ggplot(dadaux,aes(sample = resq)) +
  stat_qq_band(bandType = "pointwise", fill = "#8DA0CB", alpha = 0.4) +
  stat_qq_line(colour = "#8DA0CB") +
  stat_qq_point() +  theme_bw()
}

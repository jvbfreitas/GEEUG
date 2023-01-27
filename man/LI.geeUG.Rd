\name{LI.geeUG}
\alias{LI.geeUG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute local influence measures for geeUG model
}
\description{
This function can be used to compute the local influence measures under certain perturbations for unit gamma generalized estimating equations based in Venezuela et al. (2011) and developed in Freitas et al.(2023).
}
\usage{
LI.geeUG(model, pert = c("case-weight", "response"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an R object returned by geeUG.
}
  \item{pert}{
a character string specifying the perturbation scheme. The following are permitted: "case-weight" and "response".
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A matrix with 2 columns is returned, where the first is the local influence related to the mean model and the second for the precision model.
}
\references{
VENEZUELA, M. K.; SANDOVAL, M. C.; BOTTER, D. A. (2011). Local influence in estimating equations. Computational Statistics & Data Analysis, v. 55, n. 4, p. 1867 – 1883.

FREITAS, J.V., NOBRE, J., ESPINHEIRA, P., and RÊGO, L. (2023). Unit gamma regression models to correlated bounded data.
}
\author{
João de Freitas, jvbfreitas@ime.unicamp.br.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
LI.geeUG(aux, pert = "case-weight")
LI.geeUG(aux, pert = "response")
}

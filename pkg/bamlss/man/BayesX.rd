\name{BayesX}
\alias{BayesX}
\alias{BayesX.control}
\alias{sx}
\alias{tx}
\alias{get_BayesXsrc}

\title{
  Markov Chain Monte Carlo for BAMLSS using \pkg{BayesX}
}

\description{
  This sampler function for BAMLSS is an interface to the \pkg{BayesX} (\url{http://www.BayesX.org})
  command-line binary from \R. The sampler is based on the command line version and functions
  provided in the \pkg{BayesXsrc} package.
}

\usage{
## Sampler function:
BayesX(x, y, family, start = NULL, weights = NULL, offset = NULL,
  data = NULL, control = BayesX.control(...), ...)

## Sampler control:
BayesX.control(n.iter = 1200, thin = 1, burnin = 200,
  seed = NULL, predict = "light", model.name = "bamlss",
  data.name = "d", prg.name = NULL, dir = NULL,
  verbose = FALSE, show.prg = TRUE, ...)

## Special BayesX smooth term constructor.
sx(x, z = NULL, bs = "ps", by = NA, ...)

## Special BayesX tensor product smooth term constructor.
tx(..., k = NA, constraint = c("main", "both", "none"))

## Download the newest version of BayesXsrc.
get_BayesXsrc(dir = NULL, install = TRUE)
}

\arguments{
  \item{x}{For function \code{BayesX()} the \code{x} list, as returned from
    function \code{\link{bamlss.frame}}, holding all model matrices and other information that is
    used for fitting the model. For function \code{sx()} arguments \code{x} and \code{z} specify
    the variables the smooth should be a function of.}
  \item{y}{The model response, as returned from function \code{\link{bamlss.frame}}.}
  \item{z}{Second variable in a \code{sx()} term.}
  \item{family}{A \pkg{bamlss} family object, see \code{\link{family.bamlss}}.}
  \item{start}{A named numeric vector containing possible starting values, the names are based on
    function \code{\link{parameters}}.}
  \item{weights}{Prior weights on the data, as returned from function \code{\link{bamlss.frame}}.}
  \item{offset}{Can be used to supply model offsets for use in fitting,
    returned from function \code{\link{bamlss.frame}}.}
  \item{data}{The model frame that should be used for modeling. Note that argument \code{data} needs
    not to be specified when the \code{BayesX()} sampler function is used with \code{\link{bamlss}}.}
  \item{control}{List of control arguments to be send to \pkg{BayesX}. See below.}
  \item{n.iter}{Sets the number of MCMC iterations.}
  \item{thin}{Defines the thinning parameter for MCMC simulation. E.g., \code{thin = 10} means,
    that only every 10th sampled parameter will be stored.}
  \item{burnin}{Sets the burnin phase of the sampler, i.e., the number of starting samples that
    should be removed.}
  \item{seed}{Sets the seed.}
  \item{predict}{Not supported at the moment, do not modify!}
  \item{model.name}{The name that should be used for the model when calling \pkg{BayesX}.}
  \item{data.name}{The name that should be used for the data set when calling \pkg{BayesX}.}
  \item{prg.name}{The name that should be used for the \code{.prg} file that is send to \pkg{BayesX}.}
  \item{dir}{Specifies the directory where \pkg{BayesX} should store all output files. For function
    \code{get_BayesXsrc()}, the directory where \pkg{BayesXsrc} should be stored.}
  \item{verbose}{Print information during runtime of the algorithm.}
  \item{show.prg}{Show the \pkg{BayesX} \code{.prg} file.}
  \item{bs}{A \code{\link{character}} string, specifying the basis/type which is used for
    this model term.}
  \item{by}{A by variable for varying coefficient model terms.}
  \item{k}{The dimension(s) of the bases used to represent the \code{tx()} smooth term.}
  \item{\dots}{Not used in \code{BayesX.control}. For function \code{sx()} any extra arguments that
    should be passed to \pkg{BayesX} for this model term can be specified here. For function
    \code{tx()}, all variables the smooth should be a function of are specified here. For function
    \code{BayesX()} all arguments that should be passed to \code{BayesX.control}.}
  \item{constraint}{Specifies the type of contraints that should be applied. \code{"main"}, both
    main effects should be removed; \code{"both"}, both main effects and varying effects should
    be removed; \code{"none"}, no constraint should be applied.}
  \item{install}{Should package \pkg{BayesXsrc} be installed?}
}

\details{
  Function \code{BayesX()} writes a \pkg{BayesX} \code{.prg} file and processes the data.
  The functions calls \pkg{BayesX} via function \code{\link[BayesXrsc]{run.bayesx}}. After
  the \pkg{BayesX} sampler has finished, the function reads back in all the parameter samples
  that can then be used for further processing within \code{\link{bamlss}}, i.a.

  The smooth term constructor functions \code{\link[mgcv]{s}} and \code{\link[mgcv]{te}} can
  be used with the \code{BayesX()} sampler. When using \code{\link[mgcv]{te}} note that only
  one smoothing variance is estimated by \pkg{BayesX}.

  For anisotropic penalties use function \code{tx()}, which currently supports smooth functions
  of two variables, only!
}

\value{
  Function \code{BayesX()} returns samples of parameters. The samples are provided as a
  \code{\link[coda]{mcmc}} matrix.

  Function \code{BayesX.control()} returns a \code{list} with control arguments for
  \pkg{BayesX}.

  Function \code{sx()} a \code{list} of class \code{"xx.smooth.spec"} and \code{"no.mgcv"}, where
  \code{"xx"} is a basis/type identifying code given by the \code{bs} argument.

  Function \code{tx()} a \code{list} of class \code{tensorX.smooth.spec}.
}

\note{
  Note that this interface is still experimental and needs the newest version of the \pkg{BayesX}
  source code, which is not yet part of the \pkg{BayesXsrc} package on CRAN. The newest version
  can be installed with function \code{get_BayesXsrc}. Note that the function assumes that sh,
  subversion (svn) and \R can be run from the command line!

  Note that for setting up a new family object to be used with \code{BayesX()} additional
  information needs to be supplied. The extra information must be placed within the
  family object in an named \code{list} element named \code{"bayesx"}. For each parameter of
  the distribution a character string with the corresponding \pkg{BayesX} \code{family} name and the
  \code{equationtype} must be supplied. See, e.g., the \R code of \code{\link{gaussian.bamlss}}
  how the setup works.

  For function \code{sx()} the following basis types are currently supported:
  \itemize{
    \item \code{"ps"}: P-spline with second order difference penalty. 
    \item \code{"mrf"}: Markov random fields: Defines a Markov random field prior for a
               spatial covariate, where geographical information is provided by a map object in
               boundary or graph file format (see function \code{\link{read.bnd}}, \code{\link{read.gra}} and
               \code{\link{shp2bnd}}), as an additional argument named \code{map}.
    \item \code{"re"}: Gaussian i.i.d. Random effects of a unit or cluster identification covariate.
  }

  Function \code{tx()} currently supports smooth terms with two variables.
}



\seealso{
  \code{\link{bamlss}}, \code{\link{bamlss.frame}}
}

\examples{
\dontrun{
## Get newest version of BayesXsrc.
## Note: needs sh, svn and R build tools!
get_BayesXsrc()

if(require("BayesXsrc")) {
  ## Simulate some data
  set.seed(123)
  d <- GAMart()

  ## Estimate model with BayesX. Note
  ## that BayesX computes staring values, so
  ## these are not required by some optimizer function
  ## in bamlss()
  b1 <- bamlss(num ~ s(x1) + s(x2) + s(x3) + s(lon,lat),
    data = d, optimizer = FALSE, sampler = BayesX)

  plot(b1)

  ## Same model with anisotropic penalty.
  b2 <- bamlss(num ~ s(x1) + s(x2) + s(x3) + tx(lon,lat),
    data = d, optimizer = FALSE, sampler = BayesX)

  plot(b2)
}
}
}

\keyword{regression}


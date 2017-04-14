\name{plot.sim}
\alias{plot.sim}

\title{ Plot function for \code{S3} class \code{"sim"} }

\description{ Visualizes structure of the true graph. }
\usage{
\method{plot}{sim}( x, main = NULL, layout = layout.circle, ... )
}

\arguments{
  \item{x}{ An object of \code{S3} class \code{"sim"}, from function \code{\link{bdgraph.sim}}.  }
  
  \item{main}{ Graphical parameter (see plot).}
  
  \item{layout}{The vertex placement algorithm which is according to \pkg{igraph} package. }
  
  \item{\dots}{System reserved (no specific usage).}
}
\references{
Mohammadi, A. and E. Wit (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{arXiv:1501.05108v2} 

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C} 
}

\author{ Abdolreza Mohammadi and Ernst Wit }

\seealso{\code{\link{bdgraph}}}

\examples{
\dontrun{
# Generating synthetic multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 10, p = 15 )
  
plot( data.sim )
}
}
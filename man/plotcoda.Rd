\name{plotcoda}
\alias{plotcoda}

\title{Convergence plot}

\description{
This function according to "bdgraph" output gives us several plots. These plots can be used to check the convergence of the BD-MCMC algorithm.
}

\usage{plotcoda(output, thin = NULL, trace = TRUE, main = NULL, ...)}

\arguments{
  \item{output}{An object with S3 class "bdgraph" }
  
  \item{thin}{An option for getting fast result for a cumulative plot according to part of the iteration.}
  
  \item{trace}{Logical: if TRUE (default), tracing information is printed.}
  
  \item{main}{ Graphical parameter (see plot). }
  
  \item{\dots}{ System reserved (no specific usage). }  
}
\details{
Note that a spending time for this function depends on the graph. It should be slow for the 
high-dimensional graphs. To make it faster you can choose bigger value for 'thin'.
}
\references{
Mohammadi, A. and E. C. Wit (2013). Bayesian structure learning in sparse Gaussian 
graphical models, arXiv:1210.5371v6. \url{http://arxiv.org/abs/1210.5371v6}
}
\author{Abdolreza Mohammadi and Ernst Wit}

\seealso{\code{\link{bdgraph}}}

\examples{
\dontrun{
  # generating synthetic multivariate normal data from a 'circle' graph
  data.sim <- bdgraph.sim( n = 50, p = 6, graph = "circle", vis = TRUE )
  
  output <- bdgraph( data = data.sim, iter = 2000, save.all = TRUE )
 
  plotcoda(output)
  }
}
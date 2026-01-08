\name{shiftNA}
\alias{shiftNA}
\title{NA-Padded Lag and Lead Operator}
\description{
Shifts a vector forward or backward by a specified number of positions and
fills out-of-range values with \code{NA}. This helper is designed for 
position-based transformations where missing values should be explicitly 
propagated as \code{NA} rather than assumed to be zero.
}
\usage{
shiftNA(v, k)
}
\arguments{
  \item{v}{A vector (numeric, character, etc.).}

  \item{k}{Integer. Shift order. Positive values correspond to lags
  (shifting the series downward, i.e.\ \eqn{t-k}), while negative values
  correspond to leads (shifting the series upward, i.e.\ \eqn{t+k}).}
}
\details{
This function performs a position-based shift of the input vector \code{v}.
Unlike strict time-indexed helpers (such as \code{\link{get_lag}}), \code{shiftNA()}
does not rely on an explicit time index and does not propagate gaps based on 
timestamps. Instead, it performs a simple index shift, padding the 
resulting empty slots with \code{NA}.

Let \eqn{n} denote the length of \code{v}. The behaviour is:
\itemize{
  \item \code{k = 0}: return \code{v} unchanged,
  \item \code{k > 0}: lag by \code{k} positions; the first \code{k} elements are \code{NA},
  \item \code{k < 0}: lead by \eqn{|k|} positions; the last \eqn{|k|} elements are \code{NA}.
}

Formally, for \code{k > 0},
\deqn{
\text{out}_t =
\begin{cases}
\text{NA}, & t \le k, \\
v_{t-k}, & t > k,
\end{cases}
}
and for \code{k < 0},
\deqn{
\text{out}_t =
\begin{cases}
v_{t+|k|}, & t \le n-|k|, \\
\text{NA}, & t > n-|k|.
\end{cases}
}
}
\value{
A vector of the same length as \code{v}, containing the shifted values
with \code{NA} inserted where the shift would otherwise go out of range.
}
\section{Usage and Propagation}{
This section explains the implications of using \code{NA} padding in 
recursive or difference-based calculations.

\subsection{Why NA-padding?}{
Using \code{NA} is the standard behavior in R for out-of-bounds operations. 
It ensures that subsequent calculations (like \code{v - shiftNA(v, 1)}) 
correctly result in \code{NA} for boundary cases, preventing the 
accidental use of arbitrary values (like zero) in statistical estimations 
unless explicitly intended.
}

\subsection{Difference from get_lag()}{
Unlike \code{\link{get_lag}}, \code{shiftNA()} is position-based and ignores
time indices. Use \code{shiftNA()} when you need a fast, simple shift on 
vectors that are already correctly ordered.
}
}
\examples{
v <- c(1, 2, 3, 4, 5)

## Lag by one (k = 1)
shiftNA(v, 1)
# [1] NA 1 2 3 4

## Lead by one (k = -1)
shiftNA(v, -1)
# [1] 2 3 4 5 NA

## Larger shifts
shiftNA(v, 2)
# [1] NA NA 1 2 3

shiftNA(v, -2)
# [1] 3 4 5 NA NA
}
\seealso{
\code{\link{WesterlundBootstrap}},
\code{\link{get_lag}},
\code{\link{get_diff}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotCoverageSampleProfiles.R
\name{plotCoverageSampleProfiles}
\alias{plotCoverageSampleProfiles}
\title{Plots coverage profiles over a specific region}
\usage{
plotCoverageSampleProfiles(rrbs, CpGs, group1, group2, CpGindex,
  samples = NaN, plot_title = NaN)
}
\arguments{
\item{rrbs}{An rrbs object containing methylation and coverage data as 
created using the BiSeq pacakge}

\item{CpGs}{A GRanges object with each row being a testing region}

\item{group1}{The name of the first testing group}

\item{group2}{The name of the second testing group}

\item{CpGindex}{The index within the CpGs object of the region we are 
plotting}

\item{samples}{Names of specific samples in the groups to plot, defaults to 
all samples in both groups}

\item{plot_title}{Optional title for the plot, overrides automatic title}
}
\value{
NULL, the function plots the profiles
}
\description{
Plots the profile for each of the two testing groups. 
Within each group, the mean of 
methylation level is taken, smoothed and plotted, along with the individual 
values.
}
\examples{
# plot the 9th region in the Toy Data Set
data(rrbsDemo)
data(CpGsDemo)
plotMethProfile(rrbsDemo, CpGsDemo, 'H1-hESC', 'K562', 9)
}
\author{
Tom Mayo \email{t.mayo@ed.ac.uk}
}


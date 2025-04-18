% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Heatmap.R
\name{Heatmap}
\alias{Heatmap}
\title{Heatmap with top annotation bar}
\usage{
Heatmap(
  data,
  col = c("green", "white", "red"),
  col_lim = c(-3, 0, 3),
  reverse = F,
  distance = "euclidean",
  rownames = TRUE,
  colnames = FALSE,
  Hsize = c(12, 6),
  g_legend = "Group",
  h_legend = "Intensity",
  T_size = 10,
  R_size = 7,
  C_size = 5,
  Gcol = c(ASD = "black", HGH = "red", LAC = "blue", LUE = "grey", SDF = "yellow", WEI =
    "green"),
  Title = "Sample Heatmap",
  dend_h = 0.5,
  a_h = 0.2
)
}
\arguments{
\item{data}{data inheriting from D_tran}

\item{col}{colors for heatmap}

\item{col_lim}{color boundaries c(-3, 0, 3)}

\item{reverse}{Reverse column and rows}

\item{distance}{"pearson", "manhattan","euclidean","spearman","kendall"}

\item{rownames}{rownames}

\item{colnames}{colnames}

\item{Hsize}{Width & Height c(a,b)}

\item{g_legend}{Annotation legend title}

\item{h_legend}{Heatmap legend title}

\item{T_size}{Title text size}

\item{R_size}{row text size}

\item{C_size}{column text size}

\item{Gcol}{Color for top_annotation bar c("ASD" = "black","HGH"="red","LAC"="blue","LUE" ="grey","SDF" = "yellow","WEI"="green")}

\item{Title}{title}

\item{dend_h}{dendrite height}

\item{a_h}{top annotation height}
}
\value{
Heatmap
}
\description{
Heatmap with top annotation bar
}
\examples{
data(Data)
data <- D_tran(Data, param = "Auto", save = F)
Heatmap(data)
}

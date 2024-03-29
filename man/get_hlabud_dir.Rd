% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hlabud.R
\name{get_hlabud_dir}
\alias{get_hlabud_dir}
\title{Get the name of the folder for caching downloaded IMGTHLA files}
\usage{
get_hlabud_dir()
}
\value{
The name of the folder.
}
\description{
This function will:
\itemize{
\item Get the folder name from \code{getOption("hlabud_dir")} or else automatically choose an appropriate folder for your operating system thanks to \href{https://github.com/r-lib/rappdirs}{rappdirs}.
\item Create the folder automatically if it does not already exist.
\item Set the the \code{hlabud_dir} option to that new folder.
}
}
\details{
Here are the locations of the \code{hlabud_dir} folder on each operating system.

Linux:

\if{html}{\out{<div class="sourceCode">}}\preformatted{~/.local/share/hlabud
}\if{html}{\out{</div>}}

Mac:

\if{html}{\out{<div class="sourceCode">}}\preformatted{~/Library/Application Support/hlabud
}\if{html}{\out{</div>}}

Windows:

\if{html}{\out{<div class="sourceCode">}}\preformatted{C:\\Documents and Settings\\\{User\}\\Application Data\\slowkow\\hlabud
}\if{html}{\out{</div>}}

To set the \code{hlabud_dir} option, please use:

\if{html}{\out{<div class="sourceCode">}}\preformatted{options(hlabud_dir = "/my/favorite/path")
}\if{html}{\out{</div>}}
}
\examples{
\dontrun{
hlabud_dir <- get_hlabud_dir()
}
}

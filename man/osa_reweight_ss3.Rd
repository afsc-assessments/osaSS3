\name{osa_reweight_ss3}
\alias{osa_reweight_ss3}
\title{Experimental iterative OSA-SDNR reweighting for SS3}
\usage{
osa_reweight_ss3(model_dir, ss_exe = NULL, max_iter = 8L,
  tolerance = 0.05, damping = 0.5, step_limit = 1.6,
  weight_bounds = c(0.02, 3), seed = 10L,
  output_dir = NULL, min_years = 3L, dry_run = FALSE,
  plot_final = FALSE, show_final_plots = FALSE)
}
\arguments{
\item{model_dir}{Directory with a runnable SS3 model and starter.ss.}
\item{ss_exe}{Path to SS3 executable; auto-discovered if omitted.}
\item{max_iter}{Maximum number of refits.}
\item{tolerance}{Absolute SDNR convergence tolerance around 1.}
\item{damping}{Exponent damping strength in (0,1].}
\item{step_limit}{Maximum multiplicative change per iteration.}
\item{weight_bounds}{Lower and upper bounds on adjusted weights.}
\item{seed}{Base random seed for OSA calculations.}
\item{output_dir}{New output directory. Must not already exist.}
\item{min_years}{Minimum number of composition years in a stratum.}
\item{dry_run}{Prepare the first run folder without running SS3.}
\item{plot_final}{Produce composite afscOSA diagnostic PDFs for the final fitted iteration.}
\item{show_final_plots}{Display final plots in an interactive session.}
}
\description{An experimental method; use only with appropriate model validation.}
\value{Invisibly returns a list with output_dir, history, converged, targets, and final_plot_files.}
\details{Supports standard multinomial length and marginal age composition variance adjustments, factors 4 and 5. Does not support conditional age-at-length or Dirichlet-multinomial models. Original input directory is left unchanged.}

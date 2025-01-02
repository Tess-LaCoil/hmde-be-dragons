# Here Be Dragons: Bimodal posteriors arise from numerical error in longitudinal models

This repo is an extended investigation of a quirk that I found while building 
[hmde](https://github.com/traitecoevo/hmde). The problem arises from an interaction
between Markov Chain Monte Carlo sampling and numerical error. When building the von Bertalanffy 
(affine first order ODE) model, I observed persistent bimodality in the posterior 
parameter estimates when using numerical integration to estimate growth increments. 
The fix for the package was to implement an analytic solution for the von 
Bertalanffy model. 

Here I document and demonstrate a more thorough investigation of that interaction.
I demonstrate that posterior bimodality in is a sticky problem for numerics
in the case of

$$f\left( Y \left( t \right), \beta_0, \beta_1 \right) = \frac{dY}{dt} = \beta_0 - \beta_1 Y \left( t \right), $$

across step sizes, different methods, priors, and parameter values. 
I also demonstrate that the Canham model in [hmde](https://github.com/traitecoevo/hmde)
is robust to similar problems with a suitable numerical method.

## Set-up
This repo requires the [hmde](https://github.com/traitecoevo/hmde) package in R.

The here_be_dragons.Rmd file provides a complete walk-through to reproduce the results.

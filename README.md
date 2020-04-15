This repository uses publicly available data to estimate the case ascertainment rate (proportion of cases that get recorded in the public data) of COVID-19 in Australia.

The general approach is that described in [the LSHTM analysis](https://cmmid.github.io/topics/covid19/severity/global_cfr_estimates.html): comparing the apparent Case Fatality Rate (CFR; the ratio of deaths to cases) to an assumed 'true', or baseline CFR from a more detailed study. If the apparent CFR from reported data is too high, that suggests that cases are being missed.

Estimating something that is not observed is a difficult problem, and this approach is far from perfect. Among other assumptions, this assumes that the number of deaths is reported perfectly, and that the baseline CFR is appropriate. Nevertheless, this os probably the best that can be done with publicly available data.

The analysis here replaces the method used for the global estimates linked to above with a Bayesian analysis, which enables a number of improvements:
  - The assumed baseline CFR is represented by a prior rather than a point estimate, correctly integrating uncertainty in this key parameter
  - The ascertainment rate parameter is comstrained to be between zero and one during modelling, rather than estimates exceeding 1 being truncated post-hoc. This is why the estimates and uncertainty intervals in the LSHTM analysis are exactly 1 in some cases, but always between 0 and 1 in this analysis.
  - Temporal variation is modelled as a Gaussian process, with uncertainty in the lengthscale parameters ('wiggliness' over time) accounted for by MCMC
  - This analysis infers different ascertainment rate timeseries for each Australian state and territory, but with information shared between them via a hierarchical Gaussian process prior.

Further extensions of this analysis, including differentiating ascertainment rates for different transmission sources (overseas-acquired, local transmission with a known source, and local transmission with unknown source) are in development on other branches.


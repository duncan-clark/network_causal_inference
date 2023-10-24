# Causal Inference over Stochastic Networks

This is the software and data used for the paper:

[Causal Inference over Stochastic Networks](https://arxiv.org/abs/2106.14145), by Duncan A. Clark and
Mark S. Handcock.

It has been submitted to the *Journal of the Royal Statistical Society*, Series A.

##Summary

Claiming causal inferences in network settings necessitates careful consideration
of the often complex dependency between outcomes for actors. Of particular importance are
treatment spillover or outcome interference effects. We consider causal inference when the
actors are connected via an underlying network structure. Our key contribution is a model for
causality when the underlying network is endogenous; where the ties between actors and the
actor covariates are statistically dependent. We develop a joint model for the relational and
covariate generating process that avoids restrictive separability and fixed network assumptions,
as these rarely hold in realistic social settings. While our framework can be used with general
models, we develop the highly expressive class of Exponential-family Random Network models
(ERNM) of which Markov Random Fields (MRF) and Exponential-family Random Graph models
(ERGM) are special cases. We present potential outcome based inference within a Bayesian
framework, and propose a modification to the exchange algorithm to allow for sampling from
ERNM posteriors. We present results of a simulation study demonstrating the validity of the
approach. Finally, we demonstrate the value of the framework in a case-study of smoking in the
context of adolescent friendship networks.

*Keywords*: Causality, Social Networks, Network models, Spillover, Contagion, Interference,
Gibbs measures


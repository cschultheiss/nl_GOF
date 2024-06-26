# Causal well-specification
R Software for assessing causal well-specification of additive noise models

multi_spec.R: implementation of the FOCI assessment using multiple splits <br>
Function multi.spec() can be applied to data <br>
split.R: contains the downstream proportion test <br>
Function fisher.split() can be applied to the matrix containing information if the covariate was not selected by FOCI

To create the data for the figures in the paper:
- Figure 4: run the function stored in random_simulation.R with its default parameters
- Figure 6: run the function stored in het_simulation.R with its default parameters

The called plotting functions with explaining comments are stored in figures_fun.R

To get Table 1 (in LaTeX encoding) for the analysis of the K562 et al.: run the function stored in k562.R <br>
The data must be preprocessed following the benchmark of [Chevalley et al. (2023)](https://arxiv.org/abs/2210.17283). <br>
Our results were obtained using version 1.4.1.1 of "xgboost". With other versions, they can differ a bit. <br>
Also, "FOCI" with non-unique data is not deterministic, and the implementation ignores seeds in R. We added some checks to enforce consistent results, but there might be further edge cases not detected yet.

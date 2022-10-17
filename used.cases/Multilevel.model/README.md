# Multilevel model (or, Linear Mixed Model)

## Great intro and tutorial 
[2018-A brief introduction to mixed effects modelling and multi model inference in ecology](https://peerj.com/articles/4794/#)

**Outline**
- Understanding fixed and random effects
  - controlling for non-independence among data points.
  - improving the accuracy of parameter estimation.
  Assuming all group means are drawn from a common distribution causes the estimates of their means to drift towards the global mean mgroup. This phenomenon, known as `shrinkage` (Gelman & Hill, 2007; K ́ery, 2010), can also lead to smaller and more precise standard errors around means. Shrinkage is strongest for groups with small sample sizes, as the paucity of within-group information to estimate the mean is counteracted by the model using data from other groups to improve the precision of the estimate. This `partial pooling` of the estimates is a principal benefit of fitting something as a random effect (Gelman & Hill, 2007).   
  - estimating variance components.
  `intra-class correlation coefficient`. In particular, quantitative genetic studies rely on variance component analysis for estimating the heritability of traits such as body mass or size of secondary sexual characteristics (Wilson et al., 2010).    
  - making predictions for unmeasured groups
  
- CONSIDERATIONS WHEN FITTING RANDOM EFFECTS
  - First, they are quite ‘data hungry’; requiring `at least five ‘levels’ (groups)` for a random intercept term to achieve robust estimates of variance (Gelman & Hill, 2007; Harrison, 2015).
  - Second, models can be unstable if sample sizes across groups are `highly unbalanced`, i.e. if some groups contain very few data. These issues are especially relevant to random slope models (Grueber et al., 2011). 
  - Third, an important issue is the difficulty in deciding the ‘significance’ or ‘importance’ of variance among groups. 
  - Finally, an issue that is not often addressed is that of mis-specification of random effects. GLMMs are powerful tools, but incorrectly parameterising the random effects in the model could yield model estimates that are as unreliable as ignoring the need for random effects altogether.
 
- DECIDING MODEL STRUCTURE FOR GLMMs
  - Choosing error structures and link functions
  - Choosing random effects I: crossed or nested?
   Here, female ID is said to be nested within woodland: each woodland contains multiple females `unique` to that woodland (that never move among woodlands). The nested random effect controls for the fact that (i) clutches from the same female are not independent, and (ii) females from the same woodland may have clutch masses more similar to one another than to females from other woodlands `Clutch Mass ∼ Foraging Rate + (1|Woodland/Female ID)`
  - Choosing random effects II: random slopes
  If fitting a random slope model including correlations between intercepts and slopes, always inspect the intercept-slope correlation coefficient in the variance/covariance summary returned by packages like lme4 to look for evidence of perfect correlations, indicative of insufficient data to estimate the model.
  - Choosing fixed effect predictors and interactions
  - How complex should my global model be?
    - Assessing predictor collinearity
     First, it can cause model convergence issues as models struggle to partition variance between predictor variables. Second, positively correlated variables can have negatively correlated regression coefficients, as the marginal effect of one is estimated, given the effect of the other, leading to incorrect interpretations of the direction of effects (Fig. 2). Third, collinearity can inflate standard errors of coefficient estimates and make ‘true’ effects harder to detect (Zuur, Ieno & Elphick, 2010). Finally, collinearity can affect the accuracy of model averaged parameter estimates during multi-model inference (Freckleton, 2011; Cade, 2015)... 
      When collinearity is detected, researchers can either select one variable as representative of multiple collinear variables (Austin, 2002), ideally using biological knowledge/reasoning to select the most meaningful variable (Zuur, Ieno & Elphick, 2010); or conduct a dimension-reduction analysis (e.g. Principal Components Analysis; James & McCullugh, 1990), leaving a single variable that accounts for most of the shared variance among the correlated variables. 
           
    - `Standardising and centring predictors`
    Transformations of predictor variables are common, and can improve **model performance and interpretability** (Gelman & Hill, 2007). *center and standardising*. Rescaling the mean of predictors containing large values (e.g. rainfall measured in 1,000s of millimetre) through centring/standardising will often `solve convergence problems`, in part because the estimation of intercepts is brought into the main body of the data themselves.
    Scaling is therefore a useful tool to `improve the stability of models and likelihood of model convergence, and the accuracy of parameter estimates if variables in a model are on large` (e.g. 1,000s of millimetre of rainfall), or vastly different scales. When using scaling, care must be taken in the interpretation and graphical representation of outcomes.
    
  - Quantifying GLMM fit and performance
    - Inspection of residuals and linear model assumptions
    - Overdispersion
    - R^2
    - Stability of `variance components` and testing significance of random effects
  `REML (restricted/residual maximum likelihood)` should be used for `estimating variance components of random effects` in Gaussian models as it produces less biased estimates compared to maximum likelihood (ML) (Bolker et al., 2009). However, when `comparing two models` with the same random structure but different fixed effects, `ML estimation` cannot easily be avoided.
    - Assessing model fit through simulation
    - Dealing with missing data

- MODEL SELECTION AND MULTI-MODEL INFERENCE
  - Stepwise selection, LRTs and p values
  - Information-theory and multi-model inference
- PRACTICAL ISSUES WITH APPLYING INFORMATION THEORY TO BIOLOGICAL DATA
  - Using all-subsets selection
  - Deciding which information criterion to use
  - Choice of AIC threshold
  - Using the nesting rule to improve inference from the top model set
  - Using akaike weights to quantify variable importance
  - Model averaging when predictors are collinear
  
- CONCLUSION
1. Modern mixed effect models offer an unprecedented opportunity to explore complex biological problems by explicitly modelling non-Normal data structures and/or non- independence among observational units. However, the LMM and GLMM toolset should be used with caution.
2. Rigorous testing of both model fit (R2) and model adequacy (violation of assumptions like homogeneity of variance) must be carried out. We must recognise that satisfactory fit does not guarantee we have not violated the assumptions of LMM, and vice versa. Interpret measures of R2 for (G)LMMs with hierarchical errors cautiously, especially when OLRE are used.
3. Collinearity among predictors is difficult to deal with and can severely impair model accuracy. Be especially vigilant if data are from field surveys rather than controlled experiments, as collinearity is likely to be present.
4. When including a large number of predictors is necessary, backwards selection and NHST should be avoided, and ranking via AIC of all competing models is preferred. A critical question that remains to be addressed is whether model selection based on IT is superior to NHST even in cases of balanced experimental designs with few predictors.
5. Data simulation is a powerful but underused tool. If the analyst harbours any uncertainty regarding the fit or adequacy of the model structure, then the analysis of data simulated to recreate the perceived structure of the favoured model can provide reassurance, or justify doubt.
6. Wherever possible, provide diagnostic assessment of model adequacy, and metrics of model fit, even if in the supplemental information.  
    

## Datasets
- `flowerdata.csv` from [2018 How to analyse plant phenotypic plasticity in response to a changing climate](https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.15656) Use fitted random slope to assess the variation of plasticity across genotypes, i.e., how a phenotype (plastic trait) respond to temperature change (large slope ~ high change, zero slope ~ plasticity), and how this plasticity vary across genotypes.
- `dragons.RData` from [INTRODUCTION TO LINEAR MIXED MODELS](https://ourcodingclub.github.io/tutorials/mixed-models/)
-  `sleepstudy` [Plotting partial pooling in mixed-effects models](https://www.tjmahr.com/plotting-partial-pooling-in-mixed-effects-models/) More on `complete pooling, no pooling, and partial pooling` and visualize the `shrinkage` effect of `partial pooling` approach.
- `R package variancePartition` from [variancePartition: interpreting drivers of variation in complex gene expression studies](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1323-z)


## Note
[When Mixed Effects (Hierarchical) Models Fail: Pooling and Uncertainty](https://towardsdatascience.com/when-mixed-effects-hierarchical-models-fail-pooling-and-uncertainty-77e667823ae8)

- Complete pooling (or simple linear regression)
- No-pooling model (or separate linear regressions)
- Partial-pooling model (or linear mixed effects)

***When partial-pooling goes wrong***
In essence what we see is smaller groups “borrow information” from larger groups and get their estimates shrunk towards the population values.
1. Vanishing confidence intervals
>Simply put, this warning tells us that the model is too complex for the data to support. In our case, in counties with smaller sample sizes, the random intercept accounts for most of the variability so when we try to fit a county-specific slope there isn’t enough information left to do so.
2. Unintended Shrinkage (cautionary tale) 
>Before moving on it is important to note a crucial lesson about how partial-pooling/shrinkage might lead to unintended consequences. In 2016, George et al. showed that partial pooling can lead to artificially low predicted mortality rates for smaller hospitals. As seen in the figures below, the predicted mortality rates (y-axis) have been “artificially shrunk” for hospitals with smaller numbers of beds (x-axis). That is to say, the model is systematically down-weighting the observed (higher) mortality in smaller hospitals. In this scenario, partial-pooling from a standard mixed effects model has potentially nightmarish public health and legal consequences. This under-prediction of the observed mortality rate could affect patients down the line.

**Bayes to the rescue**
>Bayesian mixed effects method are powerful alternatives to more commonly used Frequentist approaches (as above). Particularly, they are able to account for uncertainty when estimating group-level effects and provide stable estimates for groups with smaller sample sizes with the help of weakly informative priors.

**Conclusions**
Here’s some final thoughts to consider
1. Mixed effects models are powerful methods that help us account for complex, nested structures in our data
2. Keep in mind that, depending on the context/problem, partial-pooling can hurt you rather than help you
3. Consider using a Bayesian framework when appropriate (I’ll cover this in detail in the next article)


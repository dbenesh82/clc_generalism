Look for tradeoffs driving host specificity
================

Why are complex life cycle parasites specific at one stage but not
another? Here I explore tradeoffs that might affect the evolution of
host specificity in complex life cycles.

[Elsewhere](../stage_level_analyses/make_stage_level_df.Rmd), I created
a data table at the level of stages within parasite species. I start by
looking at correlations between variables measured at this level, like
host records, taxonomic dissimilarity, development time, growth, etc.
Here’s a correlation matrix, but it might not be too useful, because
several variables need log transformation.

    ##                              study_effort num_hosts_suspicious_removed
    ## study_effort                         1.00                         0.07
    ## num_hosts_suspicious_removed         0.07                         1.00
    ## hsi_lcdb_suspcious_rem               0.10                         0.44
    ## avg_dd                               0.09                        -0.03
    ## rel_growth_len                       0.10                        -0.10
    ## rel_growth_biov                      0.12                        -0.04
    ##                              hsi_lcdb_suspcious_rem avg_dd rel_growth_len
    ## study_effort                                   0.10   0.09           0.10
    ## num_hosts_suspicious_removed                   0.44  -0.03          -0.10
    ## hsi_lcdb_suspcious_rem                         1.00  -0.03          -0.21
    ## avg_dd                                        -0.03   1.00           0.34
    ## rel_growth_len                                -0.21   0.34           1.00
    ## rel_growth_biov                               -0.19   0.47           0.81
    ##                              rel_growth_biov
    ## study_effort                            0.12
    ## num_hosts_suspicious_removed           -0.04
    ## hsi_lcdb_suspcious_rem                 -0.19
    ## avg_dd                                  0.47
    ## rel_growth_len                          0.81
    ## rel_growth_biov                         1.00

There is less development time data than growth data. And growth was
calculated with two different size variables: length and biovolume.
Biovolume had more missing data, because biovolume was calculated with
lengths and widths and there were more missing widths than lengths.

    ##          avg_dt          avg_dd  abs_growth_len abs_growth_biov 
    ##            1200            1296             644             914

Growth was less commonly reported in paratenic hosts, so I made some
assumptions to impute this data. For any paratenic host stage without
growth data, I set the growth to zero.

# Models

I want to test whether life history variables can explain variation in
generalism, which could provide evidence for a tradeoff. There are
several groups of variables to test in the model: (1) uninteresting
confounders (phylogeny, study effort), (2) life cycle variables (def vs
int, host number), and (3) life history variables (growth and
development). I first add the uninteresting confounders to the model,
simply to control for them and to be sure that variables added
afterwards explain *additional* variation. Here, I’m mainly interested
in the life history variables. I could add them to the model either
before or after the life cycle variables. If I add them after, I’m
testing whether life history explains variation beyond that explained by
the life cycle. If I add them before, I’m testing for life history
tradeoffs, but their effects are possibly dependent/confounded with life
cycle characteristics. I’ll take both approaches below.

I retained the same model structure as for the [stage-level
analyses](../stage_level_analyses/stage_level_analysis_host_range_freq.Rmd),
specifically I used a generalized linear mixed model with Poisson
errors. An observation-level random effect accounted for overdispersion.

# Host range

### Case 1

I look first at host range. The *Case 1* approach is the simpler one.
Ignoring life cycle variables, does life history matter for stage
generalism? I fit 3 linear mixed models: (1) start with taxonomic random
effects, (2) add study effort, and (3) add the life history variable, in
this first case relative growth.

#### Relative growth

Adding relative growth to a model with just parasite taxonomy and study
effort is a slight improvement.

    ## Data: filter(st_level, !is.na(rel_growth_biov))
    ## Models:
    ## reg00: num_hosts_suspicious_removed ~ 1 + (1 | obs) + (1 | Parasite.species)
    ## reg0: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg0:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg0:     (1 | parasite_class) + (1 | parasite_phylum)
    ## reg1: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg1:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg1:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort
    ## reg2: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg2:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg2:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort + 
    ## reg2:     rel_growth_biov
    ##       Df    AIC    BIC  logLik deviance    Chisq Chi Df Pr(>Chisq)    
    ## reg00  3 5466.0 5481.0 -2730.0   5460.0                               
    ## reg0   8 5451.4 5491.4 -2717.7   5435.4  24.5673      5  0.0001689 ***
    ## reg1   9 5340.0 5385.0 -2661.0   5322.0 113.4100      1  < 2.2e-16 ***
    ## reg2  10 5337.4 5387.4 -2658.7   5317.4   4.6266      1  0.0314805 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Here’s the regression parameters.

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: poisson  ( log )
    ## Formula: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) +  
    ##     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) +  
    ##     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort +  
    ##     rel_growth_biov
    ##    Data: filter(st_level, !is.na(rel_growth_biov))
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   5337.4   5387.4  -2658.7   5317.4     1089 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.3560 -0.5290 -0.1386  0.3106  1.2221 
    ## 
    ## Random effects:
    ##  Groups           Name        Variance  Std.Dev. 
    ##  obs              (Intercept) 3.661e-01 6.051e-01
    ##  Parasite.species (Intercept) 5.796e-03 7.613e-02
    ##  parasite_genus   (Intercept) 2.881e-02 1.697e-01
    ##  parasite_family  (Intercept) 8.394e-02 2.897e-01
    ##  parasite_order   (Intercept) 2.346e-06 1.532e-03
    ##  parasite_class   (Intercept) 4.526e-10 2.127e-05
    ##  parasite_phylum  (Intercept) 2.325e-08 1.525e-04
    ## Number of obs: 1099, groups:  
    ## obs, 1099; Parasite.species, 653; parasite_genus, 326; parasite_family, 110; parasite_order, 31; parasite_class, 6; parasite_phylum, 3
    ## 
    ## Fixed effects:
    ##                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)      1.21463    0.06191  19.620   <2e-16 ***
    ## zstudy_effort    0.46963    0.04089  11.486   <2e-16 ***
    ## rel_growth_biov -0.03798    0.01768  -2.148   0.0317 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) zstdy_
    ## zstudy_ffrt -0.058       
    ## rl_grwth_bv -0.636 -0.071
    ## convergence code: 0
    ## boundary (singular) fit: see ?isSingular

The R<sup>2</sup> improvement through adding parasite growth is very
small.

    ## # A tibble: 4 x 5
    ##   step                df_used marg_r2 cond_r2 rand_var_explained
    ##   <chr>                 <dbl>   <dbl>   <dbl>              <dbl>
    ## 1 within-species corr      NA   0       0.265              0.265
    ## 2 taxonomy                  0   0       0.285              0.285
    ## 3 study effort              1   0.13    0.272              0.142
    ## 4 growth                    1   0.132   0.268              0.136

The same pattern is seen if we assume worms do not grow in paratenic
hosts. Surprisingly, though, the effect is not much stronger.

    ## Data: filter(st_level, !is.na(rel_growth_paratenic))
    ## Models:
    ## reg00: num_hosts_suspicious_removed ~ 1 + (1 | obs) + (1 | Parasite.species)
    ## reg0: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg0:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg0:     (1 | parasite_class) + (1 | parasite_phylum)
    ## reg1: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg1:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg1:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort
    ## reg2: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg2:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg2:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort + 
    ## reg2:     rel_growth_paratenic
    ##       Df    AIC    BIC  logLik deviance    Chisq Chi Df Pr(>Chisq)    
    ## reg00  3 6056.2 6071.5 -3025.1   6050.2                               
    ## reg0   8 6034.4 6075.1 -3009.2   6018.4  31.8608      5   6.33e-06 ***
    ## reg1   9 5919.7 5965.5 -2950.8   5901.7 116.6614      1  < 2.2e-16 ***
    ## reg2  10 5917.4 5968.4 -2948.7   5897.4   4.2868      1    0.03841 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## # A tibble: 4 x 5
    ##   step                df_used marg_r2 cond_r2 rand_var_explained
    ##   <chr>                 <dbl>   <dbl>   <dbl>              <dbl>
    ## 1 within-species corr      NA   0       0.271              0.271
    ## 2 taxonomy                  0   0       0.276              0.276
    ## 3 study effort              1   0.124   0.265              0.141
    ## 4 growth                    1   0.128   0.262              0.134

Here’s the regression parameters.

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: poisson  ( log )
    ## Formula: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) +  
    ##     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) +  
    ##     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort +  
    ##     rel_growth_paratenic
    ##    Data: filter(st_level, !is.na(rel_growth_paratenic))
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   5917.4   5968.4  -2948.7   5897.4     1196 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.3397 -0.5268 -0.1244  0.3026  1.1681 
    ## 
    ## Random effects:
    ##  Groups           Name        Variance  Std.Dev. 
    ##  obs              (Intercept) 3.685e-01 0.6070201
    ##  Parasite.species (Intercept) 2.330e-02 0.1526515
    ##  parasite_genus   (Intercept) 2.513e-02 0.1585147
    ##  parasite_family  (Intercept) 6.730e-02 0.2594277
    ##  parasite_order   (Intercept) 5.801e-06 0.0024086
    ##  parasite_class   (Intercept) 7.088e-07 0.0008419
    ##  parasite_phylum  (Intercept) 1.538e-07 0.0003921
    ## Number of obs: 1206, groups:  
    ## obs, 1206; Parasite.species, 696; parasite_genus, 338; parasite_family, 112; parasite_order, 31; parasite_class, 6; parasite_phylum, 3
    ## 
    ## Fixed effects:
    ##                      Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)           1.20043    0.05620  21.360   <2e-16 ***
    ## zstudy_effort         0.45715    0.03950  11.575   <2e-16 ***
    ## rel_growth_paratenic -0.03433    0.01654  -2.075    0.038 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) zstdy_
    ## zstudy_ffrt -0.093       
    ## rl_grwth_pr -0.602 -0.045
    ## convergence code: 0
    ## Model failed to converge with max|grad| = 0.0418042 (tol = 0.001, component 1)

The negative correlation is mainly driven by paratenic hosts. If we
exclude them, the relationship between generalism and growth is still
negative, but it is weaker and not significant.

    ## Data: filter(st_level, !is.na(rel_growth_biov), Facultative == "no")
    ## Models:
    ## reg0: num_hosts_suspicious_removed ~ 1 + (1 | obs) + (1 | Parasite.species) + 
    ## reg0:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg0:     (1 | parasite_class) + (1 | parasite_phylum)
    ## reg1: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg1:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg1:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort
    ## reg2: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg2:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg2:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort + 
    ## reg2:     rel_growth_paratenic
    ##      Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
    ## reg0  8 5031.5 5071.0 -2507.8   5015.5                              
    ## reg1  9 4938.5 4983.0 -2460.2   4920.5 95.0035      1     <2e-16 ***
    ## reg2 10 4939.9 4989.4 -2460.0   4919.9  0.5588      1     0.4547    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#### Development time - days

Moving on to development time expressed as days to maturity. There is
not a clear relationship between host range and developmental time.

    ## Data: filter(st_level, !is.na(avg_dt))
    ## Models:
    ## reg00: num_hosts_suspicious_removed ~ 1 + (1 | obs) + (1 | Parasite.species)
    ## reg0: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg0:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg0:     (1 | parasite_class) + (1 | parasite_phylum)
    ## reg1: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg1:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg1:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort
    ## reg2: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg2:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg2:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort + 
    ## reg2:     log10(avg_dt)
    ##       Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
    ## reg00  3 4119.0 4133.1 -2056.5   4113.0                              
    ## reg0   8 4110.5 4148.2 -2047.2   4094.5 18.4970      5   0.002384 ** 
    ## reg1   9 4046.9 4089.3 -2014.5   4028.9 65.5607      1  5.635e-16 ***
    ## reg2  10 4048.4 4095.4 -2014.2   4028.4  0.5573      1   0.455367    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Here’s the model output.

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: poisson  ( log )
    ## Formula: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) +  
    ##     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) +  
    ##     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort +  
    ##     log10(avg_dt)
    ##    Data: filter(st_level, !is.na(avg_dt))
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   4048.4   4095.4  -2014.2   4028.4      807 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.3155 -0.5459 -0.1294  0.3078  1.0687 
    ## 
    ## Random effects:
    ##  Groups           Name        Variance  Std.Dev. 
    ##  obs              (Intercept) 3.799e-01 0.6163770
    ##  Parasite.species (Intercept) 1.594e-05 0.0039925
    ##  parasite_genus   (Intercept) 3.084e-02 0.1756062
    ##  parasite_family  (Intercept) 8.554e-02 0.2924678
    ##  parasite_order   (Intercept) 1.588e-03 0.0398452
    ##  parasite_class   (Intercept) 6.787e-07 0.0008238
    ##  parasite_phylum  (Intercept) 3.777e-07 0.0006146
    ## Number of obs: 817, groups:  
    ## obs, 817; Parasite.species, 609; parasite_genus, 293; parasite_family, 100; parasite_order, 26; parasite_class, 6; parasite_phylum, 3
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)    1.18187    0.15619   7.567 3.83e-14 ***
    ## zstudy_effort  0.34265    0.04291   7.986 1.39e-15 ***
    ## log10(avg_dt) -0.06789    0.09396  -0.723     0.47    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) zstdy_
    ## zstudy_ffrt -0.085       
    ## lg10(vg_dt) -0.882 -0.101
    ## convergence code: 0
    ## Model failed to converge with max|grad| = 0.0493249 (tol = 0.001, component 1)

Little variation was explained by development time.

    ## # A tibble: 4 x 5
    ##   step                df_used marg_r2 cond_r2 rand_var_explained
    ##   <chr>                 <dbl>   <dbl>   <dbl>              <dbl>
    ## 1 within-species corr      NA   0       0.143              0.143
    ## 2 taxonomy                  0   0       0.17               0.17 
    ## 3 study effort              1   0.092   0.232              0.14 
    ## 4 devo time                 1   0.093   0.232              0.139

If we assume development in paratenic hosts is short (1 day), then the
relationship is significant.

    ## Data: filter(st_level, !is.na(avg_dt_paratenic))
    ## Models:
    ## reg00: num_hosts_suspicious_removed ~ 1 + (1 | obs) + (1 | Parasite.species)
    ## reg0: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg0:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg0:     (1 | parasite_class) + (1 | parasite_phylum)
    ## reg1: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg1:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg1:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort
    ## reg2: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg2:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg2:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort + 
    ## reg2:     log10(avg_dt_paratenic)
    ##       Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
    ## reg00  3 5067.3 5082.0 -2530.7   5061.3                              
    ## reg0   8 5040.0 5079.1 -2512.0   5024.0 37.3615      5  5.069e-07 ***
    ## reg1   9 4960.9 5004.9 -2471.4   4942.9 81.0842      1  < 2.2e-16 ***
    ## reg2  10 4956.9 5005.8 -2468.5   4936.9  5.9754      1    0.01451 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Here’s the model output.

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: poisson  ( log )
    ## Formula: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) +  
    ##     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) +  
    ##     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort +  
    ##     log10(avg_dt_paratenic)
    ##    Data: filter(st_level, !is.na(avg_dt_paratenic))
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   4956.9   5005.8  -2468.5   4936.9      968 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.3202 -0.5358 -0.1213  0.2954  1.0866 
    ## 
    ## Random effects:
    ##  Groups           Name        Variance  Std.Dev. 
    ##  obs              (Intercept) 3.910e-01 6.253e-01
    ##  Parasite.species (Intercept) 1.713e-06 1.309e-03
    ##  parasite_genus   (Intercept) 3.281e-02 1.811e-01
    ##  parasite_family  (Intercept) 8.814e-02 2.969e-01
    ##  parasite_order   (Intercept) 5.288e-03 7.272e-02
    ##  parasite_class   (Intercept) 0.000e+00 0.000e+00
    ##  parasite_phylum  (Intercept) 7.853e-09 8.862e-05
    ## Number of obs: 978, groups:  
    ## obs, 978; Parasite.species, 656; parasite_genus, 307; parasite_family, 102; parasite_order, 26; parasite_class, 6; parasite_phylum, 3
    ## 
    ## Fixed effects:
    ##                         Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              1.24910    0.08972  13.922   <2e-16 ***
    ## zstudy_effort            0.36907    0.03888   9.492   <2e-16 ***
    ## log10(avg_dt_paratenic) -0.11861    0.04792  -2.475   0.0133 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) zstdy_
    ## zstudy_ffrt -0.224       
    ## lg10(vg_d_) -0.692  0.024
    ## convergence code: 0
    ## boundary (singular) fit: see ?isSingular

The effect size is still small, though, explaining less than 1% of the
variation.

    ## # A tibble: 4 x 5
    ##   step                df_used marg_r2 cond_r2 rand_var_explained
    ##   <chr>                 <dbl>   <dbl>   <dbl>              <dbl>
    ## 1 within-species corr      NA   0       0.222              0.222
    ## 2 taxonomy                  0   0       0.247              0.247
    ## 3 study effort              1   0.104   0.265              0.161
    ## 4 devo time                 1   0.108   0.253              0.145

#### Development time - degree days

Moving on to development time expressed in temp-corrected degree days.
Short development times were associated with higher generalism, which is
consistent with a tradeoff.

    ## Data: filter(st_level, !is.na(avg_dd))
    ## Models:
    ## reg00: num_hosts_suspicious_removed ~ 1 + (1 | obs) + (1 | Parasite.species)
    ## reg0: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg0:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg0:     (1 | parasite_class) + (1 | parasite_phylum)
    ## reg1: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg1:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg1:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort
    ## reg2: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg2:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg2:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort + 
    ## reg2:     log10(avg_dd)
    ##       Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
    ## reg00  3 3638.3 3652.0 -1816.2   3632.3                              
    ## reg0   8 3638.6 3675.3 -1811.3   3622.6  9.6578      5    0.08553 .  
    ## reg1   9 3583.1 3624.3 -1782.5   3565.1 57.5838      1  3.239e-14 ***
    ## reg2  10 3580.2 3626.0 -1780.1   3560.2  4.8887      1    0.02703 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Here’s the model output.

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: poisson  ( log )
    ## Formula: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) +  
    ##     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) +  
    ##     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort +  
    ##     log10(avg_dd)
    ##    Data: filter(st_level, !is.na(avg_dd))
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   3580.2   3626.0  -1780.1   3560.2      711 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.2972 -0.5525 -0.1250  0.3274  1.0325 
    ## 
    ## Random effects:
    ##  Groups           Name        Variance  Std.Dev.
    ##  obs              (Intercept) 3.799e-01 0.616345
    ##  Parasite.species (Intercept) 9.059e-05 0.009518
    ##  parasite_genus   (Intercept) 3.945e-02 0.198614
    ##  parasite_family  (Intercept) 5.263e-02 0.229412
    ##  parasite_order   (Intercept) 8.734e-05 0.009345
    ##  parasite_class   (Intercept) 7.456e-06 0.002731
    ##  parasite_phylum  (Intercept) 1.326e-05 0.003641
    ## Number of obs: 721, groups:  
    ## obs, 721; Parasite.species, 550; parasite_genus, 273; parasite_family, 94; parasite_order, 23; parasite_class, 6; parasite_phylum, 3
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)    1.63907    0.25214   6.501 8.00e-11 ***
    ## zstudy_effort  0.34468    0.04281   8.052 8.15e-16 ***
    ## log10(avg_dd) -0.19244    0.08997  -2.139   0.0324 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) zstdy_
    ## zstudy_ffrt  0.113       
    ## lg10(vg_dd) -0.979 -0.180
    ## convergence code: 0
    ## Model failed to converge with max|grad| = 0.705696 (tol = 0.001, component 1)

Despite the significant effect, little variation was explained by
development time.

    ## # A tibble: 4 x 5
    ##   step                df_used marg_r2 cond_r2 rand_var_explained
    ##   <chr>                 <dbl>   <dbl>   <dbl>              <dbl>
    ## 1 within-species corr      NA   0       0.132              0.132
    ## 2 taxonomy                  0   0       0.143              0.143
    ## 3 study effort              1   0.093   0.211              0.118
    ## 4 devo time                 1   0.096   0.209              0.113

If we assume development in paratenic hosts is short (1 day), then the
relationship is much stronger.

    ## Data: filter(st_level, !is.na(avg_dd_paratenic))
    ## Models:
    ## reg00: num_hosts_suspicious_removed ~ 1 + (1 | obs) + (1 | Parasite.species)
    ## reg0: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg0:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg0:     (1 | parasite_class) + (1 | parasite_phylum)
    ## reg1: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg1:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg1:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort
    ## reg2: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg2:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg2:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort + 
    ## reg2:     log10(avg_dd_paratenic)
    ##       Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
    ## reg00  3 4586.6 4601.0 -2290.3   4580.6                             
    ## reg0   8 4570.2 4608.5 -2277.1   4554.2 26.399      5  7.468e-05 ***
    ## reg1   9 4499.1 4542.1 -2240.5   4481.1 73.136      1  < 2.2e-16 ***
    ## reg2  10 4489.2 4537.0 -2234.6   4469.2 11.887      1  0.0005652 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Here’s the model output.

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: poisson  ( log )
    ## Formula: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) +  
    ##     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) +  
    ##     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort +  
    ##     log10(avg_dd_paratenic)
    ##    Data: filter(st_level, !is.na(avg_dd_paratenic))
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   4489.2   4537.0  -2234.6   4469.2      872 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.3384 -0.5346 -0.1113  0.3152  0.9888 
    ## 
    ## Random effects:
    ##  Groups           Name        Variance  Std.Dev. 
    ##  obs              (Intercept) 3.915e-01 6.257e-01
    ##  Parasite.species (Intercept) 1.987e-06 1.410e-03
    ##  parasite_genus   (Intercept) 4.068e-02 2.017e-01
    ##  parasite_family  (Intercept) 6.456e-02 2.541e-01
    ##  parasite_order   (Intercept) 4.578e-08 2.140e-04
    ##  parasite_class   (Intercept) 7.412e-10 2.722e-05
    ##  parasite_phylum  (Intercept) 0.000e+00 0.000e+00
    ## Number of obs: 882, groups:  
    ## obs, 882; Parasite.species, 605; parasite_genus, 292; parasite_family, 98; parasite_order, 24; parasite_class, 6; parasite_phylum, 3
    ## 
    ## Fixed effects:
    ##                         Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              1.55265    0.13041  11.906  < 2e-16 ***
    ## zstudy_effort            0.36548    0.03964   9.219  < 2e-16 ***
    ## log10(avg_dd_paratenic) -0.16507    0.04744  -3.480 0.000502 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) zstdy_
    ## zstudy_ffrt -0.086       
    ## lg10(vg_d_) -0.923 -0.033
    ## convergence code: 0
    ## boundary (singular) fit: see ?isSingular

The effect size is still small, though, explaining less than 1% of the
variation.

    ## # A tibble: 4 x 5
    ##   step                df_used marg_r2 cond_r2 rand_var_explained
    ##   <chr>                 <dbl>   <dbl>   <dbl>              <dbl>
    ## 1 within-species corr      NA   0       0.222              0.222
    ## 2 taxonomy                  0   0       0.234              0.234
    ## 3 study effort              1   0.105   0.257              0.152
    ## 4 devo time                 1   0.113   0.237              0.124

#### Relative growth rate

Now, we look at growth rate, which combines the effects of growth and
development. It is important to keep in mind that the sample size is
smaller when considering development and growth simultaneously (some
stages with growth data lack development data and vice versa). There are
different ways of measuring growth rate. Here I use this formula for
relative growth rate (rgr): (ln end size - ln starting size) / time. For
modeling, it does not matter if we use log10 or ln, because the
resulting rgrs are perfectly correlated. Interpretation of rgr is easier
with natural logs, though, as it approximates the exponential growth
rate (i.e. the % change per unit time), but only for lower growth rates.

Adding relative growth rate weakly improves the model.

    ## Data: filter(st_level, !is.na(rg_dd))
    ## Models:
    ## reg00: num_hosts_suspicious_removed ~ 1 + (1 | obs) + (1 | Parasite.species)
    ## reg0: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg0:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg0:     (1 | parasite_class) + (1 | parasite_phylum)
    ## reg1: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg1:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg1:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort
    ## reg2: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg2:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg2:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort + 
    ## reg2:     rg_dd
    ##       Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
    ## reg00  3 2560.1 2572.8 -1277.1   2554.1                              
    ## reg0   8 2557.6 2591.2 -1270.8   2541.6 12.5831      5    0.02762 *  
    ## reg1   9 2502.6 2540.4 -1242.3   2484.6 56.9940      1  4.371e-14 ***
    ## reg2  10 2500.1 2542.1 -1240.0   2480.1  4.5079      1    0.03374 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

The parameter estimate is positive, faster growth rates are associated
with more generalism, which is not what we would expect.

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: poisson  ( log )
    ## Formula: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) +  
    ##     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) +  
    ##     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort +      rg_dd
    ##    Data: filter(st_level, !is.na(rg_dd))
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   2500.1   2542.1  -1240.0   2480.1      484 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.3623 -0.5252 -0.1187  0.3014  0.9716 
    ## 
    ## Random effects:
    ##  Groups           Name        Variance  Std.Dev. 
    ##  obs              (Intercept) 3.943e-01 6.279e-01
    ##  Parasite.species (Intercept) 3.770e-08 1.942e-04
    ##  parasite_genus   (Intercept) 3.563e-02 1.887e-01
    ##  parasite_family  (Intercept) 7.750e-02 2.784e-01
    ##  parasite_order   (Intercept) 5.778e-09 7.601e-05
    ##  parasite_class   (Intercept) 1.919e-08 1.385e-04
    ##  parasite_phylum  (Intercept) 4.212e-09 6.490e-05
    ## Number of obs: 494, groups:  
    ## obs, 494; Parasite.species, 390; parasite_genus, 227; parasite_family, 89; parasite_order, 23; parasite_class, 6; parasite_phylum, 3
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)    1.00843    0.08091  12.464  < 2e-16 ***
    ## zstudy_effort  0.41216    0.05301   7.775 7.53e-15 ***
    ## rg_dd         12.94283    6.05151   2.139   0.0325 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) zstdy_
    ## zstudy_ffrt -0.211       
    ## rg_dd       -0.653 -0.026
    ## convergence code: 0
    ## boundary (singular) fit: see ?isSingular

The variance explained is still less than 1%.

    ## # A tibble: 4 x 5
    ##   step                df_used marg_r2 cond_r2 rand_var_explained
    ##   <chr>                 <dbl>   <dbl>   <dbl>              <dbl>
    ## 1 within-species corr      NA   0       0.194              0.194
    ## 2 taxonomy                  0   0       0.218              0.218
    ## 3 study effort              1   0.126   0.258              0.132
    ## 4 growth rate               1   0.131   0.259              0.128

What about if we include paratenic hosts in the analysis with relative
growth rates? This does not necessarily make sense, because we assumed
quick development times for imputation. This may lead to spuriously
large growth rates. In that direction, the relationship between relative
growth is modestly more significant.

    ## Data: filter(st_level, !is.na(rel_growth_rate_paratenic))
    ## Models:
    ## reg00: num_hosts_suspicious_removed ~ 1 + (1 | obs) + (1 | Parasite.species)
    ## reg0: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg0:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg0:     (1 | parasite_class) + (1 | parasite_phylum)
    ## reg1: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg1:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg1:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort
    ## reg2: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg2:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg2:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort + 
    ## reg2:     rel_growth_rate_paratenic
    ##       Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
    ## reg00  3 3506.0 3519.4 -1750.0   3500.0                              
    ## reg0   8 3490.3 3526.2 -1737.2   3474.3 25.6464      5  0.0001045 ***
    ## reg1   9 3420.7 3461.0 -1701.3   3402.7 71.6780      1  < 2.2e-16 ***
    ## reg2  10 3416.4 3461.3 -1698.2   3396.4  6.2383      1  0.0125018 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

But this term still does not explain much variation.

    ## # A tibble: 4 x 5
    ##   step                df_used marg_r2 cond_r2 rand_var_explained
    ##   <chr>                 <dbl>   <dbl>   <dbl>              <dbl>
    ## 1 within-species corr      NA   0       0.243              0.243
    ## 2 taxonomy                  0   0       0.255              0.255
    ## 3 study effort              1   0.131   0.258              0.127
    ## 4 growth rate               1   0.137   0.271              0.134

### Case 2

We know that growth, development, and specificity vary across life
cycles. Do life cycle variables like host number or host type
(definitive vs intermediate) explain the variation in host generalism
better? In other words, does life history explain variation in
generalism beyond that explained by the life cycle itself? Now, I
examine life history variables *after* building a model that includes
life cycle characteristics.

#### Relative growth

Here are a series of models indicating that (1) study effort affects
host counts, (2) the distinction between intermediate and definitive
host as well as between 1st, 2nd, and 3rd hosts is, and (3) adding
relative growth is not an improvement.

    ## Data: filter(st_level, !is.na(rel_growth_biov))
    ## Models:
    ## reg0: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg0:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg0:     (1 | parasite_class) + (1 | parasite_phylum)
    ## reg1: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg1:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg1:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort
    ## reg4: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg4:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg4:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort + 
    ## reg4:     Def.int + Host_no_fac + Def.int:Host_no_fac
    ## reg5: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg5:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg5:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort + 
    ## reg5:     Def.int + Host_no_fac + rel_growth_biov + Def.int:Host_no_fac
    ##      Df    AIC    BIC  logLik deviance    Chisq Chi Df Pr(>Chisq)    
    ## reg0  8 5451.4 5491.4 -2717.7   5435.4                               
    ## reg1  9 5340.0 5385.0 -2661.0   5322.0 113.4100      1  < 2.2e-16 ***
    ## reg4 16 5302.3 5382.3 -2635.2   5270.3  51.6767      7  6.761e-09 ***
    ## reg5 17 5302.7 5387.7 -2634.3   5268.7   1.6063      1      0.205    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## # A tibble: 7 x 5
    ##   step                df_used marg_r2 cond_r2 rand_var_explained
    ##   <chr>                 <dbl>   <dbl>   <dbl>              <dbl>
    ## 1 within-species corr      NA   0       0.265              0.265
    ## 2 taxonomy                  0   0       0.285              0.285
    ## 3 study effort              1   0.13    0.272              0.142
    ## 4 stage function            1   0.133   0.281              0.148
    ## 5 host number               3   0.166   0.269              0.103
    ## 6 host x stage              3   0.167   0.276              0.109
    ## 7 growth                    1   0.167   0.277              0.11

Check whether this changes after imputing data for paratenic stages. The
effect of parasite growth is actually non-significant\! Stage info must
explain that imputed variation.

    ## Data: filter(st_level, !is.na(rel_growth_paratenic))
    ## Models:
    ## reg00: num_hosts_suspicious_removed ~ 1 + (1 | obs) + (1 | Parasite.species)
    ## reg0: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg0:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg0:     (1 | parasite_class) + (1 | parasite_phylum)
    ## reg1: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg1:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg1:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort
    ## reg4: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg4:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg4:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort + 
    ## reg4:     Def.int + Host_no_fac + Def.int:Host_no_fac
    ## reg5: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg5:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg5:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort + 
    ## reg5:     Def.int + Host_no_fac + rel_growth_paratenic + Def.int:Host_no_fac
    ##       Df    AIC    BIC  logLik deviance    Chisq Chi Df Pr(>Chisq)    
    ## reg00  3 6056.2 6071.5 -3025.1   6050.2                               
    ## reg0   8 6034.4 6075.1 -3009.2   6018.4  31.8608      5  6.330e-06 ***
    ## reg1   9 5919.7 5965.5 -2950.8   5901.7 116.6614      1  < 2.2e-16 ***
    ## reg4  16 5886.8 5968.4 -2927.4   5854.8  46.8610      7  5.942e-08 ***
    ## reg5  17 5888.0 5974.6 -2927.0   5854.0   0.8415      1      0.359    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#### Development time - days

Adding development time in days to a model with life cycle
characteristics is not an improvement.

    ## Data: filter(st_level, !is.na(avg_dt))
    ## Models:
    ## reg00: num_hosts_suspicious_removed ~ 1 + (1 | obs) + (1 | Parasite.species)
    ## reg0: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg0:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg0:     (1 | parasite_class) + (1 | parasite_phylum)
    ## reg1: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg1:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg1:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort
    ## reg4: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg4:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg4:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort + 
    ## reg4:     Def.int + Host_no_fac + Def.int:Host_no_fac
    ## reg5: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg5:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg5:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort + 
    ## reg5:     Def.int + Host_no_fac + log10(avg_dt) + Def.int:Host_no_fac
    ##       Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
    ## reg00  3 4119.0 4133.1 -2056.5   4113.0                              
    ## reg0   8 4110.5 4148.2 -2047.2   4094.5 18.4970      5   0.002384 ** 
    ## reg1   9 4046.9 4089.3 -2014.5   4028.9 65.5607      1  5.635e-16 ***
    ## reg4  15 4013.6 4084.2 -1991.8   3983.6 45.3401      6  4.006e-08 ***
    ## reg5  16 4015.2 4090.5 -1991.6   3983.2  0.4386      1   0.507775    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## # A tibble: 7 x 5
    ##   step                df_used marg_r2 cond_r2 rand_var_explained
    ##   <chr>                 <dbl>   <dbl>   <dbl>              <dbl>
    ## 1 within-species corr      NA   0       0.143             0.143 
    ## 2 taxonomy                  0   0       0.17              0.17  
    ## 3 study effort              1   0.092   0.232             0.14  
    ## 4 stage function            1   0.095   0.23              0.135 
    ## 5 host number               3   0.138   0.222             0.0840
    ## 6 host x stage              2   0.14    0.227             0.087 
    ## 7 devo                      1   0.141   0.225             0.084

When we assume development does not occur in paratenic hosts, the
results are similar. Developmental time does not explain additional
variation in host range beyond that explained by the life cycle.

    ## Data: filter(st_level, !is.na(avg_dt_paratenic))
    ## Models:
    ## reg0: num_hosts_suspicious_removed ~ 1 + (1 | obs) + (1 | Parasite.species) + 
    ## reg0:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg0:     (1 | parasite_class) + (1 | parasite_phylum)
    ## reg1: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg1:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg1:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort
    ## reg4: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg4:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg4:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort + 
    ## reg4:     Def.int + Host_no_fac + Def.int:Host_no_fac
    ## reg5: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg5:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg5:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort + 
    ## reg5:     Def.int + Host_no_fac + log10(avg_dt_paratenic) + Def.int:Host_no_fac
    ##      Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
    ## reg0  8 5040.0 5079.1 -2512.0   5024.0                              
    ## reg1  9 4960.9 5004.9 -2471.4   4942.9 81.0842      1  < 2.2e-16 ***
    ## reg4 16 4922.6 5000.7 -2445.3   4890.6 52.3375      7  5.009e-09 ***
    ## reg5 17 4924.5 5007.6 -2445.2   4890.5  0.0583      1     0.8092    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#### Development time - degree days

On the other hand, adding development time in degree days to a model
with life cycle characteristics is a marginal improvement.

    ## Data: filter(st_level, !is.na(avg_dd))
    ## Models:
    ## reg00: num_hosts_suspicious_removed ~ 1 + (1 | obs) + (1 | Parasite.species)
    ## reg0: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg0:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg0:     (1 | parasite_class) + (1 | parasite_phylum)
    ## reg1: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg1:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg1:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort
    ## reg4: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg4:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg4:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort + 
    ## reg4:     Def.int + Host_no_fac + Def.int:Host_no_fac
    ## reg5: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg5:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg5:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort + 
    ## reg5:     Def.int + Host_no_fac + log10(avg_dd) + Def.int:Host_no_fac
    ##       Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
    ## reg00  3 3638.3 3652.0 -1816.2   3632.3                              
    ## reg0   8 3638.6 3675.3 -1811.3   3622.6  9.6578      5    0.08553 .  
    ## reg1   9 3583.1 3624.3 -1782.5   3565.1 57.5838      1  3.239e-14 ***
    ## reg4  15 3551.4 3620.1 -1760.7   3521.4 43.6673      6  8.604e-08 ***
    ## reg5  16 3549.8 3623.1 -1758.9   3517.8  3.5699      1    0.05884 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

But it explained little variation in host range.

    ## # A tibble: 7 x 5
    ##   step                df_used marg_r2 cond_r2 rand_var_explained
    ##   <chr>                 <dbl>   <dbl>   <dbl>              <dbl>
    ## 1 within-species corr      NA   0       0.132             0.132 
    ## 2 taxonomy                  0   0       0.143             0.143 
    ## 3 study effort              1   0.093   0.211             0.118 
    ## 4 stage function            1   0.096   0.202             0.106 
    ## 5 host number               3   0.147   0.192             0.045 
    ## 6 host x stage              2   0.149   0.197             0.048 
    ## 7 devo                      1   0.15    0.176             0.0260

When we assume development does not occur in paratenic hosts, the
results change slightly. Specifically, life cycle characteristics become
more important, i.e. being a second or third intermediate hosts, but
developmental time still does not explain additional variation in host
range beyond that explained by the life cycle.

    ## Data: filter(st_level, !is.na(avg_dd_paratenic))
    ## Models:
    ## reg0: num_hosts_suspicious_removed ~ 1 + (1 | obs) + (1 | Parasite.species) + 
    ## reg0:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg0:     (1 | parasite_class) + (1 | parasite_phylum)
    ## reg1: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg1:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg1:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort
    ## reg4: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg4:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg4:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort + 
    ## reg4:     Def.int + Host_no_fac + Def.int:Host_no_fac
    ## reg5: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg5:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg5:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort + 
    ## reg5:     Def.int + Host_no_fac + log10(avg_dd_paratenic) + Def.int:Host_no_fac
    ##      Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
    ## reg0  8 4570.2 4608.5 -2277.1   4554.2                              
    ## reg1  9 4499.1 4542.1 -2240.5   4481.1 73.1364      1  < 2.2e-16 ***
    ## reg4 16 4461.6 4538.1 -2214.8   4429.6 51.4972      7  7.334e-09 ***
    ## reg5 17 4463.0 4544.3 -2214.5   4429.0  0.5869      1     0.4436    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#### Relative growth rate

Again, the relation between relative growth rate and host range is
marginal, but not in the expected direction.

    ## Data: filter(st_level, !is.na(rg_dd))
    ## Models:
    ## reg00: num_hosts_suspicious_removed ~ 1 + (1 | obs) + (1 | Parasite.species)
    ## reg0: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg0:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg0:     (1 | parasite_class) + (1 | parasite_phylum)
    ## reg1: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg1:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg1:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort
    ## reg4: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg4:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg4:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort + 
    ## reg4:     Def.int + Host_no_fac + Def.int:Host_no_fac
    ## reg5: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg5:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg5:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort + 
    ## reg5:     Def.int + Host_no_fac + rg_dd + Def.int:Host_no_fac
    ##       Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
    ## reg00  3 2560.1 2572.8 -1277.1   2554.1                              
    ## reg0   8 2557.6 2591.2 -1270.8   2541.6 12.5831      5    0.02762 *  
    ## reg1   9 2502.6 2540.4 -1242.3   2484.6 56.9940      1  4.371e-14 ***
    ## reg4  14 2486.6 2545.4 -1229.3   2458.6 25.9770      5  9.016e-05 ***
    ## reg5  15 2483.2 2546.2 -1226.6   2453.2  5.4332      1    0.01976 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

The effect of relative growth rate is weak, explaining \<1% of the
variation

    ## # A tibble: 7 x 5
    ##   step                df_used marg_r2 cond_r2 rand_var_explained
    ##   <chr>                 <dbl>   <dbl>   <dbl>              <dbl>
    ## 1 within-species corr      NA   0       0.194             0.194 
    ## 2 taxonomy                  0   0       0.218             0.218 
    ## 3 study effort              1   0.126   0.258             0.132 
    ## 4 stage function            1   0.128   0.261             0.133 
    ## 5 host number               3   0.17    0.236             0.0660
    ## 6 host x stage              1   0.17    0.241             0.0710
    ## 7 rel growth rate           1   0.178   0.245             0.067

The marginally significant effect of growth rate remains after imputing
data for paratenic stages, though it is debatable whether this makes
sense, since imputed growth rates are usually zero.

    ## Data: filter(st_level, !is.na(rel_growth_rate_paratenic))
    ## Models:
    ## reg00: num_hosts_suspicious_removed ~ 1 + (1 | obs) + (1 | Parasite.species)
    ## reg0: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg0:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg0:     (1 | parasite_class) + (1 | parasite_phylum)
    ## reg1: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg1:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg1:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort
    ## reg4: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg4:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg4:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort + 
    ## reg4:     Def.int + Host_no_fac + Def.int:Host_no_fac
    ## reg5: num_hosts_suspicious_removed ~ (1 | obs) + (1 | Parasite.species) + 
    ## reg5:     (1 | parasite_genus) + (1 | parasite_family) + (1 | parasite_order) + 
    ## reg5:     (1 | parasite_class) + (1 | parasite_phylum) + zstudy_effort + 
    ## reg5:     Def.int + Host_no_fac + rel_growth_rate_paratenic + Def.int:Host_no_fac
    ##       Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
    ## reg00  3 3506.0 3519.4 -1750.0   3500.0                              
    ## reg0   8 3490.3 3526.2 -1737.2   3474.3 25.6464      5  0.0001045 ***
    ## reg1   9 3420.7 3461.0 -1701.3   3402.7 71.6780      1  < 2.2e-16 ***
    ## reg4  16 3400.6 3472.4 -1684.3   3368.6 34.0686      7  1.672e-05 ***
    ## reg5  17 3398.6 3474.9 -1682.3   3364.6  3.9771      1  0.0461220 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Conclusions

Parasite life cycle, life history, and patterns of host generalism are
related. In particular, parasite exhibit high generalism in second or
third intermediate hosts. These tend to be paratenic hosts, in which
parasites undergo little growth and development. But even if these hosts
are excluded, there is still a negative relationship between generalism
and worm development. Worms grow less in stages where they infect a
broader range of hosts, which is consistent with costs of generalism.

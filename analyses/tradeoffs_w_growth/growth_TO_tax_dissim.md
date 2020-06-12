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

I calculated relative and absolute growth for parasite stages, as well
as the time needed to complete development at a given stage, in days and
degree days.

There is less development time data than growth data. And growth was
calculated with two different size variables: length and biovolume.
Biovolume had more missing data, because biovolume was calculated with
lengths and widths and there were more missing widths than lengths.

    ##          avg_dt          avg_dd  abs_growth_len abs_growth_biov 
    ##            1200            1296             644             914

Therefore, for any paratenic host stage without growth data, I set the
growth to zero and then replotted the data. This makes the negative
relationship clearer, especially for worms in intermediate hosts (since
those are the paratenic hosts).

# Models

I want to test whether life history variables can explain variation in
generalism, which would be evidence for a tradeoff. There are several
groups of variables to test in the model: (1) uninteresting confounders
(phylogeny, study effort), (2) life cycle variables (def vs int, host
number), and (3) life history variables (growth and development). I
first add the uninteresting confounders to the model, simply to control
for them and to be sure that variables added afterwards explain
*additional* variation. Here, I’m mainly interested in the life history
variables. I could add them to the model either before or after the life
cycle variables. If I add them after, I’m testing whether life history
explains variation beyond that explained by the life cycle. If I add
them before, I’m testing for life history tradeoffs, but their effects
are possibly dependent/confounded with life cycle characteristics. I’ll
take both approaches below.

As for controlling for phylogeny, I added taxonomic categories into a
linear mixed model. I did not include the highest and lower taxonomic
levels (genus cuts the data quite thin while phylum is too coarse) as
this improved model convergence. I still need to look at properly
modelling phylogenetic effects

# Taxonomic dissimilarity

Moving onto the next generalism variable: taxonomic dissimilarity.

### Case 1

#### Relative growth

As above, adding relative growth to a model with just taxonomy and study
effort is an improvement. It has about the same magnitude effect as
study effort.

    ## Data: filter(st_level, !is.na(rel_growth_biov))
    ## Models:
    ## reg00: hsi_lcdb_suspcious_rem ~ 1 + (1 | Parasite.species)
    ## reg0: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg0:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg0:     (1 | parasite_phylum)
    ## reg1: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg1:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg1:     (1 | parasite_phylum) + zstudy_effort
    ## reg2: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg2:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg2:     (1 | parasite_phylum) + zstudy_effort + rel_growth_biov
    ##       Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
    ## reg00  3 3213.2 3228.2 -1603.6   3207.2                             
    ## reg0   8 3179.0 3219.0 -1581.5   3163.0 44.191      5  2.118e-08 ***
    ## reg1   9 3120.0 3165.0 -1551.0   3102.0 61.051      1  5.561e-15 ***
    ## reg2  10 3069.4 3119.5 -1524.7   3049.4 52.531      1  4.236e-13 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Here are the parameters:

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: 
    ## hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) +  
    ##     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) +  
    ##     (1 | parasite_phylum) + zstudy_effort + rel_growth_biov
    ##    Data: filter(st_level, !is.na(rel_growth_biov))
    ## 
    ## REML criterion at convergence: 3062.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.0012 -0.7934 -0.1222  0.6932  3.5410 
    ## 
    ## Random effects:
    ##  Groups           Name        Variance  Std.Dev.
    ##  Parasite.species (Intercept) 0.000e+00 0.000000
    ##  parasite_genus   (Intercept) 3.142e-02 0.177261
    ##  parasite_family  (Intercept) 5.694e-02 0.238621
    ##  parasite_order   (Intercept) 6.315e-02 0.251288
    ##  parasite_class   (Intercept) 8.959e-03 0.094653
    ##  parasite_phylum  (Intercept) 1.020e-08 0.000101
    ##  Residual                     8.657e-01 0.930444
    ## Number of obs: 1099, groups:  
    ## Parasite.species, 653; parasite_genus, 326; parasite_family, 110; parasite_order, 31; parasite_class, 6; parasite_phylum, 3
    ## 
    ## Fixed effects:
    ##                 Estimate Std. Error t value
    ## (Intercept)      2.24159    0.10257  21.855
    ## zstudy_effort    0.40332    0.04638   8.695
    ## rel_growth_biov -0.14838    0.02015  -7.364
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) zstdy_
    ## zstudy_ffrt -0.027       
    ## rl_grwth_bv -0.443 -0.073
    ## convergence code: 0
    ## boundary (singular) fit: see ?isSingular

    ## # A tibble: 4 x 6
    ##   step            df_used marg_r2 cond_r2 rand_var_explained species_var_explai~
    ##   <chr>             <dbl>   <dbl>   <dbl>              <dbl>               <dbl>
    ## 1 within-species~      NA   0       0.147              0.147               0.147
    ## 2 taxonomy              0   0       0.189              0.189               0.032
    ## 3 study effort          1   0.065   0.207              0.142               0    
    ## 4 growth                1   0.105   0.245              0.14                0

The same pattern is seen if we assume worms do not grow in paratenic
hosts, in fact the effect size is quite a bit larger.

    ## Data: filter(st_level, !is.na(rel_growth_paratenic))
    ## Models:
    ## reg00: hsi_lcdb_suspcious_rem ~ 1 + (1 | Parasite.species)
    ## reg0: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg0:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg0:     (1 | parasite_phylum)
    ## reg1: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg1:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg1:     (1 | parasite_phylum) + zstudy_effort
    ## reg2: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg2:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg2:     (1 | parasite_phylum) + zstudy_effort + rel_growth_paratenic
    ##       Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
    ## reg00  3 3682.7 3698.0 -1838.3   3676.7                             
    ## reg0   8 3611.6 3652.4 -1797.8   3595.6 81.086      5  4.973e-16 ***
    ## reg1   9 3560.5 3606.3 -1771.2   3542.5 53.137      1  3.111e-13 ***
    ## reg2  10 3482.0 3532.9 -1731.0   3462.0 80.482      1  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## # A tibble: 4 x 6
    ##   step            df_used marg_r2 cond_r2 rand_var_explained species_var_explai~
    ##   <chr>             <dbl>   <dbl>   <dbl>              <dbl>               <dbl>
    ## 1 within-species~      NA   0       0.204              0.204               0.204
    ## 2 taxonomy              0   0       0.22               0.22                0.018
    ## 3 study effort          1   0.053   0.233              0.18                0    
    ## 4 devo                  1   0.115   0.267              0.152               0

#### Developmental time - days

Now let’s fit the same models, but using the simplest metric for
developmental time, days. The addition of devo time, without imputation
in paratenic hosts is marginally significant.

    ## Data: filter(st_level, !is.na(avg_dt))
    ## Models:
    ## reg00: hsi_lcdb_suspcious_rem ~ 1 + (1 | Parasite.species)
    ## reg0: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg0:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg0:     (1 | parasite_phylum)
    ## reg1: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg1:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg1:     (1 | parasite_phylum) + zstudy_effort
    ## reg2: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg2:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg2:     (1 | parasite_phylum) + zstudy_effort + log10(avg_dt)
    ##       Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
    ## reg00  3 2272.7 2286.8 -1133.3   2266.7                             
    ## reg0   8 2262.9 2300.5 -1123.4   2246.9 19.810      5   0.001357 ** 
    ## reg1   9 2231.3 2273.7 -1106.7   2213.3 33.540      1  6.981e-09 ***
    ## reg2  10 2230.0 2277.1 -1105.0   2210.0  3.342      1   0.067531 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Here are the parameters:

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: 
    ## hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) +  
    ##     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) +  
    ##     (1 | parasite_phylum) + zstudy_effort + log10(avg_dt)
    ##    Data: filter(st_level, !is.na(avg_dt))
    ## 
    ## REML criterion at convergence: 2220.7
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.7630 -0.8461 -0.0738  0.7516  2.9351 
    ## 
    ## Random effects:
    ##  Groups           Name        Variance  Std.Dev. 
    ##  Parasite.species (Intercept) 2.598e-10 1.612e-05
    ##  parasite_genus   (Intercept) 9.974e-02 3.158e-01
    ##  parasite_family  (Intercept) 6.048e-03 7.777e-02
    ##  parasite_order   (Intercept) 4.142e-02 2.035e-01
    ##  parasite_class   (Intercept) 0.000e+00 0.000e+00
    ##  parasite_phylum  (Intercept) 1.383e-10 1.176e-05
    ##  Residual                     7.787e-01 8.825e-01
    ## Number of obs: 817, groups:  
    ## Parasite.species, 609; parasite_genus, 293; parasite_family, 100; parasite_order, 26; parasite_class, 6; parasite_phylum, 3
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error t value
    ## (Intercept)    2.21856    0.15993  13.872
    ## zstudy_effort  0.27460    0.04506   6.094
    ## log10(avg_dt) -0.18201    0.09976  -1.825
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) zstdy_
    ## zstudy_ffrt  0.066       
    ## lg10(vg_dt) -0.892 -0.144
    ## convergence code: 0
    ## boundary (singular) fit: see ?isSingular

And the R<sup>2</sup>.

    ## # A tibble: 4 x 6
    ##   step            df_used marg_r2 cond_r2 rand_var_explained species_var_explai~
    ##   <chr>             <dbl>   <dbl>   <dbl>              <dbl>               <dbl>
    ## 1 within-species~      NA   0       0.099              0.099               0.099
    ## 2 taxonomy              0   0       0.17               0.17                0    
    ## 3 study effort          1   0.048   0.204              0.156               0    
    ## 4 devo                  1   0.052   0.202              0.15                0

The negative relationship between generalism and devo time is much
clearer when we assume minimal development in paratenic hosts.

    ## Data: filter(st_level, !is.na(avg_dt_paratenic))
    ## Models:
    ## reg00: hsi_lcdb_suspcious_rem ~ 1 + (1 | Parasite.species)
    ## reg0: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg0:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg0:     (1 | parasite_phylum)
    ## reg1: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg1:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg1:     (1 | parasite_phylum) + zstudy_effort
    ## reg2: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg2:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg2:     (1 | parasite_phylum) + zstudy_effort + log10(avg_dt_paratenic)
    ##       Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
    ## reg00  3 2983.2 2997.9 -1488.6   2977.2                             
    ## reg0   8 2933.5 2972.6 -1458.8   2917.5 59.730      5  1.382e-11 ***
    ## reg1   9 2895.2 2939.2 -1438.6   2877.2 40.264      1  2.218e-10 ***
    ## reg2  10 2812.1 2861.0 -1396.1   2792.1 85.137      1  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## # A tibble: 4 x 6
    ##   step            df_used marg_r2 cond_r2 rand_var_explained species_var_explai~
    ##   <chr>             <dbl>   <dbl>   <dbl>              <dbl>               <dbl>
    ## 1 within-species~      NA   0       0.188              0.188               0.188
    ## 2 taxonomy              0   0       0.209              0.209               0.005
    ## 3 study effort          1   0.048   0.243              0.195               0    
    ## 4 devo                  1   0.134   0.272              0.138               0.005

#### Development time - degree days

Now we look at developmental time in temperature-corrected degree days.
The taxonomic diversity of hosts is negatively related to developmental
times (quick development - high generalism; long development - low
generalism).

    ## Data: filter(st_level, !is.na(avg_dd))
    ## Models:
    ## reg00: hsi_lcdb_suspcious_rem ~ 1 + (1 | Parasite.species)
    ## reg0: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg0:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg0:     (1 | parasite_phylum)
    ## reg1: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg1:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg1:     (1 | parasite_phylum) + zstudy_effort
    ## reg2: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg2:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg2:     (1 | parasite_phylum) + zstudy_effort + log10(avg_dd)
    ##       Df    AIC    BIC   logLik deviance   Chisq Chi Df Pr(>Chisq)    
    ## reg00  3 2009.7 2023.5 -1001.86   2003.7                              
    ## reg0   8 2007.4 2044.1  -995.71   1991.4 12.2983      5   0.030922 *  
    ## reg1   9 1978.7 2019.9  -980.36   1960.7 30.6978      1  3.015e-08 ***
    ## reg2  10 1972.0 2017.8  -976.02   1952.0  8.6881      1   0.003203 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Here are the parameters:

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: 
    ## hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) +  
    ##     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) +  
    ##     (1 | parasite_phylum) + zstudy_effort + log10(avg_dd)
    ##    Data: filter(st_level, !is.na(avg_dd))
    ## 
    ## REML criterion at convergence: 1962.8
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.77848 -0.84764 -0.06272  0.77666  2.90439 
    ## 
    ## Random effects:
    ##  Groups           Name        Variance  Std.Dev. 
    ##  Parasite.species (Intercept) 3.334e-09 5.774e-05
    ##  parasite_genus   (Intercept) 1.058e-01 3.253e-01
    ##  parasite_family  (Intercept) 1.126e-02 1.061e-01
    ##  parasite_order   (Intercept) 2.955e-02 1.719e-01
    ##  parasite_class   (Intercept) 0.000e+00 0.000e+00
    ##  parasite_phylum  (Intercept) 0.000e+00 0.000e+00
    ##  Residual                     7.751e-01 8.804e-01
    ## Number of obs: 721, groups:  
    ## Parasite.species, 550; parasite_genus, 273; parasite_family, 94; parasite_order, 23; parasite_class, 6; parasite_phylum, 3
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error t value
    ## (Intercept)    2.79759    0.27531  10.162
    ## zstudy_effort  0.28777    0.04739   6.072
    ## log10(avg_dd) -0.28942    0.09700  -2.984
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) zstdy_
    ## zstudy_ffrt  0.134       
    ## lg10(vg_dd) -0.967 -0.187
    ## convergence code: 0
    ## boundary (singular) fit: see ?isSingular

And the R<sup>2</sup>.

    ## # A tibble: 4 x 6
    ##   step            df_used marg_r2 cond_r2 rand_var_explained species_var_explai~
    ##   <chr>             <dbl>   <dbl>   <dbl>              <dbl>               <dbl>
    ## 1 within-species~      NA   0       0.071              0.071               0.071
    ## 2 taxonomy              0   0       0.16               0.16                0    
    ## 3 study effort          1   0.05    0.191              0.141               0    
    ## 4 devo                  1   0.058   0.208              0.15                0

This is strongly exacerbated by assuming minimal development in
paratenic hosts. Short development and generalism are associated.

    ## Data: filter(st_level, !is.na(avg_dd_paratenic))
    ## Models:
    ## reg00: hsi_lcdb_suspcious_rem ~ 1 + (1 | Parasite.species)
    ## reg0: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg0:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg0:     (1 | parasite_phylum)
    ## reg1: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg1:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg1:     (1 | parasite_phylum) + zstudy_effort
    ## reg2: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg2:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg2:     (1 | parasite_phylum) + zstudy_effort + log10(avg_dd_paratenic)
    ##       Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
    ## reg00  3 2712.4 2726.7 -1353.2   2706.4                             
    ## reg0   8 2667.5 2705.7 -1325.7   2651.5 54.907      5  1.364e-10 ***
    ## reg1   9 2635.2 2678.2 -1308.6   2617.2 34.311      1  4.698e-09 ***
    ## reg2  10 2545.3 2593.1 -1262.7   2525.3 91.852      1  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## # A tibble: 4 x 6
    ##   step            df_used marg_r2 cond_r2 rand_var_explained species_var_explai~
    ##   <chr>             <dbl>   <dbl>   <dbl>              <dbl>               <dbl>
    ## 1 within-species~      NA   0       0.19               0.19                 0.19
    ## 2 taxonomy              0   0       0.215              0.215                0   
    ## 3 study effort          1   0.045   0.255              0.21                 0   
    ## 4 devo                  1   0.145   0.293              0.148                0

#### Relative growth rate

There are different ways of measuring growth rate. Here I use this
formula for relative growth rate (rgr): (ln end size - ln starting size)
/ time. For modeling, it does not matter if we use log10 or ln, because
the resulting rgrs are perfectly correlated. Interpretation of rgr is
easier with natural logs, though, as it approximates the exponential
growth rate (i.e. the % change per unit time), but only for lower growth
rates.

Growth rate was unrelated to taxonomic dissimilarity, explaining no
variation.

    ## Data: filter(st_level, !is.na(rg_dd))
    ## Models:
    ## reg00: hsi_lcdb_suspcious_rem ~ 1 + (1 | Parasite.species)
    ## reg0: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg0:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg0:     (1 | parasite_phylum)
    ## reg1: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg1:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg1:     (1 | parasite_phylum) + zstudy_effort
    ## reg2: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg2:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg2:     (1 | parasite_phylum) + zstudy_effort + rg_dd
    ##       Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
    ## reg00  3 1384.5 1397.1 -689.23   1378.5                              
    ## reg0   8 1388.9 1422.6 -686.47   1372.9  5.5118      5     0.3567    
    ## reg1   9 1367.9 1405.8 -674.97   1349.9 23.0064      1  1.615e-06 ***
    ## reg2  10 1369.7 1411.7 -674.84   1349.7  0.2587      1     0.6110    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## # A tibble: 4 x 6
    ##   step            df_used marg_r2 cond_r2 rand_var_explained species_var_explai~
    ##   <chr>             <dbl>   <dbl>   <dbl>              <dbl>               <dbl>
    ## 1 within-species~      NA   0       0.089              0.089               0.089
    ## 2 taxonomy              0   0       0.117              0.117               0    
    ## 3 study effort          1   0.052   0.163              0.111               0    
    ## 4 rel growth rate       1   0.052   0.165              0.113               0

The same pattern is seen if we assume worms do not grow in paratenic
hosts. Growth rate does not differ for generalist vs specialist stages.

    ## Data: filter(st_level, !is.na(rel_growth_rate_paratenic))
    ## Models:
    ## reg00: hsi_lcdb_suspcious_rem ~ 1 + (1 | Parasite.species)
    ## reg0: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg0:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg0:     (1 | parasite_phylum)
    ## reg1: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg1:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg1:     (1 | parasite_phylum) + zstudy_effort
    ## reg2: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg2:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg2:     (1 | parasite_phylum) + zstudy_effort + rel_growth_rate_paratenic
    ##       Df    AIC    BIC   logLik deviance   Chisq Chi Df Pr(>Chisq)    
    ## reg00  3 2064.8 2078.3 -1029.40   2058.8                              
    ## reg0   8 2029.3 2065.2 -1006.66   2013.3 45.4839      5  1.157e-08 ***
    ## reg1   9 2006.8 2047.1  -994.38   1988.8 24.5698      1  7.166e-07 ***
    ## reg2  10 2007.4 2052.2  -993.69   1987.4  1.3825      1     0.2397    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## # A tibble: 4 x 6
    ##   step            df_used marg_r2 cond_r2 rand_var_explained species_var_explai~
    ##   <chr>             <dbl>   <dbl>   <dbl>              <dbl>               <dbl>
    ## 1 within-species~      NA   0       0.233              0.233               0.233
    ## 2 taxonomy              0   0       0.236              0.236               0    
    ## 3 study effort          1   0.043   0.268              0.225               0    
    ## 4 devo                  1   0.044   0.274              0.23                0

### Case 2

#### Relative growth

Is the variation in taxonomic dissimilarity caused by life history
variables simply due to life cycle characteristics?

In the series of models below, the distinction between intermediate vs
definitive hosts is mildly significant and the host number matters
(highest in second and third hosts). Adding relative growth is a clear
improvement.

    ## Data: filter(st_level, !is.na(rel_growth_biov))
    ## Models:
    ## reg0: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg0:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg0:     (1 | parasite_phylum)
    ## reg1: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg1:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg1:     (1 | parasite_phylum) + zstudy_effort
    ## reg4: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg4:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg4:     (1 | parasite_phylum) + zstudy_effort + Def.int + Host_no_fac + 
    ## reg4:     Def.int:Host_no_fac
    ## reg5: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg5:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg5:     (1 | parasite_phylum) + zstudy_effort + Def.int + Host_no_fac + 
    ## reg5:     rel_growth_biov + Def.int:Host_no_fac
    ##      Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
    ## reg0  8 3179.0 3219.0 -1581.5   3163.0                              
    ## reg1  9 3120.0 3165.0 -1551.0   3102.0  61.051      1  5.561e-15 ***
    ## reg4 16 3029.1 3109.1 -1498.5   2997.1 104.869      7  < 2.2e-16 ***
    ## reg5 17 3002.8 3087.8 -1484.4   2968.8  28.334      1  1.021e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## # A tibble: 7 x 6
    ##   step            df_used marg_r2 cond_r2 rand_var_explained species_var_explai~
    ##   <chr>             <dbl>   <dbl>   <dbl>              <dbl>               <dbl>
    ## 1 within-species~      NA   0       0.147              0.147               0.147
    ## 2 taxonomy              0   0       0.189              0.189               0.032
    ## 3 study effort          1   0.065   0.207              0.142               0    
    ## 4 stage function        1   0.081   0.215              0.134               0    
    ## 5 host number           3   0.138   0.292              0.154               0    
    ## 6 host x stage          3   0.145   0.288              0.143               0    
    ## 7 growth                1   0.164   0.317              0.153               0

The same results are seen under the assumption of limited growth in
paratenic hosts.

    ## Data: filter(st_level, !is.na(rel_growth_paratenic))
    ## Models:
    ## reg00: hsi_lcdb_suspcious_rem ~ 1 + (1 | Parasite.species)
    ## reg0: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg0:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg0:     (1 | parasite_phylum)
    ## reg1: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg1:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg1:     (1 | parasite_phylum) + zstudy_effort
    ## reg2: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg2:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg2:     (1 | parasite_phylum) + zstudy_effort + Def.int
    ## reg3: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg3:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg3:     (1 | parasite_phylum) + zstudy_effort + Def.int + Host_no_fac
    ## reg4: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg4:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg4:     (1 | parasite_phylum) + zstudy_effort + Def.int + Host_no_fac + 
    ## reg4:     Def.int:Host_no_fac
    ## reg5: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg5:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg5:     (1 | parasite_phylum) + zstudy_effort + Def.int + Host_no_fac + 
    ## reg5:     rel_growth_paratenic + Def.int:Host_no_fac
    ##       Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
    ## reg00  3 3682.7 3698.0 -1838.3   3676.7                              
    ## reg0   8 3611.6 3652.4 -1797.8   3595.6 81.0858      5  4.973e-16 ***
    ## reg1   9 3560.5 3606.3 -1771.2   3542.5 53.1367      1  3.111e-13 ***
    ## reg2  10 3525.6 3576.5 -1752.8   3505.6 36.8913      1  1.249e-09 ***
    ## reg3  13 3434.1 3500.3 -1704.0   3408.1 97.4944      3  < 2.2e-16 ***
    ## reg4  16 3431.2 3512.8 -1699.6   3399.2  8.8423      3    0.03146 *  
    ## reg5  17 3407.7 3494.3 -1686.8   3373.7 25.5331      1  4.349e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## # A tibble: 7 x 6
    ##   step            df_used marg_r2 cond_r2 rand_var_explained species_var_explai~
    ##   <chr>             <dbl>   <dbl>   <dbl>              <dbl>               <dbl>
    ## 1 within-species~      NA   0       0.204              0.204               0.204
    ## 2 taxonomy              0   0       0.22               0.22                0.018
    ## 3 study effort          1   0.053   0.233              0.18                0    
    ## 4 stage function        1   0.08    0.235              0.155               0    
    ## 5 host number           3   0.16    0.308              0.148               0.005
    ## 6 host x stage          3   0.166   0.299              0.133               0.005
    ## 7 growth                1   0.183   0.332              0.149               0.024

#### Development time - days

Parasite stages with a taxonomically diverse set of hosts have slightly
shorter development times, even after accounting for life cycle
characteristics.

    ## Data: filter(st_level, !is.na(avg_dt))
    ## Models:
    ## reg00: hsi_lcdb_suspcious_rem ~ 1 + (1 | Parasite.species)
    ## reg0: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg0:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg0:     (1 | parasite_phylum)
    ## reg1: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg1:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg1:     (1 | parasite_phylum) + zstudy_effort
    ## reg4: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg4:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg4:     (1 | parasite_phylum) + zstudy_effort + Def.int + Host_no_fac + 
    ## reg4:     Def.int:Host_no_fac
    ## reg5: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg5:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg5:     (1 | parasite_phylum) + zstudy_effort + Def.int + Host_no_fac + 
    ## reg5:     log10(avg_dt) + Def.int:Host_no_fac
    ##       Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
    ## reg00  3 2272.7 2286.8 -1133.3   2266.7                              
    ## reg0   8 2262.9 2300.5 -1123.4   2246.9 19.8100      5  0.0013566 ** 
    ## reg1   9 2231.3 2273.7 -1106.7   2213.3 33.5400      1  6.981e-09 ***
    ## reg4  15 2217.5 2288.1 -1093.8   2187.5 25.8488      6  0.0002376 ***
    ## reg5  16 2215.4 2290.7 -1091.7   2183.4  4.0878      1  0.0431936 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## # A tibble: 7 x 6
    ##   step            df_used marg_r2 cond_r2 rand_var_explained species_var_explai~
    ##   <chr>             <dbl>   <dbl>   <dbl>              <dbl>               <dbl>
    ## 1 within-species~      NA   0       0.099              0.099               0.099
    ## 2 taxonomy              0   0       0.17               0.17                0    
    ## 3 study effort          1   0.048   0.204              0.156               0    
    ## 4 stage function        1   0.053   0.216              0.163               0    
    ## 5 host number           3   0.077   0.221              0.144               0    
    ## 6 host x stage          2   0.08    0.225              0.145               0    
    ## 7 devo                  1   0.085   0.224              0.139               0

When we assume development does not occur in paratenic hosts, the
results appear stronger.

    ## Data: filter(st_level, !is.na(avg_dt_paratenic))
    ## Models:
    ## reg00: hsi_lcdb_suspcious_rem ~ 1 + (1 | Parasite.species)
    ## reg0: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg0:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg0:     (1 | parasite_phylum)
    ## reg1: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg1:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg1:     (1 | parasite_phylum) + zstudy_effort
    ## reg2: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg2:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg2:     (1 | parasite_phylum) + zstudy_effort + Def.int
    ## reg3: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg3:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg3:     (1 | parasite_phylum) + zstudy_effort + Def.int + Host_no_fac
    ## reg4: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg4:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg4:     (1 | parasite_phylum) + zstudy_effort + Def.int + Host_no_fac + 
    ## reg4:     Def.int:Host_no_fac
    ## reg5: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg5:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg5:     (1 | parasite_phylum) + zstudy_effort + Def.int + Host_no_fac + 
    ## reg5:     log10(avg_dt_paratenic) + Def.int:Host_no_fac
    ##       Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
    ## reg00  3 2983.2 2997.9 -1488.6   2977.2                              
    ## reg0   8 2933.5 2972.6 -1458.8   2917.5 59.7305      5  1.382e-11 ***
    ## reg1   9 2895.2 2939.2 -1438.6   2877.2 40.2645      1  2.218e-10 ***
    ## reg2  10 2865.5 2914.3 -1422.7   2845.5 31.7623      1  1.742e-08 ***
    ## reg3  13 2808.0 2871.5 -1391.0   2782.0 63.4834      3  1.059e-13 ***
    ## reg4  16 2808.0 2886.2 -1388.0   2776.0  5.9823      3  0.1124746    
    ## reg5  17 2795.8 2878.8 -1380.9   2761.8 14.2653      1  0.0001588 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## # A tibble: 7 x 6
    ##   step            df_used marg_r2 cond_r2 rand_var_explained species_var_explai~
    ##   <chr>             <dbl>   <dbl>   <dbl>              <dbl>               <dbl>
    ## 1 within-species~      NA   0       0.188              0.188               0.188
    ## 2 taxonomy              0   0       0.209              0.209               0.005
    ## 3 study effort          1   0.048   0.243              0.195               0    
    ## 4 stage function        1   0.071   0.26               0.189               0    
    ## 5 host number           3   0.151   0.297              0.146               0    
    ## 6 host x stage          3   0.152   0.286              0.134               0    
    ## 7 growth                1   0.167   0.289              0.122               0

#### Development time - degree days

Now we fit the same models, but with degree days. The results are
consistent. Even after accounting for life cycle characteristics, there
is an effect of developmental time.

    ## Data: filter(st_level, !is.na(avg_dd))
    ## Models:
    ## reg00: hsi_lcdb_suspcious_rem ~ 1 + (1 | Parasite.species)
    ## reg0: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg0:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg0:     (1 | parasite_phylum)
    ## reg1: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg1:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg1:     (1 | parasite_phylum) + zstudy_effort
    ## reg4: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg4:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg4:     (1 | parasite_phylum) + zstudy_effort + Def.int + Host_no_fac + 
    ## reg4:     Def.int:Host_no_fac
    ## reg5: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg5:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg5:     (1 | parasite_phylum) + zstudy_effort + Def.int + Host_no_fac + 
    ## reg5:     log10(avg_dd) + Def.int:Host_no_fac
    ##       Df    AIC    BIC   logLik deviance   Chisq Chi Df Pr(>Chisq)    
    ## reg00  3 2009.7 2023.5 -1001.86   2003.7                              
    ## reg0   8 2007.4 2044.1  -995.71   1991.4 12.2983      5  0.0309216 *  
    ## reg1   9 1978.7 2019.9  -980.36   1960.7 30.6978      1  3.015e-08 ***
    ## reg4  15 1966.5 2035.2  -968.23   1936.5 24.2567      6  0.0004684 ***
    ## reg5  16 1962.1 2035.4  -965.04   1930.1  6.3856      1  0.0115049 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## # A tibble: 7 x 6
    ##   step            df_used marg_r2 cond_r2 rand_var_explained species_var_explai~
    ##   <chr>             <dbl>   <dbl>   <dbl>              <dbl>               <dbl>
    ## 1 within-species~      NA   0       0.071              0.071               0.071
    ## 2 taxonomy              0   0       0.16               0.16                0    
    ## 3 study effort          1   0.05    0.191              0.141               0    
    ## 4 stage function        1   0.055   0.204              0.149               0    
    ## 5 host number           3   0.079   0.21               0.131               0    
    ## 6 host x stage          2   0.084   0.21               0.126               0    
    ## 7 devo                  1   0.091   0.225              0.134               0

When we assume development does not occur in paratenic hosts, the
results appear stronger.

    ## Data: filter(st_level, !is.na(avg_dd_paratenic))
    ## Models:
    ## reg00: hsi_lcdb_suspcious_rem ~ 1 + (1 | Parasite.species)
    ## reg0: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg0:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg0:     (1 | parasite_phylum)
    ## reg1: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg1:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg1:     (1 | parasite_phylum) + zstudy_effort
    ## reg2: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg2:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg2:     (1 | parasite_phylum) + zstudy_effort + Def.int
    ## reg3: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg3:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg3:     (1 | parasite_phylum) + zstudy_effort + Def.int + Host_no_fac
    ## reg4: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg4:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg4:     (1 | parasite_phylum) + zstudy_effort + Def.int + Host_no_fac + 
    ## reg4:     Def.int:Host_no_fac
    ## reg5: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg5:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg5:     (1 | parasite_phylum) + zstudy_effort + Def.int + Host_no_fac + 
    ## reg5:     log10(avg_dd_paratenic) + Def.int:Host_no_fac
    ##       Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
    ## reg00  3 2712.4 2726.7 -1353.2   2706.4                              
    ## reg0   8 2667.5 2705.7 -1325.7   2651.5 54.9068      5  1.364e-10 ***
    ## reg1   9 2635.2 2678.2 -1308.6   2617.2 34.3108      1  4.698e-09 ***
    ## reg2  10 2605.5 2653.4 -1292.8   2585.5 31.6151      1  1.880e-08 ***
    ## reg3  13 2556.2 2618.3 -1265.1   2530.2 55.3715      3  5.721e-12 ***
    ## reg4  16 2556.8 2633.3 -1262.4   2524.8  5.3868      3     0.1456    
    ## reg5  17 2540.9 2622.2 -1253.5   2506.9 17.8351      1  2.409e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## # A tibble: 7 x 6
    ##   step            df_used marg_r2 cond_r2 rand_var_explained species_var_explai~
    ##   <chr>             <dbl>   <dbl>   <dbl>              <dbl>               <dbl>
    ## 1 within-species~      NA   0       0.19               0.19                 0.19
    ## 2 taxonomy              0   0       0.215              0.215                0   
    ## 3 study effort          1   0.045   0.255              0.21                 0   
    ## 4 stage function        1   0.072   0.272              0.2                  0   
    ## 5 host number           3   0.153   0.301              0.148                0   
    ## 6 host x stage          3   0.153   0.288              0.135                0   
    ## 7 growth                1   0.172   0.299              0.127                0

#### Relative growth rate

Taxonomic dissimilarity is still not related to growth rate after
accounting for other life cycle characteristics.

    ## Data: filter(st_level, !is.na(rg_dd))
    ## Models:
    ## reg00: hsi_lcdb_suspcious_rem ~ 1 + (1 | Parasite.species)
    ## reg0: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg0:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg0:     (1 | parasite_phylum)
    ## reg1: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg1:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg1:     (1 | parasite_phylum) + zstudy_effort
    ## reg4: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg4:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg4:     (1 | parasite_phylum) + zstudy_effort + Def.int + Host_no_fac + 
    ## reg4:     Def.int:Host_no_fac
    ## reg5: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg5:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg5:     (1 | parasite_phylum) + zstudy_effort + Def.int + Host_no_fac + 
    ## reg5:     rg_dd + Def.int:Host_no_fac
    ##       Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
    ## reg00  3 1384.5 1397.1 -689.23   1378.5                              
    ## reg0   8 1388.9 1422.6 -686.47   1372.9  5.5118      5     0.3567    
    ## reg1   9 1367.9 1405.8 -674.97   1349.9 23.0064      1  1.615e-06 ***
    ## reg4  14 1373.5 1432.3 -672.76   1345.5  4.4235      5     0.4902    
    ## reg5  15 1375.1 1438.2 -672.56   1345.1  0.3897      1     0.5324    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## # A tibble: 7 x 6
    ##   step            df_used marg_r2 cond_r2 rand_var_explained species_var_explai~
    ##   <chr>             <dbl>   <dbl>   <dbl>              <dbl>               <dbl>
    ## 1 within-species~      NA   0       0.089              0.089               0.089
    ## 2 taxonomy              0   0       0.117              0.117               0    
    ## 3 study effort          1   0.052   0.163              0.111               0    
    ## 4 stage function        1   0.055   0.172              0.117               0    
    ## 5 host number           3   0.059   0.173              0.114               0    
    ## 6 host x stage          1   0.059   0.174              0.115               0    
    ## 7 rel growth rate       1   0.06    0.176              0.116               0

The same results are seen under the assumption of limited growth in
paratenic hosts.

    ## Data: filter(st_level, !is.na(rel_growth_rate_paratenic))
    ## Models:
    ## reg00: hsi_lcdb_suspcious_rem ~ 1 + (1 | Parasite.species)
    ## reg0: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg0:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg0:     (1 | parasite_phylum)
    ## reg1: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg1:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg1:     (1 | parasite_phylum) + zstudy_effort
    ## reg2: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg2:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg2:     (1 | parasite_phylum) + zstudy_effort + Def.int
    ## reg3: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg3:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg3:     (1 | parasite_phylum) + zstudy_effort + Def.int + Host_no_fac
    ## reg4: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg4:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg4:     (1 | parasite_phylum) + zstudy_effort + Def.int + Host_no_fac + 
    ## reg4:     Def.int:Host_no_fac
    ## reg5: hsi_lcdb_suspcious_rem ~ (1 | Parasite.species) + (1 | parasite_genus) + 
    ## reg5:     (1 | parasite_family) + (1 | parasite_order) + (1 | parasite_class) + 
    ## reg5:     (1 | parasite_phylum) + zstudy_effort + Def.int + Host_no_fac + 
    ## reg5:     rel_growth_rate_paratenic + Def.int:Host_no_fac
    ##       Df    AIC    BIC   logLik deviance   Chisq Chi Df Pr(>Chisq)    
    ## reg00  3 2064.8 2078.3 -1029.40   2058.8                              
    ## reg0   8 2029.3 2065.2 -1006.66   2013.3 45.4839      5  1.157e-08 ***
    ## reg1   9 2006.8 2047.1  -994.38   1988.8 24.5698      1  7.166e-07 ***
    ## reg2  10 1987.1 2031.9  -983.53   1967.1 21.6967      1  3.193e-06 ***
    ## reg3  13 1951.2 2009.5  -962.60   1925.2 41.8542      3  4.308e-09 ***
    ## reg4  16 1953.2 2024.9  -960.58   1921.2  4.0363      3     0.2576    
    ## reg5  17 1955.1 2031.3  -960.55   1921.1  0.0660      1     0.7973    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## # A tibble: 7 x 6
    ##   step            df_used marg_r2 cond_r2 rand_var_explained species_var_explai~
    ##   <chr>             <dbl>   <dbl>   <dbl>              <dbl>               <dbl>
    ## 1 within-species~      NA   0       0.233              0.233               0.233
    ## 2 taxonomy              0   0       0.236              0.236               0    
    ## 3 study effort          1   0.043   0.268              0.225               0    
    ## 4 stage function        1   0.07    0.281              0.211               0    
    ## 5 host number           3   0.145   0.315              0.17                0    
    ## 6 host x stage          3   0.148   0.303              0.155               0    
    ## 7 rel growth rate       1   0.148   0.304              0.156               0

# Conclusions

Parasite life cycle, life history, and patterns of host generalism are
related. In particular, parasite exhibit high generalism in second or
third intermediate hosts. These tend to be paratenic hosts, in which
parasites undergo little growth and development. But even if these hosts
are excluded, there is still a negative relationship between generalism
and worm growth. Worms grow less in stages where they infect a broader
range of hosts, which is consistent with costs of generalism.

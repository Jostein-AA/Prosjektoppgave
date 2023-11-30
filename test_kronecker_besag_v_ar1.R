###
#Want to quickly compare how besag kronecker ar1 vs ar1 kronecker besag compares.
test_besag_kronecker_ar1_formula <- deaths ~ 1 + f(county, 
                                                   model = "besag",
                                                   graph = Besag_prec,
                                                   group = year, 
                                                   control.group = list(model = "ar1"))

test_ar1_kronecker_besag_formula <- deaths ~ 1 + f(year, 
                                                   model = "ar1",
                                                   group = county, 
                                                   control.group = list(model = "besag"))

besag_kronecker_ar1_fit <- inla(test_besag_kronecker_ar1_formula, 
                                data = ohio_df,
                                family = "poisson",
                                E = pop_at_risk, 
                                control.compute = list(config = TRUE, # To see constraints later
                                                       cpo = T,   # For model selection
                                                       waic = T)) # For model selection
print(mean(-log(besag_kronecker_ar1_fit$cpo$cpo)))

ar1_kronecker_besag_fit <- inla(test_ar1_kronecker_besag_formuladeaths, 
                                data = ohio_df,
                                family = "poisson",
                                E = pop_at_risk, 
                                control.compute = list(config = TRUE, # To see constraints later
                                                       cpo = T,   # For model selection
                                                       waic = T)) # For model selection
print(mean(-log(ar1_kronecker_besag_fit$cpo$cpo)))
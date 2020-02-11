post_test <- posterior_samples(rand.brm)
head(post_test)
plot(post_test$b_Intercept)
plot(density(post_test$b_Intercept))
median(post_test$b_Intercept)
print(rand.brm)
quantile(post_test$b_Intercept, probs=c(0.025,0.975))
marginal_effects(rand.brm, method="predict")
print(rand.brm, prior=T)
pp_check(rand.brm, type="boxplot") #posterior predictive check/distribution

post_test <- posterior_samples(L.brm)
b_L <- exp(post_test$b_pca1) # exp to make interpretable
quantile(1-b_L, probs=c(0.025,0.5, 0.975)) # estimate CrI
plot(b_L)
head(b_L)
sum(b_L<1)
#For every unit increase in x, y declined by an average of 30%. There was a 0.99 probability that this decline was less than 100%.
sum(b_L<.7)/1000
print(L.brm)
# estimate no links @ pca value
b_L_neg6 <- exp(post_test$b_Intercept + post_test$b_pca1*-4)
quantile(b_L_neg6, probs = c(0.025, 0.5, 0.975))
b_L_1 <- exp(post_test$b_Intercept + post_test$b_pca1*-1)
plot(b_L_1)
b_L_6 <- exp(post_test$b_Intercept + post_test$b_pca1*6)
head(b_L_6)



# arbitrary change in response across x_values
diff <- b_L_neg61 - b_L_6
plot(density(diff))
quantile(diff, probs = c(0.025, 0.5, 0.975))

fold <- b_L_neg6/b_L_6
quantile(fold, probs = c(0.025, 0.5, 0.975))

# visually inspect chains
plot(L.brm)

# compare data and posterior predictive values
pp_check(L.brm)
pp_check(L.brm, type = "boxplot")

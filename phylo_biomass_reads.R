##### Phylogenetic signal modelling of HTS read abundance and speices biomass #####

setwd("/Users/mingxin/UTAS/MS/3MS_DNA_Metabarcoding_Correlation/Quantitative/")

##### load required packages #######
library(tidyverse)
library(ggplot2)
library(phytools)
library(lme4)
library(MuMIn)
library(picante)
library(brms)

##### violin plot of the number of mismatches for 16S and COI #####
mis.dat <- read.csv("primer_mismatch_violinplot.csv", sep=",", header = TRUE)

violin.plot <- ggplot(mis.dat, aes(x=marker, y=mismatch)) + geom_boxplot() + xlab("Markers") +
  ylab("Number of primer-template mismatches") + theme_bw() + scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7))

violin.plot

##### QUESTION 1: are there correlation between sequence abundance and species biomass?
##### load species biomass and primer mismatch data
sp_biomass <- read.csv("sp_biomass.csv", sep = ",", header = TRUE)
sp_primer_mismatch <- read.csv("sp_primer_mismatch.csv", sep = ",", header = TRUE)

##### 16S dataset #####
biomass_reads_16S <- read.csv("biomass_reads_16S.csv", sep =",", header = TRUE)
biomass_reads_16S$biomass_weight <- sp_biomass$biomass_weight[match(biomass_reads_16S$Species, sp_biomass$Species)] * biomass_reads_16S$Individual
biomass_reads_16S <- merge(biomass_reads_16S, sp_primer_mismatch, by.x = "Species")
biomass_reads_16S <- biomass_reads_16S %>% filter(Individual != 0 & Reads != 0)
biomass_reads_16S <- biomass_reads_16S %>% mutate(Prop_reads = Reads/Sample_total_reads)
biomass_reads_16S <- biomass_reads_16S %>% mutate(Prop_biomass = biomass_weight/Sample_total_biomass)

##### plot 16S linear regression
lm_model_16S <- lm(Prop_reads ~ Prop_biomass, data = biomass_reads_16S)
summary(lm_model_16S)

p_16S_lm <- ggplot(biomass_reads_16S, aes(x=Prop_biomass, y=Prop_reads)) +
  geom_smooth(formula = y ~ x, method = "lm", se=FALSE) + geom_point(aes(size=1),shape=21) + theme_bw() +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p_16S_lm

##### do AIC on all models:
lm_16S_model1 <- lm(Prop_reads ~ Prop_biomass, data = biomass_reads_16S)
summary(lm_16S_model1)

lm_16S_model2 <- lm(Prop_reads ~ mismatch_16S, data = biomass_reads_16S)
summary(lm_16S_model2)

lm_16S_model3 <- lm(Prop_reads ~ mismatch_16S_position, data = biomass_reads_16S)
summary(lm_16S_model3)

lm_16S_model4 <- lm(Prop_reads ~ mismatch_16S_type, data = biomass_reads_16S)
summary(lm_16S_model4)

lm_16S_model5 <- lm(Prop_reads ~ Prop_biomass + mismatch_16S, data = biomass_reads_16S)
summary(lm_16S_model5)

lm_16S_model6 <- lm(Prop_reads ~ Prop_biomass + mismatch_16S_position, data = biomass_reads_16S)
summary(lm_16S_model6)

lm_16S_model7 <- lm(Prop_reads ~ Prop_biomass + mismatch_16S_type, data = biomass_reads_16S)
summary(lm_16S_model7)

lm_16S_model8 <- lm(Prop_reads ~ Prop_biomass + mismatch_16S_position + mismatch_16S_type, data = biomass_reads_16S)
summary(lm_16S_model8)

AICc(lm_16S_model1, lm_16S_model2, lm_16S_model3,lm_16S_model4, lm_16S_model5, lm_16S_model6, lm_16S_model7, lm_16S_model8)
Weights(AICc(lm_16S_model1, lm_16S_model2, lm_16S_model3,lm_16S_model4, lm_16S_model5, lm_16S_model6, lm_16S_model7, lm_16S_model8))

##### add lm residuals into the dataframe
biomass_reads_16S$true_residuals <- lm_model_16S$residuals

res_16S <- as.data.frame(biomass_reads_16S[, c("Species", "true_residuals")])

res_16S <- res_16S %>% group_by(Species) %>% summarise(true_residuals = mean(true_residuals))

rownames(res_16S) <- res_16S$Species

res_16S <- setNames(res_16S$true_residuals, rownames(res_16S))

##### CO1 dataset #####
biomass_reads_CO1 <- read.csv("biomass_reads_CO1.csv", sep =",", header = TRUE)
biomass_reads_CO1$biomass_weight <- sp_biomass$biomass_weight[match(biomass_reads_CO1$Species, sp_biomass$Species)] * biomass_reads_CO1$Individual
biomass_reads_CO1 <- merge(biomass_reads_CO1, sp_primer_mismatch, by.x = "Species")
biomass_reads_CO1 <- biomass_reads_CO1 %>% filter(Individual != 0 & Reads != 0)
biomass_reads_CO1 <- biomass_reads_CO1 %>% mutate(Prop_reads = Reads/Sample_total_reads)
biomass_reads_CO1 <- biomass_reads_CO1 %>% mutate(Prop_biomass = biomass_weight/Sample_total_biomass)

##### plot CO1 linear regression
lm_model_CO1 <- lm(Prop_reads ~ Prop_biomass, data = biomass_reads_CO1)
summary(lm_model_CO1)

p_CO1_lm <- ggplot(biomass_reads_CO1, aes(x=Prop_biomass, y=Prop_reads)) + theme_bw() +
  geom_smooth(formula = y ~ x, method = "lm", se=FALSE) + geom_point(aes(size=1),shape=21) + 
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p_CO1_lm

##### a few linear models
lm_CO1_model1 <- lm(Prop_reads ~ Prop_biomass, data = biomass_reads_CO1)
summary(lm_CO1_model1)

lm_CO1_model2 <- lm(Prop_reads ~ mismatch_CO1, data = biomass_reads_CO1)
summary(lm_CO1_model2)

lm_CO1_model3 <- lm(Prop_reads ~ mismatch_CO1_position, data = biomass_reads_CO1)
summary(lm_CO1_model3)

lm_CO1_model4 <- lm(Prop_reads ~ mismatch_CO1_type, data = biomass_reads_CO1)
summary(lm_CO1_model4)

lm_CO1_model5 <- lm(Prop_reads ~ Prop_biomass + mismatch_CO1, data = biomass_reads_CO1)
summary(lm_CO1_model5)

lm_CO1_model6 <- lm(Prop_reads ~ Prop_biomass + mismatch_CO1_position, data = biomass_reads_CO1)
summary(lm_CO1_model6)

lm_CO1_model7 <- lm(Prop_reads ~ Prop_biomass + mismatch_CO1_type, data = biomass_reads_CO1)
summary(lm_CO1_model7)

lm_CO1_model8 <- lm(Prop_reads ~ Prop_biomass + mismatch_CO1_position + mismatch_CO1_type, data = biomass_reads_CO1)
summary(lm_CO1_model8)

AICc(lm_CO1_model1, lm_CO1_model2, lm_CO1_model3, lm_CO1_model4, lm_CO1_model5, lm_CO1_model6, lm_CO1_model7, lm_CO1_model8)
Weights(AICc(lm_CO1_model1, lm_CO1_model2, lm_CO1_model3, lm_CO1_model4, lm_CO1_model5, lm_CO1_model6, lm_CO1_model7, lm_CO1_model8))

##### add lm residuals into the dataframe
biomass_reads_CO1$true_residuals <- lm_model_CO1$residuals

res_CO1 <- as.data.frame(biomass_reads_CO1[, c("Species", "true_residuals")])

res_CO1 <- res_CO1 %>% group_by(Species) %>% summarise(true_residuals = mean(true_residuals))

rownames(res_CO1) <- res_CO1$Species

res_CO1 <- setNames(res_CO1$true_residuals, rownames(res_CO1))

##### plot Figure 2 #####
#### set color ramps
lims_res <- range(min(res_CO1), abs(min(res_CO1)))

lims_mis <- range(0,7)

cols_res <- colorRampPalette(c("darkblue", "white", "darkred"))

cols_mis <- colorRampPalette(c("darkblue","white", "darkred"))

cols1 <- cols_res(200)

cols2 <- cols_mis(8)

##### read the best 16S tree
tree_16S <- read.tree("beetles.tre")

##### drop species not in 16S
tree_16S <- drop.tip(tree_16S, c("Neuroptera_Apochrysa_matsumurae",
                                 "Neuroptera_Ascaloptynx_appendiculatus",
                                 "Strepsiptera_Mengenilla_moldrzyki", 
                                 "Carabidae_Pseudocensis_solicitus",
                                 "Carabidae_Rhabdotus_reflexus",
                                 "Carabidae_Percosoma_carenoides",
                                 "Carabidae_Scopodes_intermedius",
                                 "Melandryidae_Orchesia_alphabetica",
                                 "Carabidae_Chylnus_ater",
                                 "Nitidulidae_Thalycrodes_pulchrum",
                                 "Staphylinidae_Quedius_TFIC_sp02",
                                 "Staphylinidae_Quedius_TFIC_sp03",
                                 "Staphylinidae_Euconnus_clarus",
                                 "Staphylinidae_Quedius_inadequalipennis"))

plotTree(tree_16S)

##### fast estimate ancestral states for residuals
anc_res_16S <- contMap(tree_16S, res_16S, plot = FALSE, lims = lims_res, res = 1000)

anc_res_16S <- setMap(anc_res_16S, cols1)

plot(anc_res_16S)

##### read primer mismatch file
mismatch_16S <- read.csv("mismatch_16S.csv", sep = ",", row.names = 1)

mismatch_16S <- as.matrix(mismatch_16S)[,1]

##### fast estimate ancestral states for the number of mismatches
anc_mismatch_16S <- contMap(tree_16S, 7 - mismatch_16S, plot = FALSE, lims = lims_mis)

anc_mismatch_16S <- setMap(anc_mismatch_16S, cols2)

plot(anc_mismatch_16S)

##### read the best ML CO1 tree
tree_CO1 <- read.tree("beetles.tre")

##### drop the outgroup
tree_CO1 <- drop.tip(tree_CO1, c("Neuroptera_Apochrysa_matsumurae",
                                 "Neuroptera_Ascaloptynx_appendiculatus",
                                 "Strepsiptera_Mengenilla_moldrzyki", 
                                 "Staphylinidae_Zyras_TFIC_sp01",
                                 "Staphylinidae_Spanioda_carissima",
                                 "Staphylinidae_Quedius_TFIC_sp04",
                                 "Staphylinidae_Quedius_stenocephalus",
                                 "Staphylinidae_Microsilpha_TFIC_sp15",
                                 "Staphylinidae_Aleoc_R",
                                 "Leiodidae_Nargomorphus_globulus",
                                 "Leiodidae_Nargiotes_gordini",
                                 "Leiodidae_Catoposchema_tasmaniae",
                                 "Leiodidae_Austronemadus_TFIC_sp02",
                                 "Latridiidae_Cortinicara_TFIC_sp01",
                                 "Latridiidae_Aridius_nodifer",
                                 "Curculionidae_Dinichus_terreus",
                                 "Curculionidae_Decilaus_striatus"))

plotTree(tree_CO1)
nodelabels()

tree_CO1 <- rotate(tree_CO1, 25)
tree_CO1 <- rotate(tree_CO1, 29)
tree_CO1 <- rotate(tree_CO1, 30)
tree_CO1 <- rotate(tree_CO1, 34)
tree_CO1 <- rotate(tree_CO1, 32)
tree_CO1 <- rotate(tree_CO1, 37)
tree_CO1 <- rotate(tree_CO1, 40)
tree_CO1 <- rotate(tree_CO1, 43)
plotTree(tree_CO1)

##### fast estimate ancestral states for a continuous traits
anc_res_CO1 <- contMap(tree_CO1, res_CO1, plot = FALSE, lims = lims_res, res = 1000)

anc_res_CO1 <- setMap(anc_res_CO1, cols1)

plot(anc_res_CO1)

##### read CO1 mismatches
mismatch_CO1 <- read.csv("mismatch_CO1.csv", row.names = 1)

mismatch_CO1 <- as.matrix(mismatch_CO1)[,1]

##### ancestral state for the number of CO1 mismatches
anc_mismatch_CO1 <- contMap(tree_CO1, 7 - mismatch_CO1, plot = FALSE, lims = lims_mis)

anc_mismatch_CO1 <- setMap(anc_mismatch_CO1, cols2)

plot(anc_mismatch_CO1)

################# QUESTION 2: Is there any phylogenetic preference in the deviation of sequence abundance estimation?
##### brms modelling #####
# 16S dataset
C.16S <- ape::vcv.phylo(tree_16S, corr = T)
p.b.16S <- set_prior("student_t(3, 0, 2.5)", class = "b", lb = c(0,0,NA))

b.16S <- brm(Prop_reads ~ Prop_biomass + (1|gr(Species, cov=C.16S)),
             data=biomass_reads_16S, family=gaussian(),
             data2=list(C.16S=C.16S), cores=2,
             chains=2, iter=5000000, thin = 5)

# An estimate of phylogenetic signal 
# (the proportion of variation in a trait attributed to phylogenetic effects)
b.16S %>% as_tibble() %>% dplyr::select(sigma_b = sd_Species__Intercept, sigma_e = sigma) %>% 
  mutate(h2 = sigma_b^2/(sigma_b^2 + sigma_e^2)) %>% 
  pull(h2) %>% quantile(probs=c(0.025,0.5,0.975))

# CO1 dataset
C.CO1 = ape::vcv.phylo(tree_CO1, corr = T)

p.b.CO1 <- set_prior("student_t(3, 0, 2.5)", class = "b", lb = c(0,0,NA))

b.CO1 <- brm(Prop_reads ~ Prop_biomass + (1|gr(Species, cov=C.CO1)),
           data=biomass_reads_CO1, family = gaussian(),
           data2 = list(C.CO1=C.CO1),cores=2,
           chains=2, iter=5000000, thin = 5)

# An estimate of phylogenetic signal 
# (the proportion of variation in a trait attributed to phylogenetic effects)
b.CO1 %>% as_tibble() %>% 
  dplyr::select(sigma_b = sd_Species__Intercept, sigma_e = sigma) %>% 
  mutate(h2 = sigma_b^2/(sigma_b^2 + sigma_e^2)) %>% 
  pull(h2) %>% quantile(probs=c(0.025,0.5,0.975))

####### PGLS models
##### create caper object for 16S data
comp.dat.16S <- comparative.data(tree_16S, biomass_reads_16S, names.col = "Species", na.omit = FALSE, warn.dropped = TRUE)

##### Phylogenetic generalized least squares regression
model.16S.1 <- pgls(log10(reads) ~ log10(biomass_weight), data=comp.dat.16S, lambda = "ML")
summary(model.16S.1)
lambda.model.16S.1 <- pgls.profile(model.16S.1, "lambda")
plot(lambda.model.16S.1)

model.16S.2 <- pgls(log10(reads) ~ mismatch_16S, data=comp.dat.16S, lambda = "ML")
summary(model.16S.2)

model.16S.3 <- pgls(log10(reads) ~ mismatch_16S_position, data=comp.dat.16S, lambda = "ML")
summary(model.16S.3)

model.16S.4 <- pgls(log10(reads) ~ mismatch_16S_type, data=comp.dat.16S, lambda = "ML")
summary(model.16S.4)

model.16S.5 <- pgls(log10(reads) ~ log10(biomass_weight) + mismatch_16S, data=comp.dat.16S, lambda = "ML")
summary(model.16S.5)

model.16S.6 <- pgls(log10(reads) ~ log10(biomass_weight) + mismatch_16S_position, data=comp.dat.16S, lambda = "ML")
summary(model.16S.6)

model.16S.7 <- pgls(log10(reads) ~ log10(biomass_weight) + mismatch_16S_type, data=comp.dat.16S, lambda = "ML")
summary(model.16S.7)

model.16S.8 <- pgls(log10(reads) ~ log10(biomass_weight) + mismatch_16S_position + mismatch_16S_type, data=comp.dat.16S, lambda = "ML")
summary(model.16S.8)

AICc(model.16S.1, model.16S.2, model.16S.3, model.16S.4, model.16S.5, model.16S.6, model.16S.7, model.16S.8)

Weights(AICc(model.16S.1, model.16S.2, model.16S.3, model.16S.4))

# test phylogenetic signal of traits with separate variables
vec_16S_reads <- biomass_reads_16S$Prop_reads
names(vec_16S_reads) <- biomass_reads_16S$Species
# Blomberg's K
Kcalc(vec_16S_reads[tree_16S$tip.label], tree_16S)
phylosignal(vec_16S_reads[tree_16S$tip.label], tree_16S, reps=1000)
# Pagel's lambda
phylosig(tree_16S, vec_16S_reads, test=TRUE, method = "lambda")

#
vec_16S_biomass <- biomass_reads_16S$Prop_biomass
names(vec_16S_biomass) <- biomass_reads_16S$Species
Kcalc(vec_16S_biomass[tree_16S$tip.label], tree_16S)
phylosignal(vec_16S_biomass[tree_16S$tip.label], tree_16S, reps = 1000)
phylosig(tree_16S, vec_16S_biomass, test=TRUE, method = "lambda")

#
vec_16s_res <- biomass_reads_16S$true_residuals
names(vec_16s_res) <- biomass_reads_16S$Species
Kcalc(vec_16s_res[tree_16S$tip.label], tree_16S)
phylosignal(vec_16s_res[tree_16S$tip.label], tree_16S, reps=1000)
phylosig(tree_16S, vec_16s_res, test=TRUE, method="lambda")

#
vec_16S_mismatch <- biomass_reads_16S$mismatch_16S
names(vec_16S_mismatch) <- biomass_reads_16S$Species
Kcalc(vec_16S_mismatch[tree_16S$tip.label], tree_16S)
phylosignal(vec_16S_mismatch[tree_16S$tip.label], tree_16S, reps=1000)
phylosig(tree_16S, vec_16S_mismatch, test=TRUE, method="lambda")

#
vec_16S_mismatch_position <- biomass_reads_16S$mismatch_16S_position
names(vec_16S_mismatch_position) <- biomass_reads_16S$Species
Kcalc(vec_16S_mismatch_position[tree_16S$tip.label], tree_16S)
phylosignal(vec_16S_mismatch_position[tree_16S$tip.label], tree_16S, reps=1000)
phylosig(tree_16S, vec_16S_mismatch_position, test=TRUE, method="lambda")

#
vec_16S_mismatch_type <- biomass_reads_16S$mismatch_16S_type
names(vec_16S_mismatch_type) <- biomass_reads_16S$Species
Kcalc(vec_16S_mismatch_type[tree_16S$tip.label], tree_16S)
phylosignal(vec_16S_mismatch_type[tree_16S$tip.label], tree_16S, reps=1000)
phylosig(tree_16S, vec_16S_mismatch_type, test=TRUE, method="lambda")

##### create caper object for COI data
comp.dat.CO1 <- comparative.data(tree_CO1, biomass_reads_CO1, names.col = "Species", vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)

##### Phylogenetic generalized least squares regression
model.CO1.1 <- pgls(log10(reads) ~ log10(biomass_weight), data=comp.dat.CO1, lambda = "ML")
summary(model.CO1.1)
lambda.model.CO1.1 <- pgls.profile(model.CO1.1, "lambda")
plot(lambda.model.CO1.1)

model.CO1.2 <- pgls(log10(reads) ~ mismatch_CO1, data=comp.dat.CO1, lambda = "ML")
summary(model.CO1.2)

model.CO1.3 <- pgls(log10(reads) ~ mismatch_CO1_position, data=comp.dat.CO1, lambda = "ML")
summary(model.CO1.3)

model.CO1.4 <- pgls(log10(reads) ~ mismatch_CO1_type, data=comp.dat.CO1, lambda = "ML")
summary(model.CO1.4)

model.CO1.5 <- pgls(log10(reads) ~ log10(biomass_weight) + mismatch_CO1, data=comp.dat.CO1, lambda = "ML")
summary(model.CO1.5)

model.CO1.6 <- pgls(log10(reads) ~ log10(biomass_weight) + mismatch_CO1_position, data=comp.dat.CO1, lambda = "ML")
summary(model.CO1.6)

model.CO1.7 <- pgls(log10(reads) ~ log10(biomass_weight) + mismatch_CO1_type, data=comp.dat.CO1, lambda = "ML")
summary(model.CO1.7)

model.CO1.8 <- pgls(log(reads) ~ log10(biomass_weight) + mismatch_CO1_type + mismatch_CO1_position, data=comp.dat.CO1, lambda = "ML")
summary(model.CO1.8)

AICc(model.CO1.1, model.CO1.2, model.CO1.3, model.CO1.4, model.CO1.5, model.CO1.6, model.CO1.7, model.CO1.8)

Weights(AICc(model.CO1.1, model.CO1.2, model.CO1.3, model.CO1.4))

#
vec_CO1_reads <- biomass_reads_CO1$Prop_reads
names(vec_CO1_reads) <- biomass_reads_CO1$Species
Kcalc(vec_CO1_reads[tree_CO1$tip.label], tree_CO1)
phylosignal(vec_CO1_reads[tree_CO1$tip.label], tree_CO1, reps=1000)
phylosig(tree_CO1, vec_CO1_reads, test = TRUE, method = "lambda")

#
vec_CO1_biomass <- biomass_reads_CO1$Prop_biomass
names(vec_CO1_biomass) <- biomass_reads_CO1$Species
Kcalc(vec_CO1_biomass[tree_CO1$tip.label], tree_CO1)
phylosignal(vec_CO1_biomass[tree_CO1$tip.label], tree_CO1, reps=1000)
phylosig(tree_CO1, vec_CO1_biomass, test = TRUE, method = "lambda")

#
vec_CO1_res <- biomass_reads_CO1$true_residuals
names(vec_CO1_res) <- biomass_reads_CO1$Species
Kcalc(vec_CO1_res[tree_CO1$tip.label], tree_CO1)
phylosignal(vec_CO1_res[tree_CO1$tip.label], tree_CO1, reps=1000)
phylosig(tree_CO1, vec_CO1_res, test = TRUE, method = "lambda")

#
vec_CO1_mismatch <- biomass_reads_CO1$mismatch_CO1
names(vec_CO1_mismatch) <- biomass_reads_CO1$Species
Kcalc(vec_CO1_mismatch[tree_CO1$tip.label], tree_CO1)
phylosignal(vec_CO1_mismatch[tree_CO1$tip.label], tree_CO1, reps=1000)
phylosig(tree_CO1, vec_CO1_mismatch, test = TRUE, method = "lambda")

#
vec_CO1_mismatch_type <- biomass_reads_CO1$mismatch_CO1_type
names(vec_CO1_mismatch_type) <- biomass_reads_CO1$Species
Kcalc(vec_CO1_mismatch_type[tree_CO1$tip.label], tree_CO1)
phylosignal(vec_CO1_mismatch_type[tree_CO1$tip.label], tree_CO1, reps=1000)
phylosig(tree_CO1, vec_CO1_mismatch_type, test = TRUE, method = "lambda")

#
vec_CO1_mismatch_position <- biomass_reads_CO1$mismatch_CO1_position
names(vec_CO1_mismatch_position) <- biomass_reads_CO1$Species
Kcalc(vec_CO1_mismatch_position[tree_CO1$tip.label], tree_CO1)
phylosignal(vec_CO1_mismatch_position[tree_CO1$tip.label], tree_CO1, reps=1000)
phylosig(tree_CO1, vec_CO1_mismatch_position, test = TRUE, method = "lambda")

#
vec_CO1_mismatch_scores <- biomass_reads_CO1$mismatch_scores
names(vec_CO1_mismatch_scores) <- biomass_reads_CO1$Species
Kcalc(vec_CO1_mismatch_scores[tree_CO1$tip.label], tree_CO1)
phylosignal(vec_CO1_mismatch_scores[tree_CO1$tip.label], tree_CO1, reps=1000)
phylosig(tree_CO1, vec_CO1_mismatch_scores, test = TRUE, method = "lambda")


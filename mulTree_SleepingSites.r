# Run analyses on sleeping site dataset using mulTree

## Loading the package
library(mulTree)
library(geiger)
library(phytools)

#set directory
setwd("../../../MCMCglmm_analyses/")


## import dataset
data_Sleep <- read.csv("input_glmm.csv")
head(data_Sleep)
data<- data_Sleep[1:8]
dim(data)

#delete NA entries from data file
na.omit(data) -> data

## import trees
trees1000 <- read.nexus("../../Upham2019_458taxa.nex")
class(trees1000)

#select random tree
rand_tree<-trees1000$tree_9317
row.names(data)<-data$new_names
head(data)
chk<-name.check(rand_tree,data)

# prune multiPhylo object from taxa not present in dataset
pruned<-drop.tip.multiPhylo(trees1000,chk$tree_not_data)

# delete taxa not present in the tree
data<-data[-which(rownames(data) %in% chk$data_not_tree),]
name.check(pruned$tree_9317,data)

#make all tree ultrametric
pruned<-lapply(pruned,force.ultrametric)
class(pruned)<-"multiPhylo"


#add logMass to data
data<-cbind(data, log(data$BodyMass_kg))
names(data)[names(data) == 'log(data$BodyMass_kg)'] <- 'logMass'
head(data)

pruned_100<-sample(pruned,size=100)

## preparing the mulTree data: as.mulTree
## Creating the mulTree object
mulTree_data <- as.mulTree(data = data, tree = pruned, taxa = "new_names")

mulTree_data_100 <- as.mulTree(data = data, tree = pruned_100, taxa = "new_names")


## This object is classified as "mulTree" and contains different elements to be ## passed to the mulTree function
class(mulTree_data) ; names(mulTree_data)
class(mulTree_data_100) ; names(mulTree_data_100)
## Running the MCMCglmm on multiple trees: mulTree

## The glmm formula
my_formula <- tree.holes ~ activity + logMass
my_formula <- on.tree ~ activity + logMass
my_formula <- nest ~ activity + logMass

my_formula <- logMass ~ activity

## The MCMC parameters (generations, sampling, burnin)
my_parameters <- c(100000, 10, 1000)

## The MCMCglmm priors
my_priors <- list(R = list(V = 1/2, nu = 0.002), G = list(G1 = list(V = 1/2, nu = 0.002)))

## Running the MCMCglmm on multiple trees (1000)
mulTree(mulTree.data = mulTree_data, formula = my_formula, priors = my_priors, parameters = my_parameters, output = "tree_hole_example", ESS = 50, chains = 2)

## Running the MCMCglmm on multiple trees (100)
mulTree(mulTree.data = mulTree_data_100, formula = my_formula, priors = my_priors, parameters = my_parameters, output = "tree_hole_example_100", ESS = 50, chains = 2)

mulTree(mulTree.data = mulTree_data_100, formula = my_formula, priors = my_priors, parameters = my_parameters, output = "on.tree_example_100", ESS = 50, chains = 2)

mulTree(mulTree.data = mulTree_data_100, formula = my_formula, priors = my_priors, parameters = my_parameters, output = "nest_example_100", ESS = 50, chains = 2)

## Reading the models: read.mulTree

## Reading only one specific model
one_model <- read.mulTree("tree_hole_example-tree1_chain1", model = TRUE)

## This model is a normal MCMCglmm object that has been ran on one single tree 
class(one_model) ; names(one_model)

## Reading the convergence diagnosis test to see if the two chains converged for each tree
read.mulTree("tree_hole_example", convergence = TRUE)
read.mulTree("tree_hole_example_100", convergence = TRUE)
read.mulTree("on.tree_example_100", convergence = TRUE)


## As indicated here, the chains converged for both chains!
## Reading all the models to perform the MCMCglmm analysis on multiple trees
all_models <- read.mulTree("tree_hole_example")
str(all_models)

all_models_100_TR <- read.mulTree("on.tree_example_100")
str(all_models_100_TR)

all_models_100_NS <- read.mulTree("nest_example_100")
str(all_models_100_NS)

## This object contains 39600 estimations of the Intercept and the terms!
## Removing the chains from the current directory
file.remove(list.files(pattern="longevity_example"))

## Summarising the results by estimating the highest density regions ## and their associated 95 and 50 confidence intervals (default) 
summarised_results <- summary(all_models)
summarised_results

summarised_results_100 <- summary(all_models_100)
summarised_results_100

summarised_results_100_TR <- summary(all_models_100_TR)
summarised_results_100_TR

summarised_results_100_NS <- summary(all_models_100_NS)
summarised_results_100_NS

summary(all_models, use.hdr = FALSE, cent.tend = mean, prob = c(75, 25))

########################################################
### 

#best model tree holes => AP*BM file=tree_hole_example_100_M1c
all_models_100_TH <- read.mulTree("tree_hole_example_100_M1c")
str(all_models_100_TH)
summarised_results_TH <- summary(all_models_100_TH)
summarised_results_TH


# Visualising the results: plot.mulTree

## Graphical options
quartz(width = 10, height = 5) ; par(mfrow = (c(1,2)), bty = "n") 

## Plotting using the default options
plot(summarised_results_TH)

## Plotting using some more pretty options
plot(summarised_results_TH, horizontal = TRUE, ylab = "", cex.coeff = 0.8,
     main = "Posterior distributions", ylim = c(-2,2), cex.terms = 0.5,
     terms = c("Intercept", "activitynocturnal", "logMass", "activitynocturnal:logMass","Phylogeny", "Residuals"), col = "gray", cex.main = 0.8)
abline(v = 0, lty = 3)

###########
#best model nest => AP file=nest_example_100_M4
all_models_100_NS <- read.mulTree("nest_example_100_M4")
str(all_models_100_NS)
summarised_results_NS <- summary(all_models_100_NS)
summarised_results_NS


# Visualising the results: plot.mulTree

## Graphical options
quartz(width = 10, height = 5) ; par(mfrow = (c(1,2)), bty = "n") 

## Plotting using the default options
#plot(summarised_results_NS)

## Plotting using some more pretty options
plot(summarised_results_NS, horizontal = TRUE, ylab = "", cex.coeff = 0.8,
     main = "Posterior distributions", ylim = c(-2,2), cex.terms = 0.5,
     terms = c("Intercept", "activitynocturnal", "logMass", "activitynocturnal:logMass","Phylogeny", "Residuals"), col = "gray", cex.main = 0.8)
abline(v = 0, lty = 3)

###########
#best model on_trees => AP file=on.tree_example_100_M4
all_models_100_TR <- read.mulTree("on.tree_example_100_M4")
str(all_models_100_TR)
summarised_results_TR <- summary(all_models_100_TR)
summarised_results_TR


# Visualising the results: plot.mulTree

## Graphical options
quartz(width = 10, height = 5) ; par(mfrow = (c(1,2)), bty = "n") 

## Plotting using the default options
#plot(summarised_results_NS)

## Plotting using some more pretty options
plot(summarised_results_TR, horizontal = TRUE, ylab = "", cex.coeff = 0.8,
     main = "Posterior distributions", ylim = c(-2,2), cex.terms = 0.5,
     terms = c("Intercept", "activitynocturnal", "logMass", "activitynocturnal:logMass","Phylogeny", "Residuals"), col = "gray", cex.main = 0.8)
abline(v = 0, lty = 3)



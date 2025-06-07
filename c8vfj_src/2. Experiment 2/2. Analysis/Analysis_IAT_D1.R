	############################
	######    PACKAGES    ######
	############################

library(reshape2)
library(lattice)
library(latticeExtra)
library(car)
library(heplots)
library(lsr)
library(MBESS)
library(plyr)
library(dplyr)	
library(grid)
library(BayesFactor)
library(TOSTER)

	############################
	######    SETTINGS    ######
	############################

options(contrasts = c("contr.sum", "contr.poly"))

	#############################
	######    FUNCTIONS    ######
	#############################

######################################
######    READ DATA	FUNCTION    ######
######################################

make_list <- function(ppn){

	###### READ DATA ######

	if(ppn < 10){
		file <- sprintf("C:\\Users\\emiel\\OneDrive - UGent\\Studenten\\Laura\\Data\\Exp 2\\Data IAT\\dataiat0%d.txt", ppn)
	}else{
		file <- sprintf("C:\\Users\\emiel\\OneDrive - UGent\\Studenten\\Laura\\Data\\Exp 2\\Data IAT\\dataiat%d.txt", ppn)
	}

	data <- read.table(file, header = F)
	
	###### NAME VARIABLES ######

	names <- c("trialnr","ppn","age","gender","handpos","group","block","targetconcept","targetcol","attrdim","attrval","stim","color_left","loc","con","xr","r","corr","rt","re")
	names(data) <- names
	
	###### RT > 10000 MS ######
	
	data$no_resp <- ifelse(data$rt >= 2000, 0, 1)
	
	###### RETURN DATA ######
			
	return(data)
}

###############################
######    SE FUNCTION    ######
###############################

SE_Within <- function(data, colnames, ref = F){

	###### CREATE DATA FRAME WITH RELEVANT DATA ######

	data_norm <- data.frame(matrix(nrow = nrow(data), ncol = length(colnames) + 2))
	names(data_norm) <- c("subj_mean", "grand_mean", colnames)

	###### COMPLETE DATASET ######

		###### SUBJECT MEAN ######

	data_norm$subj_mean <- rowMeans(data[,c(colnames)])

		###### GRAND MEAN ######

	data_norm$grand_mean <- mean(as.matrix(data[,c(colnames)]))

		###### NORMALIZED DATA ######

	for(i in 1:length(colnames)){
		for(j in 1:nrow(data_norm)){
			data_norm[j, colnames[i]] <- data[j, colnames[i]] - data_norm$subj_mean[j] +  data_norm$grand_mean[j]
		}
	}

	###### CALCULATE STANDARD ERRORS ACCORDING TO DIFFERENT METHODS ######

		###### BETWEEN-SUBJECT SE AND CI ######

	sd.between <- as.vector(sapply(data[,c(colnames)], sd))
	se.between <- sd.between / sqrt(nrow(data_norm))
	ci95.between <- se.between * qt(0.975, nrow(data_norm) - 1)

		###### COUSINEAU SE AND CI ######

	var.within <- as.vector(sapply(data_norm[,c(colnames)], var))
	se.within <- sqrt(var.within) / sqrt(nrow(data_norm))
	ci95.within <- se.within * qt(0.975, nrow(data_norm) - 1)

		###### MOREY SE AND CI ######

	var.within.corr <- var.within * (length(colnames) / (length(colnames) - 1))	
	sd.within.corr <- sqrt(var.within.corr)	
	se.within.corr <- sd.within.corr / sqrt(nrow(data_norm))	
	ci95.within.corr <- se.within.corr * qt(0.975, nrow(data_norm) - 1)
	ci95.baguley <- ci95.within.corr * (sqrt(2)/2)

	###### RETURN STANDARD ERROR INFORMATION ######

	output <- data.frame(row.names = colnames, se.between, ci95.between, se.within, ci95.within,
				   se.within.corr, ci95.within.corr, ci95.baguley)
	
	if(ref){
		cat("references:","\n")
		within.ref <- "within: 'Cousineau, D. (2005). Confidence intervals in within-subject designs: A simpler solution to Loftus and Masson's method. Tutorials in Quantitative Methods for Psychology, 2005, 1, 42-45'"
		cat("  ",within.ref,"\n")
		within.corr.ref <- "within.corr: 'Morey, R.D. (2008). Confidence intervals from normalized data: A correction to Cousineau (2005). Tutorials in Quantitative Methods for Psychology, 4, 61-64'"
		cat("  ",within.corr.ref,"\n")
		baguley.ref <- "baguley: 'Baguley T. (2012). Calculating and graphing within-subject confidence intervals for ANOVA. Behav Res, 44, 158-175'"
		cat("  ",baguley.ref,"\n","\n")
	}

	return(output)	
}

	####################################
	######    DATA PREPARATION    ######
	####################################

#############################
######    READ DATA    ######
#############################

	###### CLEAR WORKING SPACE ######
	
rm(list = setdiff(ls(), c("make_list", "SE_Within")))	

	###### MAKE LIST ######

ppn <- c(1:39)
data.list <- lapply(ppn, make_list)
sapply(data.list, nrow)
sapply(data.list, function(x) mean(x$rt, na.rm = T))

	###### MAKE DATA FRAME ######
	
data.full <- do.call(rbind, data.list)

head(data.full)
tail(data.full)
summary(data.full)
str(data.full)

################################################
######    DELETE UNIMPORTANT VARIABLES    ######
################################################

head(data.full)
data.full <- select(data.full, -one_of("trialnr","group","targetconcept","targetcol","attrdim","attrval","stim","color_left","loc","xr","r","re"))
head(data.full)
tail(data.full)

test <- filter(data.full, ppn == 1)

#############################################
######    SPECIFY TYPE OF VARIABLES    ######
#############################################

str(data.full)

data.full$age <- as.numeric(data.full$age)
	#age
data.full$handpos <- factor(data.full$handpos)
	#position of stimulus hands
data.full$block <- factor(data.full$block + 1)
	#block number
data.full$con <- factor(data.full$con, levels = c(1,0), labels = c("C", "IC"))
	#congruency (i.e., congruent or incongruent)
data.full$corr <- factor(data.full$corr)
	#accuracy: 0 = error, 1 = correct
data.full$no_resp <- factor(data.full$no_resp)
	#is RT faster than 2000 ms?: 0 = no, 1 = yes

str(data.full)

#############################################
######    SPECIFY TYPE OF VARIABLES    ######
#############################################

	###### SELECT RELEVANT BLOCKS ######

table(filter(data.full, ppn == 1)$block)
data.full <- filter(data.full, block %in% c(5:8, 11:14))

	###### SPECIFY BLOCK TYPE ######
	
data.full$block_type <- ifelse(data.full$block %in% c(5,11), "practice", "test")
data.full$block_type <- factor(data.full$block_type)

########################################
######    ACCURACY PER SUBJECT    ######
########################################

	###### DATASET ######

data.acc.check <- dcast(ppn ~ corr, data = filter(data.full, no_resp == 1), value.var = "rt", length)
colMeans(data.acc.check)
names(data.acc.check)[c(2,3)] <- c("err","corr")

	###### PERCENTAGE ERRORS ######

data.acc.check <- mutate(data.acc.check, err.perc = (err/(err+corr))*100)
data.acc.check <- mutate(data.acc.check, corr.perc = (corr/(err+corr))*100)
data.acc.check

	###### PLOT ######

barchart(err.perc ~ factor(ppn), data = data.acc.check)
boxplot(data.acc.check$err.perc)

	###### EXCLUDE PPN WITH ERROR RATE >= 50% ######

ppn.exclude.acc.chance <- print(filter(data.acc.check, err.perc >= 40)$ppn)	
	
	###### EXCLUDE PPN WITH ERROR RATE > 3 SD ######
	
max_sd_acc <- print(mean(filter(data.acc.check, !(ppn %in% ppn.exclude.acc.chance))$err.perc)+3*sd(filter(data.acc.check, !(ppn %in% ppn.exclude.acc.chance))$err.perc))
ppn.exclude.acc <- print(filter(data.acc.check, err.perc >= max_sd_acc)$ppn)
#ppn.exclude.acc <- ppn.exclude.acc.chance

##################################
######    RT PER SUBJECT    ######
##################################

	###### DATASET ######

data.rt.check <- dcast(ppn ~ "rt", data = filter(data.full, !(ppn %in% ppn.exclude.acc.chance), no_resp == 1), value.var = "rt", mean)

	###### PLOT ######

barchart(rt ~ factor(ppn), data = data.rt.check)
boxplot(data.rt.check$rt)

	###### EXCLUDE PPN WITH RT < or > 3 SD ######

max_sd_rt <- print(mean(data.rt.check$rt)+3*sd(data.rt.check$rt))
ppn.exclude.rt <- print(filter(data.rt.check, rt >= max_sd_rt)$ppn)

####################################
######    DESCRIPTIVE DATA    ######
####################################

data.descriptive <- filter(data.full, !(ppn %in% ppn.exclude.acc),!(ppn %in% ppn.exclude.rt))
data.descriptive <- dcast(ppn + age + gender ~ "N", data = data.descriptive, value.var = "rt", length)

print(age_m <- mean(data.descriptive$age, na.rm = T))
print(age_sd <- sd(data.descriptive$age, na.rm = T))
print(gender <- ftable(data.descriptive$gender)) #male = m, female = v

################################
######    OUTLIER DATA    ######
################################

data.outlier <- filter(data.full, !(ppn %in% ppn.exclude.acc),!(ppn %in% ppn.exclude.rt))

	###### NO RESP ######

data.outlier.r <- data.outlier[data.outlier$no_resp == 1, ]
(1 - (nrow(data.outlier.r) / nrow(data.outlier)))*100

###########################
######    RT DATA    ######
###########################

	###### DATA FRAME ######

data.rt <- filter(data.full, !(ppn %in% ppn.exclude.acc),!(ppn %in% ppn.exclude.rt))
data.rt <- select(data.rt, -one_of("corr","plus100","no_resp","max3SD"))
head(data.rt)
length(unique(data.rt$ppn))

	###### CALCULATE D1 OVER BLOCKS ######

RT_mean	<- dcast(ppn ~ con, data = data.rt, value.var = "rt", mean)
RT_sd <- dcast(ppn ~ "RT_sd", data = data.rt, value.var = "rt", sd)
D1 <- data.frame("ppn" = RT_mean$ppn, "D1" = (RT_mean$IC - RT_mean$C) / RT_sd$RT_sd)
mean(D1$D1);sd(D1$D1)

	###### CALCULATE D1 PER BLOCK ######

# RT_mean	<- dcast(ppn ~ block_type + con, data = data.rt, value.var = "rt", mean)
# RT_sd <- dcast(ppn ~ block_type, data = data.rt, value.var = "rt", sd)
# D1 <- data.frame("ppn" = RT_mean$ppn, "D1a" = (RT_mean$practice_IC - RT_mean$practice_C) / RT_sd$practice, "D1b" = (RT_mean$test_IC - RT_mean$test_C) / RT_sd$test)
# D1 <- mutate(D1, D1 = (D1a+D1b)/2)
# mean(D1$D1);sd(D1$D1)

############################
######    ANALYSIS    ######
############################

my_t <- print(t.test(D1$D1))
my_d <- print(ci.sm(sm = my_t$stat / sqrt(39), N = 39))
ttest.tstat(my_t$stat, n1 = 39, simple = T)

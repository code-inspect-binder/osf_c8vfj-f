	############################
	######    PACKAGES    ######
	############################
	
library(car)
library(heplots)
library(lsr)
library(MBESS)
library(BayesFactor)
library(TOSTER)

library(reshape2)
library(plyr)
library(dplyr)

library(lattice)
library(latticeExtra)
library(ggplot2)
library(grid)
library(cowplot)
library(readr)
library(colorblindr) 

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
		file <- sprintf("MYPATH\\dataimi0%d.txt", ppn)
	}else{
		file <- sprintf("MYPATH\\dataimi%d.txt", ppn)
	}

	data <- read.table(file, header = F)
	
	###### NAME VARIABLES ######

	names <- c("trialnr","ppn","age","gender","handpos","group","block","cue","finger","stim","color_left","con","xr","r","corr","rt","re")
	names(data) <- names
	
	###### CONDITION ######

	data$cond <- ifelse(data$group == 1, # color ingroup hand: in_blue (= 1), in_green (= 2)
					# TRUE							
					ifelse(data$color_left == 1, # color left hand: blue (= 1), green (= 2)							
						#TRUE								
						ifelse(data$stim == 1, "in", ifelse(data$stim == 0, "out", "both")), # moving hand: left (= 1), right (= 0), both (= 2) 
						#FALSE
						ifelse(data$stim == 1, "out", ifelse(data$stim == 0, "in", "both"))),								
					# FALSE							
					ifelse(data$color_left == 1, 
						#TRUE
						ifelse(data$stim == 1, "out", ifelse(data$stim == 0, "in", "both")), 
						#FALSE
						ifelse(data$stim == 1, "in", ifelse(data$stim == 0, "out", "both"))))

	###### RT < 100 AND > 2000 MS ######
	
	data$plus100 <- ifelse(data$rt <= 100, 0, 1)
	data$min2000 <- ifelse(data$rt >= 2000, 0, 1)
	
	data.filter <- filter(data, corr == 1, plus100 == 1, min2000 == 1)
	min_rt <- mean(data.filter$rt) - 3*sd(data.filter$rt)
	max_rt <- mean(data.filter$rt) + 3*sd(data.filter$rt)

	###### 3 SD ######
	
	data$max3SD <- ifelse(data$rt >= min_rt, ifelse(data$rt <= max_rt, 1, 0), 0)

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
	
rm(list = setdiff(ls(), c("make_list", "SE_Within", "data.analysis.rt.exp2", "data.analysis.rt.exp3", "data.analysis.rt.dotprobe.exp3", "data.analysis.err.exp2", "data.analysis.err.exp3", "data.analysis.err.dotprobe.exp3","rb_dp_both_exp3")))	

	###### SET WD ######
	
setwd("MYPATH")

	###### RAINDCLOUD PLOT FUNCTIONS ######
	
source("R_rainclouds.R")
source("summarySE.R")

	###### MAKE LIST ######

ppn <- c(1:65)
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
data.full <- select(data.full, -one_of("trialnr","group","cue","finger","stim","color_left","xr","r","re"))
data.full <- select(data.full, ppn:block, cond, con:rt, plus100:max3SD)
head(data.full)
tail(data.full)

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
data.full$cond <- factor(data.full$cond, levels = c("in","out","both"), labels = c("in","out","both"))
	#condition: which hands moved?		
data.full$con <- factor(data.full$con, levels = c(1,0), labels = c("C", "IC"))
	#congruency (i.e., congruent or incongruent)
data.full$corr <- factor(data.full$corr)
	#accuracy: 0 = error, 1 = correct
data.full$plus100 <- factor(data.full$plus100)
	#is RT slower than 100 ms?: 0 = no, 1 = yes
data.full$min2000 <- factor(data.full$min2000)
	#is RT faster than 2000 ms?: 0 = no, 1 = yes
data.full$max3SD <- factor(data.full$max3SD)
	#is RT within the 3 SD bounds?: 0 = no, 1 = yes

str(data.full)
summary(data.full)

########################################
######    ACCURACY PER SUBJECT    ######
########################################

	###### DATASET ######

data.acc.check <- dcast(ppn ~ corr, data = filter(data.full, plus100 == 1, min2000 == 1, max3SD == 1), value.var = "rt", length)
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

ppn.exclude.acc.chance <- print(filter(data.acc.check, err.perc >= 50)$ppn)	
	
	###### EXCLUDE PPN WITH ERROR RATE > 3 SD ######
	
max_sd_acc <- print(mean(filter(data.acc.check, !(ppn %in% ppn.exclude.acc.chance))$err.perc)+3*sd(filter(data.acc.check, !(ppn %in% ppn.exclude.acc.chance))$err.perc))
ppn.exclude.acc <- print(filter(data.acc.check, err.perc >= max_sd_acc)$ppn)

##################################
######    RT PER SUBJECT    ######
##################################

	###### DATASET ######

data.rt.check <- dcast(ppn ~ "rt", data = filter(data.full, !(ppn %in% ppn.exclude.acc.chance), plus100 == 1, min2000 == 1, max3SD == 1, corr == 1), value.var = "rt", mean)

	###### PLOT ######

barchart(rt ~ factor(ppn), data = data.rt.check)
boxplot(data.rt.check$rt)

	###### EXCLUDE PPN WITH RT < or > 3 SD ######

max_sd_rt <- print(mean(data.rt.check$rt)+3*sd(data.rt.check$rt))
ppn.exclude.rt <- print(filter(data.rt.check, rt >= max_sd_rt)$ppn)

####################################
######    DESCRIPTIVE DATA    ######
####################################

data.descriptive <- data.full
data.descriptive <- filter(data.full, !(ppn %in% ppn.exclude.acc),!(ppn %in% ppn.exclude.rt))
data.descriptive <- dcast(ppn + age + gender ~ "N", data = data.descriptive, value.var = "rt", length)

print(age_m <- mean(data.descriptive$age, na.rm = T))
print(age_sd <- sd(data.descriptive$age, na.rm = T))
print(age_range <- range(data.descriptive$age, na.rm = T))
print(gender <- ftable(data.descriptive$gender)) #male = m, female = v

################################
######    OUTLIER DATA    ######
################################

data.outlier <- filter(data.full, !(ppn %in% ppn.exclude.acc),!(ppn %in% ppn.exclude.rt))

	###### TABLE ######
	
ftable(data.outlier$plus100, data.outlier$min2000, data.outlier$corr, data.outlier$max3SD, col.vars = c(1,2), dnn = c("plus100", "min2000", "corr", "max3SD"))

	###### NO RESP ######

data.outlier.r <- data.outlier[data.outlier$min2000 == 1, ]
(1 - (nrow(data.outlier.r) / nrow(data.outlier)))*100

	###### 100 MS ######

data.outlier.100 <- data.outlier.r[data.outlier.r$plus100 == 1, ]
(1 - (nrow(data.outlier.100) / nrow(data.outlier.r)))*100

	###### CORR ######

data.outlier.corr <- data.outlier.100[data.outlier.100$corr == 1, ]
(1 - (nrow(data.outlier.corr) / nrow(data.outlier.100)))*100

	###### 3 SD ######

data.outlier.sd <- data.outlier.corr[data.outlier.corr$max3SD == 1, ]
(1 - (nrow(data.outlier.sd) / nrow(data.outlier.corr)))*100

###########################
######    RT DATA    ######
###########################

data.rt <- filter(data.full, !(ppn %in% ppn.exclude.acc),!(ppn %in% ppn.exclude.rt), plus100 == 1, min2000 == 1, max3SD == 1, corr == 1)
data.rt <- select(data.rt, -one_of("corr","plus100","min2000","max3SD"))
head(data.rt)
length(unique(data.rt$ppn))

##############################
######    ERROR DATA    ######
##############################

data.corr <- filter(data.full,  !(ppn %in% ppn.exclude.acc),!(ppn %in% ppn.exclude.rt), plus100 == 1, min2000 == 1)
data.corr <- select(data.corr, -one_of("plus100","min2000","max3SD"))
head(data.corr)
length(unique(data.corr$ppn))

	###############################
	######    RT ANALYSIS    ######
	###############################

###########################################
######    VISUALISATION BAR PLOTS    ######
###########################################

	###### ERROR BAR FUNCTION ######

panel.err <- function(x, y, subscripts, groups, stderr, box.ratio, ...)
{
  d <- 1/(nlevels(groups)+nlevels(groups)/box.ratio)
  g <- (as.numeric(groups[subscripts])-1); g <- (g-median(g))*d
  panel.arrows(as.numeric(x)+g,y-stderr[subscripts], as.numeric(x)+g, y+stderr[subscripts],code=3,angle=90, length=0.3, unit = "cm")
}

	###### CALCULATE ERROR BARS ######

data.plot.rt.se <- dcast(ppn ~ cond + con, data = data.rt, value.var = "rt", mean)
SE_rt <- SE_Within(data = data.plot.rt.se, colnames = names(data.plot.rt.se)[-1])$se.within.corr

	###### MAKE PLOT ######

data.plot.rt.a <- dcast(ppn + cond + con ~ "RT", data = data.rt, value.var = "rt", mean)
data.plot.rt <- dcast(cond + con ~ "RT", data = data.plot.rt.a, value.var = "RT", mean)

mytheme = trellis.par.get()
mytheme$superpose.polygon$col = c("white", "grey")
mytheme$superpose.polygon$lwd = 2
mytheme$add.line$lwd = 2
mytheme$fontsize$text = 18
mytheme$axis.text$cex = 1
mytheme$axis.line$lwd = 2

plot_rt <- barchart(RT ~ cond, groups = con, data = data.plot.rt, box.ratio = 3,
	ylim = c(512,578), scales = list(tck = c(1,0), x = list(labels = c("In-Group","Out-Group","In+Out-Group")), y = list(at = round(seq(515, 575, 15)))),
	ylab = "Reaction Time (ms)",
	#auto.key = list(corner = c(0.975,.975), points = F, rectangles = T, height = 0.85, size = 3.5, padding.text = 2.5),
    auto.key = list(corner = c(.5,.975), text = c("Congruent","Incongruent"), points = F, rectangles = T, size = 3.5, columns = 2, between.columns = 0.75, padding.text = 2),
	par.settings = mytheme,
	panel=function(x, y, subscripts, groups, stderr, box.ratio, ...){
		panel.barchart(x, y, subscripts = subscripts, groups = groups, box.ratio = box.ratio, origin = 0,...)
		panel.err(x, y, subscripts = subscripts, groups = groups, box.ratio = box.ratio, stderr = SE_rt)
	})
print(plot_rt)

	###### CONGRUENCY EFFECT ######
	
data.plot.rt.con <- data.frame(cond = c("in","out","both"), CE = NA)
data.plot.rt.con$cond <- factor(data.plot.rt.con$cond, levels = c("in","out","both"))

for(i in 1:nrow(data.plot.rt.con)){
	data.plot.rt.con$CE[i] <- data.plot.rt$RT[i*2] - data.plot.rt$RT[(i*2)-1]
}

barchart(CE ~ cond, data = data.plot.rt.con, col = "grey",par.settings = mytheme, ylim = c(-2,37), scales = list(tck = c(1,0), y = list(at = round(seq(0, 35, 5)))))

#################################################
######    VISUALISATION RAINCLOUD PLOTS    ######
#################################################
	
data.plot.rt.ppn <- dcast(ppn + cond + con ~ "RT", data = data.rt, value.var = "rt", mean)
data.plot.rt.sum <- summarySE(data.plot.rt.ppn, measurevar = "RT", groupvars=c("cond", "con"))
  
rb_rt <- ggplot(data.plot.rt.ppn, aes(x = cond, y = RT, fill = con)) +
  geom_flat_violin(aes(fill = con),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(cond)-.15, y = RT, colour = con),position = position_jitter(width = .05), size = .25, shape = 20, show.legend = F)+
  geom_boxplot(aes(x = cond, y = RT, fill = con),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  geom_line(data = data.plot.rt.sum, aes(x = as.numeric(cond)+.1, y = RT_mean, group = con, colour = con), linetype = 3, show.legend = F)+
  geom_point(data = data.plot.rt.sum, aes(x = as.numeric(cond)+.1, y = RT_mean, group = con, colour = con), shape = 18, show.legend = F) +
  geom_errorbar(data = data.plot.rt.sum, aes(x = as.numeric(cond)+.1, y = RT_mean, group = con, colour = con, ymin = RT_mean-se, ymax = RT_mean+se), width = .05, show.legend = F)+
  scale_x_discrete(name = NULL, labels = c("In-Group","Out-Group","In+Out-Group")) + 
  scale_y_continuous(name = "Reaction Time (ms)", limits = c(350,850), breaks = seq(350,850,100))+#,expand = expand_scale(mult = c(0, .05))) + 
  scale_color_OkabeIto(order = c(5,6))+
  scale_fill_OkabeIto(order = c(5,6))+theme_cowplot()+#theme(plot.title = element_text(hjust = 0.5))+
  guides(fill= guide_legend(title=NULL))#+
  #ggtitle("Reaction Times")
rb_rt 

# ggsave('rb_rt.png', rb_rt, width = 6, height = 3)

############################
######    ANALYSIS    ######
############################

	###### ANOVA ######
	
data.analysis.rt <- dcast(ppn ~ cond + con, data = data.rt, value.var = "rt", mean)
data.analysis.rt.exp3 <- data.analysis.rt
head(data.analysis.rt)
tail(data.analysis.rt)	

cond <- factor(rep(1:3, each = 2))
con <- factor(rep(1:2, times = 3))
idata <- data.frame(cond, con)

fit.rt <- lm(cbind(in_C,in_IC,out_C,out_IC,both_C,both_IC) ~ 1, data = data.analysis.rt)
etasq(Anova(fit.rt, type = "III", test = "Wilks", idata = idata, idesign = ~cond*con), anova = T, partial = T)
etasq(Anova(fit.rt, type = "III", test = "Wilks", idata = idata, idesign = ~cond*con), anova = T, partial = T)$"approx F"
	
	###### CON X COND ######
	
###### ES ######

Lims <- conf.limits.ncf(F.value = 12.88, conf.level = 0.90, df.1 <- 2, df.2 <- 61)
Lims$Lower.Limit/(Lims$Lower.Limit + df.1 + df.2 + 1)
Lims$Upper.Limit/(Lims$Upper.Limit + df.1 + df.2 + 1)
	
###### IN vs. OUT ######

my_t_inout <- print(t.test(data.analysis.rt$in_IC - data.analysis.rt$in_C, data.analysis.rt$out_IC - data.analysis.rt$out_C, paired = T))
my_d_inout <- print(ci.sm(sm = cohensD(data.analysis.rt$in_IC - data.analysis.rt$in_C, data.analysis.rt$out_IC - data.analysis.rt$out_C, method = "paired"), N = 63))

ttest.tstat(my_t_inout$stat, n1 = 63, simple = T)
TOSTone(m = my_d_inout$Standardized.Mean, mu = 0, sd = 1, n = 63, low_eqbound_d = -0.3, high_eqbound_d = 0.3)

###### IN vs. BOTH ######

my_t_inboth <- print(t.test(data.analysis.rt$in_IC - data.analysis.rt$in_C, data.analysis.rt$both_IC - data.analysis.rt$both_C, paired = T))
my_d_inboth <- print(ci.sm(sm = cohensD(data.analysis.rt$in_IC - data.analysis.rt$in_C, data.analysis.rt$both_IC - data.analysis.rt$both_C, method = "paired"), N = 63))

ttest.tstat(my_t_inboth$stat, n1 = 63, simple = T)
TOSTone(m = my_d_inboth$Standardized.Mean, mu = 0, sd = 1, n = 63, low_eqbound_d = -0.3, high_eqbound_d = 0.3)

###### OUT vs. BOTH ######

my_t_outboth <- print(t.test(data.analysis.rt$out_IC - data.analysis.rt$out_C, data.analysis.rt$both_IC - data.analysis.rt$both_C, paired = T))
my_d_outboth <- print(ci.sm(sm = cohensD(data.analysis.rt$out_IC - data.analysis.rt$out_C, data.analysis.rt$both_IC - data.analysis.rt$both_C, method = "paired"), N = 63))

ttest.tstat(my_t_outboth$stat, n1 = 63, simple = T)
TOSTone(m = my_d_outboth$Standardized.Mean, mu = 0, sd = 1, n = 63, low_eqbound_d = -0.3, high_eqbound_d = 0.3)

	###### MAIN EFFECT CON ######

data.analysis.rt.con.a <- dcast(ppn + cond + con ~ "rt", data = data.rt, value.var = "rt", mean)
data.analysis.rt.con <- dcast(ppn  ~ con, data = data.analysis.rt.con.a, value.var = "rt", mean)

t.test(data.analysis.rt.con$IC, data.analysis.rt.con$C, paired = T)
ci.sm(sm = cohensD(data.analysis.rt.con$IC, data.analysis.rt.con$C, method = "paired"), N = 63)

Lims <- conf.limits.ncf(F.value = 46.62, conf.level = 0.90, df.1 <- 1, df.2 <- 62)
Lims$Lower.Limit/(Lims$Lower.Limit + df.1 + df.2 + 1)
Lims$Upper.Limit/(Lims$Upper.Limit + df.1 + df.2 + 1)

	###### MAIN EFFECT CONDITION ######

data.analysis.rt.cond <- dcast(ppn ~ cond, data = data.analysis.rt.con.a, value.var = "rt", mean)
colMeans(data.analysis.rt.cond)

Lims <- conf.limits.ncf(F.value = 0.52, conf.level = 0.90, df.1 <- 2, df.2 <- 61)
Lims$Lower.Limit/(Lims$Lower.Limit + df.1 + df.2 + 1)
Lims$Upper.Limit/(Lims$Upper.Limit + df.1 + df.2 + 1)

	###### ADD DOT-PROBE ######

###### DATA ######	

data.analysis.rt.dotprobe.exp3$dp_CE <- data.analysis.rt.dotprobe.exp3$IC - data.analysis.rt.dotprobe.exp3$C	
data.analysis.rt.plusdp <- join(data.analysis.rt, select(data.analysis.rt.dotprobe.exp3, ppn, dp_CE), type = "inner", by = "ppn")
head(data.analysis.rt.plusdp)

###### ANOVA ######

fit.rt <- lm(cbind(in_C,in_IC,out_C,out_IC,both_C,both_IC) ~ dp_CE, data = data.analysis.rt.plusdp)
etasq(Anova(fit.rt, type = "III", test = "Wilks", idata = idata, idesign = ~cond*con), anova = T, partial = T)
etasq(Anova(fit.rt, type = "III", test = "Wilks", idata = idata, idesign = ~cond*con), anova = T, partial = T)$"approx F"

Lims <- conf.limits.ncf(F.value = 2.82, conf.level = 0.90, df.1 <- 2, df.2 <- 56)
Lims$Lower.Limit/(Lims$Lower.Limit + df.1 + df.2 + 1)
Lims$Upper.Limit/(Lims$Upper.Limit + df.1 + df.2 + 1)

###### CORRELATION CE IN vs CE OUT ######

data.analysis.rt.plusdp$inout_CE <- (data.analysis.rt.plusdp$in_IC - data.analysis.rt.plusdp$in_C) - (data.analysis.rt.plusdp$out_IC - data.analysis.rt.plusdp$out_C)
cor.test(data.analysis.rt.plusdp$dp_CE, data.analysis.rt.plusdp$inout_CE)
plot(data.analysis.rt.plusdp$dp_CE, data.analysis.rt.plusdp$inout_CE)

###### CORRELATION CE IN vs CE IN/OUT ######

data.analysis.rt.plusdp$inboth_CE <- (data.analysis.rt.plusdp$both_IC - data.analysis.rt.plusdp$both_C) - (data.analysis.rt.plusdp$in_IC - data.analysis.rt.plusdp$in_C)

cor.test(data.analysis.rt.plusdp$dp_CE, data.analysis.rt.plusdp$inboth_CE)
plot(data.analysis.rt.plusdp$dp_CE, data.analysis.rt.plusdp$inboth_CE)

###### CORRELATION CE OUT vs CE IN/OUT ######

data.analysis.rt.plusdp$outboth_CE <- (data.analysis.rt.plusdp$both_IC - data.analysis.rt.plusdp$both_C) - (data.analysis.rt.plusdp$out_IC - data.analysis.rt.plusdp$out_C)

cor.test(data.analysis.rt.plusdp$dp_CE, data.analysis.rt.plusdp$outboth_CE)
plot(data.analysis.rt.plusdp$dp_CE, data.analysis.rt.plusdp$outboth_CE)

	###### ADD MOTIVATION AFFILIATE QUESTIONNAIRE ######

###### DATA ######	
	
motaffq <- read.csv2("MYPATH\\MotAff.csv",header = T)
motaffq <- filter(motaffq, ppn != 2)
data.analysis.rt.plusq <- join(data.analysis.rt, motaffq, type = "inner", by = "ppn")

###### AVERAGE ######

mean(motaffq$MotAff); sd(motaffq$MotAff)

###### ANOVA ######

fit.rt <- lm(cbind(in_C,in_IC,out_C,out_IC,both_C,both_IC) ~ MotAff, data = data.analysis.rt.plusq)
etasq(Anova(fit.rt, type = "III", test = "Wilks", idata = idata, idesign = ~cond*con), anova = T, partial = T)
etasq(Anova(fit.rt, type = "III", test = "Wilks", idata = idata, idesign = ~cond*con), anova = T, partial = T)$"approx F"

Lims <- conf.limits.ncf(F.value = 0.19, conf.level = 0.90, df.1 <- 2, df.2 <- 59)
Lims$Lower.Limit/(Lims$Lower.Limit + df.1 + df.2 + 1)
Lims$Upper.Limit/(Lims$Upper.Limit + df.1 + df.2 + 1)

###### CORRELATION CE IN vs CE OUT ######

data.analysis.rt.plusq$inout_CE <- (data.analysis.rt.plusq$in_IC - data.analysis.rt.plusq$in_C) - (data.analysis.rt.plusq$out_IC - data.analysis.rt.plusq$out_C)

cor.test(data.analysis.rt.plusq$MotAff, data.analysis.rt.plusq$inout_CE)
plot(data.analysis.rt.plusq$MotAff, data.analysis.rt.plusq$inout_CE)

###### CORRELATION CE IN vs CE IN/OUT ######

data.analysis.rt.plusq$inboth_CE <- (data.analysis.rt.plusq$both_IC - data.analysis.rt.plusq$both_C) - (data.analysis.rt.plusq$in_IC - data.analysis.rt.plusq$in_C)

cor.test(data.analysis.rt.plusq$MotAff, data.analysis.rt.plusq$inboth_CE)
plot(data.analysis.rt.plusq$MotAff, data.analysis.rt.plusq$inboth_CE)

###### CORRELATION CE OUT vs CE IN/OUT ######

data.analysis.rt.plusq$outboth_CE <- (data.analysis.rt.plusq$both_IC - data.analysis.rt.plusq$both_C) - (data.analysis.rt.plusq$out_IC - data.analysis.rt.plusq$out_C)

cor.test(data.analysis.rt.plusq$MotAff, data.analysis.rt.plusq$outboth_CE)
plot(data.analysis.rt.plusq$MotAff, data.analysis.rt.plusq$outboth_CE)

	################################
	######    ERR ANALYSIS    ######
	################################
				
#########################################
######    CALCULATE ERROR RATES    ######
#########################################

data.err <- dcast(ppn + cond + con ~ corr, data = data.corr, value.var = "rt", length)
names(data.err)[4:5] <- c("n_err","n_corr")
data.err <- mutate(data.err, err = (n_err/(n_err+n_corr))*100)
data.err <- select(data.err, -one_of("n_err", "n_corr"))

head(data.err)

###########################################
######    VISUALISATION BAR PLOTS    ######
###########################################

	###### ERROR BAR FUNCTION ######

panel.err <- function(x, y, subscripts, groups, stderr, box.ratio, ...)
{
  d <- 1/(nlevels(groups)+nlevels(groups)/box.ratio)
  g <- (as.numeric(groups[subscripts])-1); g <- (g-median(g))*d
  panel.arrows(as.numeric(x)+g,y-stderr[subscripts], as.numeric(x)+g, y+stderr[subscripts],code=3,angle=90, length=0.3, unit = "cm")
}

	###### CALCULATE ERROR BARS ######

data.plot.err.se <- dcast(ppn ~ cond + con, data = data.err, value.var = "err", mean)
SE_err <- SE_Within(data = data.plot.err.se, colnames = names(data.plot.err.se)[-1])$se.within.corr

	###### MAKE PLOT ######

data.plot.err.a <- dcast(ppn + cond + con ~ "ErrorRate", data = data.err, value.var = "err", mean)
data.plot.err <- dcast(cond + con ~ "ErrorRate", data = data.plot.err.a, value.var = "ErrorRate", mean)

mytheme = trellis.par.get()
mytheme$superpose.polygon$col = c("white", "grey")
mytheme$superpose.polygon$lwd = 2
mytheme$add.line$lwd = 2
mytheme$fontsize$text = 18
mytheme$axis.text$cex = 1
mytheme$axis.line$lwd = 2

plot_err <- barchart(ErrorRate ~ cond, groups = con, data = data.plot.err, box.ratio = 3,
	ylim = c(0.7,7.3), scales = list(tck = c(1,0), x = list(labels = c("In-Group","Out-Group","In+Out-Group")), y = list(at = seq(1, 7, 1.5))),
	ylab = "Error Rate (%)",
    auto.key = list(corner = c(.5,.975), text = c("Congruent","Incongruent"), points = F, rectangles = T, size = 3.5, columns = 2, between.columns = 0.75, padding.text = 2),
	par.settings = mytheme,
	panel=function(x, y, subscripts, groups, stderr, box.ratio, ...){
		panel.barchart(x, y, subscripts = subscripts, groups = groups, box.ratio = box.ratio, ...)
		panel.err(x, y, subscripts = subscripts, groups = groups, box.ratio = box.ratio, stderr = SE_err)
	})
print(plot_err)

	###### PRINT RT AND ER PLOTS ######
	
# ppi <- 300
# tiff(filename = "MYPATH\\imi_exp3.tiff", compression = "zip", width = 7*2*ppi, height = 7*1*ppi, pointsize = 7*(9/5), res = ppi)

print(plot_rt, split = c(1, 1, 2, 1), more = T)
print(plot_err, split = c(2, 1, 2, 1), more = F)

# dev.off()

	###### CONGRUENCY EFFECT ######
	
data.plot.err.con <- data.frame(cond = c("in","out","both"), CE = NA)
data.plot.err.con$cond <- factor(data.plot.err.con$cond, levels = c("in","out","both"))

for(i in 1:nrow(data.plot.err.con)){
	data.plot.err.con$CE[i] <- data.plot.err$ErrorRate[i*2] - data.plot.err$ErrorRate[(i*2)-1]
}

barchart(CE ~ cond, data = data.plot.err.con, col = "grey",par.settings = mytheme, ylim = c(-0.5,3.5), scales = list(tck = c(1,0), y = list(at = seq(0, 3, 1))))

#################################################
######    VISUALISATION RAINCLOUD PLOTS    ######
#################################################
	
data.plot.err.ppn <- dcast(ppn + cond + con ~ "ErrorRate", data = data.err, value.var = "err", mean)
data.plot.err.sum <- summarySE(data.plot.err.ppn, measurevar = "ErrorRate", groupvars=c("cond", "con"))
  
rb_err <- ggplot(data.plot.err.ppn, aes(x = cond, y = ErrorRate, fill = con)) +
  geom_flat_violin(aes(fill = con),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(cond)-.15, y = ErrorRate, colour = con),position = position_jitter(width = .05), size = .25, shape = 20, show.legend = F)+
  geom_boxplot(aes(x = cond, y = ErrorRate, fill = con),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  geom_line(data = data.plot.err.sum, aes(x = as.numeric(cond)+.1, y = ErrorRate_mean, group = con, colour = con), linetype = 3, show.legend = F)+
  geom_point(data = data.plot.err.sum, aes(x = as.numeric(cond)+.1, y = ErrorRate_mean, group = con, colour = con), shape = 18, show.legend = F) +
  geom_errorbar(data = data.plot.err.sum, aes(x = as.numeric(cond)+.1, y = ErrorRate_mean, group = con, colour = con, ymin = ErrorRate_mean-se, ymax = ErrorRate_mean+se), width = .05, show.legend = F)+
  scale_x_discrete(name = NULL, labels = c("In-Group","Out-Group","In+Out-Group")) + 
  scale_y_continuous(name = "Error Rate (%)", limits = c(-5,20), breaks = seq(-5,20,5))+#,expand = expand_scale(mult = c(0, .05))) + 
  scale_color_OkabeIto(order = c(5,6))+
  scale_fill_OkabeIto(order = c(5,6))+theme_cowplot()+#theme(plot.title = element_text(hjust = 0.5))+
  guides(fill= guide_legend(title=NULL))#+
  #ggtitle("Error Rates")
rb_err

# ggsave('rb_err.png', rb_err, width = 6, height = 3)
 
	###### PRINT RT AND ER PLOTS ######
	
rb_both_plot <- plot_grid(rb_rt, rb_err, nrow = 2)
rb_both_title <- ggdraw() + draw_label("Automatic Imitation", fontface = 'bold')
rb_both <- plot_grid(rb_both_title, rb_both_plot, ncol = 1, rel_heights = c(0.1, 1))

# ggsave('rb_rt+err.png', rb_both, width = 6, height = 6)

	###### PRINT AI AND DP PLOTS ######
	
rb_all <- plot_grid(rb_both, rb_dp_both_exp3, ncol = 1, nrow = 2, rel_heights = c(2,1))
# ggsave('rb_ai+dp.png', rb_all, width = 6, height = 9)

############################
######    ANALYSIS    ######
############################

data.analysis.err <- dcast(ppn ~ cond + con, data = data.err, value.var = "err", mean)
data.analysis.err.exp3 <- data.analysis.err
head(data.analysis.err)
tail(data.analysis.err)	

cond <- factor(rep(1:3, each = 2))
con <- factor(rep(1:2, times = 3))
idata <- data.frame(cond, con)

fit.err <- lm(cbind(in_C,in_IC,out_C,out_IC,both_C,both_IC) ~ 1, data = data.analysis.err)
etasq(Anova(fit.err, type = "III", test = "Wilks", idata = idata, idesign = ~cond*con), anova = T, partial = T)
etasq(Anova(fit.err, type = "III", test = "Wilks", idata = idata, idesign = ~cond*con), anova = T, partial = T)$"approx F"

	###### CON X COND ######
	
###### ES ######

Lims <- conf.limits.ncf(F.value = 4.10, conf.level = 0.90, df.1 <- 2, df.2 <- 61)
Lims$Lower.Limit/(Lims$Lower.Limit + df.1 + df.2 + 1)
Lims$Upper.Limit/(Lims$Upper.Limit + df.1 + df.2 + 1)	
	
###### IN vs. OUT ######

my_t_inout <- print(t.test(data.analysis.err$in_IC - data.analysis.err$in_C, data.analysis.err$out_IC - data.analysis.err$out_C, paired = T))
my_d_inout <- print(ci.sm(sm = cohensD(data.analysis.err$in_IC - data.analysis.err$in_C, data.analysis.err$out_IC - data.analysis.err$out_C, method = "paired"), N = 63))

ttest.tstat(my_t_inout$stat, n1 = 63, simple = T)
TOSTone(m = my_d_inout$Standardized.Mean, mu = 0, sd = 1, n = 63, low_eqbound_d = -0.3, high_eqbound_d = 0.3)

###### IN vs. BOTH ######

my_t_inboth <- print(t.test(data.analysis.err$in_IC - data.analysis.err$in_C, data.analysis.err$both_IC - data.analysis.err$both_C, paired = T))
my_d_inboth <- print(ci.sm(sm = cohensD(data.analysis.err$in_IC - data.analysis.err$in_C, data.analysis.err$both_IC - data.analysis.err$both_C, method = "paired"), N = 63))

ttest.tstat(my_t_inboth$stat, n1 = 63, simple = T)
TOSTone(m = my_d_inboth$Standardized.Mean, mu = 0, sd = 1, n = 63, low_eqbound_d = -0.3, high_eqbound_d = 0.3)

###### OUT vs. BOTH ######

my_t_outboth <- print(t.test(data.analysis.err$out_IC - data.analysis.err$out_C, data.analysis.err$both_IC - data.analysis.err$both_C, paired = T))
my_d_outboth <- print(ci.sm(sm = cohensD(data.analysis.err$out_IC - data.analysis.err$out_C, data.analysis.err$both_IC - data.analysis.err$both_C, method = "paired"), N = 63))

ttest.tstat(my_t_outboth$stat, n1 = 63, simple = T)
TOSTone(m = my_d_outboth$Standardized.Mean, mu = 0, sd = 1, n = 63, low_eqbound_d = -0.3, high_eqbound_d = 0.3)

	###### MAIN EFFECT CON ######

data.analysis.err.con.a <- dcast(ppn + cond + con ~ "err", data = data.err, value.var = "err", mean)
data.analysis.err.con <- dcast(ppn  ~ con, data = data.analysis.err.con.a, value.var = "err", mean)

t.test(data.analysis.err.con$IC, data.analysis.err.con$C, paired = T)
ci.sm(sm = cohensD(data.analysis.err.con$IC, data.analysis.err.con$C, method = "paired"), N = 63)

Lims <- conf.limits.ncf(F.value = 35.14, conf.level = 0.90, df.1 <- 1, df.2 <- 62)
Lims$Lower.Limit/(Lims$Lower.Limit + df.1 + df.2 + 1)
Lims$Upper.Limit/(Lims$Upper.Limit + df.1 + df.2 + 1)	

	###### MAIN EFFECT CONDITION NUMBER ######

data.analysis.err.cond <- dcast(ppn ~ cond, data = data.analysis.err.con.a, value.var = "err", mean)
colMeans(data.analysis.err.cond)

Lims <- conf.limits.ncf(F.value = 7.00, conf.level = 0.90, df.1 <- 2, df.2 <- 61)
Lims$Lower.Limit/(Lims$Lower.Limit + df.1 + df.2 + 1)
Lims$Upper.Limit/(Lims$Upper.Limit + df.1 + df.2 + 1)	

	#####################################
	######    COMBINED ANALYSIS    ######
	#####################################
		
###############################
######    RT ANALYSIS    ######
###############################

	###### DATA ######
	
data.analysis.rt.exp2$exp <- factor(2)
data.analysis.rt.exp3$exp <- factor(3)
data.analysis.rt.combo <- rbind(data.analysis.rt.exp2, data.analysis.rt.exp3)
data.analysis.rt.combo$ppn <- 1:nrow(data.analysis.rt.combo)

	###### VISUALIZE ######
	
###### CALCULATE ERROR BARS ######

SE_rt_combo <- SE_Within(data = data.analysis.rt.combo, colnames = names(data.analysis.rt.combo)[-c(1,8)])$se.within.corr

###### MAKE PLOT ######
head(data.plot.rt.combo)
data.plot.rt.combo <- data.plot.rt
data.plot.rt.combo$RT <- as.numeric(colMeans(data.analysis.rt.combo[-c(1,8)]))

mytheme = trellis.par.get()
mytheme$superpose.polygon$col = c("white", "grey")
mytheme$fontsize$text = 22
mytheme$axis.text$cex = 1

plot_rt_combo <- barchart(RT ~ cond, groups = con, data = data.plot.rt.combo, box.ratio = 3,
	ylim = c(515,575), scales = list(tck = c(1,0), y = list(at = round(seq(520, 570, 10)))),
	xlab = "Condition", ylab = "RT",
	auto.key = list(corner = c(0.975,.975), points = F, rectangles = T, height = 0.85, size = 3.5, padding.text = 2.5),
	par.settings = mytheme,
	panel=function(x, y, subscripts, groups, stderr, box.ratio, ...){
		panel.barchart(x, y, subscripts = subscripts, groups = groups, box.ratio = box.ratio, origin = 0,...)
		panel.err(x, y, subscripts = subscripts, groups = groups, box.ratio = box.ratio, stderr = SE_rt_combo)
		grid.text("(a)", x = unit(.05, "npc"), y = unit(.95, "npc"))
	})
print(plot_rt_combo)

###### CONGRUENCY EFFECT ######

data.plot.rt.combo.con <- data.frame(cond = c("in","out","both"), CE = NA)
data.plot.rt.combo.con$cond <- factor(data.plot.rt.combo.con$cond, levels = c("in","out","both"))

for(i in 1:nrow(data.plot.rt.combo.con)){
	data.plot.rt.combo.con$CE[i] <- data.plot.rt.combo$RT[i*2] - data.plot.rt.combo$RT[(i*2)-1]
}

barchart(CE ~ cond, data = data.plot.rt.combo.con, col = "grey",par.settings = mytheme, ylim = c(-2,37), scales = list(tck = c(1,0), y = list(at = round(seq(0, 35, 5)))))

###### RAINCLOUD PLOT ######

data.plot.rt.combo.ppn <- melt(id.vars = c("ppn","exp"), value.name = "RT", data = data.analysis.rt.combo)
data.plot.rt.combo.ppn <- cbind(data.plot.rt.combo.ppn, colsplit(data.plot.rt.combo.ppn$variable, "_", c("cond","con"))) 
data.plot.rt.combo.ppn <- select(data.plot.rt.combo.ppn, ppn, cond, con, RT)
data.plot.rt.combo.ppn$cond <- factor(data.plot.rt.combo.ppn$cond, levels = c("in", "out", "both"))
data.plot.rt.combo.ppn$con <- factor(data.plot.rt.combo.ppn$con, levels = c("C", "IC"))

data.plot.rt.combo.sum <- summarySE(data.plot.rt.combo.ppn, measurevar = "RT", groupvars=c("cond", "con"))

rb_rt_combo <- ggplot(data.plot.rt.combo.ppn, aes(x = cond, y = RT, fill = con)) +
  geom_flat_violin(aes(fill = con),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(cond)-.15, y = RT, colour = con),position = position_jitter(width = .05), size = .25, shape = 20, show.legend = F)+
  geom_boxplot(aes(x = cond, y = RT, fill = con),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  geom_line(data = data.plot.rt.combo.sum, aes(x = as.numeric(cond)+.1, y = RT_mean, group = con, colour = con), linetype = 3, show.legend = F)+
  geom_point(data = data.plot.rt.combo.sum, aes(x = as.numeric(cond)+.1, y = RT_mean, group = con, colour = con), shape = 18, show.legend = F) +
  geom_errorbar(data = data.plot.rt.combo.sum, aes(x = as.numeric(cond)+.1, y = RT_mean, group = con, colour = con, ymin = RT_mean-se, ymax = RT_mean+se), width = .05, show.legend = F)+
  scale_x_discrete(name = NULL, labels = c("In-Group","Out-Group","In+Out-Group")) + 
  scale_y_continuous(name = "Reaction Time (ms)", limits = c(350,850), breaks = seq(350,850,100))+#,expand = expand_scale(mult = c(0, .05))) + 
  scale_color_OkabeIto(order = c(5,6))+
  scale_fill_OkabeIto(order = c(5,6))+theme_cowplot()+#theme(plot.title = element_text(hjust = 0.5))+
  guides(fill= guide_legend(title=NULL))#+
  #ggtitle("Reaction Times")
rb_rt_combo

# ggsave('rb_rt_combo.png', rb_rt_combo, width = 6, height = 3)

	###### ANOVA ######

cond <- factor(rep(1:3, each = 2))
con <- factor(rep(1:2, times = 3))
idata <- data.frame(cond, con)

fit.rt.combo <- lm(cbind(in_C,in_IC,out_C,out_IC,both_C,both_IC) ~ exp, data = data.analysis.rt.combo)
etasq(Anova(fit.rt.combo, type = "III", test = "Wilks", idata = idata, idesign = ~cond*con), anova = T, partial = T)
etasq(Anova(fit.rt.combo, type = "III", test = "Wilks", idata = idata, idesign = ~cond*con), anova = T, partial = T)$"approx F"

	###### CON X COND ######
		
###### ES ######

Lims <- conf.limits.ncf(F.value = 12.50, conf.level = 0.90, df.1 <- 2, df.2 <- 99)
Lims$Lower.Limit/(Lims$Lower.Limit + df.1 + df.2 + 1)
Lims$Upper.Limit/(Lims$Upper.Limit + df.1 + df.2 + 1)
	
###### IN vs. OUT ######

my_t_inout <- print(t.test(data.analysis.rt.combo$in_IC - data.analysis.rt.combo$in_C, data.analysis.rt.combo$out_IC - data.analysis.rt.combo$out_C, paired = T))
my_d_inout <- print(ci.sm(sm = cohensD(data.analysis.rt.combo$in_IC - data.analysis.rt.combo$in_C, data.analysis.rt.combo$out_IC - data.analysis.rt.combo$out_C, method = "paired"), N = 102))

ttest.tstat(my_t_inout$stat, n1 = 102, simple = T)
TOSTone(m = my_d_inout$Standardized.Mean, mu = 0, sd = 1, n = 102, low_eqbound_d = -0.3, high_eqbound_d = 0.3)

###### IN vs. BOTH ######

my_t_inboth <- print(t.test(data.analysis.rt.combo$in_IC - data.analysis.rt.combo$in_C, data.analysis.rt.combo$both_IC - data.analysis.rt.combo$both_C, paired = T))
my_d_inboth <- print(ci.sm(sm = cohensD(data.analysis.rt.combo$in_IC - data.analysis.rt.combo$in_C, data.analysis.rt.combo$both_IC - data.analysis.rt.combo$both_C, method = "paired"), N = 102))

ttest.tstat(my_t_inboth$stat, n1 = 102, simple = T)
TOSTone(m = my_d_inboth$Standardized.Mean, mu = 0, sd = 1, n = 102, low_eqbound_d = -0.3, high_eqbound_d = 0.3)

###### OUT vs. BOTH ######

my_t_outboth <- print(t.test(data.analysis.rt.combo$out_IC - data.analysis.rt.combo$out_C, data.analysis.rt.combo$both_IC - data.analysis.rt.combo$both_C, paired = T))
my_d_outboth <- print(ci.sm(sm = cohensD(data.analysis.rt.combo$out_IC - data.analysis.rt.combo$out_C, data.analysis.rt.combo$both_IC - data.analysis.rt.combo$both_C, method = "paired"), N = 102))

ttest.tstat(my_t_outboth$stat, n1 = 102, simple = T)
TOSTone(m = my_d_outboth$Standardized.Mean, mu = 0, sd = 1, n = 102, low_eqbound_d = -0.3, high_eqbound_d = 0.3)

	###### MAIN EFFECT CON ######

data.analysis.rt.combo.con.a <- dcast(ppn + cond + con ~ "rt", data = data.rt, value.var = "rt", mean)
data.analysis.rt.combo.con <- dcast(ppn  ~ con, data = data.analysis.rt.combo.con.a, value.var = "rt", mean)

t.test(data.analysis.rt.combo.con$IC, data.analysis.rt.combo.con$C, paired = T)
ci.sm(sm = cohensD(data.analysis.rt.combo.con$IC, data.analysis.rt.combo.con$C, method = "paired"), N = 102)

Lims <- conf.limits.ncf(F.value = 108.61, conf.level = 0.90, df.1 <- 1, df.2 <- 100)
Lims$Lower.Limit/(Lims$Lower.Limit + df.1 + df.2 + 1)
Lims$Upper.Limit/(Lims$Upper.Limit + df.1 + df.2 + 1)

	###### MAIN EFFECT CONDITION ######

data.analysis.rt.combo.cond <- dcast(ppn ~ cond, data = data.analysis.rt.combo.con.a, value.var = "rt", mean)
colMeans(data.analysis.rt.combo.cond)

Lims <- conf.limits.ncf(F.value = 1.49, conf.level = 0.90, df.1 <- 2, df.2 <- 99)
Lims$Lower.Limit/(Lims$Lower.Limit + df.1 + df.2 + 1)
Lims$Upper.Limit/(Lims$Upper.Limit + df.1 + df.2 + 1)

##################################
######    ERROR ANALYSIS    ######
##################################

	###### DATA ######
	
data.analysis.err.exp2$exp <- factor(2)
data.analysis.err.exp3$exp <- factor(3)
data.analysis.err.combo <- rbind(data.analysis.err.exp2, data.analysis.err.exp3)
data.analysis.err.combo$ppn <- 1:nrow(data.analysis.err.combo)

	###### VISUALIZE ######
	
###### CALCULATE ERROR BARS ######

data.plot.err.se <- dcast(ppn ~ cond + con, data = data.err, value.var = "err", mean)
SE_err_combo <- SE_Within(data = data.analysis.err.combo, colnames = names(data.analysis.err.combo)[-c(1,8)])$se.within.corr

###### MAKE PLOT ######

data.plot.err.combo <- data.plot.err
data.plot.err.combo$RT <- as.numeric(colMeans(data.analysis.err.combo[-c(1,8)]))

mytheme = trellis.par.get()
mytheme$superpose.polygon$col = c("white", "grey")
mytheme$fontsize$text = 22
mytheme$axis.text$cex = 1

plot_err_combo <- barchart(ErrorRate ~ cond, groups = con, data = data.plot.err.combo, box.ratio = 3,
	ylim = c(-0.5,8.5), scales = list(tck = c(1,0), y = list(at = seq(0, 8, 1))),
	xlab = "Condition", ylab = "ErrorRate",
	auto.key = list(corner = c(0.975,.975), points = F, rectangles = T, height = 0.85, size = 3.5, padding.text = 2.5),
	par.settings = mytheme,
	panel=function(x, y, subscripts, groups, stderr, box.ratio, ...){
		panel.barchart(x, y, subscripts = subscripts, groups = groups, box.ratio = box.ratio, ...)
		panel.err(x, y, subscripts = subscripts, groups = groups, box.ratio = box.ratio, stderr = SE_err)
		grid.text("(a)", x = unit(.05, "npc"), y = unit(.95, "npc"))
	})
print(plot_err_combo)

###### CONGRUENCY EFFECT ######

data.plot.err.combo.con <- data.frame(cond = c("in","out","both"), CE = NA)
data.plot.err.combo.con$cond <- factor(data.plot.err.combo.con$cond, levels = c("in","out","both"))

for(i in 1:nrow(data.plot.err.combo.con)){
	data.plot.err.combo.con$CE[i] <- data.plot.err.combo$RT[i*2] - data.plot.err.combo$RT[(i*2)-1]
}

barchart(CE ~ cond, data = data.plot.err.combo.con, col = "grey",par.settings = mytheme, ylim = c(-0.5,5), scales = list(tck = c(1,0), y = list(at = seq(0, 5, 1))))

###### RAINCLOUD PLOT ######

data.plot.err.combo.ppn <- melt(id.vars = c("ppn","exp"), value.name = "ErrorRate", data = data.analysis.err.combo)
data.plot.err.combo.ppn <- cbind(data.plot.err.combo.ppn, colsplit(data.plot.err.combo.ppn$variable, "_", c("cond","con"))) 
data.plot.err.combo.ppn <- select(data.plot.err.combo.ppn, ppn, cond, con, ErrorRate)
data.plot.err.combo.ppn$cond <- factor(data.plot.err.combo.ppn$cond, levels = c("in", "out", "both"))
data.plot.err.combo.ppn$con <- factor(data.plot.err.combo.ppn$con, levels = c("C", "IC"))

data.plot.err.combo.sum <- summarySE(data.plot.err.combo.ppn, measurevar = "ErrorRate", groupvars=c("cond", "con"))

rb_err_combo <- ggplot(data.plot.err.combo.ppn, aes(x = cond, y = ErrorRate, fill = con)) +
  geom_flat_violin(aes(fill = con),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(cond)-.15, y = ErrorRate, colour = con),position = position_jitter(width = .05), size = .25, shape = 20, show.legend = F)+
  geom_boxplot(aes(x = cond, y = ErrorRate, fill = con),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  geom_line(data = data.plot.err.combo.sum, aes(x = as.numeric(cond)+.1, y = ErrorRate_mean, group = con, colour = con), linetype = 3, show.legend = F)+
  geom_point(data = data.plot.err.combo.sum, aes(x = as.numeric(cond)+.1, y = ErrorRate_mean, group = con, colour = con), shape = 18, show.legend = F) +
  geom_errorbar(data = data.plot.err.combo.sum, aes(x = as.numeric(cond)+.1, y = ErrorRate_mean, group = con, colour = con, ymin = ErrorRate_mean-se, ymax = ErrorRate_mean+se), width = .05, show.legend = F)+
  scale_x_discrete(name = NULL, labels = c("In-Group","Out-Group","In+Out-Group")) + 
  scale_y_continuous(name = "Error Rate (%)", limits = c(-5,20), breaks = seq(-5,20,5))+#,expand = expand_scale(mult = c(0, .05))) + 
  scale_color_OkabeIto(order = c(5,6))+
  scale_fill_OkabeIto(order = c(5,6))+theme_cowplot()+#theme(plot.title = element_text(hjust = 0.5))+
  guides(fill= guide_legend(title=NULL))#+
  #ggtitle("Error Rates")
rb_err_combo

# ggsave('rb_err_combo.png', rb_err_combo, width = 6, height = 3)
 
	###### PRINT RT AND ER PLOTS ######
	
rb_both_combo_plot <- plot_grid(rb_rt_combo, rb_err_combo, nrow = 2)
rb_both_combo_title <- ggdraw() + draw_label("Automatic Imitation", fontface = 'bold')
rb_both_combo <- plot_grid(rb_both_combo_title, rb_both_combo_plot, ncol = 1, rel_heights = c(0.1, 1))

# ggsave('rb_rt+err_combo.png', rb_both_combo, width = 6, height = 6)

	###### ANOVA ######

cond <- factor(rep(1:3, each = 2))
con <- factor(rep(1:2, times = 3))
idata <- data.frame(cond, con)

fit.err.combo <- lm(cbind(in_C,in_IC,out_C,out_IC,both_C,both_IC) ~ exp, data = data.analysis.err.combo)
etasq(Anova(fit.err.combo, type = "III", test = "Wilks", idata = idata, idesign = ~cond*con), anova = T, partial = T)
etasq(Anova(fit.err.combo, type = "III", test = "Wilks", idata = idata, idesign = ~cond*con), anova = T, partial = T)$"approx F"

	###### CON X COND ######
	
###### ES ######

Lims <- conf.limits.ncf(F.value = 3.15, conf.level = 0.90, df.1 <- 2, df.2 <- 99)
Lims$Lower.Limit/(Lims$Lower.Limit + df.1 + df.2 + 1)
Lims$Upper.Limit/(Lims$Upper.Limit + df.1 + df.2 + 1)
	
###### IN vs. OUT ######

my_t_inout <- print(t.test(data.analysis.err.combo$in_IC - data.analysis.err.combo$in_C, data.analysis.err.combo$out_IC - data.analysis.err.combo$out_C, paired = T))
my_d_inout <- print(ci.sm(sm = cohensD(data.analysis.err.combo$in_IC - data.analysis.err.combo$in_C, data.analysis.err.combo$out_IC - data.analysis.err.combo$out_C, method = "paired"), N = 102))

ttest.tstat(my_t_inout$stat, n1 = 102, simple = T)
TOSTone(m = my_d_inout$Standardized.Mean, mu = 0, sd = 1, n = 102, low_eqbound_d = -0.3, high_eqbound_d = 0.3)

###### IN vs. BOTH ######

my_t_inboth <- print(t.test(data.analysis.err.combo$in_IC - data.analysis.err.combo$in_C, data.analysis.err.combo$both_IC - data.analysis.err.combo$both_C, paired = T))
my_d_inboth <- print(ci.sm(sm = cohensD(data.analysis.err.combo$in_IC - data.analysis.err.combo$in_C, data.analysis.err.combo$both_IC - data.analysis.err.combo$both_C, method = "paired"), N = 102))

ttest.tstat(my_t_inboth$stat, n1 = 102, simple = T)
TOSTone(m = my_d_inboth$Standardized.Mean, mu = 0, sd = 1, n = 102, low_eqbound_d = -0.3, high_eqbound_d = 0.3)

###### OUT vs. BOTH ######

my_t_outboth <- print(t.test(data.analysis.err.combo$out_IC - data.analysis.err.combo$out_C, data.analysis.err.combo$both_IC - data.analysis.err.combo$both_C, paired = T))
my_d_outboth <- print(ci.sm(sm = cohensD(data.analysis.err.combo$out_IC - data.analysis.err.combo$out_C, data.analysis.err.combo$both_IC - data.analysis.err.combo$both_C, method = "paired"), N = 102))

ttest.tstat(my_t_outboth$stat, n1 = 102, simple = T)
TOSTone(m = my_d_outboth$Standardized.Mean, mu = 0, sd = 1, n = 102, low_eqbound_d = -0.3, high_eqbound_d = 0.3)

	###### MAIN EFFECT CON ######

data.analysis.err.combo.con.a <- dcast(ppn + cond + con ~ "err", data = data.err, value.var = "err", mean)
data.analysis.err.combo.con <- dcast(ppn  ~ con, data = data.analysis.err.combo.con.a, value.var = "err", mean)

t.test(data.analysis.err.combo.con$IC, data.analysis.err.combo.con$C, paired = T)
ci.sm(sm = cohensD(data.analysis.err.combo.con$IC, data.analysis.err.combo.con$C, method = "paired"), N = 102)

Lims <- conf.limits.ncf(F.value = 40.16, conf.level = 0.90, df.1 <- 2, df.2 <- 100)
Lims$Lower.Limit/(Lims$Lower.Limit + df.1 + df.2 + 1)
Lims$Upper.Limit/(Lims$Upper.Limit + df.1 + df.2 + 1)

	###### MAIN EFFECT CONDITION ######

data.analysis.err.combo.cond <- dcast(ppn ~ cond, data = data.analysis.err.combo.con.a, value.var = "err", mean)
colMeans(data.analysis.err.combo.cond)

Lims <- conf.limits.ncf(F.value = 0.34, conf.level = 0.90, df.1 <- 2, df.2 <- 99)
Lims$Lower.Limit/(Lims$Lower.Limit + df.1 + df.2 + 1)
Lims$Upper.Limit/(Lims$Upper.Limit + df.1 + df.2 + 1)

	###### EXP x CONDITION ######
	
Lims <- conf.limits.ncf(F.value = 8.02, conf.level = 0.90, df.1 <- 2, df.2 <- 99)
Lims$Lower.Limit/(Lims$Lower.Limit + df.1 + df.2 + 1)
Lims$Upper.Limit/(Lims$Upper.Limit + df.1 + df.2 + 1)
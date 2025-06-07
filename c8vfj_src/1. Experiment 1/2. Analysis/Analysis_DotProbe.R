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
		file <- sprintf("MYPATH\\datadotprobe0%d.txt", ppn)
	}else{
		file <- sprintf("MYPATH\\datadotprobe%d.txt", ppn)
	}

	data <- read.table(file, header = F)
	
	###### NAME VARIABLES ######

	names <- c("trialnr","ppn","age","gender","handpos","group","block","stim","targetloc","ingrouploc","color_left","con","xr","r","corr","rt","re")
	names(data) <- names
	
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
	
rm(list = setdiff(ls(), c("make_list", "SE_Within", "data.analysis.rt.dotprobe.exp1", "data.analysis.err.dotprobe.exp1", "rb_dp_both_exp1")))	

	###### SET WD ######
	
setwd("MYPATH")

	###### RAINDCLOUD PLOT FUNCTIONS ######
	
source("R_rainclouds.R")
source("summarySE.R")

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
data.full <- select(data.full, -one_of("trialnr","group","stim","targetloc","ingrouploc","color_left","xr","r","re"))
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
data.full$con <- factor(data.full$con, levels = c(1,0), labels = c("C", "IC"))
	#congruency (i.e., congruent or incongruent) --> congruent = in-group side, incongruent = out-group side
data.full$corr <- factor(data.full$corr)
	#accuracy: 0 = error, 1 = correct
data.full$plus100 <- factor(data.full$plus100)
	#is RT slower than 100 ms?: 0 = no, 1 = yes
data.full$min2000 <- factor(data.full$min2000)
	#is RT faster than 2000 ms?: 0 = no, 1 = yes
data.full$max3SD <- factor(data.full$max3SD)
	#is RT within the 3 SD bounds?: 0 = no, 1 = yes

str(data.full)

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

data.descriptive <- filter(data.full, !(ppn %in% ppn.exclude.acc),!(ppn %in% ppn.exclude.rt))
data.descriptive <- dcast(ppn + age + gender ~ "N", data = data.descriptive, value.var = "rt", length)

print(age_m <- mean(data.descriptive$age, na.rm = T))
print(age_sd <- sd(data.descriptive$age, na.rm = T))
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

panel.err <- function(x, y, subscripts, stderr, box.ratio, ...)
{
  d <- 1/(nlevels(y)+nlevels(y)/box.ratio)
  panel.arrows(as.numeric(x),y-stderr[subscripts], as.numeric(x), y+stderr[subscripts],code=3,angle=90, length=0.3, unit = "cm")
}

	###### CALCULATE ERROR BARS ######

data.plot.rt.se <- dcast(ppn ~ con, data = data.rt, value.var = "rt", mean)
SE_rt <- SE_Within(data = data.plot.rt.se, colnames = names(data.plot.rt.se)[-1])$se.within.corr

	###### MAKE PLOT ######

data.plot.rt.a <- dcast(ppn + con ~ "RT", data = data.rt, value.var = "rt", mean)
data.plot.rt <- dcast(con ~ "RT", data = data.plot.rt.a, value.var = "RT", mean)

mytheme = trellis.par.get()
mytheme$plot.polygon$col = "grey"
mytheme$plot.polygon$lwd = 2
mytheme$add.line$lwd = 2
mytheme$fontsize$text = 18
mytheme$axis.text$cex = 1
mytheme$axis.line$lwd = 2

plot_rt <- barchart(RT ~ con, data = data.plot.rt, box.ratio = 3,
	ylim = c(404,436), scales = list(tck = c(1,0), x = list(labels = c("In-Group","Out-Group")), y = list(at = round(seq(405, 435, 10)))),
	ylab = "Reaction Time (ms)",
	auto.key = list(corner = c(0.975,.975), points = F, rectangles = T, height = 0.85, size = 3.5, padding.text = 2.5),
	par.settings = mytheme,
	panel=function(x, y, subscripts, stderr, box.ratio, ...)
	{
		panel.barchart(x, y, subscripts = subscripts, box.ratio = box.ratio, ...)
		panel.err(x, y, subscripts = subscripts, box.ratio = box.ratio, stderr = SE_rt)
	})
print(plot_rt)

#################################################
######    VISUALISATION RAINCLOUD PLOTS    ######
#################################################

data.plot.rt.ppn <- dcast(ppn + con ~ "RT", data = data.rt, value.var = "rt", mean)
data.plot.rt.sum <- summarySE(data.plot.rt.ppn, measurevar = "RT", groupvars=c("con"))
  
rb_dp_rt_exp1 <- ggplot(data.plot.rt.ppn, aes(x = con, y = RT, fill = con)) +
  geom_flat_violin(aes(fill = con),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(con)-.15, y = RT, colour = con),position = position_jitter(width = .05), size = .25, shape = 20, show.legend = F)+
  geom_boxplot(aes(x = con, y = RT, fill = con),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  geom_line(data = data.plot.rt.sum, aes(x = as.numeric(con)+.1, y = RT_mean, group = 1), linetype = 3, show.legend = F)+
  geom_point(data = data.plot.rt.sum, aes(x = as.numeric(con)+.1, y = RT_mean, colour = con), shape = 18, show.legend = F) +
  geom_errorbar(data = data.plot.rt.sum, aes(x = as.numeric(con)+.1, y = RT_mean, colour = con, ymin = RT_mean-se, ymax = RT_mean+se), width = .05, show.legend = F)+
  scale_x_discrete(name = NULL, labels = c("In-Group","Out-Group")) + 
  scale_y_continuous(name = "Reaction Time (ms)", limits = c(250,750), breaks = seq(250,750,100))+#,expand = expand_scale(mult = c(0, .05))) + 
  scale_color_OkabeIto(order = c(5,6))+
  scale_fill_OkabeIto(order = c(5,6))+theme_cowplot()+#theme(plot.title = element_text(hjust = 0.5))+
  guides(fill=FALSE)#+
  #ggtitle("Reaction Times")
rb_dp_rt_exp1
 
# ggsave('rb_dp_rt.png', rb_dp_rt_exp1, width = 3, height = 3, scale = 1.25)

############################
######    ANALYSIS    ######
############################

data.analysis.rt <- dcast(ppn ~ con, data = data.rt, value.var = "rt", mean)
data.analysis.rt.dotprobe.exp1 <- data.analysis.rt

my_t <- print(t.test(data.analysis.rt$IC, data.analysis.rt$C, paired = T))
my_d <- print(ci.sm(sm = cohensD(data.analysis.rt$IC, data.analysis.rt$C, method = "paired"), N = 38))

ttest.tstat(my_t$stat, n1 = 38, simple = T)
TOSTone(m = my_d$Standardized.Mean, mu = 0, sd = 1, n = 38, low_eqbound_d = -0.3, high_eqbound_d = 0.3)

	################################
	######    ERR ANALYSIS    ######
	################################
				
#########################################
######    CALCULATE ERROR RATES    ######
#########################################

data.err <- dcast(ppn + con ~ corr, data = data.corr, value.var = "rt", length)
names(data.err)[3:4] <- c("n_err","n_corr")
data.err <- mutate(data.err, err = (n_err/(n_err+n_corr))*100)
data.err <- select(data.err, -one_of("n_err", "n_corr"))

head(data.err)

###########################################
######    VISUALISATION BAR PLOTS    ######
###########################################

	###### ERROR BAR FUNCTION ######

panel.err <- function(x, y, subscripts, stderr, box.ratio, ...)
{
  d <- 1/(nlevels(y)+nlevels(y)/box.ratio)
  panel.arrows(as.numeric(x),y-stderr[subscripts], as.numeric(x), y+stderr[subscripts],code=3,angle=90, length=0.3, unit = "cm")
}

	###### CALCULATE ERROR BARS ######

data.plot.err.se <- dcast(ppn ~ con, data = data.err, value.var = "err", mean)
SE_err <- SE_Within(data = data.plot.err.se, colnames = names(data.plot.err.se)[-1])$se.within.corr

	###### MAKE PLOT ######

data.plot.err.a <- dcast(ppn + con ~ "ErrorRate", data = data.err, value.var = "err", mean)
data.plot.err <- dcast(con ~ "ErrorRate", data = data.plot.err.a, value.var = "ErrorRate", mean)

mytheme = trellis.par.get()
mytheme$plot.polygon$col = "grey"
mytheme$plot.polygon$lwd = 2
mytheme$add.line$lwd = 2
mytheme$fontsize$text = 18
mytheme$axis.text$cex = 1
mytheme$axis.line$lwd = 2

plot_err <- barchart(ErrorRate ~ con, data = data.plot.err, box.ratio = 3,
	ylim = c(-0.1,3.1), scales = list(tck = c(1,0), x = list(labels = c("In-Group","Out-Group")), y = list(at = format(seq(0, 3, 1),nsmall = 1))),
	ylab = "Error Rate (%)",
	auto.key = list(corner = c(0.975,.975), points = F, rectangles = T, height = 0.85, size = 3.5, padding.text = 2.5),
	par.settings = mytheme,
	panel=function(x, y, subscripts, stderr, box.ratio, ...)
	{
		panel.barchart(x, y, subscripts = subscripts, box.ratio = box.ratio, ...)
		panel.err(x, y, subscripts = subscripts, box.ratio = box.ratio, stderr = SE_err)
	})
print(plot_err)

	###### PRINT RT AND ER PLOTS ######
	
# ppi <- 300
# tiff(filename = "MYPATH\\dotprobe_exp1.tiff", compression = "zip", width = 7*2*ppi, height = 7*1*ppi, pointsize = 7*(9/5), res = ppi)

print(plot_rt, split = c(1, 1, 2, 1), more = T)
print(plot_err, split = c(2, 1, 2, 1), more = F)

# dev.off()

#################################################
######    VISUALISATION RAINCLOUD PLOTS    ######
#################################################

data.plot.err.ppn <- dcast(ppn + con ~ "ErrorRate", data = data.err, value.var = "err", mean)
data.plot.err.sum <- summarySE(data.plot.err.ppn, measurevar = "ErrorRate", groupvars=c("con"))
  
rb_dp_err_exp1 <- ggplot(data.plot.err.ppn, aes(x = con, y = ErrorRate, fill = con)) +
  geom_flat_violin(aes(fill = con),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(con)-.15, y = ErrorRate, colour = con),position = position_jitter(width = .05), size = .25, shape = 20, show.legend = F)+
  geom_boxplot(aes(x = con, y = ErrorRate, fill = con),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  geom_line(data = data.plot.err.sum, aes(x = as.numeric(con)+.1, y = ErrorRate_mean, group = 1), linetype = 3, show.legend = F)+
  geom_point(data = data.plot.err.sum, aes(x = as.numeric(con)+.1, y = ErrorRate_mean, colour = con), shape = 18, show.legend = F) +
  geom_errorbar(data = data.plot.err.sum, aes(x = as.numeric(con)+.1, y = ErrorRate_mean, colour = con, ymin = ErrorRate_mean-se, ymax = ErrorRate_mean+se), width = .05, show.legend = F)+
  scale_x_discrete(name = NULL, labels = c("In-Group","Out-Group")) + 
  scale_y_continuous(name = "Error Rate (%)", limits = c(-2,8), breaks = seq(-2,8,2))+#,expand = expand_scale(mult = c(0, .05))) + 
  scale_color_OkabeIto(order = c(5,6))+
  scale_fill_OkabeIto(order = c(5,6))+theme_cowplot()+#theme(plot.title = element_text(hjust = 0.5))+
  guides(fill=FALSE)#+
  #ggtitle("Reaction Times")
rb_dp_err_exp1
 
# ggsave('rb_dp_err.png', rb_dp_err_exp1, width = 3, height = 3, scale = 1.25)
	
	###### PRINT RT AND ER PLOTS ######
	
rb_dp_both_plot <- plot_grid(rb_dp_rt_exp1, rb_dp_err_exp1, ncol = 2)
rb_dp_both_title <- ggdraw() + draw_label("Dot-Probe", fontface = 'bold')
rb_dp_both_exp1 <- plot_grid(rb_dp_both_title, rb_dp_both_plot, ncol = 1, rel_heights = c(0.1, 1))

# ggsave('rb_dp_rt+err.png', rb_dp_both_exp1, width = 6, height = 3, scale = 1.25)	

############################
######    ANALYSIS    ######
############################

data.analysis.err <- dcast(ppn ~ con, data = data.err, value.var = "err", mean)
data.analysis.err.dotprobe.exp1 <- data.analysis.rt

my_t <- print(t.test(data.analysis.err$IC, data.analysis.err$C, paired = T))
my_d <- print(ci.sm(sm = cohensD(data.analysis.err$IC, data.analysis.err$C, method = "paired"), N = 38))

ttest.tstat(my_t$stat, n1 = 38, simple = T)
TOSTone(m = my_d$Standardized.Mean, mu = 0, sd = 1, n = 38, low_eqbound_d = -0.3, high_eqbound_d = 0.3)

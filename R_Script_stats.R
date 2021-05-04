
#restructuring and aggregating data formats
install.packages("reshape2")
#fits parametric distributions of data
install.packages("fitdistrplus")

#fitting GAM models
install.packages("gam")
#GAM model extensions
install.packages("mgcv")

#methods for auto-determiniation of axis breaks and scales
install.packages("scales")
#creates plots
install.packages("ggplot2")
#accessory themes for ggplot2
install.packages("ggthemes")
#customise ggplot2 outputs
install.packages("ggpubr")
#creating high quality figure outputs from ggplot2
install.packages("cowplot")
#grid graphic formats for plots
install.packages("gridExtra")

#models to determin egrowth rates from experimental data
install.packages("growthrates")
#fit simple non-linear models
install.packages("easynls")
#fit and compare linear and nl mixed effect models
install.packages("nlme")
#modification of nls algorithm
install.packages("minpack.lm")
#analysis of deose-response curves
install.packages("drc")
#tests and confidence intervals for general linear hypotheses
install.packages("multcomp")
#pairwise multiple comparisons of mean ranked sums
install.packages("PMCMRplus")
#time series analysis
install.packagaboutes("xts")
#fitting models to perform sensitivity analysis
install.packages("FME")
#breaks down data/models into smaller pieces
install.packages("plyr")
#visualises fits of regression models
install.packages("visreg")


###############################

library(mgcv)
library(gam)
library(cowplot)
library(fitdistrplus)
library(scales)
library(multcomp)
library(ggplot2)
library(ggthemes)
library(reshape2)
library(gridExtra)
library(PMCMR)
library(drc)
library(FSA)
library(minpack.lm)
library(easynls)
library(FME)
library(growthrates)
library(nlme)
library(plyr)
library(visreg)
library(ggpubr)


############################################# CALCIUM AND BICARBONATE GRAPHS

data = read.table("Calcium_juvenile_calc_rates.txt", header = TRUE)
data
data2 = read.table("Bicarbonate_juvenile_calc_rates.txt", header = TRUE)
data2


###VBGM


data2 = read.table("Bicarbonate_juvenile_calc_rates.txt", header = TRUE)
data2$salinity = as.factor(data2$salinity)
data2$salinity = as.character(data2$salinity)
data2$salinity = factor(data2$salinity, levels = c("16", "11", "6"))


F2 = formula(SM~Linf*(1-exp(-K*(bicarbonate_cont-380)))^3|salinity)
VBGM_all = nlsList(F2, data=data2,
              start=list(Linf=20.33,K=0.01))
summary(VBGM_all)
AIC(VBGM_all)

VBGM = nlsLM(SM~Linf*(1-exp(-K*(bicarbonate_cont-t0))), data=data2,
              start=list(Linf=20.33,K=1.27,t0=354))
summary(VBGM)
AIC(VBGM)

VBGM6 = nlsLM(SM~Linf*(1-exp(-K*(bicarbonate_cont-380)))^3, data=data26,
              start=list(Linf=4.26,K=0.27))
VBGM11 = nlsLM(SM~Linf*(1-exp(-K*(bicarbonate_cont-380)))^3, data=data211,
              start=list(Linf=30.8,K=0.36))
VBGM16 = nlsLM(SM~Linf*(1-exp(-K*(bicarbonate_cont-380)))^3, data=data216,
              start=list(Linf=31.8,K=0.2))

summary(VBGM16)
AIC(VBGM6, VBGM11, VBGM16)
AIC(VBGM6, GP6, LC6, ed6, mm6)
AIC(VBGM11, GP11, LC11, ed11, mm11)
AIC(VBGM16, GP16, LC16, ed16, mm16)
summary(VBGM11)
AIC(VBGM11)
summary(VBGM16)
AIC(VBGM16)

### nls logistic curve

LC6 = nls(SM~a*(1+b*(exp(-c*bicarbonate_cont)))^-1,
          data=data26,
          start=list(a= 5, b=380, c=0.01))
LC11 = nls(SM~a*(1+b*(exp(-c*bicarbonate_cont)))^-1,
          data=data211,
          start=list(a= 30, b=380, c=0.01))
LC16 = nls(SM~a*(1+b*(exp(-c*bicarbonate_cont)))^-1,
           data=data216,
           start=list(a= 30, b=380, c=0.01))
summary(LC6)
summary(LC11)
summary(LC16)

F3 = formula(SM~a*(1+380*(exp(-c*bicarbonate_cont)))^-1|salinity)
LC_all = nlsList(F3, data=data2,
                   start=list(a=30,c=0.01))
summary(LC_all)

### Gompertz model

data2m6 = data2m[c(1:5), c(1:5)]
data2m6
data2m11 = data2m[c(6:10), c(1:5)]
data2m11
data2m16 = data2m[c(11:15), c(1:5)]
data2m16

GP6 = nls(SM~a*exp(-380*exp(-c*bc)),
               start=list(a=5, c=0.01),
               data=data26)
GP11 = nls(SM~a*exp(-380*exp(-c*bicarbonate_cont)),
               start=list(a=40, c=0.01), data=data211)
GP16 = nls(SM~a*exp(-380*exp(-c*bicarbonate_cont)),
                start=list(a=35, c=0.01), data=data216)
summary(GP6)
summary(GP11)
summary(GP16)

#data26$VBGM = predict(fitTypical) # create new vector for mode predictions

plot(SM ~ bicarbonate_cont, data=data2, col = data2$salinity,
     xlim=c(0, 2500),
     ylim=c(0, 60))
lines(data26$bicarbonate_cont, predict(fitTypical6), col = "green")
lines(data26$bicarbonate_cont, predict(fitTypical11), col = "red")
lines(data26$bicarbonate_cont, predict(fitTypical16), col = "black")

#### michaelis menten model

mm6 = nls(SM ~ Vm * bicarbonate_cont / (K + bicarbonate_cont),
               start=list(Vm=4, K=500), data=data26)
mm11 = nls(SM ~ Vm * bicarbonate_cont / (K + bicarbonate_cont),
               start=list(Vm=30, K=500), data=data211)
mm16 = nls(SM ~ Vm * bicarbonate_cont / (K + bicarbonate_cont),
                start=list(Vm=30, K=500), data=data216)
summary(mm6)
summary(mm11)
summary(mm16)

######Yield-loss density

yld6 = nls(SM ~ (Vm/K) * bicarbonate_cont / (1+ ((Vm/K)*bicarbonate_cont)/Vm),
           start=list(Vm=6, K=400), data=data26)
summary(yld6)

#exponential decay model

ed6 = nls(SM ~ Vm * (1-exp(-bicarbonate_cont/K)),
              start=list(Vm = 5, K = 300), data=data26)
ed11 = nls(SM ~ Vm * (1-exp(-bicarbonate_cont/K)),
               start=list(Vm = 30, K = 300), data=data211)
ed16 = nls(SM ~ Vm * (1-exp(-bicarbonate_cont/K)),
               start=list(Vm = 30, K = 300), data=data216)

summary(ed6)
summary(ed11)
summary(ed16)

F2 = formula(log(SM+3) ~ Vm * (1-exp(-bicarbonate_cont/K))|salinity)
edm_all = nlsList(F2, data=data2,
                   start=list(Vm=10.33,K=300))
summary(edm_all)
AIC(VBGM_all)

############



data2_2$salinity = as.factor(data2_2$salinity)
data2_2$salinity = as.character(data2_2$salinity)
data2_2$salinity = factor(data2_2$salinity, levels = c("6", "11", "16"))


data2$salinity = as.factor(data2$salinity)
data2$salinity = as.character(data2$salinity)
data2$salinity = factor(data2$salinity, levels = c("16", "11", "6"))

data$salinity = as.factor(data$salinity)
data$salinity = as.character(data$salinity)
data$salinity = factor(data$salinity, levels = c("16", "11", "6"))


SM_CA = 
ggplot(data, aes(x = Ca, y = SM, colour = salinity)) + 
geom_point(size = 2.5, aes(x = Ca, y = SM, shape = salinity)) +
stat_smooth(method = "lm", formula = y ~ x, size = 0.5, se = F, linetype = "dashed") +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
xlab(expression(Calcium ~ conc. ~ (mmol ~ kg^{-1}  ))) +
ylab(expression(Calcification ~ rate ~ (?g ~ day^{-1}  ))) + 
theme_bw() + 
theme(legend.background = element_rect(fill="white", size=.5, 
linetype="dotted"), 
legend.justification=c(1,0), 
legend.position=c("none"), 
legend.title=element_blank(),
legend.text = element_text(size = 7),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_y_continuous(limits=c(-5, 60), expand = c(0, 0)) +
scale_x_continuous(limits=c(0, 4), expand = c(0, 0))

SM_BI = 
ggplot(data2, aes(x = bicarbonate_cont, y = SM, colour = salinity)) +
geom_point(size=2.5, aes(x = bicarbonate_cont, y = SM, shape = salinity)) +
geom_line(data = data2, aes(x=bicarbonate_cont, y=VBGM), linetype = "dashed", size = 0.5) +
#geom_smooth(method="nls",
#             formula=y~Vmax*(1-exp(-x/K)),
#             method.args = list(start=c(K=1000,Vmax=2)),
#             se=FALSE, colour="gray") +
#stat_smooth(method = "lm", formula = y ~ log(x), size = 1.2, se = F) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
xlab(expression(Bicarbonate ~ conc. ~ (?mol ~ kg^{-1}  ))) +
ylab(expression(Calcification ~ rate ~ (?g ~ day^{-1}  ))) + 
theme_bw() + 
theme(legend.background = element_rect(fill="white", size=.5, 
linetype="dotted"), 
legend.justification=c(1,0), 
legend.position=c(0.13, 0.45), 
legend.title=element_blank(),
legend.text = element_text(size = 7),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_y_continuous(limits=c(-5, 60), expand = c(0, 0)) +
scale_x_continuous(limits=c(0, 2500), expand = c(0, 0))

Linf_BI = 
ggplot(data2_2, aes(x = salinity, y = Linf, colour = salinity)) +
geom_point(size=4) +
geom_errorbar(aes(ymin=Linf-X95_CI, ymax=Linf+X95_CI), width=.1) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
xlab("Salinity") +
ylab("Linf") + 
theme_bw() + 
theme(legend.background = element_rect(fill="white", size=.5), legend.position="none", 
        panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
        panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(),
        axis.text=element_text(size=10, colour = "black")) +
  scale_y_continuous(limits=c(0, 40), expand = c(0, 0)) +
  scale_x_discrete()

#plot all graphs onto one output file

ion_plots_SM = ggdraw(xlim = c(0, 1.15), ylim = c(0, 1.15)) +
draw_plot(SM_BI, x = 0, y = 0.55, width = 1.1, height = 0.55) +
draw_plot(SM_CA, x = 0, y = 0, width = 0.55, height = 0.55) +
draw_plot(Linf_BI, x = 0.55, y = 0.001, width = 0.55, height = 0.55) 
pdf("ion_growth_plots_SM.pdf", width = 6, height = 5)
ion_plots_SM
dev.off()

######################################## Graph 2
######################################## Field SL growth rates

###############SL-SM relationships

dataSLSM = read.table("SL_SM_relationship.txt", header=T)
dataSLSM

dataSLSM$Population = as.factor(dataSLSM$Population)

shapiro.test(dataSLSM$CaCO3_ug)
leveneTest(CaCO3_ug ~ Population, data=dataSLSM)

lmLM1 = lm(CaCO3_ug ~ (SL_mm)^2 * as.factor(Population), data=dataSLSM)
summary(lmLM1)
anova(lmLM1)

###################### subsett SL-SM data by population

dataSLkie = dataSLSM[c(1:124), c(2:3)]
dataSLkie
dataSLahp = dataSLSM[c(125:216), c(2:3)]
dataSLahp
dataSLuse = dataSLSM[c(217:296), c(2:3)]
dataSLuse


power_kie = nlsLM(CaCO3_ug~a*SL_mm^b, data=dataSLkie,
              start=list(a=1, b=1))
power_ahp = nlsLM(CaCO3_ug~a*SL_mm^b, data=dataSLahp,
              start=list(a=1, b=1))
power_use = nlsLM(CaCO3_ug~a*SL_mm^b, data=dataSLuse,
              start=list(a=1, b=1))

summary(power_use)
anova(power_kie, power_ahp, power_use)

plot(CaCO3_ug ~ (SL_mm)^2, data=dataSLSM)

dataSLSM$Population = as.factor(dataSLSM$Population)
dataSLSM$Population = factor(dataSLSM$Population, levels = c("Kiel", "Ahp", "Use"))

pdf("SL_SM_relationship.pdf", width = 5, height = 5)

ggplot(dataSLSM, aes(x = SL_mm, y = CaCO3_ug, colour = Population)) + 
geom_point(size = 1, position = position_dodge(width = .2)) +
stat_smooth(method = 'nls', formula = 'y~a*x^b', start = list(a = 1,b=1),se=F) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
xlab("Shell length (mm)") +
ylab("Shell mass (?g)") +
theme_bw() +
theme(legend.background = element_rect(fill="white", size=.5, linetype="dotted"), 
legend.justification=c(1,0), 
legend.position=c(0.3,0.77), 
legend.title=element_blank(),
legend.text = element_text(size = 10),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_x_continuous(limits=c(0, 20), expand = c(0, 0)) +
scale_y_continuous(limits=c(0, 250))


dev.off()




dataSL = read.table("field_shell_growth_2016-2017.txt", header = T)
dataSL

dataSL6 = subset(dataSL, population=="Usedom", drop=T)
dataSL11 = subset(dataSL, population=="Ahrenshoop", drop=T)
dataSL16 = subset(dataSL, population=="Kiel", drop=T)

dataSL$month = as.numeric(dataSL$month)
dataSLm$month = as.numeric(dataSLm$month)
dataSL$population = as.factor(dataSL$population)

cdata = ddply(dataSL, c("month", "population"), summarise,
              N = length(SM),
              mean = mean(SM),
              sd = sd(SM),
              se = sd/sqrt(N))
cdata
dataSMm = data.frame(cdata)
dataSMm


write.table(dataSMm,"SM_growth_means.txt",sep="\t",row.names=T)
dataSMm = read.table("SM_growth_means.txt", header = T)
dataSMm
dataSMm17$SM_ug = (dataSMm17$mean*1000)
dataSMm17$day = (dataSMm17$month*30)

dataSMm17 = dataSMm[-c(10, 11, 12, 13, 14, 15, 16, 17, 18), ]
dataSMm17 

lm17 = lm(log(SM_ug) ~ day * population, data=dataSMm17)
summary(lm17)


shapiro.test(log(dataSL$SM_pow))

plot1 = plot(log(SM_ug) ~ day, colour = population, data = dataSMm17)
lm1 = aov(mean ~ I(month^2) * as.factor(population), data=dataSMm)
anova(lm1)
summary(lm1)
coef(lm1)

lm2 = aov(mean ~ I(month^2) * as.factor(population), data=dataSMm)
anova(lm2)
TukeyHSD(lm1, "as.factor(population)", conf.level = 0.95)

z = dataSL$population = as.factor(dataSL$population)


plotdist(dataSLm$mean, histo = TRUE, demp = TRUE)
descdist(dataSL$SM, discrete = F, boot=1000)
qqnorm(dataSL$SM);qqline(log(dataSL$SM), col = 2)
hist(log(dataSL$SM))

TukeyHSD(aov1, "as.factor(population)", conf.level = 0.95)

dataSMm

plot(mean ~ month, colour = population, data = dataSMm19)

#plot mean environmental conditions with calcification rates and compare with experiments

dataSMm$population = as.factor(dataSMm$population)
dataSMm$population = factor(dataSMm$population, levels = c("kiel", "ahp", "use"))


pdf("SM_field.pdf", width = 5, height = 5)

ggplot(dataSMm, aes(x = month, y = mean, colour = population)) + 
geom_point(size = 3.5, position = position_dodge(width = .2)) +
#geom_smooth(method = "nls", formula=y~(a*x)^b, start=c(a=1, b=1), se=F) +
stat_smooth(method='lm', formula = y~I(x^2), se = F,
#            method.args = list(family = "exp"),
             linetype = "dashed", size = 0.55) +
geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(.2)) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
xlab("Month") +
ylab("Shell mass (mg)") +
theme_bw() +
theme(legend.background = element_rect(fill="white", size=.5, linetype="dotted"), 
legend.justification=c(1,0), 
legend.position=c(0.3,0.77), 
legend.title=element_blank(),
legend.text = element_text(size = 10),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_x_continuous(limits=c(0, 16), expand = c(0, 0)) +
scale_y_continuous(limits=c(0, 2500))


dev.off()


###################################### Graph 3
###################################### Field Temperature

install.packages("dunn.test")


descdist(datac$pH, discrete = F)
fit.norm = fitdist(log(datac$HCO3), "norm")
plot(fit.norm)

leveneTest(datac$pH ~ as.factor(site), data=datac)


shapiro.test(datac$HCO3)

lmph = aov(pH ~ as.factor(site), data = datac)
anova(lmph)

TukeyHSD(lmph, "as.factor(site)", conf.level = 0.95)


datac = read.table("Carb_chem_all.txt", header = TRUE)
datac

datac$date = as.Date(datac$date, "%d/%m/%Y")
#datac$date = as.Date(datac$date, "%Y-%m-%d")

datac$site = as.character(datac$site)
datac$site = factor(datac$site, levels = c("Kiel", "Ahrenshoop", "Usedom"))

pH_all_sp = 
ggplot(datac, aes(x = date, y = pH, colour = site)) + 
geom_point(size = 2) +
#geom_smooth(method = "lm", se = F) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
xlab("Date") + 
#ylab(expression(italic(C)[T]  ~ (?mol ~ kg^{-1}))) +
theme_bw() + 
theme(legend.background = element_rect(fill="white", size=.5, linetype="dotted"), 
legend.justification=c(1,0), 
legend.position="none", 
legend.title=element_blank(),
legend.text = element_text(size = 12),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_x_date(breaks=seq(from=as.Date("2014-12-01"), to=as.Date("2017-12-01"), by="6 month"), 
labels=date_format("%b. '%y"), limits=as.Date(c("2014-12-01", "2017-12-01"), expand = c(0, 0))) +
scale_y_continuous(limits=c(7.4, 8.4), expand = c(0, 0))

pH_all_bp = 
ggplot(datac, aes(x = site, y = pH, colour = site)) + 
geom_boxplot(size = .5, outlier.size = .05) +
#geom_smooth(method = "lm", se = F) +
scale_color_manual(values=c("deepskyblue1", "seagreen4", "coral1")) +
xlab("Site") + 
ylab("") +
theme_bw() + 
theme(legend.background = element_rect(fill="white", size=.5, linetype="dotted"), 
legend.justification=c(1,0), 
legend.position="none", 
legend.title=element_blank(),
legend.text = element_text(size = 12),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_y_continuous(limits=c(7.4, 8.4), expand = c(0, 0))

field_pH_plots = ggdraw(xlim = c(0, 3), ylim = c(0, 1)) +
draw_plot(pH_all_sp, x = 0, y = 0, width = 2, height = 1) +
draw_plot(pH_all_bp, x = 2, y = 0, width = 1, height = 1) 

pdf("field_pH_plots.pdf", width = 9, height = 3)

field_pH_plots

dev.off()


####################################### Graph 4
####################################### Field Salinity plots

descdist(log(data_S$value[1:5000]), discrete = F)
fit.norm = fitdist(log(data_S$value[1:5000]), "norm")
plot(fit.norm)


lm2 = aov(log(value+3) ~ site, data = data_S)
anova(lm2)
TukeyHSD(lm2, "site", conf.level = 0.95)

kruskal.test(value ~ site, data = data_S)
dunnTest(data_S$value, g = data_S$site, kw=TRUE, label=TRUE, wrap=FALSE, alpha=0.05)




data_S = read.table("Salinity_all_2015_2018.txt", sep = "\t", header = TRUE)
data_S

data_S$date = as.Date(data_S$date, "%d.%m.%Y")
data_S$date = as.Date(data_S$date, "%Y-%m-%d")

data_S$site = as.character(data_S$site)
data_S$site = factor(data_S$site, levels = c("Kiel", "Ahrenshoop", "Usedom"))


sal_all_sp = ggplot(data_S, aes(x = date, y = value, colour = site)) + 
geom_point(size = .1) +
#stat_smooth(method = "gam", formula = y ~ s(x, bs = "ps", k = 12), size = 1, se = F) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
xlab("Date") + 
ylab("Salinity (psu)") +
theme_bw() + 
theme(legend.background = element_rect(fill="white", size=.5, 
linetype="dotted"), 
legend.justification=c(1,0), 
legend.position="none", 
legend.title=element_blank(),
legend.text = element_text(size = 8),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_x_date(breaks=seq(from=as.Date("2015-8-01"), to=as.Date("2017-12-01"), by="4 month"), 
labels=date_format("%b. '%y"), limits=c(), expand = c(0, 0)) +
scale_y_continuous(limits=c(0, 25), expand = c(0, 0))


sal_all_bp = ggplot(data_S, aes(x = site, y = value, colour = site)) + 
geom_boxplot(size = .5, outlier.size = .05) +
stat_smooth(method = "gam", formula = y ~ s(x, bs = "ps", k = 12), size = 1.5, se = F) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
xlab("Site") + 
ylab("") +
theme_bw() + 
theme(legend.background = element_rect(fill="white", size=.5, 
linetype="dotted"), 
legend.justification=c(1,0), 
legend.position="none", 
legend.title=element_blank(),
legend.text = element_text(size = 8),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_y_continuous(limits=c(0, 25), expand = c(0, 0))

S_plots = ggdraw(xlim = c(0, 3), ylim = c(0, 1)) +
draw_plot(sal_all_sp, x = 0, y = 0, width = 2, height = 1) +
draw_plot(sal_all_bp, x = 2, y = 0, width = 1, height = 1) 

pdf("S_plots.pdf", width = 9, height = 3)

S_plots

dev.off()

############################################## Graph 5
############################################## Chl_a monitoring data

descdist(log(data_chla$chla), discrete = F)
fit.norm = fitdist(log(data_chla$chla), "norm")
plot(fit.norm)


lm3 = aov(log(chla) ~ Site, data = data_chla)
anova(lm3)
TukeyHSD(lm3, "Site", conf.level = 0.95)

kruskal.test(value ~ site, data = data_S)
dunnTest(data_S$value, g = data_S$site, kw=TRUE, label=TRUE, wrap=FALSE, alpha=0.05)


data_chla = read.table("chla_monitored_2015-2017.txt", sep="\t", header = T)
data_chla

mean(data_chla$chla[1:25])

data_chla$date = as.Date(data_chla$date, "%d/%m/%Y")

data_chla$Site = as.character(data_chla$Site)
data_chla$Site = factor(data_chla$Site, levels = c("Kiel", "Ahrenshoop", "Usedom"))

data_chla$chla = as.numeric(data_chla$chla)

chl_monitor_sp = 
ggplot(data_chla, aes(x = date, y = chla, colour = Site, shape = Site)) + 
geom_point(size = 4) +
#stat_smooth(method = "loess", formula = y ~ x, size = .6, se = F) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
xlab("") + 
ylab(expression(Chl-italic(a) ~ (?g ~ L^{-1} ))) +
theme_bw() + 
theme(legend.background = element_rect(fill="white", size=.5, linetype="dotted"), 
legend.justification=c(1,0), 
legend.position="none", 
legend.title=element_blank(),
legend.text = element_text(size = 8),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_x_date(breaks=seq(from=as.Date("2014-12-01"), to=as.Date("2017-12-01"), by="6 month"), 
labels=date_format("%b. '%y"), limits=as.Date(c("2014-12-01", "2017-12-01"), expand = c(0, 0))) +
scale_y_continuous(limits=c(0, 20), expand = c(0, 0))

chl_monitor_bp = 
ggplot(data_chla, aes(x = Site, y = chla, colour = Site)) + 
geom_boxplot(size = .5, outlier.size = .05) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
xlab("") + 
ylab("") +
theme_bw() + 
theme(legend.background = element_rect(fill="white", size=.5, linetype="dotted"), 
legend.justification=c(1,0), 
legend.position="none", 
legend.title=element_text(size=8),
legend.text = element_text(size = 8),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_y_continuous(limits=c(0, 20), expand = c(0, 0))


chla_plots = ggdraw(xlim = c(0, 3.03), ylim = c(0, 1)) +
draw_plot(chl_monitor_sp, x = 0, y = 0, width = 2, height = 1) +
draw_plot(chl_monitor_bp, x = 2.03, y = 0, width = 1, height = 1) 

pdf("chl_monitor.pdf", width = 9, height = 3)

chla_plots

dev.off()

######################################## modelled chla data

descdist(log(data_chla$chla), discrete = F)
fit.norm = fitdist(log(data_chla$chla), "norm")
plot(fit.norm)


lm3 = aov(log(chla) ~ Site, data = data_chla)
anova(lm3)
TukeyHSD(lm3, "Site", conf.level = 0.95)

kruskal.test(value ~ site, data = data_S)
dunnTest(data_S$value, g = data_S$site, kw=TRUE, label=TRUE, wrap=FALSE, alpha=0.05)


data_chla2 = read.table("chla_modelled_2015-2017.txt", sep="\t", header = T)
data_chla2

data_chla3 = read.table("mean_chla_modelled_2015-2017.txt", sep="\t", header = T)
data_chla3

data_chla3$date = as.Date(data_chla3$date, "%d/%m/%Y")

data_chla3$site = as.character(data_chla3$site)
data_chla3$site = factor(data_chla3$site, levels = c("Kiel", "Ahrenshoop", "Usedom"))

data_chla3$chla = as.numeric(data_chla3$chl)

chl2_model_sp = 
ggplot(data_chla3, aes(x = date, y = chl, colour = site)) + 
geom_line(size = 1) +
#stat_smooth(method = "gam", formula = y ~ s(x), size = .6, se = F) +
scale_color_manual(values=c("deepskyblue1", "seagreen4", "coral1")) +
xlab("Date") + 
ylab(expression(Chl-italic(a) ~ (?g ~ L^{-1} ))) +
theme_bw() + 
theme(legend.background = element_rect(fill="white", size=.5, linetype="dotted"), 
legend.justification=c(1,0), 
legend.position="none", 
legend.title=element_blank(),
legend.text = element_text(size = 8),
panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_x_date(breaks=seq(from=as.Date("2016-07-01"),
	 to=as.Date("2017-12-01"),
	 by="3 month"), 
	 labels=date_format("%b. %y")) +
scale_y_continuous(limits=c(0, 25), expand = c(0, 0))

chl2_model_bp = 
ggplot(data_chla3, aes(x = site, y = chl, colour = site)) + 
geom_boxplot(size = .5, outlier.size = .05) +
scale_color_manual(values=c("deepskyblue1", "seagreen4", "coral1")) +
xlab("Site") + 
ylab("") +
theme_bw() + 
theme(legend.background = element_rect(fill="white", size=.5, linetype="dotted"), 
legend.justification=c(1,0), 
legend.position="none", 
legend.title=element_text(size=8),
legend.text = element_text(size = 8),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_y_continuous(limits=c(0, 25), expand = c(0, 0))


chla2_plots = ggdraw(xlim = c(0, 3.03), ylim = c(0, 1)) +
draw_plot(chl2_model_sp, x = 0, y = 0, width = 2, height = 1) +
draw_plot(chl2_model_bp, x = 2.03, y = 0, width = 1, height = 1) 

pdf("chl_modelled.pdf", width = 9, height = 3)

chla2_plots

dev.off()


####################################### Field HCO3

descdist(log(datac$HCO3), discrete = F)
fit.norm = fitdist(log(datac$HCO3), "norm")
plot(fit.norm)

leveneTest(datac$HCO3 ~ as.factor(site), data=datac)


shapiro.test(datac$HCO3)

lm5 = aov(HCO3 ~ as.factor(site), data = datac)
anova(lm5)

TukeyHSD(lm5, "as.factor(site)", conf.level = 0.95)


datac = read.table("Carb_chem_all.txt", header = TRUE)
datac

datac$date = as.Date(datac$date, "%d/%m/%Y")
#datac$date = as.Date(datac$date, "%Y-%m-%d")

datac$site = as.character(datac$site)
datac$site = factor(datac$site, levels = c("Kiel", "Ahrenshoop", "Usedom"))

cc_all_sp = 
ggplot(datac, aes(x = date, y = HCO3, colour = site, shape = site)) + 
geom_point(size = 4) +
#geom_smooth(method = "lm", se = F) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
xlab("Date") + 
#ylab(expression(italic(C)[T]  ~ (?mol ~ kg^{-1}))) +
theme_bw() + 
theme(legend.background = element_rect(fill="white", size=.5, linetype="dotted"), 
legend.justification=c(1,0), 
legend.position="none", 
legend.title=element_blank(),
legend.text = element_text(size = 12),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_x_date(breaks=seq(from=as.Date("2014-12-01"), to=as.Date("2017-12-01"), by="6 month"), 
labels=date_format("%b. '%y"), limits=as.Date(c("2014-12-01", "2017-12-01"), expand = c(0, 0))) +
scale_y_continuous(limits=c(1500, 2100), expand = c(0, 0))

cc_all_bp = 
ggplot(datac, aes(x = site, y = HCO3, colour = site)) + 
geom_boxplot(size = .5, outlier.size = .05) +
#geom_smooth(method = "lm", se = F) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
xlab("Site") + 
ylab("") +
theme_bw() + 
theme(legend.background = element_rect(fill="white", size=.5, linetype="dotted"), 
legend.justification=c(1,0), 
legend.position="none", 
legend.title=element_blank(),
legend.text = element_text(size = 12),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_y_continuous(limits=c(1500, 2100), expand = c(0, 0))

field_hco3_plots = ggdraw(xlim = c(0, 3), ylim = c(0, 1)) +
draw_plot(cc_all_sp, x = 0, y = 0, width = 2, height = 1) +
draw_plot(cc_all_bp, x = 2, y = 0, width = 1, height = 1) 

pdf("field_hco3_plots.pdf", width = 9, height = 3)

field_hco3_plots

dev.off()

############################################ Field omega values

descdist(log(datac$omega_ara), discrete = F)
fit.norm = fitdist(log(datac$omega_ara), "norm")
plot(fit.norm)

leveneTest(datac$omega_ara ~ as.factor(site), data=datac)


shapiro.test(log(datac$omega_ara))

lm6 = aov(log(omega_ara) ~ as.factor(site), data = datac)
anova(lm6)

TukeyHSD(lm6, "as.factor(site)", conf.level = 0.95)


datac = read.table("Carb_chem_all.txt", header = TRUE)
datac

datac$date = as.Date(datac$date, "%d/%m/%Y")
#datac$date = as.Date(datac$date, "%Y-%m-%d")

datac$site = as.character(datac$site)
datac$site = factor(datac$site, levels = c("Kiel", "Ahrenshoop", "Usedom"))

cc_omega_sp = 
ggplot(datac, aes(x = date, y = omega_ara, colour = site, shape = site)) + 
geom_point(size = 4) +
#geom_smooth(method = "lm", se = F) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
xlab("Date") + 
#ylab(expression(italic(C)[T]  ~ (?mol ~ kg^{-1}))) +
theme_bw() + 
theme(legend.background = element_rect(fill="white", size=.5, linetype="dotted"), 
legend.justification=c(1,0), 
legend.position="none", 
legend.title=element_blank(),
legend.text = element_text(size = 12),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_x_date(breaks=seq(from=as.Date("2014-12-01"), to=as.Date("2017-12-01"), by="6 month"), 
labels=date_format("%b. '%y"), limits=as.Date(c("2014-12-01", "2017-12-01"), expand = c(0, 0))) +
scale_y_continuous(limits=c(0,2), expand = c(0, 0))

cc_omega_bp = 
ggplot(datac, aes(x = site, y = omega_ara, colour = site)) + 
geom_boxplot(size = .5, outlier.size = .05) +
#geom_smooth(method = "lm", se = F) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
xlab("Site") + 
ylab("") +
theme_bw() + 
theme(legend.background = element_rect(fill="white", size=.5, linetype="dotted"), 
legend.justification=c(1,0), 
legend.position="none", 
legend.title=element_blank(),
legend.text = element_text(size = 12),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_y_continuous(limits=c(0,2), expand = c(0, 0))

field_omega_plots = ggdraw(xlim = c(0, 3), ylim = c(0, 1)) +
draw_plot(cc_omega_sp, x = 0, y = 0, width = 2, height = 1) +
draw_plot(cc_omega_bp, x = 2, y = 0, width = 1, height = 1) 

pdf("field_omega_plots.pdf", width = 9, height = 3)

field_omega_plots

dev.off()

############################################ Field ESIR values

descdist(log(datac$ex_SIR), discrete = F)
fit.norm = fitdist(log(datac$ex_SIR), "norm")
plot(fit.norm)

leveneTest(datac$ex_SIR ~ as.factor(site), data=datac)


shapiro.test(log(datac$ex_SIR))

lm7 = aov(log(ex_SIR) ~ as.factor(site), data = datac)
anova(lm7)

TukeyHSD(lm7, "as.factor(site)", conf.level = 0.95)


datac = read.table("Carb_chem_all.txt", header = TRUE)
datac

datac$date = as.Date(datac$date, "%d/%m/%Y")
#datac$date = as.Date(datac$date, "%Y-%m-%d")

datac$site = as.character(datac$site)
datac$site = factor(datac$site, levels = c("Kiel", "Ahrenshoop", "Usedom"))

cc_ESIR_sp = 
ggplot(datac, aes(x = date, y = ex_SIR, colour = site, shape = site)) + 
geom_point(size = 4) +
#geom_smooth(method = "lm", se = F) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
xlab("Date") + 
#ylab(expression(italic(C)[T]  ~ (?mol ~ kg^{-1}))) +
theme_bw() + 
theme(legend.background = element_rect(fill="white", size=.5, linetype="dotted"), 
legend.justification=c(1,0), 
legend.position="none", 
legend.title=element_blank(),
legend.text = element_text(size = 12),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_x_date(breaks=seq(from=as.Date("2014-12-01"), to=as.Date("2017-12-01"), by="6 month"), 
labels=date_format("%b. '%y"), limits=as.Date(c("2014-12-01", "2017-12-01"), expand = c(0, 0))) +
scale_y_continuous(limits=c(0,1.6), breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2), expand = c(0, 0))

cc_ESIR_bp = 
ggplot(datac, aes(x = site, y = ex_SIR, colour = site)) + 
geom_boxplot(size = .5, outlier.size = .05) +
#geom_smooth(method = "lm", se = F) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
xlab("Site") + 
ylab("") +
theme_bw() + 
theme(legend.background = element_rect(fill="white", size=.5, linetype="dotted"), 
legend.justification=c(1,0), 
legend.position="none", 
legend.title=element_blank(),
legend.text = element_text(size = 12),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_y_continuous(limits=c(0,1.6), breaks=c(0, 0.4, 0.8, 1.2,  1.6),
expand = c(0, 0))


field_ESIR_plots = ggdraw(xlim = c(0, 3), ylim = c(0, 1)) +
draw_plot(cc_ESIR_sp, x = 0, y = 0, width = 2, height = 1) +
draw_plot(cc_ESIR_bp, x = 2, y = 0, width = 1, height = 1) 

pdf("field_ESIR_plots.pdf", width = 9, height = 3)

field_ESIR_plots

dev.off()

###################################### Field pH

datac = read.table("Carb_chem_all.txt", header = TRUE)
datac

shapiro.test(datac$pH)

lm7 = aov(pH ~ as.factor(site), data = datac)
anova(lm7)

TukeyHSD(lm7, "as.factor(site)", conf.level = 0.95)



datac$date = as.Date(datac$date, "%d/%m/%Y")
#datac$date = as.Date(datac$date, "%Y-%m-%d")

datac$site = as.character(datac$site)
datac$site = factor(datac$site, levels = c("Kiel", "Ahrenshoop", "Usedom"))

cc_pH_sp = 
ggplot(datac, aes(x = date, y = pH, colour = site, shape = site)) + 
geom_point(size = 4) +
#geom_smooth(method = "lm", se = F) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
xlab("Date") + 
#ylab(expression(italic(C)[T]  ~ (?mol ~ kg^{-1}))) +
theme_bw() + 
theme(legend.background = element_rect(fill="white", size=.5, linetype="dotted"), 
legend.justification=c(1,0), 
legend.position="none", 
legend.title=element_blank(),
legend.text = element_text(size = 12),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_x_date(breaks=seq(from=as.Date("2014-12-01"), to=as.Date("2017-12-01"), by="6 month"), 
labels=date_format("%b. '%y"), limits=as.Date(c("2014-12-01", "2017-12-01"), expand = c(0, 0))) +
scale_y_continuous(limits=c(7.4, 8.4), expand = c(0, 0))

cc_pH_bp = 
ggplot(datac, aes(x = site, y = pH, colour = site)) + 
geom_boxplot(size = .5, outlier.size = .05) +
#geom_smooth(method = "lm", se = F) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
xlab("Site") + 
ylab("") +
theme_bw() + 
theme(legend.background = element_rect(fill="white", size=.5, linetype="dotted"), 
legend.justification=c(1,0), 
legend.position="none", 
legend.title=element_blank(),
legend.text = element_text(size = 12),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_y_continuous(limits=c(7.4, 8.4),
expand = c(0, 0))


field_pH_plots = ggdraw(xlim = c(0, 3), ylim = c(0, 1)) +
draw_plot(cc_pH_sp, x = 0, y = 0, width = 2, height = 1) +
draw_plot(cc_pH_bp, x = 2, y = 0, width = 1, height = 1) 

pdf("field_pH_plots.pdf", width = 9, height = 3)

field_pH_plots

dev.off()


###################################### SIR Omega

descdist(log(data_S_O$ESIR), discrete = F)
descdist(data_S_O$shell_growth, discrete = F)
fit.norm = fitdist(log(data_S_O$ESIR), "norm")
plot(fit.norm)

leveneTest(log(data_S_O$ESIR) ~ as.factor(salinity), data=data_S_O)

shapiro.test(log(data_S_O$ESIR))

lm9 = lm(log(shell_growth+3) ~ ESIR, data = data_S_O)
anova(lm9)

coef(lm8)

TukeyHSD(lm7, "as.factor(site)", conf.level = 0.95)


data_S_O = read.table("SIR_omega_calcification.txt", sep = "\t", header = T)
data_S_O

data_S_O$salinity = as.character(data_S_O$salinity)
data_S_O$salinity = factor(data_S_O$salinity, levels = c("16", "11", "6"))

# saturation curve / MM constant

MM_SIR = nlsLM(shell_growth ~ Vm * SIR / (K + SIR),
          start=list(Vm=0.7, K=35), data=data_S_O)

EM_SIR = nlsLM(shell_growth ~ 1 + Vm * (1-exp(-SIR/K)),
          start=list(Vm=0.7, K=32), data=data_S_O)

MM_ESIR = nlsLM(shell_growth ~ Vm * ESIR / (K + ESIR),
          start=list(Vm=0.7, K=35), data=data_S_O)

EM_ESIR = nlsLM(shell_growth ~ 1 + Vm * (1-exp(-ESIR/K)),
          start=list(Vm=0.7, K=32), data=data_S_O)

MM_omega = nlsLM(shell_growth ~ Vm * omega_ara / (K + omega_ara),
                start=list(Vm=0.7, K=35), data=data_S_O)

EM_omega = nlsLM(shell_growth ~ 1 + Vm * (1-exp(-omega_ara/K)),
                start=list(Vm=0.7, K=32), data=data_S_O)

AIC(MM_SIR,EM_SIR)
summary(EM_SIR)
 	
AIC(MM_omega,EM_omega)
summary(EM_omega)

anova(EM_ESIR, EM_omega)
summary()

######################## Graph 11 calcification vs field parameters PANEL FIGURE

data_FCL = read.table("Field_calc_linear.txt", sep = "\t", header = T)
data_FCL

lm_lc = lm(SM_log ~ days * as.factor(population), data = data_FCL)
summary(lm_lc)

data_CvE = read.table("Field_calc_vs_environment.txt", sep = "\t", header = T)
data_CvE

mod_CvE = lm(ESIR ~ calc_rate_log, data=data_CvE)
summary(mod_CvE)

data_CvE$site = as.character(data_CvE$site)
data_CvE$site = factor(data_CvE$site, levels = c("Kiel", "Ahrenshoop", "Usedom"))


sal_plot = 
ggplot(data_CvE, aes(x = salinity, y = calc_rate_log, colour = site)) + 
geom_point(size = 3, aes(shape=site)) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
#stat_smooth(method = "lm", formula = y ~ x, size = 1, se = F) +
geom_errorbar(aes(ymin=calc_rate_log-calc_log_se, 
	ymax=calc_rate_log+calc_log_se), 
	width=.3, 
	position=position_dodge(0.05)) +
geom_errorbarh(aes(xmin = salinity-sal_se, xmax = salinity+sal_se)) +
ylab("") +
xlab("salinity") +
theme_bw() +
theme(,legend.position="none", panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_x_continuous(limits=c(0, 20), expand = c(0, 0))

temp_plot = 
ggplot(data_CvE, aes(x = temp, y = calc_rate_log, colour = site)) + 
geom_point(size = 3, aes(shape=site)) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
geom_errorbar(aes(ymin=calc_rate_log-calc_log_se, 
	ymax=calc_rate_log+calc_log_se), 
	width=.3, 
	position=position_dodge(0.05)) +
geom_errorbarh(aes(xmin = temp-temp_se, xmax = temp+temp_se)) +
ylab("") +
xlab("temp") +
theme_bw() +
theme(,legend.position="none", panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_x_continuous(limits=c(0, 20), expand = c(0, 0))

chla_plot = 
ggplot(data_CvE, aes(x = chla_mon, y = calc_rate_log, colour = site)) + 
geom_point(size = 3, aes(shape=site)) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
geom_errorbar(aes(ymin=calc_rate_log-calc_log_se, 
	ymax=calc_rate_log+calc_log_se), 
	width=.3) +
geom_errorbarh(aes(xmin = chla_mon-chlmon_se, xmax = chla_mon+chlmon_se)) +
ylab("") +
xlab("Chl_a") +
theme_bw() +
theme(,legend.position="none", panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_x_continuous(limits=c(0, 10), expand = c(0, 0)) 

ph_plot = 
ggplot(data_CvE, aes(x = ph, y = calc_rate_log, colour=site)) + 
geom_point(size = 3, aes(shape=site)) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
geom_errorbar(aes(ymin=calc_rate_log-calc_log_se, 
	ymax=calc_rate_log+calc_log_se), 
	width=.02) +
geom_errorbarh(aes(xmin = ph-ph_se, xmax = ph+ph_se)) +
ylab("") +
xlab("pH") +
theme_bw() +
theme(,legend.position="none", panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) 

ca_plot = 
ggplot(data_CvE, aes(x = ca, y = calc_rate_log, colour = site)) + 
geom_point(size = 3, aes(shape=site)) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
geom_errorbar(aes(ymin=calc_rate_log-calc_log_se, 
	ymax=calc_rate_log+calc_log_se), 
	width=.2) +
geom_errorbarh(aes(xmin = ca-ca_Se, xmax = ca+ca_Se)) +
ylab("") +
xlab("Calcium") +
theme_bw() +
theme(,legend.position="none", panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_x_continuous(limits=c(0, 8), expand = c(0, 0))

co3_plot = 
ggplot(data_CvE, aes(x = co3, y = calc_rate_log, colour=site)) + 
geom_point(size = 3, aes(shape=site)) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
geom_errorbar(aes(ymin=calc_rate_log-calc_log_se, 
	ymax=calc_rate_log+calc_log_se), 
	width=2) +
geom_errorbarh(aes(xmin = co3-co3_se, xmax = co3+co3_se)) +
ylab("") +
xlab("CO3") +
theme_bw() +
theme(,legend.position="none", panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black"))+ 
scale_x_continuous(limits=c(-5, 120), expand = c(0, 0))

hco3_plot = 
ggplot(data_CvE, aes(x = hco3, y = calc_rate_log, colour=site)) + 
geom_point(size = 3, aes(shape=site)) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
geom_errorbar(aes(ymin=calc_rate_log-calc_log_se, 
	ymax=calc_rate_log+calc_log_se), 
	width=50) +
geom_errorbarh(aes(xmin = hco3-hco3_se, xmax = hco3+hco3_se)) +
ylab("") +
xlab("HCO3") +
theme_bw() +
theme(,legend.position="none", panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_x_continuous(limits=c(0, 2500), expand = c(0, 0))

omega_plot = 
ggplot(data_CvE, aes(x = omega, y = calc_rate_log, colour=site)) + 
geom_point(size = 3, aes(shape=site)) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
geom_errorbar(aes(ymin=calc_rate_log-calc_log_se, 
	ymax=calc_rate_log+calc_log_se), 
	width=0.05) +
geom_errorbarh(aes(xmin = omega-omega_se, xmax = omega+omega_se)) +
ylab("") +
xlab("omega") +
theme_bw() +
theme(,legend.position="none", panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_x_continuous(limits=c(0, 1.8), expand = c(0, 0)) 

ESIR_plot = 
ggplot(data_CvE, aes(x = ESIR, y = calc_rate_log, colour=site)) + 
geom_point(size = 3, aes(shape=site)) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
geom_errorbar(aes(ymin=calc_rate_log-calc_log_se, 
	ymax=calc_rate_log+calc_log_se), 
	width=0.05) +
geom_errorbarh(aes(xmin = ESIR-ESIR_se, xmax = ESIR+ESIR_se)) +
ylab("") +
xlab("ESIR") +
theme_bw() +
theme(,legend.position="none", panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_x_continuous(limits=c(0, 1.8), expand = c(0, 0)) 

Calc_v_env = ggdraw(xlim = c(0, 3.1), ylim = c(0, 3)) +
draw_plot(sal_plot, x = 0, y = 2, width = 1, height = 1) +
draw_plot(temp_plot, x = 1, y = 2, width = 1, height = 1) +
draw_plot(chla_plot, x = 2, y = 2, width = 1, height = 1) +
draw_plot(ph_plot, x = 0, y = 1, width = 1, height = 1) +
draw_plot(ca_plot, x = 1, y = 1, width = 1, height = 1) +
draw_plot(co3_plot, x = 2, y = 1, width = 1, height = 1) +
draw_plot(hco3_plot, x = 0, y = 0, width = 1, height = 1) +
draw_plot(omega_plot, x = 1, y = 0, width = 1, height = 1) +
draw_plot(ESIR_plot, x = 2, y = 0, width = 1, height = 1) 

pdf("calc_v_enviro.pdf", width = 9,height = 9)

Calc_v_env

dev.off()


########### Expontential decay model hco3 and ca2+

EM_calcium = nlsLM(shell_growth ~ 1 + Vm * (1-exp(-calcium/K)),
          start=list(Vm=0., K=32), data=data_S_O)
summary(EM_calcium)

EM_bicarb = nlsLM(shell_growth ~ 1 + Vm * (1-exp(-bicarbonate/K)),
          start=list(Vm=5, K=32), data=data_S_O)
summary(EM_bicarb)

############ plot model over data

data_S_O = read.table("SIR_omega_calcification.txt", sep = "\t", header = T)
data_S_O

data_S_O$salinity = as.character(data_S_O$salinity)
data_S_O$salinity = factor(data_S_O$salinity, levels = c("16", "11", "6"))

hco3 = 
ggplot(data_S_O, aes(x = bicarbonate, y = shell_growth, colour = salinity)) + 
geom_point(size = 2, aes(shape = exp)) +
#geom_smooth(method="nls", fullrange=TRUE, 
#            formula=y~1+Vmax*(1-exp(-x/tau)),
#            method.args = list(start=c(tau=1.0,Vmax=30)),
#            se=FALSE, colour="gray") +
labs(y=expression(Calcification ~ rate ~ (?g ~ day^{-1}  ))) +
xlab("aragonite") +
theme_bw() +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
theme(legend.background = element_rect(fill="white", size=.5, linetype="dotted"), 
legend.justification=c(1,0), 
legend.position="none", 
legend.title=element_blank(),
legend.text = element_text(size = 8),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_x_continuous(limits=c(0, 2200), expand = c(0, 0)) +
scale_y_continuous(limits=c(-1, 50), expand = c(0, 0)) 

ca = 
ggplot(data_S_O, aes(x = calcium, y = shell_growth, colour = salinity)) + 
geom_point(size = 2, aes(shape = exp)) +
#geom_smooth(method="nls", fullrange=TRUE, 
#            formula=y~1+Vmax*(1-exp(-x/tau)),
#            method.args = list(start=c(tau=1.0,Vmax=30)),
#            se=FALSE, colour="gray") +
labs(y=expression(Calcification ~ rate ~ (?g ~ day^{-1}  ))) +
xlab("ESIR") +
theme_bw() +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
theme(legend.background = element_rect(fill="white", size=.5, linetype="dotted"), 
legend.justification=c(1,0), 
legend.position="none", 
legend.title=element_blank(),
legend.text = element_text(size = 8),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_x_continuous(limits=c(0, 6), expand = c(0, 0)) +
scale_y_continuous(limits=c(-1, 50), expand = c(0, 0)) 

SIR = 
ggplot(data_S_O, aes(x = SIR, y = shell_growth, colour = salinity)) + 
geom_point(size = 2, aes(shape = exp)) +
#geom_smooth(method="nls", fullrange=TRUE, 
#            formula=y~1+Vmax*(1-exp(-x/tau)),
#            method.args = list(start=c(tau=1.0,Vmax=30)),
#            se=FALSE, colour="gray") +
labs(y=expression(Calcification ~ rate ~ (?g ~ day^{-1}  ))) +
xlab("SIR") +
theme_bw() +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
theme(legend.background = element_rect(fill="white", size=.5, linetype="dotted"), 
legend.justification=c(1,0), 
legend.position="none", 
legend.title=element_blank(),
legend.text = element_text(size = 8),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_x_continuous(limits=c(0, 0.35), expand = c(0, 0)) +
scale_y_continuous(limits=c(-1, 50), expand = c(0, 0)) 


sub = ggdraw(xlim = c(0, 3), ylim = c(0, 3)) +
draw_plot(SIR, x = 0, y = 0, width = 1.5, height = 3) +
draw_plot(ca, x = 1.5, y = 0, width = 1.5, height = 3) 

pdf("SIR_plot.pdf", width = 5.5, height = 5)

sub

dev.off()

Omega = 
ggplot(data_S_O, aes(x = omega_ara, y = shell_growth, colour = salinity)) + 
geom_point(size = 2, aes(shape = exp)) +
geom_smooth(method="nls", fullrange=TRUE, 
            formula=y~1+Vmax*(1-exp(-x/tau)),
            method.args = list(start=c(tau=1.0,Vmax=30)),
            se=FALSE, colour="gray") +
labs(y=expression(Calcification ~ rate ~ (?g ~ day^{-1}  ))) +
xlab("aragonite") +
theme_bw() +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
theme(legend.background = element_rect(fill="white", size=.5, linetype="dotted"), 
legend.justification=c(1,0), 
legend.position="none", 
legend.title=element_blank(),
legend.text = element_text(size = 8),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_x_continuous(limits=c(0, 2), expand = c(0, 0)) +
scale_y_continuous(limits=c(-1, 50), expand = c(0, 0)) 

ESIR = 
ggplot(data_S_O, aes(x = ESIR, y = shell_growth, colour = salinity)) + 
geom_point(size = 2, aes(shape = exp)) +
geom_smooth(method="nls", fullrange=TRUE, 
            formula=y~1+Vmax*(1-exp(-x/tau)),
            method.args = list(start=c(tau=1.0,Vmax=30)),
            se=FALSE, colour="gray") +
labs(y=expression(Calcification ~ rate ~ (?g ~ day^{-1}  ))) +
xlab("ESIR") +
theme_bw() +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
theme(legend.background = element_rect(fill="white", size=.5, linetype="dotted"), 
legend.justification=c(1,0), 
legend.position="none", 
legend.title=element_blank(),
legend.text = element_text(size = 8),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_x_continuous(limits=c(0, 1.5), expand = c(0, 0)) +
scale_y_continuous(limits=c(-1, 50), expand = c(0, 0)) 

SIR = ggdraw(xlim = c(0, 3), ylim = c(0, 3)) +
draw_plot(Omega, x = 0, y = 0, width = 1.5, height = 3) +
draw_plot(ESIR, x = 1.5, y = 0, width = 1.5, height = 3) 

pdf("SIR_plots.pdf", width = 5.5, height = 5)

SIR

dev.off()


######################## Cmax and K paramaters confidence intervals

data_mp = read.table("mp1.txt", header = TRUE)
data_mp

cmax = 
ggplot(data_mp, aes(x = treat, y = cmax)) +
geom_point(size=4) +
geom_errorbar(aes(ymin=cmax-cmax_ci, ymax=cmax+cmax_ci), width=.1, colour = "gray61") +
theme_bw() + 
theme(legend.background = element_rect(fill="white", size=.5), legend.position="none", 
        panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
        panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(),
        axis.text=element_text(size=10, colour = "black")) +
  scale_y_continuous(limits=c(0, 50), expand = c(0, 0)) +
  scale_x_discrete() +
  ylab("") +
  xlab("")

K1 = 
ggplot(data_mp, aes(x = treat, y = K)) +
geom_point(size=4) +
geom_errorbar(aes(ymin=K-K_ci, ymax=K+K_ci), width=.1, colour = "gray61") +
theme_bw() + 
theme(legend.background = element_rect(fill="white", size=.5), legend.position="none", 
        panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), 
        panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(),
        axis.text=element_text(size=10, colour = "black")) +
  scale_y_continuous(limits=c(0, 0.8), expand = c(0, 0)) +
  scale_x_discrete() +
  ylab("") +
  xlab("")


cmax_K = ggdraw(xlim = c(0, 1.2), ylim = c(0, 0.55)) +
draw_plot(cmax, x = 0, y = 0, width = 0.55, height = 0.55) +
draw_plot(K1, x = 0.6, y = 0, width = 0.55, height = 0.55) 

pdf("cmax_K_ESIR_omega.pdf", width = 6, height = 4)

cmax_K

dev.off()


#################################### ESIR ~ Omega

datac = read.table("Carb_chem_all.txt", header = TRUE)
datac

datac$site = as.factor(datac$site)
datac$site = as.character(datac$site)
datac$site = factor(datac$site, levels = c("Kiel", "Ahrenshoop", "Usedom"))


lm10 = lm(omega_ara ~ ex_SIR, data=datac)
summary(lm10)

Omega_ESIR_plot = 
ggplot(datac, aes(x = omega_ara, y =ex_SIR, colour = site)) + 
geom_point(size = 3, aes(x = omega_ara, y =ex_SIR, shape = site)) +
geom_smooth(method = "lm", formula = y ~ x, se = T,  size = .5, colour = "black") +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
xlab("Omega") + 
ylab("ESIR") +
theme_bw() + 
theme(legend.background = element_rect(fill="white", size=.5, linetype="dotted"), 
legend.justification=c(1,0), 
legend.position=c(0.4, 0.7), 
legend.title=element_blank(),
legend.text = element_text(size = 11),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_x_continuous(limits=c(0, 2), expand = c(0, 0)) +
scale_y_continuous(limits=c(0, 2), expand = c(0, 0))

pdf("Omega_ESIR_plot.pdf", width = 4, height = 4)

Omega_ESIR_plot

dev.off()



#################################### AT S relationship

descdist(data_S_TA$AT_meas, discrete = F)
fit.norm = fitdist(data_S_TA$AT_meas, "norm")
plot(fit.norm)
shapiro.test(log(data_S_TA$AT_meas))

bartlett.test(AT_meas ~ site, data = data_S_TA)


leveneTest(data_S_TA$AT_meas ~ as.factor(site), data=data_S_TA)
lm_STA = lm(AT_meas ~ salinity, data = data_S_TA)
anova(lm_STA)
summary(lm_STA)

data_S_TA = read.table("Carb_chem_all.txt", header = TRUE)
data_S_TA
dataSTA_m = read.table("Mueller_Basin_data.txt", header = TRUE)
dataSTA_m

dataSTA_m$Basin = as.character(dataSTA_m$Basin)
# dataSTA_m$Basin = factor(dataSTA_m$site, 
# levels = c("Kat", "Both", "Fin", "Riga"))

data_S_TA$site = as.character(data_S_TA$site)
data_S_TA$site = factor(data_S_TA$site, levels = c("Kiel", "Ahrenshoop", "Usedom"))

AT_S_field = 
ggplot(data_S_TA, aes(x = salinity, y =AT_meas, colour = site)) + 
geom_point(size = 4, aes(x = salinity, y =AT_meas, shape = site)) +
geom_line(data=dataSTA_m, aes(x = salinity, y = AT_2014, colour = Basin), size = 1, linetype = "dashed") +
geom_smooth(method = "lm", se = T, linetype = "dashed", size = .5, 
            colour = "black", (aes(group = 1))) +
scale_color_manual(values=c("aquamarine3", "black", "black", "black", "yellow3", "black", "steelblue4")) +
xlab("Salinity") + 
#ylab(expression(italic(C)[T]  ~ (?mol ~ kg^{-1}))) +
theme_bw() + 
guides(shape = FALSE) +
theme(legend.background = element_rect(fill="white", size=.5, linetype="dotted"), 
legend.justification=c(1,0), 
legend.position=c(0.9, 0.01), 
legend.title=element_blank(),
legend.text = element_text(size = 6),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_x_continuous(limits=c(0, 30), expand = c(0, 0)) +
scale_y_continuous(limits=c(0, 3000), expand = c(0, 0))



pdf("field_salinity_TA_plot.pdf", width = 6, height = 4)

AT_S_field

dev.off()

################## Bicarbonate experiment mortality rates

data_mort = read.table("bicarb_mortality.txt", header=TRUE)
data_mort

######stats

shapiro.test(data_mort$percent_mort)
bartlett.test(percent_mort ~ as.factor(salinity), data=data_mort)

aov_mort = aov(percent_mort ~ bicarb * as.factor(salinity), data=data_mort)
anova(aov_mort)

TukeyHSD(aov_mort, "as.factor(salinity)", conf.level = 0.95)

######graph

data_mort$salinity = as.factor(data_mort$salinity)
data_mort$salinity = as.character(data_mort$salinity)
data_mort$salinity = factor(data_mort$salinity, levels = c("16", "11", "6"))

exp_mort = 
ggplot(data_mort, aes(x = bicarb, y =percent_mort, colour = salinity)) + 
geom_point(size = 2) +
geom_smooth(method = "lm", se = F, linetype = "dashed", size = .5) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
xlab("") + 
ylab("") +
theme_bw() + 
theme(legend.background = element_rect(fill="white", size=.5, linetype="dotted"), 
legend.justification=c(1,0), 
legend.position=c(0.2, 0.1), 
legend.title=element_blank(),
legend.text = element_text(size = 8),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_x_continuous(limits=c(0, 2400), expand = c(0, 0)) +
scale_y_continuous(limits=c(0, 100), expand = c(0, 0))

pdf("bicarbonate_experiment_mortality.pdf", width = 5, height = 5)

exp_mort

dev.off()

###################################### SIR Omega

data_S_O = read.table("SIR_omega_calcification.txt", sep = "\t", header = T)
data_S_O

data_S_O$exp = as.character(data_S_O$exp)
data_SI$salinity = factor(data_SI$salinity, levels = c("16", "11", "6"))

omega = ggplot(data_S_O, aes(x = omega_ara, y = shell_growth, colour = salinity)) + 
geom_point(size = 2) +
stat_smooth(method='lm', formula = y ~ log2(x), se = T, size = 1.2, colour = "black") +
scale_color_manual(values=c("gray0", "gray47", "gray75")) +
labs(y=expression(Calcification ~ rate ~ (?g ~ day^{-1}  ))) +
xlab("? aragonite") +
theme_bw() +
theme(legend.background = element_rect(fill="white", 
   size=.5, linetype="solid", colour = "black"), 
   legend.justification=c(1,0), 
   legend.position="none", 
   legend.title=element_blank(),
   legend.text = element_text(size = 9),
   panel.grid.major.x=element_blank(), 
   panel.grid.minor.x=element_blank(), 
   panel.grid.minor.y=element_blank(), 
   panel.grid.major.y=element_blank(),
   axis.text=element_text(size=10, colour = "black"),
   plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm") 
   )+
scale_x_continuous(limits=c(0, 1.5), expand = c(0, 0)) +
scale_y_continuous(limits=c(-5, 50), expand = c(0, 0)) 


SIR = ggplot(data_SI, aes(x = sub_inh, y = SM, colour = salinity)) + 
geom_point(size = 2) +
stat_smooth(method='lm', formula = y ~ log2(x), se = T, size = 1.2, colour = "black") +
scale_color_manual(values=c("cornflowerblue", "forestgreen", "firebrick3")) +
#labs(y=expression(Calcification ~ rate ~ (?g ~ day^{-1}  ))) +
xlab("SIR") +
theme_bw() +
theme(legend.background = element_rect(fill="white", 
   size=.5, linetype="solid", colour = "black"), 
   legend.justification=c(1,0), 
   legend.position=c(0.9, 0.1), 
   #legend.title=element_blank(),
   legend.text = element_text(size = 10),
   panel.grid.major.x=element_blank(), 
   panel.grid.minor.x=element_blank(), 
   panel.grid.minor.y=element_blank(), 
   panel.grid.major.y=element_blank(),
   axis.title.y=element_blank(),
   axis.text=element_text(size=10, colour = "black"),
   plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm") 
   )+
scale_x_continuous(limits=c(0, 1.5), expand = c(0, 0)) +
scale_y_continuous(limits=c(-5, 50), expand = c(0, 0)) 


O_SIR = ggdraw(xlim = c(0, 2), ylim = c(0, 1)) +
draw_plot(omega, x = 0, y = 0, width = 1, height = 1) +
draw_plot(SIR, x = 1, y = 0, width = 1, height = 1) 

pdf("O_SIR_plots_0527.pdf", width = 9, height = 3)

O_SIR

dev.off()



dev.off()



########################## SUPPLEMENTARY: salinity-calcium relationships

data_SCa = read.table("salinity_calcium_relationships.txt", sep = "", header = T)
data_SCa

data_SCa$type = as.character(data_SCa$type)

pdf("salinity_calcium_relationships.pdf", width = 4, height = 4)

ggplot(data_SCa, aes(x = salinity, y = Ca, colour = type)) + 
  geom_point(size = 3) +
  stat_smooth(method='lm', formula = y ~ x, se = F, size = 1, linetype = "dashed") +
  scale_color_manual(values=c("black", "firebrick3")) +
  ylab(expression(Calcium ~ conc. ~ (mmol ~ kg^{-1}  ))) +
  xlab("salinity") +
  theme_bw() +
  theme(legend.background = element_rect(fill="white", 
                                         size=.5, linetype="solid", colour = "black"), 
        legend.justification=c(1,0), 
        legend.position=c(0.9, 0.1), 
        #legend.title=element_blank(),
        legend.text = element_text(size = 10),
        panel.grid.major.x=element_blank(), 
        panel.grid.minor.x=element_blank(), 
        panel.grid.minor.y=element_blank(), 
        panel.grid.major.y=element_blank(),
#        axis.title.y=element_blank(),
        axis.text=element_text(size=10, colour = "black"),
        plot.margin=unit(c(0.4,0.4,0.4,0.4),"cm") 
  )+
  scale_x_continuous(limits=c(0, 20), expand = c(0, 0)) +
  scale_y_continuous(limits=c(0, 6), expand = c(0, 0)) 

dev.off()

###################################### Field Temperature

install.packages("dunn.test")

descdist(log(data_T$value[1:5000]), discrete = F)
fit.norm = fitdist(data_T$value, "norm")
plot(fit.norm)

kruskal.test(value ~ site, data = data_T)
dunnTest(data_T$value, g = data_T$site, kw=TRUE, label=TRUE, wrap=FALSE, alpha=0.05)


data_T = read.table("Temperature_all_2015_2018.txt", sep = "\t", header = TRUE)
data_T

data_T$date = as.Date(data_T$date, "%d.%m.%Y")
data_T$date = as.Date(data_T$date, "%Y-%m-%d")

data_T$site = as.character(data_T$site)
data_T$site = factor(data_T$site, levels = c("Kiel", "Ahrenshoop", "Usedom"))

#pdf("Temperature_time_series_2015-2018.pdf", width = 6, height = 3)

temp_all_sp = 
ggplot(data_T, aes(x = date, y = value, colour = site)) + 
geom_point(size = .1) +
#stat_smooth(method = "gam", formula = y ~ s(x, bs = "ps", k = 12), size = 1, se = F) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
xlab("Date") + 
ylab("Temperature (?C)") +
theme_bw() + 
theme(legend.background = element_rect(fill="white", size=.5, 
linetype="dotted"), 
legend.justification=c(1,0), 
legend.position="none", 
legend.title=element_blank(),
legend.text = element_text(size = 8),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_x_date(breaks=seq(from=as.Date("2015-8-01"), to=as.Date("2017-12-01"), by="4 month"), 
labels=date_format("%b. '%y"), limits=c(), expand = c(0, 0)) +
scale_y_continuous(limits=c(-5, 25), expand = c(0, 0))


temp_all_bp = 
ggplot(data_T, aes(x = site, y = value, colour = site)) + 
geom_boxplot(size = .5, outlier.size = .05) +
stat_smooth(method = "gam", formula = y ~ s(x, bs = "ps", k = 12), size = 1.5, se = F) +
scale_color_manual(values=c("yellow3", "aquamarine3", "steelblue4")) +
xlab("Site") + 
ylab("") +
theme_bw() + 
theme(legend.background = element_rect(fill="white", size=.5, 
linetype="dotted"), 
legend.justification=c(1,0), 
legend.position="none", 
legend.title=element_blank(),
legend.text = element_text(size = 8),
panel.grid.major.x=element_blank(), 
panel.grid.minor.x=element_blank(), 
panel.grid.minor.y=element_blank(), 
panel.grid.major.y=element_blank(),
axis.text=element_text(size=10, colour = "black")) +
scale_y_continuous(limits=c(-5, 25), expand = c(0, 0))

T_plots = ggdraw(xlim = c(0, 3), ylim = c(0, 1)) +
draw_plot(temp_all_sp, x = 0, y = 0, width = 2, height = 1) +
draw_plot(temp_all_bp, x = 2, y = 0, width = 1, height = 1) 

pdf("T_plots.pdf", width = 9, height = 3)

T_plots

dev.off()


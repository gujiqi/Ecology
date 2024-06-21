#1.glmm.hp
library(lmerTest)
library(readxl)
library(lme4)
library(vegan)
library(Matrix)
library(ade4)
library(nlme)
library(lme4)
library(GGally)
library(lmerTest)
library(effectsize)
library(MuMIn)
library(MuMIn)
library(glmm.hp)
setwd("C:/Users/Google/Documents/R/TibetPlot08_19-ENV2024")
data2 <- read.csv("env结构方程画图-土生 - 20240417.CSV", header = TRUE)
head(data2)
data <- read.csv("env结构方程画图-土生 - 森林.CSV", header = TRUE)
glmm1 <- lm(Richness ~ Annual.Mean.Temperature + Mean.Diurnal.Range + Annual.Precipitation + SR + 
              Sand.content + pH + Soil.water,data = data)
data3 <- read.csv("env结构方程画图-土生 - 草原.CSV", header = TRUE)
glmm3 <- lm(Richness ~ Annual.Mean.Temperature + Mean.Diurnal.Range + Annual.Precipitation + SR + 
              Sand.content + pH + Soil.water,data = data3)
data4 <- read.csv("env结构方程画图-土生 - 灌丛.CSV", header = TRUE)
glmm4 <- lm(Richness ~ Annual.Mean.Temperature + Mean.Diurnal.Range + Annual.Precipitation + SR + 
              Sand.content + pH + Soil.water,data = data4)
b<- glmm.hp(glmm1,type="R2")
d<- glmm.hp(glmm3,type="R2")
e<- glmm.hp(glmm4,type="R2")
summary(glmm1)
r.squaredGLMM(glmm1) 
summary(glmm3)
r.squaredGLMM(glmm3) 
summary(glmm4)
r.squaredGLMM(glmm4) 



#2.EcoSimR co-occurrence patterns
library(EcoSimR)
setwd("C:/Users/Google/Documents/R/TibetPlot08_19-ENV")
comm_ENV1 <- read.csv("forests.csv",row.names = 1)
comm_ENV1[comm_ENV1>0] <- 1
comm_ENV1 <- as.data.frame(t(comm_ENV1))
comm_ENV2 <- read.csv("Thicket.csv",row.names = 1)
comm_ENV2[comm_ENV2>0] <- 1
comm_ENV2 <- as.data.frame(t(comm_ENV2))
comm_ENV3 <- read.csv("Herbaceous vegetation.csv",row.names = 1)
comm_ENV3[comm_ENV3>0] <- 1
comm_ENV3 <- as.data.frame(t(comm_ENV3))
#3.RDA
library(vegan)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(ggsci)
library(plyr)
library(cowplot)
library(RColorBrewer)
display.brewer.all()
library(ggsci)
library(pacman)
library(vegan)
library(ggplot2)
library(ggforce)
library(ggrepel)

p_load(ggplot2,patchwork,vegan,geosphere,psych,corrplot,
       permute,lattice,ggpubr,RColorBrewer,tidyverse,graphics)
setwd("C:/Users/Google/Documents/R/TibetPlot08_19-ENV2024")
env <- read.csv("env.CSV", header = TRUE)
head(env)
data=read.csv("spe0819soil.csv",row.names = 1)
head(data,n=3)
head(env,n=3)
data <- decostand(data, method = 'hellinger')
print(decorana(data))
B.rda=rda(data,env,scale = T)#RDA分析，如果没有环境因子参数，就是PCA分析
env1 <- read.csv("envsoil.CSV", header = TRUE)
head(env1,n=3)
Elevations <- env1$Vegetation.types
B.rda.data=data.frame(B.rda$CCA$u[,1:2],Elevations)#为了仿师兄的图添加的
colnames(B.rda.data)=c("RDA1","RDA2","Elevations")
colnames(B.rda.data)
head(B.rda.data,n=3)
B.rda.spe=data.frame(B.rda$CCA$v[,1:2])
B.rda.spe=as.data.frame(B.rda.spe)
B.rda.spe$Species<-rownames(B.rda.spe)
head(B.rda.spe,n=3)
B.rda.env <- B.rda$CCA$biplot[,1:2]
B.rda.env <- as.data.frame(B.rda.env)
head(B.rda.env,n=9)
B.rda$CCA$eig[1]/sum(B.rda$CCA$eig)*100
B.rda$CCA$eig[2]/sum(B.rda$CCA$eig)*100
gg3 <- ggplot(data = B.rda.data, aes(RDA1, RDA2)) +
  geom_point(aes(colour = Elevations), size = 3, alpha = 1) + 
  scale_shape_manual(values = c(19:24)) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_jco() +
  stat_ellipse(aes(color = Elevations), level = 0.95, linetype = 1, show.legend = FALSE)+
  geom_segment(aes(x = 0, y = 0, xend = B.rda.env[,1]/4, yend = B.rda.env[,2]/4), 
               data = B.rda.env, size = 1, alpha = 0.5, colour = "grey30",
               arrow = arrow(angle = 35, length = unit(0.3, "cm"))) +
  labs(
    x = paste("RDA1", round(B.rda$CCA$eig[1] / sum(B.rda$CCA$eig) * 100, 2), "%"),
    y = paste("RDA2", round(B.rda$CCA$eig[2] / sum(B.rda$CCA$eig) * 100, 2), "%"),
    colour = "Elevations"
  ) +
  geom_text(data = B.rda.env, aes(x = B.rda.env[,1]/4, y = B.rda.env[,2]/4,
                                  label = rownames(B.rda.env)), size = 5, colour = "grey30", 
            fontface = "bold", hjust = (1 - sign(B.rda.env[,1])) / 2,
            angle = (180 / pi) * atan(B.rda.env[,2] / B.rda.env[,1])) +
  coord_cartesian(xlim = c(-0.3, 0.3), ylim = c(-0.3, 0.3)) +
  theme(axis.title = element_text(size = 15, face = "bold", colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"), 
        axis.ticks.length = unit(-0.1, "lines"), legend.key = element_blank(), axis.text = element_text(size = 15),
        legend.title = element_text(size = 15, face = "bold", colour = "black"), legend.position = c(2, 3), 
        legend.text = element_text(size = 15, colour = "black"))

gg3


# 4.plspm 
library(plspm)
setwd("C:/Users/Google/Documents/R/TibetPlot08_19-ENV")
dat <- read.csv("群落稳定性结构方程.CSV", header = TRUE)
###########气候土壤水
#您可以直接指定列名称，或者指定列的下标都可以，我个人习惯指定列名称
dat_blocks <- list(
  climate = c('Annual.Mean.Temperature', 'Annual.Precipitation','SR','Mean.Diurnal.Range'), 
  soil = c('pH', 'Soil.water', 'Sand.content'), 
  community = 'icv'
)
dat_blocks
climate <- c(0, 0, 0)
soil <- c(1, 0, 0)
community <- c(1, 1, 0)

dat_path <- rbind(climate, soil, community)
colnames(dat_path) <- rownames(dat_path)
dat_path
dat_modes <- rep('A', 3)
dat_modes
dat_pls <- plspm(dat, dat_path, dat_blocks, modes = dat_modes)
dat_pls
summary(dat_pls)

###############5.mantel text
library(dplyr) 
library(ggplot2) 
library(linkET)  
library(devtools)
setwd("C:/Users/Google/Documents/R/TibetPlot08_19-ENV2024")
varespec=read.csv("spe0819soil.csv",row.names = 1)#读入物种(以Phylum水平为例)矩阵表
varechem=read.csv("env结构方程manteltext.CSV",row.names = 1)#读入环境因子数据(示例为随机数)
mantel <- mantel_test(varespec, varechem,
                      spec_select = list(Community = 1:324)) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.1, 0.2, Inf),
                  labels = c("< 0.1", "0.1-0.2", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
qcorrplot(correlate(varechem), type = "upper",fixed = TRUE,    parse = FALSE, drop = TRUE,diag = TRUE) +
  geom_square() +
  geom_couple(aes(colour = pd, size = rd), 
              data = mantel, 
              curvature = nice_curvature()) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))

#6.neutral model
#6.1.neutral model森林
library(Hmisc)
library(minpack.lm)
library(stats4)
setwd("C:/Users/Google/Documents/R/TibetPlot08_19-ENV")
spp=read.csv("样方矩阵-森林.csv",row.names = 1)
N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  #获取 m 值
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr 
bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'
library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 
draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N,3)), x=x[j], y=y[i], just=just)
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)


#6.2.neutral model灌木
library(Hmisc)
library(minpack.lm)
library(stats4)
setwd("C:/Users/Google/Documents/R/TibetPlot08_19-ENV")
spp=read.csv("样方矩阵-灌木.csv",row.names = 1)

N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  
bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'
library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 
draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N,3)), x=x[j], y=y[i], just=just)
  
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)
#6.3.neutral model草地
library(Hmisc)
library(minpack.lm)
library(stats4)
setwd("C:/Users/Google/Documents/R/TibetPlot08_19-ENV")
spp=read.csv("样方矩阵-草地.csv",row.names = 1)
N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit 
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  
bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'
library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 
draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N,3)), x=x[j], y=y[i], just=just)
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)
#7.mst
library(NST)
library(vegan)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(ggsci)
library(plyr)
setwd("C:/Users/Google/Documents/R/TibetPlot08_19")
comm=read.csv("spe0819soil.csv",header = TRUE)
comm <- comm[,-c(1)]
group <- tda$group
env <- read.csv("env08-19allcluster.CSV", header = TRUE)
head(env)
group <- env[,15]
group <- as.data.frame(group)
head(env)
comm.b=comm
comm.b[comm.b>0]=1
samp.ab=rowSums(comm)
prob.ab=matrix(colSums(comm),nrow=nrow(comm),ncol=ncol(comm),byrow=TRUE)
comm.rand=ab.assign(comm.b,samp.ab,prob.ab)
set.seed(123)
nst.boot
tnst1 <- tNST(comm = comm, group = group, dist.method = 'bray', null.model = 'PP', 
              abundance.weighted=TRUE, rand = 1000, nworker = 6)
nst_group1 <- tnst1$index.pair.grp
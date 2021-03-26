##################################################################################################################
# Name: plot_trajectories_20210122
#
# Purpose: Generate Figure 1 for MS
# Study: Cognition EWAS
#
# Created by Ida Karlsson (IK)
# Institute of Gerontology, Jönköping University &
# Department of Medical Epidemiology and Biostatistics (MEB) Karolinska Institutet
#
# Created:  20210122 by Ida Karlsson (IK) (based on script for BMI)
# Updated:	
#
##################################################################################################################
.libPaths("Z:/Programs/R/Packages")
library(ggplot2)
library(gridExtra)

age = seq(55, 90, len = 36)

growth_curve <- function(icpt, linear, quadratic){
  icpt + linear * ((age-65)/10) + quadratic * ((age-65)/10)^2
}


##### Processing speed
null_speed <- growth_curve(52.61, -3.43, -1.40)
meth_speed <- growth_curve(52.61+(-1.35), -3.43+(-0.38), -1.40+0.03)
speed <- as.data.frame(cbind(age, null_speed, meth_speed))
colnames(speed) <- c("Age", "Null", "Meth")

##### Spatial ability
null_spat <- growth_curve(54.70, -2.22, -0.58)
meth_spat1 <- growth_curve(54.70+1.37, -2.22+0.24, -0.58+(-0.15))
meth_spat2 <- growth_curve(54.70+(-1.71), -2.22+(-0.35), -0.58+0.10)
spat <- as.data.frame(cbind(age, null_spat, meth_spat1, meth_spat2))
colnames(spat) <- c("Age", "Null", "Meth1", "Meth2")

##### Working memory
null_wmem <- growth_curve(51.88, -1.44, 0)
meth_wmem1 <- growth_curve(51.88+(-1.61), -1.44+0.30, 0)
meth_wmem2 <- growth_curve(51.88+1.83, -1.44+(-0.85), 0)
meth_wmem3 <- growth_curve(51.88+1.41, -1.44+(-0.11), 0)

wmem <- as.data.frame(cbind(age, null_wmem, meth_wmem1, meth_wmem2, meth_wmem3))
colnames(wmem) <- c("Age", "Null", "Meth1", "Meth2", "Meth3")


########## Plots

##### Spatial ability
p1 <- ggplot() +
  geom_line(data=spat, aes(x=Age, y=Null, color="black"), linetype="longdash", size=1.1) +
  geom_line(data=spat, aes(x=Age, y=Meth1, color="blue"), size=1.1) +
  geom_line(data=spat, aes(x=Age, y=Meth2, color="magenta3"), size=1.1) +
  coord_cartesian(ylim = c(30, 60), xlim=c(55, 90)) +
  labs(x = "Age", y= "Cognitive ability", title = "Spatial ability") +
  theme(axis.title = element_text(size=20, face="bold"),
        axis.text = element_text(size=20),
        title = element_text(size=20, face="bold"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20),
        legend.position="top") +
  scale_color_identity(name = "",
                       breaks = c("black", "blue", "magenta3"),
                       labels = c("Null model", "cg04549090", "cg18064256"),
                       guide = "legend")

#ggsave("P:/Dementia_IK/Cognition/EWAS/Output/Plot_spat_20200517.png",width = 10, height = 9)

##### Processing speed
p2 <- ggplot() +
  geom_line(data=speed, aes(x=Age, y=Null, color="black"), linetype="longdash", size=1.1) +
  geom_line(data=speed, aes(x=Age, y=Meth, color="magenta3"), size=1.1) +
  coord_cartesian(ylim = c(30, 60), xlim=c(55, 90)) +
  labs(x = "Age", y= "Cognitive ability", title = "Processing speed") +
  theme(axis.title = element_text(size=20, face="bold"),
        axis.text = element_text(size=20),
        title = element_text(size=20, face="bold"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20),
        legend.position="top") +
  scale_color_identity(name = "",
                       breaks = c("black", "magenta3"),
                       labels = c("Null model", "cg18064256"),
                       guide = "legend")

#ggsave("P:/Dementia_IK/Cognition/EWAS/Output/Plot_speed_20200517.png",width = 10, height = 9)

##### Working memory
p3 <- ggplot() +
  geom_line(data=wmem, aes(x=Age, y=Null, color="black"), linetype="longdash", size=1.1) +
  geom_line(data=wmem, aes(x=Age, y=Meth1, color="blue"), size=1.1) +
  geom_line(data=wmem, aes(x=Age, y=Meth2, color="magenta3"), size=1.1) +
  geom_line(data=wmem, aes(x=Age, y=Meth3, color="cyan3"), size=1.1) +
  coord_cartesian(ylim = c(30, 60), xlim=c(55, 90)) +
  labs(x = "Age", y= "Cognitive ability", title = "Working memory") +
  theme(axis.title = element_text(size=20, face="bold"),
        axis.text = element_text(size=20),
        title = element_text(size=20, face="bold"),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20),
        legend.position="top") +
  scale_color_identity(name = "",
                       breaks = c("black", "blue", "magenta3", "cyan3"),
                       labels = c("Null model", "cg09988380", "cg25651129", "cg08011941"),
                       guide = "legend")

#ggsave("P:/Dementia_IK/Cognition/EWAS/Output/Plot_wmem_20200517.png",width = 10, height = 9)

### Combine and save
pdf(file="P:/Dementia_IK/Cognition/EWAS/Output/Figure_1.pdf",width = 25, height = 8)
plotgrid <- grid.arrange(p1, p2, p3, ncol=3, nrow=1, respect=T)
dev.off()

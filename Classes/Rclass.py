from rpy2.situation import r_home_from_registry
import os

os.environ["R_HOME"] = r_home_from_registry()

from rpy2.robjects.vectors import StrVector
import rpy2.robjects as robjects

import rpy2.robjects.packages as rpackages
import Classes.config as conf

def calculate_EC50_SE_plots(table_filename, compound, line, sheet):
    utils = rpackages.importr('utils')
    utils.chooseCRANmirror(ind=1)
    packnames = conf.packnames_for_r

    names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
    if len(names_to_install) > 0:
        utils.install_packages(StrVector(names_to_install))

    robjects.globalenv["xlsname"] = table_filename
    robjects.globalenv["compound"] = compound
    robjects.globalenv["line"] = line
    robjects.globalenv["sheet"] = sheet

    robjects.r('''
            
    library(openxlsx)
    library(magrittr)
    library(drc)
    library(ggplot2)
    library(Rmisc)
    library(data.table)
    library(sjPlot)
    
    
    
    A172_AV98 <- read.xlsx(xlsname, sheet = sheet)
    
    A172_AV98 <- rbind(data.frame(variable = rep(0,3), value = A172_AV98[,1]),
                       melt(as.data.table(A172_AV98[,-1]), measure.vars = colnames(A172_AV98)[-1]))
    
    A172_AV98$variable <- as.numeric(as.character(A172_AV98$variable))
    
    first <- ggplot(A172_AV98)+
      aes(as.factor(variable), value)+
      geom_point(position = position_dodge(.4))
    
    second <- ggplot(A172_AV98)+
      aes(as.factor(variable), value)+
      geom_boxplot()
    
    
    summarySE(A172_AV98, measurevar = 'value', groupvars = c('variable'))
    
    A172_AV98$vitality <- A172_AV98$value * 100 / 0.86350000
    
    summarySE(A172_AV98, measurevar = 'vitality', groupvars = c('variable'))
    
    dr <- drm(vitality ~ variable, data = A172_AV98,
        fct = LL.4(fixed=c(NA, 0, NA, NA),
                   names=c("Slope","Lower Limit","Upper Limit", "EC50"))) %>%
      ED(. , 50)
    
    A172_AV98_curve <- drm(vitality ~ variable, data = A172_AV98,
                           fct = LL.4(fixed=c(NA, 0, NA, NA),
                                      names=c("Slope","Lower Limit","Upper Limit", "EC50"))) %>%
      plot(.) %>%
      as.data.frame()
    
    A172_AV98_sd <- summarySE(A172_AV98, measurevar = 'vitality', groupvars = 'variable') %>%
      subset(. , variable > 0)
    
    third <- ggplot(subset(A172_AV98, variable > 1))+
      aes(variable, vitality)+
      scale_x_log10(breaks = c(1,2,4,8,16,32,64,128), limits = c(0.5,172))+
      scale_y_continuous(breaks = c(0,25,50,75,100))+
      geom_errorbar(data=A172_AV98_sd, aes(ymin=vitality-sd, ymax=vitality+sd), 
                    width=.2) +
      stat_summary(geom = 'point', size = 1)+
      geom_line(aes(variable, `1`), lty = 2,
                data=A172_AV98_curve)+
      ylab('Cell vitality, %')+
      xlab("Concentration, \\u03bcM")+
      ggtitle(paste("Cells ", line, " - compound ", compound))
      
    save_plot("./external/resources/plots/plot1.svg", fig = second, width=10, height=8)
    save_plot("./external/resources/plots/plot2.svg", fig = third, width=10, height=8)
    dev.off()
    
            ''')

    r_f = robjects.globalenv['dr']
    for i in list(r_f):
        print(i)
    return r_f


if __name__ == "__main__":
    calculate_EC50_SE_plots("C:/Users/vofan/PycharmProjects/Mag/A172_AV94_AV98.xlsx", "jj", "b52", 2)

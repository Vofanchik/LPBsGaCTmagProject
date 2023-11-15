import time

from rpy2.situation import r_home_from_registry
import os

from Classes.Testin import timeit

try:
    os.environ["R_HOME"] = r_home_from_registry()

except:
    os.environ["R_HOME"] = 'C:\\Program Files\\R\\R-4.3.1'

from rpy2.robjects.vectors import StrVector
import rpy2.robjects as robjects

import rpy2.robjects.packages as rpackages
import Classes.config as conf


@timeit
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
    
    
    
    compound_line <- read.xlsx(xlsname, sheet = sheet)
    
    compound_line <- rbind(data.frame(variable = rep(0,3), value = compound_line[,1]),
                       melt(as.data.table(compound_line[,-1]), measure.vars = colnames(compound_line)[-1]))
    
    compound_line$variable <- as.numeric(as.character(compound_line$variable))
        
    second <- ggplot(compound_line)+
      aes(as.factor(variable), value)+
      geom_boxplot()    
    
    neutrilize <- summarySE(compound_line, measurevar = 'value', groupvars = c('variable'))$value[1]
    
    compound_line$vitality <- compound_line$value * 100 / neutrilize
    
    summarySE(compound_line, measurevar = 'vitality', groupvars = c('variable'))
    
    dr <- drm(vitality ~ variable, data = compound_line,
        fct = LL.4(fixed=c(NA, 0, NA, NA),
                   names=c("Slope","Lower Limit","Upper Limit", "EC50"))) %>%
      ED(. , 50)
    
    compound_line_curve <- drm(vitality ~ variable, data = compound_line,
                           fct = LL.4(fixed=c(NA, 0, NA, NA),
                                      names=c("Slope","Lower Limit","Upper Limit", "EC50"))) %>%
      plot(.) %>%
      as.data.frame()      
    
    compound_line_sd <- summarySE(compound_line, measurevar = 'vitality', groupvars = 'variable') %>%
      subset(. , variable > 0)
    
    third <- ggplot(subset(compound_line, variable > 1))+
      aes(variable, vitality)+
      scale_x_log10(breaks = c(1,2,4,8,16,32,64,128), limits = c(0.5,172))+
      scale_y_continuous(breaks = c(0,25,50,75,100))+
      geom_errorbar(data=compound_line_sd, aes(ymin=vitality-sd, ymax=vitality+sd), 
                    width=.2) +
      stat_summary(geom = 'point', size = 1)+
      geom_line(aes(variable, `1`), lty = 2,
                data=compound_line_curve)+
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

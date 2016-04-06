# HMMTutorial
all .Rmd files should be loaded from R using R Markdown
Pandoc needs to be installed in your computer, check 
http://pandoc.org/installing.html

install.packages("rmarkdown")
setwd("~/your local directory/")
rmarkdown::render("1_BiparentalSim.Rmd")
rmarkdown::render("2_BiparentalSim.Rmd")
rmarkdown::render("3_HMMPathFinder.Rmd")
rmarkdown::render("4_Viterbi.Rmd")
rmarkdown::render("5_Simulations.Rmd") #this one takes ~10 minutes to render
rmarkdown::render("6_MoreSimulations.Rmd") #this one takes ~10 minutes to render


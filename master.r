rm(list = ls())
WD.PATH = paste(getwd(),"/Data-and-Code/",sep="")

# Re-produce Figure 1
source(paste(WD.PATH, 'code/fig1.r', sep=''))

# Re-produce Figure 2
source(paste(WD.PATH, 'code/fig2.r', sep=''))

# Re-generate simulation results: adjusted for double truncation
ubd = 1
a.low=rep(-ubd, 3); a.upp=rep(ubd, 3)                         # doubly
WD.PATH1 = paste(WD.PATH, 'results/doubly/', sep='')
source(paste(WD.PATH, 'code/simulation.r', sep=''))

# Re-generate simulation results: adjusted for right truncation
ubd = 1
a.low=rep(-Inf, 3); a.upp=rep(ubd, 3)                         # right
WD.PATH1 = paste(WD.PATH, 'results/right/', sep='')
source(paste(WD.PATH, 'code/simulation.r', sep=''))

# Re-generate simulation results: adjusted for left truncation
ubd = 1
a.low=rep(-ubd, 3); a.upp=rep(Inf, 3)                         # left
WD.PATH1 = paste(WD.PATH, 'results/left/', sep='')
source(paste(WD.PATH, 'code/simulation.r', sep=''))

# Re-produce simulation results: Table 2 & Table S1
source(paste(WD.PATH, 'code/SIMtab2tabS1.r', sep=''))

# Re-produce simulation results: fig3, fig4, figS1, figS2, figS3, figS4
source(paste(WD.PATH, 'code/SIMfig.r', sep=''))

library(matlab)
library(chron)
library(PBSmodelling)


setwd("~/Dropbox/2023 segundo semestre UNAL/Tomgro model/Codigo/")
rm(list = ls(all = TRUE)) # Delete all info

climfile = 'C99AF1B.txt'
RefIntClim = read.table(climfile) # Data frame

# Values for many variables
source('cropparm.R') 
source('mgt.R')

# Initial conditions
results = matrix(0, nrow = NDAYS, ncol = 33)
source('init2.R')
RCDRW = 0
SOSIR = 0 # Source-sink ratio
GRESP = 0
TRESP = 0 # Daily total respiration

# Daily loop
for(JDAY in 1:NDAYS){
  TIME = (JDAY - 1) * DELT
  DATE = mod(TIME + NSTART, 365)
  
  # Initialize variables that accumulate during the day
  GP = 0 # Gross photosynthesis
  TEMFAC = 0
  RMAINT = 0
  RDVLV = 0
  RDVFR = 0
  TEMEXP = 0
  SLAMIN = 0
  
  # Fast loop (if time-step = 1 hour, NFAST = 24)
  for(JF in 1:NFAST){
    TFAST = (JF - 1) * 24/NFAST # Hour of the day (h)
    
    # ***************************************************************
    # Input of hourly climate data 
    # --- source('IPW_AL94.R) ---
    
    d=find(RefIntClim[,1] == DATE & RefIntClim[ ,2] == (TFAST+1))
    if(length(d) > 1){print('more than 1 found!!' )}
    if(isempty(d)){print(' invalid index of climate data')}
    JUL = RefIntClim[d, 1]
    IHE = RefIntClim[d, 2]
    solrad = RefIntClim[d, 3] # Solar radiation
    TMPA = RefIntClim[d, 4] # Daily average temperature
    VPD = RefIntClim[d, 5] # Air vapor pressure deficit
    CO2L = RefIntClim[d, 6]
    PAR = solrad * TRGH # Photosynthetic active radiation
    PPFD = PAR * 4.57 # Photosynthetic photon flux density
    
    
    # ***************************************************************
    # Development rate
    # --- source('DEVFAST.R') ---
    
    # Temperature effects on:
    TEMFCF = approx(XTEM,GENTEM,TMPA)$y # Development rate
    TEMEXP1 = approx(TEMPLE,PLETEM,TMPA)$y # Expansion rate
    SLAMIN1 = approx(TEMSLA,SLAM,TMPA,yleft=0)$y # Dry matter distribution
    
    RDVLVF = approx(XLV,RDVLVT,TMPA)$y # Leaf aging
    RDVFRF = approx(XFR,RDVFRT,TMPA)$y # Fruit aging

    
    # ***************************************************************
    # Gross photosynthesis 
    # --- source('PHOTO2.R') ---
    
    GPF = 0
    TAU2 = 0.03 # Effect of CO2 on PMAX
    PMAX = TAU1 * CO2L # Maximum leaf photosynthesis
    if(CO2L>1500){PMAX = (TAU1*1500) + (TAU2*(CO2L-1500))}
    
    PMAX = PMAX * approx(TMPG,PGRED,TMPA)$y # Reduction at extreme temperatures
    
    if(PPFD >= 0.001){
      if(VPD>=VPDL){PMAX = PMAX * exp(CK*(VPD-VPDL))} # Reduction of PMAX by VPD

      # Acock's model
      TOP = (1-XM)*PMAX + QE*XK*PPFD
      
      # PLM2 Plant density
      BOT = (1-XM)*PMAX + QE*XK*PPFD*exp(-XK*PLAR2*PLM2)
      GPF = (PMAX/XK)*log(TOP/BOT)
      
      # Convert GPF from CO2 to CH2O 
      # 30/44 = 0.682
      GPF = GPF * 0.682
      
      # Convert GPF from ?M/m2-s to g/m2-day
      # ?M/m2-s x 0.000044g/?M x 3600s/h x 24h/d = 3.8016 g/m2-day
      GPF = GPF * 3.8016/PLM2
    }
    
    
    # ****************************************************************
    # MRES Maintenance respiration 
    # --- source('resp.R') ---
    
    # Effect of temperature on MRES
    # RMRL Respiration coefficient for stem and leaf tissues
    # RMRF Respiration coefficient for growing fruits
    RMAINTF = (RMRL*(TDML2+TDMS)+RMRF*DMGF) * (Q10^(0.1*TMPA-2.0))
    
    # Integration of variables on 24 hours
    TEMFAC = TEMFAC + TEMFCF*DTFAST
    TEMEXP = TEMEXP + TEMEXP1*DTFAST
    SLAMIN = SLAMIN + SLAMIN1*DTFAST
    RMAINT = RMAINT + RMAINTF*DTFAST
    RDVLV = RDVLV + RDVLVF*DTFAST
    RDVFR = RDVFR + RDVFRF*DTFAST
    GP = GP + GPF*DTFAST
  }

  # ***********************************************
  # Development and sink strength
  source('DEVRATE6.R')
  
  # Daily production of biomass
  # GREF Growth efficiency coefficient
  RCDRW = GREF*(GP+CPOOL-RMAINT)
  if(RCDRW<0){
    RCDRW = 0
    RMAINT = GP+CPOOL}
  
  RCDRW = RCDRW * (1 - approx(XROOT,PROOT,PLSTN,yright=min(PROOT))$y)
  SOSIR = min(1,RCDRW/(PNGP+EPS))
  
  # GRESP Daily growth respiration
  # aangepast door ALEX zodat SOSIR wordt in rekening gebracht
  
  # convert GREF from gDM/gCH2O to gCH2O/gCH2O
  # (30/12)*0.45 = 1.125
  GRESP = ((PNGP*SOSIR)/GREF) * (1 - GREF*1.125);
  TRESP = RMAINT + GRESP;
  
  # *********************************************** 
  # fruit setting --- source('LOSRATE5.R') ---
  TABNF = 0
  if(TDMF<EPS){dummy=0}else
  {FABOR = min(1,(0.67-ABORMX*SOSIR))
  FABOR = max(0,FABOR)
  TABNF =FABOR*TRCNF
  }
  
  TABF = TABF + TABNF
  TNSF = TNF - TABF
  
  #C location of "aborted" fruits: distal position
  B=0
  for(I in 1:NBRUP){
    # if no new frt or less than 2 frt on truss I or no more abortion...
    if(RCNF[I]==0 || XNFT[I]<=2 || B>=TABNF){dummy=0}else
    {ABNF[I] = min(4,RCNF[I],TABNF-B,XNFT[I]-2-ABOR[I])
    ABNF[I] = max(0,ABNF[I])
    B = B + ABNF[I]
    ABOR[I] = ABOR[I] + ABNF[I]
    NSF[I] = floor(XNFT[I] - ABOR[I]) 
    # attention ceci conduit ? compter 2 fruits avort?s pour ABOR = 1,1 !         
    if(ABOR[I]<1){dummy=0}else
    {for(J in c((NSF[I]+1):floor(XNFT[I]))){
      DWF[I,J] = -EPS
    }}}}
  
  # *********************************************************
  # dry matter partitioning
  source('DMRATE6.R')
  
  # *********************************************************
  # outputs
  results[JDAY,]<-c(DATE,PLSTN,NBLV,NBRU,TDML+CPOOL,TDMS+WSTMI,TDMF,
                    TDML+TDMS+TDMF+CPOOL+WSTMI,PLAR,TDML,TDML2,TDML2+CPOOL,PLAR2,
                    PLAR2*10000/((TDML2+CPOOL)/(1.0+FRPT)),
                    DMMF,DMGF,TNSF,TNF-TNSF,TNMF,TDMF/(TDML+TDMS+TDMF+CPOOL),
                    GP,RMAINT,GRESP,SOSIR,CPOOL,PTNLVS,PTNSTM,PTNFRT, NBRU, NBRUP,TABNF,TNF,TNSF)
  
  # info<-paste(DATE,'days simulated')	
  # setWinProgressBar(pb, DATE,paste('TOMGRO V.2.0 -',NDAYS,'days to be simulated'),info)
  
}
# end of daily loop
# *********************************************************************

dimnames(results)[[2]]<-c('DATE','PLSTN','NBLV','NBRU','TDML+CPOOL','TDMS','TDMF', 'TDML+TDMS+TDMF+CPOOL','PLAR','TDML','TDML2','TDML2+CPOOL','PLAR2',
                          'PLAR2*10000./((TDML2+CPOOL)/(1.0+FRPT))',
                          'DMMF','DMGF','TNSF','TNF-TNSF','TNMF','TDMF/(TDML+TDMS+TDMF+CPOOL)',
                          'GP','RMAINT','GRESP','SOSIR','CPOOL','PTNLVS','PTNSTM','PTNFRT','NBRU', 'NBRUP','TABNF','TNF','TNSF')
write.csv(results,paste(getwd(),'/output_',climfile,'.csv',sep=''),row.names=F)

results2 <- as.data.frame(results)

plot(results2$DATE, results2$PLAR2)


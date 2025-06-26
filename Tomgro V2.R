# ......................... T O M G R O  V 2.0..........................
#    C. Gary, N. Bertin, M. Tchamitchian, L. Charasse and J.F. Barczi
#                              1993-1994
# ......................................................................
#                                 from
# --------------------------T O M G R O  V 1.0--------------------------
#       MAIN PROGRAM FOR THE TOMATO MODEL
#       DEVELOPED BY E. DAYAN, J. W. JONES, H. VAN KEULEN, AND H. CHALLA
#       JULY, 1988
# ----------------------------------------------------------------------
#       with individual organ simulations,
#            a carbohydrate pool,
#            independant simulation of leaf d.wt and area increase.
#       All variable are calculated on a per plant basis.
# ......................................................................
library(matlab)
library(chron)
library(PBSmodelling)
#setwdGUI()

# Adicionar la ruta de acuerdo con el computador donde se esté trabajando
# NO Modificar esta
setwd("~/Dropbox/2023 segundo semestre UNAL/Tomgro model/Codigo/")
rm(list = ls(all = TRUE))

#RUND = format(Sys.time(), "%Y%m%dT%H%M%S")
#START <- as.chron(Sys.time())

climfile = 'C99AF1B.txt'
RefIntClim = read.table(climfile)
#name = paste(strsplit(climfile, '\\.')[[1]][1], RUND, sep = '_')

source('cropparm.R')
source('mgt.R')

#Windows progress bar
# pb <- winProgressBar("TOMGRO V.2.0 - Simulation progress", paste("Simulation will start at day ",NSTART),
#                      0, NSTART+NDAYS)
# Sys.sleep(5)

# initial conditions
results = matrix(0, nrow = NDAYS, ncol = 33)
source('init2.R')
RCDRW = 0
SOSIR = 0
GRESP = 0
TRESP = 0

# **********************************************************************
# daily loop
# **********
for(JDAY in 1:NDAYS){
	TIME = (JDAY - 1) * DELT
	DATE = mod(TIME + NSTART, 365)

# initialize variables that accumulate during the day
      GP = 0
      TEMFAC = 0
      RMAINT = 0
      RDVLV = 0
      RDVFR = 0
      TEMEXP = 0
      SLAMIN = 0

  # *********************************************************************
    # fast loop (if time-step = 1 hour, NFAST = 24)
    # *********************************************
    for(JF in 1:NFAST){
        
        # TFAST is the hour of the day (h)
        TFAST = (JF - 1) * 24/NFAST
        
        # ***************************************************************
	  # input of hourly climate data --- source('IPW_AL94.R) ---
      d=find(RefIntClim[,1] == DATE & RefIntClim[ ,2] == (TFAST+1))
		if(length(d)>1){print('more than 1 found!!' )}
		if(isempty(d)){print(' invalid index of climate data')}
		JUL = RefIntClim[d, 1]
		IHE = RefIntClim[d, 2]
		solrad = RefIntClim[d, 3]
		TMPA=RefIntClim[d, 4]
		VPD=RefIntClim[d, 5]
		CO2L=RefIntClim[d, 6]
		PAR = solrad * TRGH
	      PPFD = PAR * 4.57

	  # ***************************************************************
        # calculation of development rate --- source('DEVFAST.R') ---
        	# temperature effect on development rate
      		TEMFCF = approx(XTEM,GENTEM,TMPA)$y
		# temperature effect on expansion rate
      		TEMEXP1 = approx(TEMPLE,PLETEM,TMPA)$y
		# temperature effect on dry matter distribution
      		SLAMIN1 = approx(TEMSLA,SLAM,TMPA,yleft=0)$y
		# leaf aging
      		RDVLVF = approx(XLV,RDVLVT,TMPA)$y
		# fruit aging
      		RDVFRF = approx(XFR,RDVFRT,TMPA)$y
        
	  # ***************************************************************
	  # calculation of gross photosynthesis --- source('PHOTO2.R') ---
		      GPF = 0
			# effect of CO2 on Pmax (from Gainesville)
			TAU2 = 0.03
			PMAX = TAU1 * CO2L
			if(CO2L>1500){PMAX = (TAU1*1500) + (TAU2*(CO2L-1500))}
		
			# reduction of Pmax at extreme temperatures
			PMAX = PMAX * approx(TMPG,PGRED,TMPA)$y
	
			if(PPFD >= 0.001){
			      # reduction of Pmax by VPD
			     if(VPD>=VPDL){PMAX = PMAX * exp(CK*(VPD-VPDL))}
			    # Acock's model
			    TOP = (1-XM)*PMAX + QE*XK*PPFD
			    BOT = (1-XM)*PMAX + QE*XK*PPFD*exp(-XK*PLAR2*PLM2)
			    GPF = (PMAX/XK)*log(TOP/BOT)
  
			    # conversion from CO2 to CH2O (30/44 = 0.682)
			    GPF = GPF * 0.682
			    # conversion of GPF from ?M/m2-s into g/m2-day
			    #?M/m2-s x 0.000044g/?M x 3600s/h x 24h/d = 3.8016 g/m2-day
			    GPF = GPF * 3.8016/PLM2
			}

	  # ****************************************************************
        # calculation of maintenance respiration --- source('resp.R') ---
        	# Effect of temperature on maintenance respiration
      		RMAINTF = (RMRL*(TDML2+TDMS)+RMRF*DMGF) * (Q10^(0.1*TMPA-2.0))
        
        #integration of variables on 24 hours
        TEMFAC = TEMFAC + TEMFCF*DTFAST
        TEMEXP = TEMEXP + TEMEXP1*DTFAST
        SLAMIN = SLAMIN + SLAMIN1*DTFAST
        RMAINT = RMAINT + RMAINTF*DTFAST
        RDVLV = RDVLV + RDVLVF*DTFAST
        RDVFR = RDVFR + RDVFRF*DTFAST
        GP = GP + GPF*DTFAST
    }
    # end of fast loop
    # **********************************************************************

# ***********************************************
# development and sink strength
    source('DEVRATE6.R')

    # daily production of biomass
    RCDRW = GREF*(GP+CPOOL-RMAINT)
    if(RCDRW<0){
        RCDRW = 0
        RMAINT = GP+CPOOL}

    RCDRW = RCDRW * (1 - approx(XROOT,PROOT,PLSTN,yright=min(PROOT))$y)
    
    # calculation of source/sink ratio
    SOSIR = min(1,RCDRW/(PNGP+EPS))
    
    # daily growth respiration
    # aangepast door ALEX zodat SOSIR wordt in rekening gebracht
    	GRESP = ((PNGP*SOSIR)/GREF) * (1 - GREF*1.125);
    # 1.125 is (30/12)*0.45 to convert GREF from gDM/gCH2O to gCH2O/gCH2O
    # daily total respiration
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

# END<-as.chron(Sys.time())
# info2<-paste('Execution time was...',END-START,sep='')
# setWinProgressBar(pb, DATE,'TOMGRO V.2.0 - Simulation finished',info2)
# Sys.sleep(5)
# close(pb) 

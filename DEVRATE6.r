#       SUBROUTINE DEVRATE6
# C
#       INCLUDE'COMMON.INC'
# c version discute avec L. Charasse (hiver 94-95)
# c et L. Gauthier (3/7/95)
# C *********************************************************************

# C node number
GENR = TEMFAC*approx(XGEN,GENRAT,PLSTN)$y
PLSTN = PLSTN + GENR

# truss number
NBRU = floor((PLSTN-FTRUSN+(1+TPL)/TPL)*TPL/(1+TPL))
NBRU = max(0,NBRU)

# number of trusses bearing fruit
NBRUP=floor((PLSTN-FTRUSN-FRLG+(1+TPL)/TPL)*TPL/(1+TPL))
NBRUP=max(0,NBRUP)

# leaf number
NBLV = floor(PLSTN) - NBRU

# number of fruits (diameter>20mm) on each truss
# Note: FPN is now the ratio of fruit initiation on a truss per node initiation
TRCNF = 0

if(NBRUP==0){NBRUP_length<-1}else{NBRUP_length<-c(1:NBRUP)}
for(I in NBRUP_length){
    XX=XNFT[I]
    # FPN is now replaced by FPNPT=0.8
    XNFT[I] = min(NBFPT[I],XNFT[I]+GENR*FPNPT)
    RCNF[I]=XNFT[I]-XX
    TRCNF=TRCNF+RCNF[I]
    # rncf(i) et trcnf sont utilises dans LOSRATE
}
TNF=TNF+TRCNF
#
# *********************************************************************
# aging and sink strength (determined for leaves and stems
# and for each fruit on each reproductive unit)
PTNLVS = 0
PTNSTM = 0
PTNFRT = 0

# leaves and stems
# new leaves start aging when a new truss appears (2 leaves below
# and 1 leaf above)
for(I in 1:NBRU){
    if(AGLS[I]==-EPS){
        PLE[I] = 0
        PGL[I] = 0
        PGS[I] = 0}else
        # quand la feuille est taille, on n'affecte plus de biomasse
        # ni aux feuilles ni aux tiges de la mme UR
    	 {AGLS[I] = min(1,AGLS[I]+RDVLV)
        XBOX = 100*AGLS[I]
        # PLE(I) is potential leaf expansion rate
        # 1st RU has FTRUSN leaves, others (1.+TPL)/TPL
        # 5/10/95 : remplacement de (1.+TPL)/TPL par 1/TPL !!!!!
        if(I==1){
		PLE[I] = TEMFAC*max(0,approx(BOX,POL,XBOX)$y)*FTRUSN
	  } else if(I<3){
	  PLE[I] = TEMFAC*max(0,approx(BOX,POL,XBOX)$y)/TPL
	  } else {
	  PLE[I]=TEMEXP*max(0,approx(BOX,POL,XBOX)$y)/TPL
	  }
        PGL[I] = PLE[I]*(1+FRPT)/SLAMIN
        PGS[I] = PGL[I]*FRST
    }
    PTNLVS = PTNLVS + PGL[I]
    PTNSTM = PTNSTM + PGS[I]
} #loop

#  fruit
for(I in 1:NBRU){
    # c nouvelle version 22/11/95
    for(J in 1:floor(XNFT[I])){
        AGF[I,J] = min(1,AGF[I,J]+RDVFR)
        XBOX = 100*AGF[I,J]
        if(length(DWF[I,J])==0){PGF[I,J] = 0}else{
	  if(DWF[I,J]<0){PGF[I,J] = 0}else{
            PGF[I,J]= TEMFAC*max(0,approx(BOX,POF,XBOX)$y)
            #c        pgf(i,j) = pgf(i,j) * position(j)
            #c  position = coeff de 1.0  valeurs d?croissantes, ? ajouter dans CROPPARM.TOM
            PTNFRT = PTNFRT + PGF[I,J]
}}}}

# total sink demand
	PNGP = PTNLVS + PTNSTM + PTNFRT

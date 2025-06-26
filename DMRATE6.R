#       SUBROUTINE DMRATE6
#       INCLUDE'COMMON.INC'
# C
# C *********************************************************************
# C
ASCSP = 0
TDML = 0
TDML2 = 0
DMGL = 0
TDMS = 0
TDMF = 0
DMGF = 0
TNMF = 0
DMMF = 0
PLAR = 0
PLAR2 = 0

#C leaf dry weight
for(I in 1:NBRU){
    J=I+1
    ASC = min(PGL[I],PGL[I]*SOSIR)
    ASCSP = ASCSP + ASC
    DWL[I] = DWL[I] + ASC
    
    # dry weight of growing leaves
    if(AGLS[I]<1){DMGL = DMGL + DWL[I]}

    # Now leaf area expansion is either potential or limited by a max SLA
    # It is available for each unit
    XLA[I] = XLA[I] + min(PLE[I],(ASC*SLAMX/(1+FRPT)))
    TDML = TDML + DWL[I]
    # total leaf area
    PLAR = PLAR + XLA[I]
    # total less pruned leaves (leaves pruned when truss harvested)
    if(AGLS[I]>0){ 
        PLAR2 = PLAR2 + XLA[I]
        TDML2 = TDML2 + DWL[I]
    }
    
    # when truss harvested, also half of leaves above  are pruned
               if(AGLS[I]<0){ 
                if(AGLS[J]>0){ 
                   PLAR2 = PLAR2 - XLA[J] - XLA[J]
                   TDML2 = TDML2 - DWL[J] - DWL[J]
               }
           }
    
} #   10 CONTINUE

# stem dry weight
for(I in 1:NBRU+1){
    ASC = min(PGS[I],PGS[I]*SOSIR)
    ASCSP = ASCSP + ASC
    DWS[I] = DWS[I] + ASC
    TDMS = TDMS + DWS[I]
}   #   20 CONT

# fruit dry weight
TDMF = 0
ASC = 0
for(I in 1:NBRU){
    DWTR[I] = 0
    NFT=floor(XNFT[I])
    for(J in 1:NFT){
        # for "aborted" fruit, DWF = -EPS
	 if(J==0){dummy=0}else
	 {if(AGF[I,J]<1 & DWF[I,J]>=0){
            ASC = min(PGF[I,J],PGF[I,J]*SOSIR)
            ASCSP = ASCSP + ASC
            DWF[I,J] = DWF[I,J] + ASC
        }
        if(DWF[I,J]>0){
            DWTR[I] = DWTR[I] + DWF[I,J]
        }
        #    C mature fruit
        if(AGF[I,J]==1 & DWF[I,J]>0){
            DMMF = DMMF + DWF[I,J]
            TNMF = TNMF + 1
        }
        if(AGF[I,J]<1 & DWF[I,J]>0){
            DMGF = DMGF + DWF[I,J]
        }}
    } #    40     CONTINUE
    #C leaf pruning (when the truss is harvested)
    # the age of pruned leaves is set to -EPS to be identified
    
    if(AGF[I,NBFPT[I]]>=1){AGLS[I] = -EPS}
    TDMF = TDMF + DWTR[I]
} #50 CONTINUE

#     C carbohydrate pool
CPOOL = max(0,(RCDRW-ASCSP)/GREF)

#     C when CPOOL is higher than a threshold value then GP is limited
CPOOLMX = 0.06*TDML/((1+FRPT)*GREF)

if(CPOOL > CPOOLMX){
    #C Statement changed by ALEX, ik zie niet in wat dit veranderd, enkel ouput.
    #GP = GP - (CPOOL - CPOOLMX)
    CPOOL = CPOOLMX
}

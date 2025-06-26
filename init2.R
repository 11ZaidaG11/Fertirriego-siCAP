# init 2
DTFAST = 1.0/NFAST

# number of "aborted" fruits, set fruits and total fruit number
TABF = 0
TNSF = 0
TNF = 0

DWL=zeros(100,1)
XLA=zeros(100,1)
DWS=zeros(100,1)
AGLS=zeros(100,1)
ABOR=zeros(100,1)
XNFT=zeros(100,1)
DWTR=zeros(100,1)
PLE=zeros(100,1)
PGL=zeros(100,1)
PGS=zeros(100,1)
ABNF=zeros(100,1)
NSF=zeros(100,1)

# forspeed
RCNF=zeros(100,1)

# C all variables at planting stage (1st truss flowering) FIRST SECTION
      DWL[1] = WLVSI*0.83
      AGLS[1] = 1.0
      DWS[1] = WSTMI*0.93
      DWTR[1] = 0
      XLA[1] = PLARI*0.85
      ABOR[1] = 0
      XNFT[1] = 0
      DWTR[1] = 0
# all variables at planting stage (1st truss flowering) SECOND SECTION
      DWL[2] = WLVSI*0.17
      AGLS[2] = 0.9
      DWS[2] = WSTMI*0.07
      DWTR[2] = 0.
      XLA[2] = PLARI*0.15
      ABOR[2] = 0
      XNFT[2] = 0
      DWTR[2] = 0
# !!! the organ age at starting date is very important for the
# final wxlaeight of first organs !!!

DWF=zeros(100,7)
AGF=zeros(100,7)
PGF=zeros(100,7)

# set initial conditions on a per plant basis
      PLSTN = PLSTNI
      PLAR = PLARI
	PLAR2 = PLARI
# CPOOL is supposed to start at 10% leaf d.wt
# aangepast ALEX op 5%
      CPOOL = 0.05 * WLVSI
      TDML = WLVSI
	TDML2 = WLVSI
      TDMS = WSTMI
      TDMF = 0
      DMGF = 0
      DMGL = WLVSI
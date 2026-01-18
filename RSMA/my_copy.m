A=R_RSMA_iE(:,6)
B=E_RSMA_iE(:,6)
C=C_RSMA_iE(:,6)
Er=E_range

E_range=Er
R_RSMA_iE(:,end)=A
E_RSMA_iE(:,end)=B
C_RSMA_iE(:,end)=C

R_RSMA_iE = [R_RSMA_iE,A]
E_RSMA_iE = [E_RSMA_iE,B]
C_RSMA_iE = [C_RSMA_iE,C]
E_range   = [E_range,Er]

A=[R_RSMA_iE(:,1:3),A]
B=[E_RSMA_iE(:,1:3),B]
C=[C_RSMA_iE(:,1:3),C]
Er=[E_range(:,1:3),Er]

R_RSMA_iE = [A,R_RSMA_iE]
E_RSMA_iE = [B,E_RSMA_iE]
C_RSMA_iE = [C,C_RSMA_iE]
E_range   = [1e-7,E_range]

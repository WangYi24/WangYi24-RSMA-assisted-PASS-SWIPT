A=R_RSMA_iE
B=E_RSMA_iE
C=C_RSMA_iE
Er=E_range

E_range=Er
R_RSMA_iE(:,2)= A
E_RSMA_iE(:,2)= B
C_RSMA_iE(:,2)= C

R_RSMA_iE = [R_RSMA_iE(:,1:4),A,R_RSMA_iE(:,7:end)]
E_RSMA_iE = [E_RSMA_iE(:,1:4),B,E_RSMA_iE(:,7:end)]
C_RSMA_iE = [C_RSMA_iE(:,1:4),C,C_RSMA_iE(:,7:end)]
E_range   = [E_range(:,1:4),Er,E_range(:,7:end)]

A=[R_RSMA_iE(:,1:3),A]
B=[E_RSMA_iE(:,1:3),B]
C=[C_RSMA_iE(:,1:3),C]
Er=[E_range(:,1:3),Er]

R_RSMA_iE = [A,R_RSMA_iE]
E_RSMA_iE = [B,E_RSMA_iE]
C_RSMA_iE = [C,C_RSMA_iE]
E_range   = [Er,E_range]

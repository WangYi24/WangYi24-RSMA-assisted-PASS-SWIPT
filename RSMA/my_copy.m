A=R_RSMA_iP
B=E_RSMA_iP
C=C_RSMA_iP
Er=P_range

P_range=Er
R_RSMA_iP(:,1)=A
E_RSMA_iP(:,1)=B
C_RSMA_iP(:,1)=C

R_RSMA_iP = [R_RSMA_iP(:,1:2),A]
E_RSMA_iP = [E_RSMA_iP(:,1:2),B]
C_RSMA_iP = [C_RSMA_iP(:,1:2),C]
P_range   = [P_range(:,1:2),Er]

A=[R_RSMA_iP(:,1:3),A]
B=[E_RSMA_iP(:,1:3),B]
C=[C_RSMA_iP(:,1:3),C]
Er=[P_range(:,1:3),Er]

R_RSMA_iP = [A,R_RSMA_iP]
E_RSMA_iP = [B,E_RSMA_iP]
C_RSMA_iP = [C,C_RSMA_iP]
P_range   = [1e-7,P_range]

A=R_SDMA_iE(6)
B=E_SDMA_iE(6)

Er= E_range

E_range(2)=Er
R_SDMA_iE(2)=A
E_SDMA_iE(2)=B
C_SDMA_iE(:,2)=C

R_SDMA_iE = [R_SDMA_iE(:,1:4),A]
E_SDMA_iE = [E_SDMA_iE(:,1:4),B]
C_SDMA_iE = [C_SDMA_iE(:,1),C]
E_range = [E_range(:,1:2),Er]

A = R_ct(:,9:15);
B = EHR_ct(:,9:15);
C = c_ct(:,9:15);

R_ct(:,9:15) = A;
EHR_ct(:,9:15) = B;
c_ct(:,9:15) = C;

R_SDMA_iE(:,3) = sum(R_ct,"all")/K / (ct-num_badpoint);
E_SDMA_iE(:,3) = sum(EHR_ct,2) / (ct-num_badpoint);
C_SDMA_iE(:,1) = sum(c_ct,2) / (ct-num_badpoint);

R_SDMA_iE = [A*3,R_SDMA_iE]
E_SDMA_iE = [B,E_SDMA_iE]
E_range  = [1e-7,E_range]

A = 12.2203
B =2.0565e-07
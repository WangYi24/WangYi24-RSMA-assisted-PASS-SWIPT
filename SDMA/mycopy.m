A=R_SDMA_iE(:,1:11)
B=E_SDMA_iE(:,1:11)
C=C_SDMA_iE(:,1:11)
Er= E_range(:,1:11)

E_range(2)=Er
R_SDMA_iE(2)=A
E_SDMA_iE(2)=B
C_SDMA_iE(:,2)=C

R_SDMA_iE = [R_SDMA_iE,A]
E_SDMA_iE = [E_SDMA_iE,B]
C_SDMA_iE = [C_SDMA_iE,C]
E_range = [E_range,Er]
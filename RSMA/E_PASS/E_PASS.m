
figure

load('E[0.80]221_D20_P38.mat')
% plot(E_range,R_RSMA_iE*2,'ko-'),hold on
plot(E_range,R_RSMA_iE,'ko-'),hold on
% plot(E_range,sum(C_RSMA_iE),'rs-'),hold on

load('E[0.1-1]221_D20_P40.mat')
% plot(E_range,R_RSMA_iE*3,'ko--'),hold on
plot(E_range,R_RSMA_iE,'ro-'),hold on
% plot(E_range,sum(C_RSMA_iE),'rs--'),hold on



xlabel('E(uW)')
ylabel('WSR (bps/Hz)')
% ylim([0 20]);
% xlim([32 44])
% title('RSMA vs SDMA')
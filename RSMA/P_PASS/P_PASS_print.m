
figure
% load('P[30-44]221_D10.mat')
% % plot(P_range,R_RSMA_iP*2,'ko-'),hold on
% plot(P_range,R_RSMA_iP,'ko-'),hold on
% % plot(P_range,sum(C_RSMA_iP),'rs-'),hold on
% 
% load('P[30-44]221_D20.mat')
% % plot(P_range,R_RSMA_iP*3,'ko--'),hold on
% plot(P_range,R_RSMA_iP,'ro-'),hold on
% % plot(P_range,sum(C_RSMA_iP),'rs--'),hold on
% 
% load('P[30-44]221_D30.mat')
% % plot(P_range,R_RSMA_iP*2,'ko-'),hold on
% plot(P_range,R_RSMA_iP,'go-'),hold on
% % plot(P_range,sum(C_RSMA_iP),'rs-'),hold on

load('P[48]221_D10.mat')
% plot(P_range,R_RSMA_iP*2,'ko-'),hold on
plot(P_range,R_RSMA_iP,'ko-'),hold on
% plot(P_range,sum(C_RSMA_iP),'rs-'),hold on

load('P[48]221_D20.mat')
% plot(P_range,R_RSMA_iP*3,'ko--'),hold on
plot(P_range,R_RSMA_iP,'ro-'),hold on
% plot(P_range,sum(C_RSMA_iP),'rs--'),hold on

load('P[48]221_D30.mat')
% plot(P_range,R_RSMA_iP*2,'ko-'),hold on
plot(P_range,R_RSMA_iP,'go-'),hold on
% plot(P_range,sum(C_RSMA_iP),'rs-'),hold on

% load('P[44]221_D40.mat')
% % plot(P_range,R_RSMA_iP*2,'ko-'),hold on
% plot(P_range,R_RSMA_iP,'bo-'),hold on
% % plot(P_range,sum(C_RSMA_iP),'rs-'),hold on


xlabel('P(dBm)')
ylabel('WSR (bps/Hz)')
% ylim([0 20]);
% xlim([32 44])
% title('RSMA vs SDMA')
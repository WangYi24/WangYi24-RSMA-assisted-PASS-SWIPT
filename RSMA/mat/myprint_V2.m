%% E
figure
load('E[0.1-1]221_D20_P40.mat')
plot(E_range*1e6,R_RSMA_iE,'ko-'),hold on
% plot(E_range*1e6,sum(C_RSMA_iE),'ro-'),hold on

load('E[0.2-3]221_D20_P45.mat')
plot(E_range*1e6,R_RSMA_iE,'ko--'),hold on
% plot(E_range*1e6,sum(C_RSMA_iE),'ro--'),hold on


load('E_Conv[0.1-0.5]221_D20_P40.mat')
plot(E_range*1e6,R_RSMA_iE/K,'bs-'),hold on
% plot(E_range*1e6,sum(C_RSMA_iE),'rs-'),hold on

load('E_Conv[0.1-1.5]221_D20_P45.mat')
plot(E_range*1e6,R_RSMA_iE,'bs--'),hold on
% plot(E_range*1e6,sum(C_RSMA_iE),'rs--'),hold on

legend('WSR','common rate')
xlabel('E_{min}(uW)')
ylabel('WSR (bps/Hz)')
ylim([0 20]);
xlim([1e-1 30e-1])
line([E_range(end-1)*1e6, E_range(end-1)*1e6], [0, R_RSMA_iE(end-1)],'Color', 'k')
title('RSMA')

%% P
figure
load('P[30-44]221_D10.mat')
% plot(P_range,R_RSMA_iP*2,'ko-'),hold on
plot(P_range,R_RSMA_iP,'ko-'),hold on
% plot(P_range,sum(C_RSMA_iP),'rs-'),hold on

load('P[30-44]221_D20.mat')
% plot(P_range,R_RSMA_iP*3,'ko--'),hold on
plot(P_range,R_RSMA_iP,'ko--'),hold on
% plot(P_range,sum(C_RSMA_iP),'rs--'),hold on

load('P[30-44]221_D30.mat')
% plot(P_range,R_RSMA_iP*2,'ko-'),hold on
plot(P_range,R_RSMA_iP,'b>-'),hold on
% plot(P_range,sum(C_RSMA_iP),'rs-'),hold on

load('P_Conv[44]221_D10.mat')
% plot(P_range,R_RSMA_iP*3,'ko--'),hold on
plot(P_range,R_RSMA_iP,'b>--'),hold on

load('P_Conv[44]221_D20.mat')
% plot(P_range,R_RSMA_iP*3,'ko--'),hold on
plot(P_range,R_RSMA_iP,'b>--'),hold on

load('P_Conv[44]221_D30.mat')
% plot(P_range,R_RSMA_iP*3,'ko--'),hold on
plot(P_range,R_RSMA_iP,'b>--'),hold on


legend('SR,K=2','WSR,K=2','common rate,K=2','SR,K=3','WSR,K=3','common rate,K=3')
xlabel('P(dBm)')
ylabel('WSR (bps/Hz)')
% ylim([0 20]);
xlim([32 44])
title('RSMA vs SDMA')
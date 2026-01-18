figure
load('E[0.4-1]221.mat')
figure
plot(E_range*1e6,R_RSMA_iE,'ko-'),hold on
plot(E_range*1e6,sum(C_RSMA_iE),'ro-'),hold on

load('E[0.2-1]231.mat')
plot(E_range*1e6,R_RSMA_iE/K,'ks--'),hold on
plot(E_range*1e6,sum(C_RSMA_iE),'rs--'),hold on


load('E_Conv[0.1-0.5]221.mat')
plot(E_range*1e6,R_RSMA_iE/K,'gs-'),hold on
% plot(E_range*1e6,sum(C_RSMA_iE),'rs-'),hold on

load('E_Conv[0.1-0.5]231.mat')
plot(E_range*1e6,R_RSMA_iE/K,'gs--'),hold on
% plot(E_range*1e6,sum(C_RSMA_iE),'rs--'),hold on

legend('WSR','common rate')
xlabel('E_{min}(W)')
ylabel('Rate (bps/Hz)')
ylim([0 20]);
xlim([1e-1 10e-1])
line([E_range(end-1)*1e6, E_range(end-1)*1e6], [0, R_RSMA_iE(end-1)],'Color', 'k')
title('RSMA')
figure
load('E[0.2-1]221.mat')
idx = [1:7, 10:17];  % 手动列出要保留的索引

plot(E_range(idx)*1e6,R_SDMA_iE(idx)/K,'rs-'),hold on
% plot(E_range*1e6,sum(C_SDMA_iE),'ro-'),hold on

load('E[0.4-1]231.mat')
plot(E_range*1e6,R_SDMA_iE/K,'rs--'),hold on
% plot(E_range*1e6,sum(C_SDMA_iE),'ro--'),hold on


legend('WSR','common rate')
xlabel('E_{min}(W)')
ylabel('WSR (bps/Hz)')
ylim([0 20]);
xlim([4e-1 10e-1])
line([E_range(end-1)*1e6, E_range(end-1)*1e6], [0, R_SDMA_iE(end-1)],'Color', 'k')

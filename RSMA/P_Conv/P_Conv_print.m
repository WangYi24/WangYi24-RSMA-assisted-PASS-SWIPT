%% P

load('P_Conv[48]221_D10.mat')
% plot(P_range,R_RSMA_iP*3,'ko--'),hold on
plot(P_range,R_RSMA_iP,'r>--'),hold on

load('P_Conv[48]221_D20.mat')
% plot(P_range,R_RSMA_iP*3,'ko--'),hold on
plot(P_range,R_RSMA_iP,'g>--'),hold on

load('P_Conv[48]221_D30.mat')
% plot(P_range,R_RSMA_iP*3,'ko--'),hold on
plot(P_range,R_RSMA_iP,'b>--'),hold on

xlabel('P(dBm)')
ylabel('WSR (bps/Hz)')
% ylim([0 20]);
% xlim([32 44])
% title('RSMA vs SDMA')
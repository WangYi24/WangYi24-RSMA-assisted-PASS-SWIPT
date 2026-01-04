clear
cvx_clear

%% P
clear
cvx_clear

M = 4;
N = 2;
P_range = -15:5:25;

R_MW_iP = zeros(M,length(P_range)); 
R_MW_ZF_iP = zeros(M,length(P_range)); 
R_Conv_ZF_iP = zeros(M,length(P_range)); 
for iP = 1:length(P_range)
    P = P_range(iP)
    [R_MW_iP(:,iP),R_MW_ZF_iP(:,iP),R_Conv_ZF_iP(:,iP)] = SDMA_fix(M,N,P);
end
figure
plot(P_range,min(R_MW_iP),'ro-'),hold on
plot(P_range,min(R_MW_ZF_iP),'bs-'),hold on
plot(P_range,min(R_Conv_ZF_iP),'k>-'),hold on

legend('MW','MW-ZF')
xlabel('P_{max}(dBm)')
ylabel('Minimum data rate (bps)')

% %% N
% clear
% cvx_clear
% 
% M = 4;
% N_range = 4:2:12;
% P = 0;
% 
% R_MW_iN = zeros(M,length(N_range)); 
% R_MW_ZF_iN = zeros(M,length(N_range)); 
% R_Conv_ZF_iN = zeros(M,length(N_range)); 
% for iN = 1:length(N_range)
%     N = N_range(iN)
%     [R_MW_iN(:,iN),R_MW_ZF_iN(:,iN),R_Conv_ZF_iN(:,iN)] = test_func_v4(M,N,P);
% end
% figure
% plot(N_range,min(R_MW_iN),'ro-'),hold on
% plot(N_range,min(R_MW_ZF_iN),'bs-'),hold on
% plot(N_range,min(R_Conv_ZF_iN),'k>-'),hold on
% 
% legend('MW','MW-ZF')
% xlabel('N')
% ylabel('Minimum data rate (bps)')
% % 
% %% M
% clear
% cvx_clear
% 
% M_range = 2:2:8;
% N = 8;
% P = 10;
% 
% R_MW_iM = 100*ones(8,length(M_range)); 
% R_MW_ZF_iM = 100*ones(8,length(M_range)); 
% R_Conv_ZF_iM = 100*ones(8,length(M_range)); 
% for iM = 1:length(M_range)
%     M = M_range(iM)
%     [R_MW_iM(1:M,iM),R_MW_ZF_iM(1:M,iM),R_Conv_ZF_iM(1:M,iM)] = test_func_v5(M,N,P);
% end
% figure
% plot(M_range,min(R_MW_iM),'ro-'),hold on
% plot(M_range,min(R_MW_ZF_iM),'bs-'),hold on
% plot(M_range,min(R_Conv_ZF_iM),'k>-'),hold on
% 
% legend('MW','MW-ZF')
% xlabel('M')
% ylabel('Minimum data rate (bps)')
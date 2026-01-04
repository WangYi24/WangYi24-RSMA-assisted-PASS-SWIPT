clear
cvx_clear

%% Emin
clear
cvx_clear

K = 2;
J = 1;
N = 2;
E_range = 2e-7:1e-7:1e-6;

P=40;

R_RSMA_iE = zeros(1,length(E_range)); 
E_RSMA_iE = zeros(J,length(E_range));
C_RSMA_iE = zeros(K,length(E_range));
R_RSMA_ct = cell(1,length(E_range));
EHR_ct = cell(1,length(E_range));
c_ct = cell(1,length(E_range));

for iE = 1:length(E_range)
    E = E_range(iE)
    [R_RSMA_ct{iE},EHR_ct{iE},c_ct{iE},R_RSMA_iE(:,iE),E_RSMA_iE(:,iE),C_RSMA_iE(:,iE)]= RSMA_search(K,J,N,P,E);
    % [R_RSMA_iE(:,iE),E_RSMA_iE(:,iE),C_RSMA_iE(:,iE)] = RSMA_search(K,J,N,P,E);
    save("E.mat")
end

figure

plot(E_range*1e6,R_RSMA_iE,'ko-'),hold on

legend('RSMA')
xlabel('E_{min}(W)')
ylabel('sum rate (bps/Hz)')
ylim([30 46]);
xlim([0 7e-1])
line([E_range(5)*1e6, E_range(5)*1e6], [0, R_RSMA_iE(5)],'Color', 'k')


% %% P
% clear
% cvx_clear
% 
% K = 3;
% J = 2;
% N = 4;
% P_range = -15:5:25;
% 
% R_RSMA_iP = zeros(1,length(P_range)); 
% E_RSMA_iP = zeros(J,length(P_range));
% 
% for iP = 6:length(P_range)
%     P = P_range(iP)
%     [R_RSMA_iP(:,iP),E_RSMA_iP(:,iP)] = RSMA_search(K,J,N,P);
% end
% figure
% plot(P_range,R_RSMA_iP,'k--'),hold on
% 
% legend('RSMA')
% xlabel('P_{max}(dBm)')
% ylabel('Minimum data rate (bps)')

% %% N
% clear
% cvx_clear
% 
% M = 4;
% N_range = 4:2:12;
% P = 0;
% 
% R_MW_iN = zeros(M,length(N_range)); 

% for iN = 1:length(N_range)
%     N = N_range(iN)
%     [R_MW_iN(:,iN)] = test_func_v4(M,N,P);
% end
% figure
% plot(N_range,min(R_MW_iN),'ro-'),hold on
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

% for iM = 1:length(M_range)
%     M = M_range(iM)
%     [R_MW_iM(1:M,iM)] = test_func_v5(M,N,P);
% end
% figure
% plot(M_range,min(R_MW_iM),'ro-'),hold on

% 
% legend('MW','MW-ZF')
% xlabel('M')
% ylabel('Minimum data rate (bps)')
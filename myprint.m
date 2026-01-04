%% P-OFDM
clear
load('R_SW_OFDM_iP.mat')
load('OFDM_PdBm[-10_30]_CP.mat')
figure
plot(P_range,min(R_MW_OFDMA_MRT_iP)/1e9,'ko-'),hold on
plot(P_range,min(R_MW_OFDM_ZF_iP)/1e9,'ks--'),hold on
plot(P_range,min(R_SW_OFDM_iP)/1e9,'k>-'),hold on

legend('MW-OFDM-MRT','MW-OFDM-ZF','SW-OFDM')
xlabel('P_{T}(dBm)')
ylabel('Minimum data rate (Gbps)')

%% P
clear
load('PdBm[-15_25]_D[30_10]_search.mat')
figure
plot(P_range,min(R_MW_iP),'ro-'),hold on
plot(P_range,min(R_MW_ZF_iP),'bs-'),hold on
plot(P_range,min(R_Conv_ZF_iP),'k>-'),hold on

legend('MW','MW-ZF','Conv-ZF')
xlabel('P_{T}(dBm)')
ylabel('Minimum data rate (bps)')

%% N
clear
load('N[4_12]_D[30_10]_search.mat')
figure
plot(N_range,min(R_MW_iN),'ro-'),hold on
plot(N_range,min(R_MW_ZF_iN),'bs-'),hold on
plot(N_range,min(R_Conv_ZF_iN),'k>-'),hold on

legend('MW','MW-ZF','Conv-ZF')
xlabel('N')
ylabel('Minimum data rate (bps)')

%% M
clear
load('N[4_12]_D[30_10]_search.mat')
figure
plot(M_range,min(R_MW_iM),'ro-'),hold on
plot(M_range,min(R_MW_ZF_iM),'bs-'),hold on
plot(M_range,min(R_Conv_ZF_iM),'k>-'),hold on

legend('MW','MW-ZF','Conv-ZF')
xlabel('N')
ylabel('Minimum data rate (bps)')
%% D
clear
load('D[10_40].mat')
figure
plot(D_range,min(R_MW_iD),'r-o'),hold on
plot(D_range,min(R_SW_iD),'g-o'),hold on
plot(D_range,min(R_ZF_iD),'b-o')
legend('MW-OFDMA','SW-OFDMA','Conv-ZF')
xlabel('P_{T}(dBm)')
ylabel('Minimum data rate (Gbps)')

plot(P_range,min(R_MW_OFDMA_MRT_iP)/1e9,'bo-'),hold on
plot(P_range,min(R_MW_OFDM_ZF_iP)/1e9,'bs--'),hold on

legend('MW-OFDMA-MRT','MW-OFDM-ZF','SW-OFDMA')
xlabel('P_{T}(dBm)')
ylabel('Minimum data rate (Gbps)')

%% 迭代图
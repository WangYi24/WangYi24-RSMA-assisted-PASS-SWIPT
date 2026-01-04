function [R_RSMA_ct,EHR_ct] = RSMA_fix(K,J,N,P_dBm_max,Emin)
%解广播情况下的多波导多用户最大化最小速率问题
%天线位置固定到最近用户的位置
%% parameter
f_c = 28*1e9;%载波频率Hz
c = 3e8;%光速m/s
n_eff = 1.4;%波导折射率
lambda_c = c/f_c;
lambda_g = lambda_c / n_eff;
eta = lambda_c/4/pi;%自由空间损失因子

sigma2_dBm = -90;  % 噪声功率 (dBm)，公式给的是带宽乘以噪声功率谱密度，文中直接给大小有问题吧？
sigma2 = 10.^((sigma2_dBm-30)./10);%噪声功率s
% P_dBm_max = 0;
P_max = 10.^((P_dBm_max-30)./10);

D_x = 30;     %波导长度
D_y = 10;   %区域宽度
d = 3;      %波导高 (k)

ct = 50;%用户实现次数
num_badpoint = 0;%坏点

%% 
Loc_PA = zeros(3,N);%PA位置
%初始化PA位置，均为0
for i=1:N
    Loc_PA(:,i) = [0;(i-1)*D_y/(N-1);d];
end

R_RSMA_ct = zeros(1,ct); 
EHR_ct = zeros(J,ct); 
%固定EHR的位置
deltax = D_x/2/(J-1);deltay = D_y/2/(J-1);
Loc_EHR = [zeros(1,J)+ deltax*(-J/2 + 0.5 : J/2 - 0.5); D_y*(ones(1,J)-0.5) + deltay*(-J/2 + 0.5 : J/2 - 0.5); zeros(1,J)];
for it = 1:ct %Monte Carlo simulations     
    it
    %生成IDR位置
    Loc_IDR = [D_x*(rand(1,K)-0.5); D_y*rand(1,K); zeros(1,K)];
    %Loc_EHR = [D_x*(rand(1,J)-0.5); D_y*rand(1,J); zeros(1,J)];

    %初始化PA位置，距离EHR近
    for n = 1:N
        % 计算所有用户到第 i 个波导的y轴距离
        dist_y = Loc_EHR(2,:) - Loc_PA(2,n);      % K*1 - 1*1 = K*1
        [~, idx] = min(abs(dist_y));            % 获得最小距离点索引
        Loc_PA(1,n) = Loc_EHR(1,idx);            
    end

    %PAn到用户的信道：h_IDR[N,K],h_EHR[N,J]
    h_IDR = channel(Loc_PA, Loc_IDR, eta, lambda_c, lambda_g);
    h_EHR = channel(Loc_PA, Loc_EHR, eta, lambda_c, lambda_g);
    [W_p, wc, c, SumRate] = rsma_wmmse(h_IDR,h_EHR,sigma2, P_max,Emin);
    if isnan(SumRate) || SumRate==-inf
        num_badpoint = num_badpoint + 1;
        R_RSMA_ct(:,it) = 0;
        continue
    end
    R_RSMA_ct(:,it) = SumRate;
    term_c = abs(h_EHR' * wc).^2;     % [J, 1]
    h_W = h_EHR' * W_p;              % [J, K]
    term_p = sum(abs(h_W).^2, 2);  % 对每一行求和，得到 [J, 1]
    EHR_ct(:,it)= term_p + term_c;
end

R_RSMA_ct = sum(R_RSMA_ct,2) / (ct-num_badpoint);
EHR_ct = sum(EHR_ct,2) / (ct-num_badpoint);
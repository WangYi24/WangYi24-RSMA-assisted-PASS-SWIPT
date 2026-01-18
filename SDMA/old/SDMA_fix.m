function [R_MW_ct,R_MW_ZF_ct,R_Conv_ZF_ct] = SDMA_fix(M,N,P_dBm_max)
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
d = 3;      %波导高 (m)

ct = 10;%用户实现次数
ct_ZF = 1e4;
num_badpoint = 0;%坏点

%% 
Loc_PA = zeros(3,N);%PA位置
Loc_U = zeros(3,M);%用户位置
%初始化PA位置，均为0
for i=1:N
    Loc_PA(:,i) = [0;(i-1)*D_y/(N-1);d];
end
h = zeros(M,N);%信道，子载波k从波导n上的PA到用户m的信道
tau = zeros(M,N);%时延

R_MW_ct = zeros(M,ct); 
R_MW_ZF_ct = zeros(M,ct); 
R_Conv_ZF_ct = zeros(M,ct); 
%% MW_SDMA
for it = 1:ct %Monte Carlo simulations     
    %生成用户位置
    for i=1:M
        Loc_U(:,i) = [D_x*rand(1)-D_x/2; D_y*rand(1); 0];
    end

     %PA位置，先简单按照Downlink Beamforming with Pinching-Antenna Assisted MIMO Systems公式40，做了文中I的初始化 
    for n = 1:N
        % 计算所有用户到第 i 个波导的y轴距离
        dist_y = Loc_U(2,:) - Loc_PA(2,n);      % K*1 - 1*1 = K*1
        [~, idx] = min(abs(dist_y));            % 获得最小距离点索引
        Loc_PA(1,n) = Loc_U(1,idx);            
    end

    %PAn到用户m的信道：h[M,N]
    h = channel(Loc_PA, Loc_U, eta, lambda_c, lambda_g);
    %[t_opt, W_opt, out] = QT_maxmin_beamforming(h, sigma2, P_dBm_max);
    % out.history
    [W_opt, SINR_hist] = beamforming_L_v2(h', sigma2, P_max);
    if any(sum(abs(W_opt))<1e-8)
        num_badpoint = num_badpoint + 1;
        R_MW_ct(:,it) = 0;
        continue
    end
    G = abs((h*W_opt)).^2;%[M,M]
    R_test = log2(1+SINR_hist(end));
    for m = 1:M
        SINR = G(m,m)/(sum(G(m,:))-G(m,m)+sigma2);
        R_MW_ct(m,it) = log2(1+SINR);
    end
end

%% ZF
Loc_Conv = [zeros(1,N); D_y/2 + lambda_c/2 * [-N/2:N/2-1]; d * ones(1,N)]; 
for it = 1:ct_ZF %Monte Carlo simulations     
    %生成用户位置
    for i=1:M
        Loc_U(:,i) = [D_x*rand(1)-D_x/2; D_y*rand(1); 0];
    end
    %PA位置，先简单按照Downlink Beamforming with Pinching-Antenna Assisted MIMO Systems公式40，做了文中I的初始化
 
    for n = 1:N
        % 计算所有用户到第 i 个波导的y轴距离
        dist_y = Loc_U(2,:) - Loc_PA(2,n);      % K*1 - 1*1 = K*1
        [~, idx] = min(abs(dist_y));            % 获得最小距离点索引
        Loc_PA(1,n) = Loc_U(1,idx);            
    end
    R_MW_ZF_ct(:,it) = ZF(M,N,lambda_c,lambda_g,sigma2,P_max,Loc_PA,Loc_U);
        
    R_Conv_ZF_ct(:,it) = Conv_ZF(M,N,lambda_c,sigma2,P_max,Loc_Conv,Loc_U);
  
end

R_MW_ct = sum(R_MW_ct,2)/(ct-num_badpoint);
R_MW_ZF_ct = sum(R_MW_ZF_ct,2)/ct_ZF;
R_Conv_ZF_ct = sum(R_Conv_ZF_ct,2)/ct_ZF;
% SINR_hist
num_badpoint


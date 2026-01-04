K = 3;
J = 2;
N = 4;
P_dBm_max = 10;
%解广播情况下的多波导多用户最大化最小速率问题
%搜索天线位置
%% parameter
f_c = 28*1e9;%载波频率Hz
c = 3e8;%光速m/s
n_eff = 1.4;%波导折射率
lambda_c = c/f_c;
lambda_g = lambda_c / n_eff;
eta = lambda_c/4/pi;%自由空间损失因子
Emin = 5e-11;

sigma2_dBm = -90;  % 噪声功率 (dBm)，公式给的是带宽乘以噪声功率谱密度，文中直接给大小有问题吧？
sigma2 = 10.^((sigma2_dBm-30)./10);%噪声功率s
% P_dBm_max = 0;
P_max = 10.^((P_dBm_max-30)./10);

D_x = 30;     %波导长度
D_y = 10;   %区域宽度
d = 3;      %波导高 (k)


num_badpoint = 0;%坏点
num_point = 1000;%波导上的离散位置个数
candidates = linspace(-D_x/2,D_x/2,num_point);
maxiter = 50;
tolerance = 1e-3;

%% 
%初始化PA位置，均为0
Loc_PA=[zeros(1,N); (0:N-1)*D_y/(N-1); zeros(1,N)];

%生成设备位置，列为该设备坐标
Loc_IDR = [D_x*(rand(1,K)-0.5); D_y*rand(1,K); zeros(1,K)];
Loc_EHR = [D_x*(rand(1,J)-0.5); D_y*rand(1,J); zeros(1,J)];

%PA位置，先简单按照Downlink Beamforming with Pinching-Antenna Assisted MIMO Systems公式40，做了文中I的初始化 
for n = 1:N
    % 计算所有用户到第 i 个波导的y轴距离
    dist_y = Loc_IDR(2,:) - Loc_PA(2,n);      % K*1 - 1*1 = K*1
    [~, idx] = min(abs(dist_y));            % 获得最小距离点索引
    Loc_PA(1,n) = Loc_IDR(1,idx);            
end
h_IDR = channel(Loc_PA, Loc_IDR, eta, lambda_c, lambda_g);
h_EHR = channel(Loc_PA, Loc_EHR, eta, lambda_c, lambda_g);
% 内层固定PA位置，优化W                  
[W_p, wc, c, SumRate] = rsma_wmmse(h_IDR,h_EHR,sigma2, P_max,Emin); 

Loc_PA_temp = Loc_PA;
for n = 1:N
    xn_opt = Loc_PA_temp(1,n);
    %PAn到用户的信道：h_IDR[N,K],h_EHR[N,J]

    sumrate = zeros(1,length(candidates)); 
    EHR = zeros(J,length(candidates));
    for sam = 1:length(candidates)
        xn = candidates(sam);
        
        Loc_PA_temp(1,n) = xn;
        %PAn到用户的信道：h_IDR[N,K],h_EHR[N,J]
        h_IDR = channel(Loc_PA_temp, Loc_IDR, eta, lambda_c, lambda_g);
        h_EHR = channel(Loc_PA_temp, Loc_EHR, eta, lambda_c, lambda_g);
        
        G = abs(h_IDR' * W_p).^2;
        for k = 1:K
            sumrate(:,sam) = sumrate(:,sam) + log2(1+G(k,k)/(sum(G(k,:))-G(k,k)+sigma2));
        end
       
        term_c = abs(h_EHR' * wc).^2;     % [J, 1]
        h_W = h_EHR' * W_p;              % [J, K]
        term_p = sum(abs(h_W).^2, 2);  % 对每一行求和，得到 [J, 1]
        EHR(:,sam)= term_p + term_c;
       
    end
    plot(candidates,sumrate)
    xlabel('x_{1}^{p}')
    ylabel('sum rate')
    Loc_PA_temp(1,n) = xn_opt;
end
Loc_PA = Loc_PA_temp;
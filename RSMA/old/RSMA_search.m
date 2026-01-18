function [R_RSMA_ct,EHR_ct,c_ct,R,EHR,c] = RSMA_search(K,J,N,P_dBm_max,Emin)

%搜索天线位置，每个位置做一次优化
%% parameter
sigma2_dBm = -90;  % 噪声功率 (dBm)，公式给的是带宽乘以噪声功率谱密度，文中直接给大小有问题吧？
sigma2 = 10.^((sigma2_dBm-30)./10);%噪声功率s

P_max = 10.^((P_dBm_max-30)./10);

D_x = 20;     %波导长度
D_y = 10;   %区域宽度
d = 3;      %波导高 (k)

ct = 40;%用户实现次数
num_badpoint = 0;%坏点
Ns = 50;%波导上的离散位置个数
candidates = linspace(-D_x/2,D_x/2,Ns+1);
maxiter = 5;%目标函数是和速率，感觉不需要外层迭代了
tolerance = 1e-1;

%% 
%初始化PA位置，均为0
Loc_PA=[zeros(1,N); ((0:N-1)+0.5)*D_y/N; d * ones(1,N)];
R_RSMA_ct = zeros(1,ct); 
EHR_ct = zeros(J,ct); 
c_ct = zeros(K,ct); 

%固定EHR的位置
Loc_EHR = [zeros(1,J) - D_x/2 + ((0:J-1)+0.5)*D_x/J; ((0:J-1)+0.5)*D_y/J; zeros(1,J)];
%PA位置，初始化到EHR附近
for n = 1:N
    % 计算所有IDR到第 i 个波导的y轴距离
    dist_y = Loc_EHR(2,:) - Loc_PA(2,n);      % K*1 - 1*1 = K*1
    [~, idx] = min(abs(dist_y));            % 获得最小距离点索引
    Loc_PA(1,n) = Loc_EHR(1,idx);            
end

for it = 1:ct %Monte Carlo simulations     
    it
    %生成设备位置，列为该设备坐标
    Loc_IDR = [D_x*(rand(1,K)-0.5); D_y*rand(1,K); zeros(1,K)];
    SumRate_hist=[];%记录迭代的目标函数
    xn_hist = [];
    SumRate_best = -inf;
    %best记录外层迭代的结果
    for iter = 1:maxiter
        % ==== 外层搜索PA位置 ====
        Loc_PA_temp = Loc_PA;
        SumRate_opt = -inf;%opt记录搜索天线的结果
        W_opt = zeros(N,K+1);
        c_opt = zeros(K,1);
        for n = 1:N
            n
            xn_opt = Loc_PA_temp(1,n);     
            
            for xn = candidates
                Loc_PA_temp(1,n) = xn;
                %PAn到用户的信道：h_IDR[N,K],h_EHR[N,J]
                h_IDR = channel(Loc_PA_temp, Loc_IDR);
                h_EHR = channel(Loc_PA_temp, Loc_EHR);
                % 内层固定PA位置，优化W                  
                [W_p, wc, SumRate] = rsma_wmmse(h_IDR,h_EHR,sigma2, P_max,Emin);
                if isnan(SumRate) || SumRate==-inf
                    continue
                end
                if SumRate > SumRate_opt
                    SumRate_opt = SumRate;
                    xn_opt = xn;
                    W_opt = [wc,W_p];
                end
            end
            Loc_PA_temp(1,n) = xn_opt;
        end
        Loc_PA = Loc_PA_temp;

        % 检查目标函数是否收敛
        if isnan(SumRate) || SumRate==-inf
            break %此次用户的分布范围无法满足约束
        end
        % xn_hist = [xn_hist;Loc_PA(1,:)];
        % SumRate_hist = [SumRate_hist,SumRate_opt];
        
        if abs(SumRate_opt-SumRate_best)<tolerance
            break;
        end
        if SumRate_opt >= SumRate_best
            SumRate_best = SumRate_opt;
            Loc_PA_best = Loc_PA;
            W_best = W_opt;
        end
    end

    if SumRate_best==-inf
        num_badpoint = num_badpoint + 1;
        R_RSMA_ct(:,it) = 0;
        continue %进行下一次用户随机采样
    end
    
    
    wc = W_best(:,1);
    W_p = W_best(:,2:end);
    h_IDR = channel(Loc_PA_best, Loc_IDR);
    term_c = abs(h_IDR' * wc).^2;   % [K, 1]
    h_W = h_EHR' * W_p;              % [K, K]
    term_p = sum(abs(h_W).^2, 2);  % 对每一行求和，得到 [K, 1]
    c_ct(:,it) = log2(1+term_c./term_p);

    h_EHR = channel(Loc_PA_best, Loc_EHR);
    term_c = abs(h_EHR' * wc).^2;     % [J, 1]
    h_W = h_EHR' * W_p;              % [J, K]
    term_p = sum(abs(h_W).^2, 2);  % 对每一行求和，得到 [J, 1]
    EHR_ct(:,it)= term_p + term_c;
    R_RSMA_ct(:,it) = SumRate_best;
end

R = sum(R_RSMA_ct,2) / (ct-num_badpoint);
EHR = sum(EHR_ct,2) / (ct-num_badpoint);
c = sum(c_ct,2) / (ct-num_badpoint);
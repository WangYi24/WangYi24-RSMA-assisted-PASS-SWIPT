function [R_RSMA_ct,EHR_ct,c_ct] = RSMA_search_AO(K,J,N,P_dBm_max,Emin)

%搜索天线位置，AO
%% parameter


sigma2_dBm = -90;  % 噪声功率 (dBm)，公式给的是带宽乘以噪声功率谱密度，文中直接给大小有问题吧？
sigma2 = 10.^((sigma2_dBm-30)./10);%噪声功率s

P_max = 10.^((P_dBm_max-30)./10);

D_x = 30;     %波导长度
D_y = 10;   %区域宽度
d = 3;      %波导高 (k)

ct = 20;%用户实现次数
num_badpoint = 0;%坏点
Ns = 100;%波导上的离散位置个数
candidates = linspace(-D_x/2,D_x/2,Ns+1);
maxiter = 100;%目标函数是和速率，感觉不需要外层迭代了
tolerance = 1e-1;

%% 
%初始化PA位置，均为0
Loc_PA=[zeros(1,N); (0:N-1)*D_y/(N-1); zeros(1,N)];
R_RSMA_ct = zeros(1,ct); 
EHR_ct = zeros(J,ct); 
c_ct = zeros(K,ct); 

%固定EHR的位置
deltax = D_x/2/(J-1);deltay = D_y/2/(J-1);
Loc_EHR = [zeros(1,J)+ deltax*(-J/2 + 0.5 : J/2 - 0.5); D_y*(ones(1,J)-0.5) + deltay*(-J/2 + 0.5 : J/2 - 0.5); zeros(1,J)];
%PA位置，初始化到EHR附近
for n = 1:N
    % 计算所有EHR到第 i 个波导的y轴距离
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
      
        % ==== 优化W,c ====
        h_IDR = channel(Loc_PA, Loc_IDR);
        h_EHR = channel(Loc_PA, Loc_EHR);
        [W_p, wc, SumRate] = rsma_wmmse(h_IDR,h_EHR,sigma2, P_max,Emin);
        if isnan(SumRate) || SumRate==-inf
            continue
        end
        % ==== 搜索PA位置 ====
        [Loc_PA(1,:), SumRate_opt, convergence, EHR] = PSO_PA([wc,W_p], Loc_PA,Loc_IDR,Loc_EHR, D_x, sigma2,Emin);
        
        % 检查目标函数是否收敛
        if isnan(SumRate_opt) || SumRate_opt==-inf
            break %此次用户的分布范围无法满足约束
        end
        xn_hist = [xn_hist;Loc_PA(1,:)];
        SumRate_hist = [SumRate_hist,SumRate_opt];
        
        if abs(SumRate_opt-SumRate_best)<tolerance
            break;
        end
        if SumRate_opt >= SumRate_best
            SumRate_best = SumRate_opt;
            Loc_PA_best = Loc_PA;
            W_best = [wc,W_p];
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
    h_EHR = channel(Loc_PA_best, Loc_EHR);
    
    [Rate,c_ct(:,it),EHR_ct(:,it)] = getObj(W_p, wc,h_IDR,h_EHR,sigma2);
    R_RSMA_ct(:,it) = sum(Rate);
end

R_RSMA_ct = sum(R_RSMA_ct,2) / (ct-num_badpoint);
EHR_ct = sum(EHR_ct,2) / (ct-num_badpoint);
c_ct = sum(c_ct,2) / (ct-num_badpoint);
end


function [Rate,c,EHR] = getObj(W_p, wc,h_IDR,h_EHR,sigma2)
    
    [N,K] = size(W_p);
    term_c = abs(h_IDR' * wc).^2;   % [K, 1]
    G = abs((h_IDR'*W_p)).^2;%得到 [K, K]
    c = log2(1+term_c./sum(G, 2));
    R_p = zeros(K,1);
    for k = 1:K
        SINR = G(k,k)/(sum(G(k,:))-G(k,k)+sigma2);
        %R_c(k)= log2(1+abs(h_IDR(:,k)'*wc).^2/(sum(G(k,:))+sigma2));
        R_p(k) = log2(1+SINR);
    end
    Rate = R_p+c;

    term_c = abs(h_EHR' * wc).^2;     % [J, 1]
    h_W = h_EHR' * W_p;              % [J, K]
    term_p = sum(abs(h_W).^2, 2);  % 对每一行求和，得到 [J, 1]
    EHR= term_p + term_c;
end

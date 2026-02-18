function [SumRate_best,Loc_PA_best,W_best] = RSMA_search_v2(K,J,N,P_dBm_max,Emin,Loc_IDR,Loc_EHR,D_x,D_y,d,sigma2)
%搜索天线位置，每个位置做一次优化
%% parameter

P_max = 10.^((P_dBm_max-30)./10);

step = 0.5;%波导上的离散位置间距

maxiter = 5;%目标函数是和速率，感觉不需要外层迭代了
tolerance = 1e-1;

%% 
%可选位置
candidates = -D_x/2:step:D_x/2;
minPA = min([Loc_IDR(1,:),Loc_EHR(1,:)]) - step;
maxPA = max([Loc_IDR(1,:),Loc_EHR(1,:)]) + step;
candidates = candidates(candidates >=minPA  & candidates <= maxPA);
if isempty(candidates)
    candidates = (minPA + maxPA)/2;
end

%初始化PA位置，均为0
Loc_PA=[zeros(1,N); ((0:N-1)+0.5)*D_y/N; d * ones(1,N)];
%PA位置，初始化到EHR附近
for n = 1:N
    % 计算所有EHR到第 i 个波导的y轴距离
    dist_y = Loc_EHR(2,:) - Loc_PA(2,n);      % K*1 - 1*1 = K*1
    [~, idx] = min(abs(dist_y));            % 获得最小距离点索引
    Loc_PA(1,n) = Loc_EHR(1,idx);            
end
%%
SumRate_hist=[];%记录迭代的目标函数
xn_hist = [];
SumRate_best = -inf;
Loc_PA_best = Loc_PA;
W_best = NaN;
%先计算能否满足能量收集约束
%best记录外层迭代的结果
for iter = 1:maxiter
    % ==== 外层搜索PA位置 ====
    Loc_PA_temp = Loc_PA;
    SumRate_opt = -inf;%opt记录搜索天线的结果
    W_opt = zeros(N,K+1);
    c_opt = zeros(K,1);
    for n = 1:N
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
    if isnan(SumRate_opt) || SumRate_opt==-inf
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
end

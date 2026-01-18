function [R_MW_ct,R_MW_ZF_ct,R_Conv_ZF_ct] = SDMA_search(M,N,P_dBm_max)
% 解广播情况下的多波导多用户最大化最小速率问题
% 搜索天线位置
%% parameter
f_c = 28*1e9;%载波频率Hz
c = 3e8;%光速m/s
n_eff = 1.4;%波导折射率
lambda_c = c/f_c;
lambda_g = lambda_c / n_eff;
eta = lambda_c/4/pi;%自由空间损失因子

sigma2_dBm = -90;  % 噪声功率 (dBm)，公式给的是带宽乘以噪声功率谱密度，文中直接给大小有问题吧？
sigma2 = 10.^((sigma2_dBm-30)./10);%噪声功率s
P_max = 10.^((P_dBm_max-30)./10);

D_x = 30;     %波导长度
D_y = 10;   %区域宽度
d = 3;      %波导高 (m)

ct = 15;%用户实现次数15
ct_ZF = 1e3;
num_point = 1000;%波导上的离散位置个数
candidates = linspace(-D_x/2,D_x/2,num_point);
maxiter = 50;


h = zeros(M,N);%信道，子载波k从波导n上的PA到用户m的信道

R_MW_ct = zeros(M,ct); 
R_MW_ZF_ct = zeros(M,ct); 
R_Conv_ZF_ct = zeros(M,ct); 
epsilon = 1e-3;
%% MW
for it = 1:ct %Monte Carlo simulations 
    %初始化PA位置，均为0
    Loc_PA = zeros(3,N);
    for i=1:N
        Loc_PA(:,i) = [0;(i-1)*D_y/(N-1);d];
    end
    %生成用户位置
    Loc_U = zeros(3,M);
    for i=1:M
        Loc_U(:,i) = [D_x*rand(1)-D_x/2; D_y*rand(1); 0];
    end

    obj_hist=[];
    xn_hist = [];
    obj_best = zeros(M,1);%记录最优迭代的目标函数
    for iter = 1:maxiter
        % ==== 外层搜索PA位置 ====
        Loc_PA_temp = Loc_PA;
        for n = 1:N
            xn_best = Loc_PA_temp(1,n);
            SINR_opt = zeros(M,1);%记录最优位置的目标函数
            %这里应该可以把范围缩小一些，在离这两个波导最近的两个位置中间搜，能减少迭代次数
            % 计算所有用户到第 i 个波导的y轴距离
            dist_y = Loc_U(2,:) - Loc_PA(2,n);      % K*1 - 1*1 = K*1
            [~, idx] = sort(abs(dist_y));            % 获得最小距离点索引
            x1 = min(Loc_U(1,idx(1)),Loc_U(1,idx(2)));
            x2 = max(Loc_U(1,idx(1)),Loc_U(1,idx(2)));
            for xn = candidates(candidates>=x1 & candidates<=x2)
                Loc_PA_temp(1,n) = xn;
                %PAn到用户m的信道：h[M,N] 
                h = channel(Loc_PA_temp, Loc_U, eta, lambda_c, lambda_g);
                % 内层固定PA位置，优化W                  
                [W, SINR_hist] = beamforming_L_v2(h', sigma2, P_max);
                if any(sum(abs(W))<1e-8)
                    %这里暂时不知道如何处理，之后调试一下
                    continue;
                end
                SINR = SINR_hist(:,end);
                if min(SINR) > min(SINR_opt)
                    SINR_opt = SINR;
                    xn_best = xn;
                    W_opt = W;
                end
            end
            Loc_PA_temp(1,n) = xn_best;
        end
        Loc_PA = Loc_PA_temp;

        % 检查目标函数是否收敛
        xn_hist = [xn_hist;Loc_PA(1,:)];
        obj_hist = [obj_hist,SINR_opt];
        obj_current = SINR_opt;
        if min(obj_current) >= min(obj_best)
            if min(obj_current) - min(obj_best) < epsilon
                obj_best = obj_current;
                Loc_PA_best = Loc_PA;
                W_best = W_opt;
                break;
            end
            obj_best = obj_current;
            Loc_PA_best = Loc_PA;
            W_best = W_opt;
        end
    end
    h = channel(Loc_PA_best, Loc_U, eta, lambda_c, lambda_g);
    G = abs((h*W_best)).^2;%[M,M]
    R_test = log2(1+min(obj_best));
    for m = 1:M
        SINR = G(m,m)/(sum(G(m,:))-G(m,m)+sigma2);
        R_MW_ct(m,it) = log2(1+SINR);
    end
end

%% MW_ZF
for it = 1:ct_ZF %Monte Carlo simulations 
    %初始化PA位置，均为0
    Loc_PA = zeros(3,N);
    for i=1:N
        Loc_PA(:,i) = [0;(i-1)*D_y/(N-1);d];
    end
    %生成用户位置
    Loc_U = zeros(3,M);
    for i=1:M
        Loc_U(:,i) = [D_x*rand(1)-D_x/2; D_y*rand(1); 0];
    end
    %这里把用户按y坐标从小到大排序，仅方便观察对应与PA的位置，实际运行时可以删掉
    [~, idx] = sort(Loc_U(2,:));   % 获取按第2行排序的列索引   
    Loc_U = Loc_U(:, idx);             % 按索引重排列

    obj_hist=[];
    xn_hist = [];
    obj_best = zeros(M,1);%记录最优迭代的目标函数
    for iter = 1:maxiter
        % ==== 外层搜索PA位置 ====
        Loc_PA_temp = Loc_PA;
        for n = 1:N
            xn_best = Loc_PA_temp(1,n);
            SINR_opt = zeros(M,1);%记录最优位置的目标函数
            %这里应该可以把范围缩小一些，在离这两个波导最近的两个位置中间搜，能减少迭代次数
            % 计算所有用户到第 i 个波导的y轴距离
            dist_y = Loc_U(2,:) - Loc_PA(2,n);      % K*1 - 1*1 = K*1
            [~, idx] = sort(abs(dist_y));            % 获得最小距离点索引
            x1 = min(Loc_U(1,idx(1)),Loc_U(1,idx(2)));
            x2 = max(Loc_U(1,idx(1)),Loc_U(1,idx(2)));
            for xn = candidates(candidates>=x1 & candidates<=x2)
                Loc_PA_temp(1,n) = xn;
                %PAn到用户m的信道：h[M,N] 
                h = channel(Loc_PA_temp, Loc_U, eta, lambda_c, lambda_g);
                % 内层固定PA位置，优化W                  
                W = pinv(h);
                W = W  * sqrt(P_max/trace(W'*W));
                SINR = zeros(M,1);
                G = abs((h*W)).^2;%[M,M]
                for m = 1:M
                    SINR(m) = G(m,m)/(sum(G(m,:))-G(m,m)+sigma2);
                end

                if min(SINR) > min(SINR_opt)
                    SINR_opt = SINR;
                    xn_best = xn;
                    W_opt = W;
                end
            end
            Loc_PA_temp(1,n) = xn_best;
        end
        Loc_PA = Loc_PA_temp;

        % 检查目标函数是否收敛
        xn_hist = [xn_hist;Loc_PA(1,:)];
        obj_hist = [obj_hist,SINR_opt];
        obj_current = SINR_opt;
        if min(obj_current) >= min(obj_best)
            if min(obj_current) - min(obj_best) < epsilon
                obj_best = obj_current;
                Loc_PA_best = Loc_PA;
                W_best = W_opt;
                break;
            end
            obj_best = obj_current;
            Loc_PA_best = Loc_PA;
            W_best = W_opt;
        end
    end
    h = channel(Loc_PA_best, Loc_U, eta, lambda_c, lambda_g);
    G = abs((h*W_best)).^2;%[M,M]
    R_test = log2(1+min(obj_best));
    for m = 1:M
        SINR = G(m,m)/(sum(G(m,:))-G(m,m)+sigma2);
        R_MW_ZF_ct(m,it) = log2(1+SINR);
    end
end

%% Conv_ZF
Loc_Conv = [zeros(1,N); D_y/2 + lambda_c/2 * [-N/2:N/2-1]; d * ones(1,N)]; 
for it = 1:ct_ZF %Monte Carlo simulations     
    %生成用户位置
    for i=1:M
        Loc_U(:,i) = [D_x*rand(1)-D_x/2; D_y*rand(1); 0];
    end
    R_Conv_ZF_ct(:,it) = ZF(M,N,f_c,lambda_c,lambda_g,sigma2,P_max,Loc_Conv,Loc_U);
  
end

R_MW_ct = sum(R_MW_ct,2)/(ct);
R_MW_ZF_ct = sum(R_MW_ZF_ct,2)/ct_ZF;
R_Conv_ZF_ct = sum(R_Conv_ZF_ct,2)/ct_ZF;
% SINR_hist
end

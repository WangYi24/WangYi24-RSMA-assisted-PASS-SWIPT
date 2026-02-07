function [Wp_best, wc_best, SumRate_best] = rsma_wmmse(h_IDR,h_EHR,sigma2, P_max,Emin)
% Rate-Splitting for Multi-User Multi-Antenna Wireless Information and Power Transfer
[Nt, K] = size(h_IDR);
[~, J] = size(h_EHR);

P_c = P_max * 0.25;
P_p = (P_max-P_c) / K;

% 生成列索引：1,2,...,M,1,2,... 直到 N 个（当 M>=N 时自动截断为 1:N）
col_idx = mod(0:K-1, J) + 1;
% 对 B 按列索引重排，得到 R x N 矩阵
h_EHR_resized = h_EHR(:, col_idx);

% [U, ~, ~] = svd(h_EHR_resized);
% Wp = h_EHR_resized./vecnorm(h_EHR_resized)  * sqrt(P_p);
% wc = U(:, 1) * sqrt(P_c);
% maxIter = 1000;
[U, ~, ~] = svd(h_IDR);
Wp = h_IDR./vecnorm(h_IDR)  * sqrt(P_p);
wc = U(:, 1) * sqrt(P_c);
maxIter = 100;

Wp_best = Wp;wc_best = wc;
SumRate_hist=[];

tolerance = 1e-1;
SumRate_best = -inf;
maxIter = 100;
flag_max=5;%flag_max次没出现更大的目标函数就break
flag = 0;

% 归一化
h_IDR = h_IDR / sqrt(sigma2);
h_EHR = h_EHR / sqrt(sigma2);
Emin = Emin / sigma2;

for n = 1:maxIter
    T_c = zeros(K,1); % T_c(k)为第k个协作用户的T_c值
    T = zeros(K,1); 
    g_p = zeros(K,1);
    g_c = zeros(K,1);
    u_p = zeros(K,1);
    u_c = zeros(K,1);
    for k = 1:K
        % 计算T值
        T(k) = sum_square_abs( Wp'*h_IDR(:, k)) + 1;
        T_c(k) = square_abs(h_IDR(:, k)'*wc) + T(k);

        %公式（8）
        g_p(k) = Wp(:, k)' * h_IDR(:, k) / T(k);
        g_c(k) = wc' * h_IDR(:, k) / T_c(k);

        %公式12
        u_p(k) = T(k) / (T(k) - square_abs(h_IDR(:, k)'*Wp(:, k))) / log(2);
        u_c(k) = T_c(k) / T(k) / log(2);
    end

    [Wp, wc, SumRate] = rsma_SCA_cvx(h_IDR,h_EHR, P_max,Emin,wc,Wp, g_c, g_p, u_c, u_p);
    %无解的情况下，使用对准EHR的初始化
    if n ==1 && (isnan(SumRate) || SumRate==-inf)
        [U, ~, ~] = svd(h_EHR_resized);
        Wp = h_EHR_resized./vecnorm(h_EHR_resized)  * sqrt(P_p);
        wc = U(:, 1) * sqrt(P_c);
        continue
    end
    if isnan(SumRate) || SumRate==-inf
        break;%
    end
    % SumRate_hist = [SumRate_hist,SumRate];
    % ehr_str = mat2str(EHR, 2);
    % fprintf("RSMA | %3d | obj = %.3f | EHR =%s | |obj - obj_last| = %.3f  \n", n, SumRate,ehr_str, abs(SumRate - SumRate_last));
    %发现这里会震荡
    if SumRate >= SumRate_best
        if SumRate - SumRate_best < tolerance
            SumRate_best = SumRate;
            Wp_best = Wp;wc_best=wc;
            % c_best = c;EHR_best=EHR;
            break;
        end
        SumRate_best = SumRate;
        Wp_best = Wp;wc_best=wc;
        % c_best = c;EHR_best=EHR;
        flag = 0;
    else
        flag=flag+1;
        if flag==flag_max
            break;
        end
    end
end
% SumRate_hist
end
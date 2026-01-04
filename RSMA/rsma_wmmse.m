function [W_p_best, wc_best, SumRate_best] = rsma_wmmse(h_IDR,h_EHR,sigma2, P_max,Emin)
% Rate-Splitting for Multi-User Multi-Antenna Wireless Information and Power Transfer
[Nt, K] = size(h_IDR);

P_c = P_max * 0.25;
P_p = (P_max-P_c) / K;

[U, ~, ~] = svd(h_IDR);
wc = U(:, 1) * sqrt(P_c);
W_p = h_IDR./vecnorm(h_IDR)  * sqrt(P_p);
W_p_best = W_p;wc_best = wc;
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
        T(k) = sum_square_abs( W_p'*h_IDR(:, k)) + 1;
        T_c(k) = square_abs(h_IDR(:, k)'*wc) + T(k);

        %公式（8）
        g_p(k) = W_p(:, k)' * h_IDR(:, k) / T(k);
        g_c(k) = wc' * h_IDR(:, k) / T_c(k);

        %公式12
        u_p(k) = T(k) / (T(k) - square_abs(h_IDR(:, k)'*W_p(:, k))) / log(2);
        u_c(k) = T_c(k) / T(k) / log(2);
    end

    [W_p, wc, SumRate] = rsma_SCA(h_IDR,h_EHR, P_max,Emin,wc,W_p, g_c, g_p, u_c, u_p);
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
            W_p_best = W_p;wc_best=wc;
            % c_best = c;EHR_best=EHR;
            break;
        end
        SumRate_best = SumRate;
        W_p_best = W_p;wc_best=wc;
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


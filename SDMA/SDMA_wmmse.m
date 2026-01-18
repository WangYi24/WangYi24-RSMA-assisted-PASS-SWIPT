function [W_p_best, SumRate_best] = SDMA_wmmse(h_IDR,h_EHR,sigma2, P_max,Emin)
% Rate-Splitting for Multi-User Multi-Antenna Wireless Information and Power Transfer
[Nt, K] = size(h_IDR);

W_p = h_IDR./vecnorm(h_IDR)  * sqrt(P_max);
W_p_best = W_p;
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
    T = zeros(K,1); 
    g_p = zeros(K,1);
    u_p = zeros(K,1);
    for k = 1:K
        % 计算T值
        T(k) = sum_square_abs( W_p'*h_IDR(:, k)) + 1;
        %公式（8）
        g_p(k) = W_p(:, k)' * h_IDR(:, k) / T(k);
        %公式12
        u_p(k) = T(k) / (T(k) - square_abs(h_IDR(:, k)'*W_p(:, k))) / log(2);
    end

    [W_p, SumRate] = SDMA_SCA_cvx(h_IDR,h_EHR, P_max,Emin,W_p, g_p, u_p);
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
            W_p_best = W_p;
            %EHR_best=EHR;
            break;
        end
        SumRate_best = SumRate;
        W_p_best = W_p;
        % EHR_best=EHR;
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


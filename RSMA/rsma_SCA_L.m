function [W_p_best, wc_best, SumRate_best] = rsma_SCA_L(h_IDR,h_EHR,P_max,Emin,wc,W_p, g_c, g_p, u_c, u_p)
[~, K] = size(h_IDR);
[Nt, J] = size(h_EHR);


maxIter = 10;
tolerance = 1e-1;
SumRate_best = -inf;
SumRate_hist = [];
W_p_best = zeros(Nt,K);wc_best= zeros(Nt,1);
% c_best = 0;EHR_best=0;

for it = 1:maxIter
    W_pre = W_p;
    wc_pre = wc;

    abs2_g_c = abs(g_c).^2;   
    abs2_g_p  = abs(g_p).^2;
    % A[J,1]公式15最后平方项,对k求和
    A = sum(abs(h_EHR'*W_pre).^2,2) + abs(h_EHR'*wc_pre).^2;
        
    for k = 1:K
        % 计算T值
        T(k) = sum_square_abs( h_IDR(:, k)'*W_pre) + 1;
        T_c(k) = square_abs(h_IDR(:, k)'*wc_pre) + T(k);
                
        % 公式14,计算epsilon值（MMSE）
        epsilon_c(k) = abs2_g_c(k)*T_c(k) - 2*real(g_c(k)*h_IDR(:, k)'*wc_pre) + 1;
        epsilon_p(k) = abs2_g_p(k)*T(k) - 2*real(g_p(k)*h_IDR(:, k)'*W_pre(:,k)) + 1;
        
        % 公式16,计算xi值（WMMSE）
        xi_c(k) = u_c(k)*epsilon_c(k) - log2(u_c(k));
        xi_p(k) = u_p(k)*epsilon_p(k) - log2(u_p(k));
    end
    %更新mu
    [~, k_star] = max(xi_c);
    mu = 1;
    %更新lambda
    
    %更新gamma
    for j = 1:J
        EHR(j) = 2*real(wc_pre'*h_EHR(:,j)*h_EHR(:,j)'*wc)- A(j);
        for k =1:K
            wk_hh = W_pre(:,k)'*h_EHR(:,j)*h_EHR(:,j)';
            EHR(j) = EHR(j) + 2*real(wk_hh*W_p(:,k));
        end
    end

    t = sum(x) + sum(xi_p- 1/log(2) - log2(log(2)));
    

    SumRate = -t;
    c = -x;
    if isnan(t)
        break;
    end
    % SumRate_hist = [SumRate_hist,SumRate];
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
    end    
end
end

    
    % h = h * sqrt(sigma2);
    % G = abs((h'*W)).^2;%[M,M]
    % R_c = zeros(K,1);
    % R_p = zeros(K,1);
    % for k = 1:K
    %     SINR = G(k,k)/(sum(G(k,:))-G(k,k)+sigma2);
    %     R_c(k)= log2(1+abs(h(:,k)'*wc).^2/(sum(G(k,:))+sigma2));
    %     R_p(k) = log2(1+SINR);
    % end
    % R = R_p+c;
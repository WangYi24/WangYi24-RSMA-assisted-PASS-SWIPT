function [W_p_best, wc_best, SumRate_best] = rsma_SCA_L(h_IDR,h_EHR,Pmax,Emin,wc,Wp, g_c, g_p, u_c, u_p)
[~, K] = size(h_IDR);
[Nt, J] = size(h_EHR);


maxIter = 10;
tolerance = 1e-1;
SumRate_best = -inf;
SumRate_hist = [];
W_p_best = zeros(Nt,K);wc_best= zeros(Nt,1);
% c_best = 0;EHR_best=0;
gamma = 1/trace(h_EHR'*h_EHR) * ones(1,J);

for it = 1:maxIter
    Wp_pre = Wp;
    wc_pre = wc;

    abs2_g_c = abs(g_c).^2;   
    abs2_g_p  = abs(g_p).^2;
    % A[J,1]公式15最后平方项,对k求和
    A = sum(abs(h_EHR'*Wp_pre).^2,2) + abs(h_EHR'*wc_pre).^2;
        
    for k = 1:K
        % 计算T值
        T(k) = sum_square_abs( h_IDR(:, k)'*Wp_pre) + 1;
        T_c(k) = square_abs(h_IDR(:, k)'*wc_pre) + T(k);
                
        % 公式14,计算epsilon值（MMSE）
        epsilon_c(k) = abs2_g_c(k)*T_c(k) - 2*real(g_c(k)*h_IDR(:, k)'*wc_pre) + 1;
        epsilon_p(k) = abs2_g_p(k)*T(k) - 2*real(g_p(k)*h_IDR(:, k)'*Wp_pre(:,k)) + 1;
        
        % 公式16,计算xi值（WMMSE）
        xi_c(k) = u_c(k)*epsilon_c(k) - log2(u_c(k));
        xi_p(k) = u_p(k)*epsilon_p(k) - log2(u_p(k));
    end
    %更新mu
    [~, k_star] = max(xi_c);
    %更新lambda
    [lambda,W] = bisect_lambdas(h_IDR,h_EHR,  wc_pre,Wp_pre, g_c, g_p, u_c, u_p, gamma, k_star, Pmax);
    %更新gamma
    wc = W(:,1);Wp = W(:,2:end);
    for j = 1:J
        EHR(j) = 2*real(wc_pre'*h_EHR(:,j)*h_EHR(:,j)'*wc)- A(j);
        for k =1:K
            wk_hh = Wp_pre(:,k)'*h_EHR(:,j)*h_EHR(:,j)';
            EHR(j) = EHR(j) + 2*real(wk_hh*Wp(:,k));
        end
    end
    grand = EHR - Emin;



    % SumRate_hist = [SumRate_hist,SumRate];
    %发现这里会震荡
    if SumRate >= SumRate_best
        if SumRate - SumRate_best < tolerance
            SumRate_best = SumRate;
            W_p_best = Wp;wc_best=wc;
            % c_best = c;EHR_best=EHR;
            break;
        end
        SumRate_best = SumRate;
        W_p_best = Wp;wc_best=wc;
        % c_best = c;EHR_best=EHR;
    end    
end
end

function [lambda,W] = bisect_lambdas(h_IDR,h_EHR,  wc_pre, Wp_pre, g_c, g_p, u_c, u_p, gamma, k_star, Pmax)
   lambda_low = 0; lambda_high = 1e6;
    for iter = 1:50
        lambda = (lambda_low + lambda_high)/2;
        W = compute_w(h_IDR,h_EHR, wc_pre,Wp_pre, g_c, g_p, u_c, u_p, lambda, gamma, k_star);
        P = trace(W'*W);
        
        if P > Pmax
            lambda_low = lambda;
        elseif (Pmax-P)/Pmax < 1e-3
            break;
        else
            lambda_high = lambda;
        end
    end
end

function W = compute_w(h_IDR,h_EHR, wc_pre,Wp_pre, g_c, g_p, u_c, u_p, lambda, gamma, k_star)

    [N,K]=size(h_IDR);
    [~,J]=size(h_EHR);
    Wp = zeros(N,K);

    %wc
    A_c_k_star = u_c(k_star) * abs(g_c(k_star))^2 * h_IDR(:,k_star) * h_IDR(:,k_star)';
    b_c_k_star = u_c(k_star)' * g_c(k_star)' * h_IDR(:,k_star);
    c = trace(diag(gamma) * (h_EHR'*h_EHR));
    wc = (A_c_k_star + lambda * eye(N)) \ (b_c_k_star + c * wc_pre);

    %Wp
    sumA = trace(diag(u_p .* abs(g_p).^2) * (h_IDR'*h_IDR));
    for k = 1:K
        b_p_k = u_p(k)' * g_p(k)' * h_IDR(:,k);
        Wp(:,k) = (sumA + A_c_k_star + lambda * eye(N)) \ (b_p_k + c * Wp_pre(:,k));
    end
    W = [wc,Wp];
end
function [W_p_best, SumRate_best] = SDMA_SCA_cvx(h_IDR,h_EHR,P_max,Emin,W_p, g_p, u_p)
[~, K] = size(h_IDR);
[Nt, J] = size(h_EHR);

maxIter = 10;
tolerance = 1e-1;
SumRate_best = -inf;
SumRate_hist = [];
W_p_best = zeros(Nt,K);
% EHR_best=0;
for it = 1:maxIter
    W_pre = W_p;
 
    abs2_g_p  = abs(g_p).^2;
    % A[J,1]SCA最后平方项,对k求和
    A = sum(abs(h_EHR'*W_pre).^2,2);
    cvx_solver mosek
    cvx_begin quiet
        variable W_p(Nt,K) complex
       
        expression t
        expression T(K)
        expression epsilon_p(K)
        expression xi_p(K)  
        expression EHR(J) 
        
        for k = 1:K
            % 计算T值
            T(k) = sum_square_abs( h_IDR(:, k)'*W_p) + 1;           
            % 公式14,计算epsilon值（MMSE）
            epsilon_p(k) = abs2_g_p(k)*T(k) - 2*real(g_p(k)*h_IDR(:, k)'*W_p(:,k)) + 1;        
            % 公式16,计算xi值（WMMSE）
            xi_p(k) = u_p(k)*epsilon_p(k) - log2(u_p(k));
        end

        for j = 1:J
            EHR(j) = - A(j);
            for k =1:K
                wk_hh = W_pre(:,k)'*h_EHR(:,j)*h_EHR(:,j)';
                EHR(j) = EHR(j) + 2*real(wk_hh*W_p(:,k));
            end
        end

        t = sum(xi_p- 1/log(2) - log2(log(2)));
    
        minimise(t)
        subject to            
            % 功率约束
            sum(sum_square_abs(W_p)) <= P_max  % 所有用户的总功率不超过P_max
            % 能量收集约束
            EHR >= Emin
    cvx_end

    SumRate = -t;
    if isnan(t)
        break;
    end
    % SumRate_hist = [SumRate_hist,SumRate];
    %发现这里会震荡
    if SumRate >= SumRate_best
        if SumRate - SumRate_best < tolerance
            SumRate_best = SumRate;
            W_p_best = W_p;
            % EHR_best=EHR;
            break;
        end
        SumRate_best = SumRate;
        W_p_best = W_p;
        % EHR_best=EHR;
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
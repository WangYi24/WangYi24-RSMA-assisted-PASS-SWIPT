function [P,t]=getP_alloc(A,b,Pmax)
% A[M,K]：abs(h(m,:,k)*W(:,k)).^2/(deltaf*N0)
% b[M,K]:子载波分配
% df:deltaf
% Pmax

% P[M,1]:功率分配
% t:最小速率
[M,K] = size(A);
km = sum(b,2);          % [M,1]每个用户的子载波数量

cvx_begin quiet
    variables P(M) t
    maximize( t )
    subject to
        sum(P) <= Pmax;
        P >= 0;
        for m = 1:M
            % 计算 Rm(P_m)
            idx = find(b(m,:)==1);
            % Rm = df/log(2) * sum( log(1 + A(m,idx) * (P(m)/km(m))) )
            Rm = 0;
            for kk = idx
                % x = P(m)/km(m) * A(m,kk);
                % Rm = Rm - rel_entr(1,x);
                Rm = Rm + log( 1 +  P(m)/km(m) * A(m,kk) );                
            end
            Rm >= t;
        end
cvx_end

% 输出：P (用户总功率), t (最大化得到的 min-rate)

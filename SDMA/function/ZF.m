function [R]=ZF(M,N,lambda_c,lambda_g,sigma2,P_max,Loc_ant,Loc_U)
% M:用户数
% N:天线数
% d:天线高度
% K:子载波数量
% f_c:载波频率
% lambda_c:波长
% N0：噪声功率谱密度
% P_max:发射功率
% Loc_U:用户位置

    eta = lambda_c/4/pi;
    R = zeros(M,1);

    % 信道h[M,N]
    
    h = channel(Loc_ant, Loc_U, eta, lambda_c, lambda_g);
    %ZF预编码
    W_ZF = pinv(h);%想要得到右逆必须行满秩也就是列大于行
    W_ZF = sqrt(P_max / trace(W_ZF'*W_ZF)) * W_ZF;
    G = abs(h * W_ZF).^2;
    % W_ZF = sqrt(P_max/K) * W_ZF./vecnorm(W_ZF);
    for m = 1:M
        R(m) = R(m) + log2(1 + G(m,m)/(sum(G(m,:))-G(m,m)+sigma2));
    end

end

function [R]=MW_OFDM_ZF(M,N,B,K,f_c,lambda_c,lambda_g,N0,P_max,Loc_PA,Loc_U)
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
    n_eff = lambda_c/lambda_g;
    eta_B=1;
    delta_f = B/K;
    [eta_B,K,delta_f] = get_OFDMAframe(Loc_PA,Loc_U,B,n_eff);

    % 信道h[M,N,K]
    h = channel_OFDM(Loc_PA, Loc_U, eta, lambda_c, lambda_g, f_c, K, B, 1);
    %遍历每个子载波
    for k = 1:K
        %ZF预编码
        H_ZF = h(:,:,k);
        W_ZF = pinv(H_ZF);%想要得到右逆必须行满秩也就是列大于行
        W_ZF = sqrt(P_max/K / trace(W_ZF'*W_ZF)) * W_ZF;
        H_ZF * W_ZF;
        % W_ZF = sqrt(P_max/K) * W_ZF./vecnorm(W_ZF);
        for m = 1:M
            R(m) = R(m) + delta_f*log2(1 + abs(H_ZF(m, :) * W_ZF(:, m)).^2/N0/delta_f);
        end
    end
    R = R * eta_B;
end

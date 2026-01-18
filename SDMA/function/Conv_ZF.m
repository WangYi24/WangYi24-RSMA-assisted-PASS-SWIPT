function [R]=Conv_ZF(M,N,lambda_c,sigma2,P_max,Loc_ant,Loc_U)
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

    % 信道h[M,N,K]
    h = channel_Conv(Loc_ant, Loc_U, eta, lambda_c);
    %ZF预编码
    W_ZF = pinv(h);%想要得到右逆必须行满秩也就是列大于行
    W_ZF = sqrt(P_max / trace(W_ZF'*W_ZF)) * W_ZF;
    G = abs(h * W_ZF).^2;
    % W_ZF = sqrt(P_max/K) * W_ZF./vecnorm(W_ZF);
    for m = 1:M
        R(m) = R(m) + log2(1 + G(m,m)/(sum(G(m,:))-G(m,m)+sigma2));
    end

end

%function [R_ct]=Conv_ZF(M,N,D,d,B,K,f_c,lambda_c,lambda_g,sigma2,P_max)
% % M:用户数
% % N:天线数
% % d:天线高度
% % K:子载波数量
% % f_c:载波频率
% % lambda_c:波长
% % sigma2：噪声
% % P_max:发射功率
% % Loc_U:用户位置
%     ct = 1e2;%重复次数
%     c = 3e8;%光速
%     eta = lambda_c/4/pi;
%     D_x = D;
%     D_y = D;   
%     R_ct = zeros(M,ct);
%     N0 = sigma2/B;%噪声功率谱密度
%     deltaf = B/K;
%     %传统天线位置,[3,N]
%     Loc_Conv = [zeros(1,N); D/2 + lambda_c/2 * (-N/2:N/2-1); d * ones(1,N)];
%     for it = 1:ct       
%         %生成用户位置
%         for m=1:M
%             Loc_U(:,m) = [D_x*rand(1)-D_x/2; D_y*rand(1); 0];
%         end
%         % 信道h[M,N,K]
%         h = channel_OFDM(Loc_Conv, Loc_U, eta, lambda_c, lambda_g, f_c, K, B, 0);
%         %遍历每个子载波
%         for k = 1:K
%             %ZF预编码
%             H_ZF = h(:,:,k);
%             W_ZF = pinv(H_ZF);%想要得到右逆必须行满秩也就是列大于行
%             W_ZF = sqrt(P_max/K / trace(W_ZF'*W_ZF)) * W_ZF;
%             H_ZF * W_ZF;
%             % W_ZF = sqrt(P_max/K) * W_ZF./vecnorm(W_ZF);
%             for m = 1:M
%                 R_ct(m,it) = R_ct(m,it) + deltaf*log2(1 + abs(H_ZF(m, :) * W_ZF(:, m)).^2/N0/deltaf);
%             end
%         end
% 
%     end
%     R_ct = sum(R_ct,2)/ct;
% end

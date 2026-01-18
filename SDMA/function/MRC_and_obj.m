function [W,P_alloc,R] = MRC_and_obj(h,b,delta_f,P_max,N0)  
%不怎么用了，把W，P，R合在一起不易调试
%h[M,N,K]:信道
%b[M,K]:子载波分配策略
%N0：噪声功率谱密度

%W[N.K]:每一列是子载波k的预编码
%R[M,1]:M个用户的速率
    [M,N,K] = size(h);
    W = zeros(N,K);%每一列是子载波k的预编码
    g = zeros(M,K);%信道增益g(m,k)=abs(h(m,:,k)*W(:,k)).^2
    R = zeros(M,1);
    R_1 = zeros(M,1);%均分功率
    for m = 1:M
        for k = 1:K
            if b(m,k)==1
                W(:,k) = h(m,:,k)'/norm(h(m,:,k));
                g(m,k) = abs(h(m,:,k)*W(:,k)).^2;
            end
        end
    end

    [P_alloc,t]=getP_alloc(g/(delta_f*N0),b,P_max);
    for m=1:M
         R(m) = sum( delta_f * log2( 1 + (g(m,:) .* P_alloc(m)/sum(b(m,:))) / (delta_f*N0) ) );
         R_1(m) = sum( delta_f * log2( 1 + (g(m,:) .* P_max/K) / (delta_f*N0) ) );
    end
   
end




 % %子载波功率注水
    % Pm = P_max/M;
    % 
    % %对每个用户的子载波进行注水算法
    % p_alloc = zeros(1,K);
    % R = zeros(M,1);
    % R_1 = zeros(M,1);%子载波等功率速率
    % for m = 1:M
    %     idx = find(b(m,:)==1);      % 该用户的子载波索引
    %     L = length(idx);%返回用户m子载波个数
    %     if L==0
    %         continue; 
    %     end
    % 
    %     gm = g(m, idx).';           % L，1
    %     % 注水算法：对 1/SNR 升序搜索水位
    %     snr_inv = (delta_f*N0) ./ gm;      % 1/SNR
    %     [snr_inv_sorted, ord] = sort(snr_inv, 'ascend');%升序排列
    % 
    %     % 水位mu搜索
    %     water = 0;
    %     for t = 1:L
    %         mu = (Pm + sum(snr_inv_sorted(1:t))) / t;
    %         if t==L || mu <= snr_inv_sorted(t+1)
    %             water = mu; break;
    %         end
    %     end
    %     % 分配功率
    %     p_mk_sorted = max(0, water - snr_inv_sorted);
    %     p_mk = zeros(1,L);
    %     p_mk(ord) = p_mk_sorted;
    %     p_alloc(idx) = p_mk;
    %     %计算速率        
    %     R(m) = sum( delta_f * log2( 1 + (gm .* p_mk') / (delta_f*N0) ) );
    % 
    %     % % 子载波功率均分
    %     % R_1(m) = sum( delta_f * log2( 1 + (gm .* P_max/K) / (delta_f*N0) ) );
    %     % p_alloc = P_max/K;
    % end   
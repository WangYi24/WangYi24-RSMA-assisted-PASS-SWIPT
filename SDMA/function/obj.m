function [R] = obj(h,b,delta_f,P_M,N0)  
%h[M,N,K]:信道
%b[M,K]:子载波分配策略
%N0：噪声功率谱密度

%W[N.K]:每一列是子载波k的预编码
%R[M,1]:M个用户的速率
    [M,N,K] = size(h);
    W = zeros(N,K);%每一列是子载波k的预编码
    [M,N,K] = size(h);
    R = zeros(M,1);
    for m = 1:M
        for k = 1:K
            if b(m,k)==1
                W(:,k) = h(m,:,k)'/norm(h(m,:,k));
                g = (abs(h(m,:,k)* W(:,k))^2 .* P_M(m)/sum(b(m,:))) / (delta_f*N0);
                R(m) = R(m) + delta_f * log2( 1 + (abs(h(m,:,k)* W(:,k))^2 .* P_M(m)/sum(b(m,:))) / (delta_f*N0) ) ;
            end
        end
    end
end
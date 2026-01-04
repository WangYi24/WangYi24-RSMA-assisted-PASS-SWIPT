function [W_p, w_c, c, mmf_rate, rate_p] = gpi_method(h, sigma2, P_max)
% 输出：
% P[N,K]:私有信号的预编码
% p_c:共有信号的预编码
% c:
% mmf_rate:

[Nt, M] = size(h);
% 噪声归一化
h = h / sqrt(sigma2);
noise = ones(M, 1);

%初始化功率分配
P_c = P_max * 0.5;
P_p = (P_max-P_c) / M;

[U, ~, ~] = svd(h);
w_c = U(:,1) * sqrt(P_c);
W_p = h ./ vecnorm(h) * sqrt(P_p);%MRT

gamma_min = 0;
gamma_max = 100;%100
gamma_num = 1001;

alpha = 0.1;
tolerance = 1e-3;

T_m = sum(square_abs(h'*W_p), 2) + noise;
T_cm = T_m + square_abs(h'*w_c);

rate_c = log2(T_cm ./ T_m);
rate_p = log2(T_m ./ (T_m-square_abs(diag(h'*W_p))));

rate_p_sorted = sort(rate_p, 'ascend');
x = min(rate_c);
y = min((x+cumsum(rate_p_sorted)) ./ (1:M)');

max_mmf_rate = y;
corresponding_f = [w_c; reshape(W_p, M*Nt, 1)];

f = corresponding_f;
c = max(y - rate_p, 0);
f_last = f;

c_hist = [];
y_hist = [];
rate_p_hist = [];

for gamma = linspace(gamma_min, gamma_max, gamma_num)
    % 一阶段，迭代波束赋形
    for n = 1:1000
        % 公式13~16
        A = zeros(Nt*(M+1), Nt*(M+1), M);
        B = zeros(Nt*(M+1));
        C = zeros(Nt*(M+1));
        D = zeros(Nt*(M+1));

        for k = 1:M
            A_k = zeros(Nt);
            for m = 1:M
                A_k = blkdiag(A_k, h(:, k)*h(:, k)');
            end
            A(:, :, k) = A_k + eye(Nt*(M+1))/P_max;
            B(:, :, k) = A(:, :, k) - blkdiag(zeros(Nt*k), h(:, k)*h(:, k)', zeros(Nt*(M-k)));
            C(:, :, k) = A(:, :, k) + blkdiag(h(:, k)*h(:, k)', zeros(Nt*M));
            D(:, :, k) = A(:, :, k);
        end

        T_m = sum(square_abs(h'*W_p), 2) + noise;
        T_cm = T_m + square_abs(h'*w_c);
        
        rate_c = log2(T_cm ./ T_m);
        rate_p = log2(T_m ./ (T_m-square_abs(diag(h'*W_p))));
        % 公式19~20的一部分
        ratio_c = exp(-1/alpha*(rate_c)) / sum(exp(-1/alpha*(rate_c)));
        ratio_t = exp(-1/alpha*(c+rate_p)) / sum(exp(-1/alpha*(c+rate_p)));

        E = zeros(Nt*(M+1));
        F = zeros(Nt*(M+1));
        for k = 1:M
            E = E + ratio_t(k)*A(:, :, k)/(f'*A(:, :, k)*f) + ratio_c(k)*gamma*C(:, :, k)/(f'*C(:, :, k)*f);
            F = F + ratio_t(k)*B(:, :, k)/(f'*B(:, :, k)*f) + ratio_c(k)*gamma*D(:, :, k)/(f'*D(:, :, k)*f);
        end
        % 算法中更新f的公式
        f = F \ (E*f);
        f = f / norm(f) * sqrt(P_max);

        e = norm(f - f_last, 'inf');
        if e < tolerance
            break;
        end

        f_last = f;
    end

    %二阶段通用速率分配，注水算法
    w_c = f(1:Nt);
    W_p = reshape(f(Nt+1:end), Nt, M);

    T_m = sum(square_abs(h'*W_p), 2) + noise;
    T_cm = T_m + square_abs(h'*w_c);
    
    rate_c = log2(T_cm ./ T_m);
    rate_p = log2(T_m ./ (T_m-square_abs(diag(h'*W_p))));

    rate_p_sorted = sort(rate_p, 'ascend');
    x = min(rate_c);
    %确定水位y，感觉怪怪的
    y = min((x+cumsum(rate_p_sorted)) ./ (1:M)');

    c = max(y - rate_p, 0);
    c_hist = [c_hist,c];
    y_hist = [y_hist,y];
    rate_p_hist = [rate_p_hist,rate_p];

    if y > max_mmf_rate
        max_mmf_rate = y;
        corresponding_f = f;
    end

end


f = corresponding_f;
w_c = f(1:Nt);

W_p = reshape(f(Nt+1:end), Nt, M);

T_m = sum(square_abs(h'*W_p), 2) + noise;
T_cm = T_m + square_abs(h'*w_c);

rate_c = log2(T_cm ./ T_m);
rate_p = log2(T_m ./ (T_m-square_abs(diag(h'*W_p))));

rate_p_sorted = sort(rate_p, 'ascend');
x = min(rate_c);
y = min((x+cumsum(rate_p_sorted)) ./ (1:M)');
c = max(y - rate_p, 0);
mmf_rate = y;

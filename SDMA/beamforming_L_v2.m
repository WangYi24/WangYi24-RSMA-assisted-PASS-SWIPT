function [W, SINR_hist] = beamforming_L_v2(h, sigma2, Pmax, opts)
% 每次迭代先更新λ

% h: [N, M] 信道向量 (第 m 列是 h_m)
% sigma2: 噪声功率
% Pmax: 总功率约束
% opts: 可选结构体 (maxIter, tol, verbose)

if nargin < 4, opts = struct(); end
if ~isfield(opts, 'maxIter'), opts.maxIter = 50; end
if ~isfield(opts, 'tol'), opts.tol = 1e-2; end
if ~isfield(opts, 'verbose'), opts.verbose = true; end

[N, M] = size(h);

% === 初始化 ===
% W = pinv(h);
% W = W' * sqrt(Pmax/trace(W'*W));
W = h* sqrt(Pmax/trace(h'*h));
SINR_hist = [];

y_hist=[];
for iter = 1:opts.maxIter
% ==== Step 1: 固定 W，更新 y ====    
    y = compute_y(h, W, sigma2);
    y_hist=[y_hist,abs(y)];
    %有时候会出现ym=0的情况，视为坏点直接扔了
    if any(abs(y)<1e-6)
        break;
    end

% ==== Step 2: 固定 y，更新 W ====
    % 初始化
    lambda = ones(M,1) / M;   % 初值: 满足 sum λ = 1
    maxlambda = 1000;
    tollambda = 1e-3;    
    % --- 迭代 λ ---
    lambda_hist=[];
    u_hist=[];
    eta = 0.5;  % 初始学习率
    decay = 0.9; %弄一个衰减，让步长减少减少震荡
    min_urange = inf;
    for lambda_iter = 1:maxlambda
        % 对给定 λ,二分法求mu
        [mu,W] = bisect_mu(h, y, lambda, Pmax, sigma2);
        if trace(W'*W)>Pmax
            warning('二分失败');
            keyboard
        end
        % 计算各约束 u
        u_temp = compute_u(y,h, W, sigma2);
        urange = range(u_temp);
        u_hist=[u_hist,u_temp];
        
        % --- 更新 λ（保证 sum λ = 1） ---
        g = mean(u_temp) - u_temp;
        g = g / max(abs(g));
        eta = eta * decay;
        lambda = lambda .* exp( eta * g ); % 元素乘
        lambda = lambda / sum(lambda);           % 归一化
        if isnan(lambda)
            keyboard;
        end
        lambda_hist = [lambda_hist,lambda];
    
        if urange<min_urange            
            if (min_urange-urange)<1e-2
                break;
            end
            min_urange = urange;
        end
    end
% 收敛判定
    SINR = compute_SINR(h, W, sigma2);
    SINR_hist = [SINR_hist,SINR];
    %发现最后会震荡，怎么处理这样的震荡呢？
    if range(SINR) < opts.tol
        break;
    end
end
    % plot(1:20,log2(min(SINR_hist(:,1:20))),'bo-')
    % xlabel('Number of iterations')
    % ylabel('Minimum date rate(bps/Hz)')
end


%% 
function [mu,W] = bisect_mu(h, y, lambda, Pmax, sigma2)
    mu_low = 0; mu_high = 1e6;
    for iter = 1:50
        mu = (mu_low + mu_high)/2;
        W = compute_w(h, y, lambda,mu);
        P = trace(W'*W);
        
        if P > Pmax
            mu_low = mu;
        elseif (Pmax-P)/Pmax < 1e-3
            break;
        else
            mu_high = mu;
        end
    end
end

function W = compute_w(h, y, lambda,mu)
    [N,M]=size(h);
    W = zeros(N,M);
    for m = 1:M
        Hsum = zeros(N,N);
        for i = 1:M
            if i ~= m
                Hsum = Hsum + lambda(i)*abs(y(i))^2*(h(:,i)*h(:,i)');
            end
        end
        A = mu*eye(N) + Hsum;

        W(:,m) = lambda(m)*y(m)*(A\h(:,m));
    end
end

function y = compute_y(h, W, sigma2)
    [N,M]=size(h);
    y = zeros(M,1);
    G = abs((h'*W)).^2;%[M,M]
    for m = 1:M
        interf = 0;
        for i = 1:M
            if i ~= m
                interf = interf + abs(h(:,m)' * W(:,i))^2;
            end            
        end
        interf = sum(G(m,:))-G(m,m);
        y(m) = (h(:,m)' * W(:,m)) / (interf + sigma2);
    end
end

function u_temp = compute_u(y,h, W, sigma2)
    [N,M]=size(h);
    u_temp = zeros(M,1);
    G = abs((h'*W)).^2;%[M,M]
    for m = 1:M
        interf = sum(G(m,:))-G(m,m);
        % u_temp(m) = h(:,m)' * W(:,m)/(interf + sigma2);
        u_temp(m) = 2*real(y(m)' * (h(:,m)' * W(:,m))) ...
             - abs(y(m))^2 * (interf + sigma2);
    end
end

function SINR = compute_SINR(h, W, sigma2)
    [N,M]=size(h);
    SINR = zeros(M,1);
    G = abs((h'*W)).^2;%[M,M]
    for m = 1:M
        SINR(m) = G(m,m)/(sum(G(m,:))-G(m,m)+sigma2);
    end
end

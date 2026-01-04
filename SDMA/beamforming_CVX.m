function [t_opt, W_opt, out] = beamforming_CVX(h, sigma2, P_dBm_max, opts)
% h: [M,N]
% sigma2: 
% Ptot: scalar
% opts: struct with fields maxIter, tol, verbose
%
% returns t_opt (scalar), W_opt (N x M), out.history

if nargin < 4, opts = struct(); end
t_old = -inf;
if ~isfield(opts,'tol'), opts.tol = 1e-3; end
if ~isfield(opts,'maxIter'), opts.maxIter = 50; end
if ~isfield(opts,'verbose'), opts.verbose = false; end


[M, N] = size(h);
% init W (MRT)
P_max = 10.^((P_dBm_max-30)./10);

% 归一化
h_scale = trace(h*h');
h = h / sqrt(h_scale);
sigma2 = sigma2 / h_scale;
W = h';
W = W *  sqrt(P_max / trace(W*W'));

out.history = [];

for iter = 1:opts.maxIter
    % 1) update u
    G = h * W;               % M x M, G(m,l) = h(m,:) * W(:,l)
    S = abs(G).^2;           
    rowSum = sum(S,2);       % Mx1
    u = zeros(1,M);
    for m=1:M
        denom = rowSum(m) - S(m,m)+ sigma2;
        u(m) = G(m,m) / denom;   % complex scalar
    end

    % 2）update W,t
    cvx_solver mosek
    cvx_begin quiet
    variable t nonnegative              % 目标值
    variable Wvar(N, M) complex         % 所有波束，w_m 是第 m 列

    % 构造 M 个二次变换约束
    expression constr(M)
    for m = 1:M
        hw = h(m,:)*Wvar;               % 1 x M 向量，hw(l) = h_m * w_l
        % 2*R{u_m*h_m w_m}
        num = 2*real(u(m)*hw(m));    
        % 分母项：sum_{l≠m} |h_m w_l|^2 + sigma2
        den = 0;
        for l = 1:M
            if l ~= m
                den = den + square_pos(abs(hw(l)));
            end            
        end
        den = den + sigma2;
        constr(m) = num - abs(u(m))^2 * den;
    end

    % 目标：最大化 t
    maximize t
    subject to
        constr >= t;   % M 个约束
        square_pos(norm(Wvar,'fro')) <= P_max;
cvx_end
    % get solution

    if strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved')
        W = Wvar;
        t_new = t;
    else
        warning('CVX not solved: %s. Stopping.', cvx_status);

        break;
    end

    out.history = [out.history, double(t_new)];
    if opts.verbose
        fprintf('iter %d: t = %.3f\n', iter, t_new);
    end

    if abs(double(t_new) - double(t_old)) < opts.tol
        break;
    end
    t_old = double(t_new);
end

t_opt = double(t_new);
W_opt = W;
end
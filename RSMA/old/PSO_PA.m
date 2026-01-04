function [x_best, rate_best, convergence,EHR_gbest] = PSO_PA(W, Loc_PA,Loc_IDR,Loc_EHR, Dx, sigma2,Emin)
% PSO to optimize antenna x-positions on N parallel waveguides
% Inputs:
%   W          : N x K complex precoding matrix (fixed)
%   Loc_IDR   : 3*K matrix, [x, y] positions of K users
%   Dx         : total length of x-movement range (antenna in [-Dx/2, Dx/2])
%   sigma2     : noise power

options = struct();

% Default PSO parameters
num_particles = getOption(options, 'num_particles', 30);
max_iter      = getOption(options, 'max_iter', 100);
w             = getOption(options, 'w', 0.5);        % inertia weight
c1            = getOption(options, 'c1', 0.25);       % cognitive
c2            = getOption(options, 'c2', 0.25);       % social
lb            = -Dx/2;
ub            = Dx/2;

[~,K]=size(Loc_IDR);
[~,J]=size(Loc_EHR);
[N,~]=size(W);
% Initialize particles: each particle is an N-dimensional x-position vector

xn = lb + (ub - lb) * rand(num_particles-1, N);   % positions
xn = [Loc_PA(1,:);xn];
v = zeros(num_particles, N);                   % velocities
xn_best = xn;                                 % personal best positions
val_best = -inf(num_particles, 1);              % personal best values
rho_penalty = 20/Emin;

% Evaluate initial fitness
for i = 1:num_particles
    Loc_PA(1,:) = xn(i,:);
    [rate,~] = compute_fit(Loc_PA, W, Loc_IDR,Loc_EHR, sigma2,Emin,rho_penalty);    
    val_best(i) = rate;
end

[val_gbest, val_gbest_idx] = max(val_best);
xn_gbest = xn_best(val_gbest_idx, :);

convergence = zeros(max_iter, 1);

% Main PSO loop
for iter = 1:max_iter
    for i = 1:num_particles
        % Update velocity
        r1 = rand(1, N);
        r2 = rand(1, N);
        v(i,:) = w * v(i,:) ...
                 + c1 * r1 .* (xn_best(i,:) - xn(i,:)) ...
                 + c2 * r2 .* (xn_gbest - xn(i,:));
        
        % Update position
        xn(i,:) = xn(i,:) + v(i,:);
        
        % Enforce bounds
        xn(i,:) = min(max(xn(i,:), lb), ub);
        
        % Evaluate fitness
        Loc_PA(1,:) = xn(i,:);
        [rate,EHR] = compute_fit(Loc_PA, W, Loc_IDR,Loc_EHR, sigma2,Emin,rho_penalty);
        
        % Update personal best
        if rate > val_best(i)
            val_best(i) = rate;
            xn_best(i,:) = xn(i,:);
        end
        
        % Update global best
        if val_best(i) >= val_gbest
            val_gbest = val_best(i);
            xn_gbest = xn_best(i,:);
            EHR_gbest = EHR;
        end
    end
    
    convergence(iter) = val_gbest;
    
    if mod(iter, 20) == 0
        fprintf('Iter %d: Best sum rate = %.3f bps/Hz\n', iter, val_gbest);
    end
end

x_best = xn_gbest;
rate_best = val_gbest;

end

% -------------------------------------------------------------------------
function [fitness,EHR] = compute_fit(Loc_PA, W, Loc_IDR,Loc_EHR, sigma2,Emin,rho)
% x_pos: 1 x N
% y_pos: N x 1
% W: N x K
% Loc_IDR: 3*K
    wc = W(:,1);
    W_p = W(:,2:end);
    [N,K]=size(W_p);
    
    h_IDR = channel(Loc_PA, Loc_IDR);
    h_EHR = channel(Loc_PA, Loc_EHR);
    
    
    term_c = abs(h_IDR' * wc).^2;   % [K, 1]
    G = abs((h_IDR'*W_p)).^2;%得到 [K, K]
    c = log2(1+term_c./sum(G, 2));
    R_p = zeros(K,1);
    for k = 1:K
        SINR = G(k,k)/(sum(G(k,:))-G(k,k)+sigma2);
        %R_c(k)= log2(1+abs(h_IDR(:,k)'*wc).^2/(sum(G(k,:))+sigma2));
        R_p(k) = log2(1+SINR);
    end
    sum_rate = sum(R_p+c);

    term_c = abs(h_EHR' * wc).^2;     % [J, 1]
    h_W = h_EHR' * W_p;              % [J, K]
    term_p = sum(abs(h_W).^2, 2);  % 对每一行求和，得到 [J, 1]
    EHR= term_p + term_c;
    
    penalty = sum(Emin - EHR(EHR<Emin));
    if penalty > 0
        fitness = sum_rate - rho* penalty;
    else
        fitness = sum_rate;  % feasible solution
    end
end



% -------------------------------------------------------------------------
function val = getOption(s, field, default)
if isfield(s, field)
    val = s.(field);
else
    val = default;
end
end
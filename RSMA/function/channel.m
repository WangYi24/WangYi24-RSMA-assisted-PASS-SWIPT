function [h] = channel(Loc_ant, Loc_U)
% Loc_ant：3*N，PA位置
% Loc_U：3*M，用户位置
% lambda_c：自由空间波长
% lambda_p：波导内波长

%h:[M,N,K]
    
    f_c = 28*1e9;%载波频率Hz
    c = 3e8;%光速m/s
    n_eff = 1.4;%波导折射率
    lambda_c = c/f_c;
    lambda_g = lambda_c / n_eff;
    eta = lambda_c/4/pi;%自由空间损失因子
    [~,N] = size(Loc_ant);
    [~,M] = size(Loc_U);
    h = zeros(N,M);
        %公式7
    for m=1:M
        for n=1:N
            dist = norm(Loc_U(:,m)-Loc_ant(:,n));
            h(n,m) = eta /dist * exp(-1j*2*pi*dist/lambda_c) * exp(-1j*2*pi*Loc_ant(1,n)/lambda_g);
        end
    end
end
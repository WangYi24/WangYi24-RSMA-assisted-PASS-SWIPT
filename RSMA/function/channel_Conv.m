function [h] = channel_Conv(Loc_ant, Loc_U, eta, lambda_c)
% Loc_ant：3*N，PA位置
% Loc_U：3*M，用户位置
% lambda_c：自由空间波长

%h:[M,N,K]

    [~,N] = size(Loc_ant);
    [~,M] = size(Loc_U);
    h = zeros(M,N);
        %公式7
    for m=1:M
        for n=1:N
            dist = norm(Loc_U(:,m)-Loc_ant(:,n));
            h(m,n) = eta /dist * exp(-1j*2*pi*dist/lambda_c);
        end
    end
end
clear
cvx_clear
%% Emin
% clear
% cvx_clear
% 
% D_x = 20;     %波导长度
% D_y = 10;   %区域宽度
% d = 3;      %波导高 (k)
% sigma2_dBm = -90;  % 噪声功率 (dBm)
% sigma2 = 10.^((sigma2_dBm-30)./10);%噪声功率s
% 
% N = 1;
% K = 2;
% J = 1;
% 
% E_range = 8e-7:0.2e-7:10e-7;
% 
% P=40;
% ct = 20;%用户实现次数
% 
% R_SDMA_iE = zeros(1,length(E_range)); 
% E_SDMA_iE = zeros(J,length(E_range));
% C_SDMA_iE = zeros(K,length(E_range));
% 
% load('Loc_IDR_ct[60].mat')
% Loc_IDR_ct = Loc_IDR_ct(:,1:K*ct);
% Loc_IDR_cell = mat2cell(Loc_IDR_ct, 3, repmat(K, 1, ct));
% 
% for iE = 1:length(E_range)
%     E = E_range(iE)
% 
%     R_ct = zeros(K,ct); 
%     EHR_ct = zeros(J,ct);
%     c_ct = zeros(K,ct);
%     num_badpoint = 0;%坏点
%     %固定EHR的位置
%     % Loc_EHR = [zeros(1,J) - D_x/2 + ((0:J-1)+0.5)*D_x/J; ((0:J-1)+0.5)*D_y/J; zeros(1,J)];
%     Loc_EHR = [D_x/4; D_y/2; 0];
%     parfor it = 1:ct %Monte Carlo simulations   
%         fprintf('Worker processed it = %d\n', it);
%         %生成设备位置，列为该设备坐标
%         Loc_IDR = Loc_IDR_cell{it};
%         SumRate_hist=[];%记录迭代的目标函数
%         xn_hist = [];
%         SumRate_best = -inf;
%         [SumRate_best,Loc_PA_best,W_best]= SDMA_search(K,J,N,P,E,Loc_IDR,Loc_EHR,D_x,D_y,d,sigma2);
%         if SumRate_best==-inf
%             num_badpoint = num_badpoint + 1;
%         else
%             W_p = W_best;
%             h_IDR = channel(Loc_PA_best, Loc_IDR);
%             G = abs((h_IDR'*W_p)).^2;%得到 [K, K]
%             R_p = zeros(K,1);
%             for k = 1:K
%                 SINR = G(k,k)/(sum(G(k,:))-G(k,k)+sigma2);
%                 %R_c(k)= log2(1+abs(h_IDR(:,k)'*wc).^2/(sum(G(k,:))+sigma2));
%                 R_p(k) = log2(1+SINR);
%             end
%             R_ct(:,it) = R_p;
% 
%             h_EHR = channel(Loc_PA_best, Loc_EHR);
%             h_W = h_EHR' * W_p;              % [J, K]
%             term_p = sum(abs(h_W).^2, 2);  % 对每一行求和，得到 [J, 1]
%             EHR_ct(:,it)= term_p;
%         end        
%     end
% 
%     R_SDMA_iE(:,iE) = sum(R_ct,"all") / (ct-num_badpoint);
%     E_SDMA_iE(:,iE) = sum(EHR_ct,2) / (ct-num_badpoint);
%     C_SDMA_iE(:,iE) = sum(c_ct,2) / (ct-num_badpoint);
%     save("E.mat")
% end
% 
% figure
% 
% plot(E_range*1e6,R_SDMA_iE,'ko-'),hold on
% 
% legend('sum rate')
% xlabel('E_{min}(W)')
% ylabel('rate (bps/Hz)')
% ylim([0 40]);
% xlim([2e-1 10e-1])
% line([E_range(end-1)*1e6, E_range(end-1)*1e6], [0, R_SDMA_iE(end-1)],'Color', 'k')

%% P
clear
cvx_clear

D_x = 20;     %波导长度
D_y = 10;   %区域宽度
d = 3;      %波导高 (k)
sigma2_dBm = -90;  % 噪声功率 (dBm)
sigma2 = 10.^((sigma2_dBm-30)./10);%噪声功率s

N = 2;
K = 3;
J = 1;

E = 1e-7;

P_range = 40:2:42;
ct = 15;%用户实现次数

R_SDMA_iP = zeros(1,length(P_range)); 
E_SDMA_iP = zeros(J,length(P_range));
C_SDMA_iP = zeros(K,length(P_range));

% Loc_IDR_ct = [D_x*(rand(1,60)-0.5); D_y*rand(1,60); zeros(1,60)];
load('Loc_IDR_ct[60].mat')
Loc_IDR_ct = Loc_IDR_ct(:,1:K*ct);
Loc_IDR_cell = mat2cell(Loc_IDR_ct, 3, repmat(K, 1, ct));

for iP = 1:length(P_range)
    P = P_range(iP)

    R_ct = zeros(K,ct); 
    EHR_ct = zeros(J,ct);
    c_ct = zeros(K,ct);
    num_badpoint = 0;%坏点
    %固定EHR的位置
    % Loc_EHR = [zeros(1,J) - D_x/2 + ((0:J-1)+0.5)*D_x/J; ((0:J-1)+0.5)*D_y/J; zeros(1,J)];
    Loc_EHR = [D_x/4; D_y/2; 0];
    for it = 1:ct %Monte Carlo simulations   
        fprintf('Worker processed it = %d\n', it);
        %生成设备位置，列为该设备坐标
        Loc_IDR = Loc_IDR_cell{it};
        SumRate_hist=[];%记录迭代的目标函数
        xn_hist = [];
        SumRate_best = -inf;
        [SumRate_best,Loc_PA_best,W_best]= SDMA_search(K,J,N,P,E,Loc_IDR,Loc_EHR,D_x,D_y,d,sigma2);
        % [SumRate_best,Loc_PA_best,W_best]= SDMA_fix_v2(K,J,N,P,E,Loc_IDR,Loc_EHR,D_x,D_y,d,sigma2);
        if SumRate_best==-inf
            num_badpoint = num_badpoint + 1;
        else
            W_p = W_best;
            h_IDR = channel(Loc_PA_best, Loc_IDR);
            G = abs((h_IDR'*W_p)).^2;%得到 [K, K]
            R_p = zeros(K,1);
            for k = 1:K
                SINR = G(k,k)/(sum(G(k,:))-G(k,k)+sigma2);
                %R_c(k)= log2(1+abs(h_IDR(:,k)'*wc).^2/(sum(G(k,:))+sigma2));
                R_p(k) = log2(1+SINR);
            end
            R_ct(:,it) = R_p;
            
            h_EHR = channel(Loc_PA_best, Loc_EHR);
            h_W = h_EHR' * W_p;              % [J, K]
            term_p = sum(abs(h_W).^2, 2);  % 对每一行求和，得到 [J, 1]
            EHR_ct(:,it)= term_p;
        end        
    end

    R_SDMA_iP(:,iP) = sum(R_ct,"all")/K / (ct-num_badpoint);
    E_SDMA_iP(:,iP) = sum(EHR_ct,2) / (ct-num_badpoint);
    C_SDMA_iP(:,iP) = sum(c_ct,2) / (ct-num_badpoint);
    filename = sprintf('P[%d]231.mat', P);
    save(filename);
end

figure

plot(P_range,R_SDMA_iP,'ro-'),hold on

legend('WSR')
xlabel('P(dBm)')
ylabel('rate (bps/Hz)')
% ylim([30 46]);
% line([P_range(end-1)*1e6, P_range(end-1)*1e6], [0, R_SDMA_iP(end-1)],'Color', 'k')

% %% N
% clear
% cvx_clear
% 
% M = 4;
% N_range = 4:2:12;
% P = 0;
% 
% R_MW_iN = zeros(M,length(N_range)); 

% for iN = 1:length(N_range)
%     N = N_range(iN)
%     [R_MW_iN(:,iN)] = test_func_v4(M,N,P);
% end
% figure
% plot(N_range,min(R_MW_iN),'ro-'),hold on
% 
% legend('MW','MW-ZF')
% xlabel('N')
% ylabel('Minimum data rate (bps)')
% % 
% %% M
% clear
% cvx_clear
% 
% M_range = 2:2:8;
% N = 8;
% P = 10;
% 
% R_MW_iM = 100*ones(8,length(M_range)); 

% for iM = 1:length(M_range)
%     M = M_range(iM)
%     [R_MW_iM(1:M,iM)] = test_func_v5(M,N,P);
% end
% figure
% plot(M_range,min(R_MW_iM),'ro-'),hold on

% 
% legend('MW','MW-ZF')
% xlabel('M')
% ylabel('Minimum data rate (bps)')
    val_best = zeros(2,1);%记录每次iter的最优值
    val_prev = zeros(2,1);%用在迭代终止条件
    for iter = 1:max_iter
        %遍历每根波导上的PA位置
        for n = 1:N
            best_x = Loc_PA(1,n);
            candidates = linspace(-D_x/2, D_x/2, num_point);
            %ln遍历每一个位置
            for xn = candidates
                Loc_PA(1,n) = xn;
                %子载波k从PAn到用户m的信道：h[M,N,K]
                h = channel(Loc_PA, Loc_U, eta, lambda_c, lambda_g, f_c, K, B); 
                %根据OFDMA for Pinching Antenna Systems算法1stage1分配子载波,此时无预编码，子载波功率均分
                b = Subcarrier_assignment(h, P_max, sigma2, delta_f);                
                %完成了子载波分配，得到分配策略向量b。根据分配情况为不同子载波确定MRC预编码
                [W,val] = MRC_and_obj(h,b,delta_f,P_max,sigma2);%W[N.K]每一列是子载波k的预编码，val[M,1]M个用户的速率

                % %为了验证子载波分配的作用，随机生成均分子载波的b
                % b_test = generate_binary_matrix(m,k);
                % [W,val_test] = MRC_and_obj(h,b_test,delta_f,P_max,sigma2);

                if min(val) > min(val_best)
                    val_best = val;
                    best_x = xn;
                end
            end   
            Loc_PA(1,n) = best_x;  % 更新lm最优位置
        end
        % 检查目标函数是否收敛

        val_current = val_best;
        R_iter(:,iter) = val_best;
        if abs(min(val_current) - min(val_prev)) < epsilon
            break;
        end
        val_prev = val_current;        
    end
function b = generate_binary_matrix(m,k)
%m：矩阵行数
%n：矩阵列数
    b = zeros(m, k);

    % randperm返回行向量，其中包含从 1 到 n 没有重复元素的整数随机排列。
    idx = randperm(k);

    % 分配前k/2个给第一行，后k个给第二行
    b(1, idx(1:k/2)) = 1;
    b(2, idx(k/2+1:end)) = 1;
end

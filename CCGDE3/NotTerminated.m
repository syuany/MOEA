function b=NotTerminated(pop, nfe)
    problem = 'ZDT1';
    ref_point = [1.1, 1.1]; 
    hv_threshold_ratio = 0.95;
    
    % 使用 persistent 变量缓存计算量大的数据
    persistent true_hv hv_threshold
    if isempty(true_hv)
        n_true_points = 2000; 
        true_front = generate_true_front(problem, n_true_points);
        true_hv = calculate_hv(true_front, ref_point);
        hv_threshold = true_hv * hv_threshold_ratio;
    end

    b = true;
    
    if nfe >= 1e7 
        b = false;
        return;
    end
    
    % 降低 HV 计算频率：每 500 次评估才检查一次 HV (大幅提升速度)
    % 因为 calculate_hv 也是 O(N*logN)
    if mod(nfe, 500) ~= 0
        return;
    end

    % 计算当前 HV
    costs = [pop.Cost]';
    % 过滤非支配解（简单过滤，不需要完整排序）
    [~, idx] = sort(costs(:,1));
    current_hv = calculate_hv(costs, ref_point);
    
    % 进度条显示
    % fprintf('NFE: %d | Target HV: %.4f | Current HV: %.4f\n', nfe, hv_threshold, current_hv);

    if current_hv >= hv_threshold - 1e-8
        b = false;
    end
end

function hv = calculate_hv(points, ref_point)
    % 计算2目标最小化问题的超体积（参考点在右上）
    % 输入：
    %   points: [m×2] 目标函数值矩阵
    %   ref_point: [1×2] 参考点
    
    if isempty(points)
        hv = 0;
        return;
    end

    % 按 f1 降序排序 (从最右边的点开始向左扫描)
    [~, idx] = sort(points(:, 1), 'descend');
    sorted_points = points(idx, :);
    m = size(sorted_points, 1);

    % 初始化
    hv = 0;
    prev_f1 = ref_point(1);  % 上一个切片的右边界（初始为参考点x）
    
    % 扫描线算法
    for i = 1:m
        current_f1 = sorted_points(i, 1);
        current_f2 = sorted_points(i, 2);

        % 计算宽度：从上一个点的 f1 到当前点的 f1
        width = prev_f1 - current_f1;
        
        % 计算高度：从当前点的 f2 到参考点的 f2 (修正点在此)
        % 注意：这里假设输入是非支配解集。如果是从右向左扫，
        % 当前点的 f2 构成了该区间内的下界。
        height = ref_point(2) - current_f2;

        % 累加体积
        if width > 0 && height > 0
            hv = hv + width * height;
        end

        % 更新右边界
        prev_f1 = current_f1;
    end
end

function true_front = generate_true_front(problem, n_points)
    % 生成ZDT系列问题的真实帕累托前沿（离散样本点）
    % input:
    %   problem: 问题名称，如'ZDT1','ZDT2','ZDT3','ZDT6'
    %   n_points: 离散样本点数量（建议≥1000，平衡精度与效率）
    % output:
    %   true_front: [n_points × 2] 矩阵，每行是一个真实前沿点(f1,f2)

    f1 = [];
    switch problem
        case 'ZDT1'
            f1 = linspace(0, 1, n_points);  % f1∈[0,1]，均匀采样
            f2 = 1 - sqrt(f1);
        case 'ZDT2'
            f1 = linspace(0, 1, n_points);
            f2 = 1 - f1.^2;
        case 'ZDT3'
            f1 = linspace(0, 1, n_points);
            f2 = 1 - sqrt(f1) - f1.*sin(10*pi*f1);
        case 'ZDT6'
            f1 = linspace(0.2, 1, n_points);  % ZDT6的f1∈[0.2,1]
            f2 = 1 - (f1/0.5).^0.5;
        otherwise
            error('不支持该问题，请输入ZDT1-ZDT6');
    end
    true_front = [f1', f2'];  % 转换为矩阵（每行一个点）
end


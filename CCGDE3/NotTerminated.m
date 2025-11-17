function b=NotTerminated(pop, nfe)
    problem = 'ZDT1';          % 测试问题（可替换为ZDT2-ZDT6）
    n_true_points = 2000;      % 真实前沿样本点数量（精度足够）
    ref_point = [1.1, 1.1];    % 参考点（ZDT系列通用）
    hv_threshold_ratio = 0.95; % 超体积阈值比例（95%）
    max_evaluations = 1e7;     % 兜底条件：最大函数评估次数（论文设1e7）

    true_front = generate_true_front(problem, n_true_points);
    true_hv = calculate_hv(true_front, ref_point);
    hv_threshold = true_hv * hv_threshold_ratio;

    current_hv=calculate_hv(pop.Cost, ref_point);
    if current_hv >= hv_threshold - 1e-8
        b=false;
    elseif nfe>=1e7
        b=false;
    end
    b=true;
end

function hv = calculate_hv(points, ref_point)
% 计算2目标问题中，非支配解集合相对于参考点的超体积
% input:
%   points: [m × 2] 矩阵，每行是一个目标空间点(f1,f2)（需先筛选为非支配解）
%   ref_point: [1 × 2] 向量，参考点（如[1.1,1.1]）
% output:
%   hv: 超体积值（若points为空，返回0）

if isempty(points)
    hv = 0;
    return;
end

% 步骤1：确保points是非支配解（若输入未筛选，先执行非支配排序）
points = non_dominated_sort(points);

% 步骤2：按f1升序排序（Sweep Line算法前提）
[~, idx] = sort(points(:,1));
sorted_points = points(idx, :);
m = size(sorted_points, 1);

% 步骤3：计算超体积（基于参考点ref_point = [r1, r2]）
hv = 0;
prev_f1 = ref_point(1);  % 初始位置：参考点的f1
prev_f2 = ref_point(2);  % 初始位置：参考点的f2

for i = 1:m
    current_f1 = sorted_points(i, 1);
    current_f2 = sorted_points(i, 2);
    
    % 计算当前区间的体积（矩形面积：Δf1 × Δf2）
    delta_f1 = prev_f1 - current_f1;
    delta_f2 = prev_f2 - current_f2;
    if delta_f1 > 1e-10 && delta_f2 > 1e-10  % 避免数值误差导致的负体积
        hv = hv + delta_f1 * delta_f2;
    end
    
    % 更新prev_f2（Sweep Line核心：沿f1减小方向，f2取最小值）
    prev_f2 = min(prev_f2, current_f2);
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


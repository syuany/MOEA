% function f = ZDT1(x)
%     m = length(x);  
    
%     f1 = x(1);
    
%     sum_x = sum(x(2:m));  
%     g = 1 + 9 * (sum_x / (m - 1));  
    
%     f2 = g * (1 - sqrt(f1 / g));
    
%     f = [f1; f2];
% end

function f = ZDT1(x)
    
    [N, D] = size(x);
    
    f1 = x(:, 1);
    
    % 向量化计算 g (对第2列到最后一列求和)
    % sum(..., 2) 表示按行求和
    sum_x = sum(x(:, 2:end), 2);
    g = 1 + 9 * (sum_x / (D - 1));
    
    f2 = g .* (1 - sqrt(f1 ./ g));
    
    f = [f1, f2]'; 
end
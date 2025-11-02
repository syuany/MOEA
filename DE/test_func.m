% function y = test_func(xx, o, M)
%     D = length(xx);
%     f_bias=0;
%     z = (xx - o) * M; 
%     sum_term = 0;
%     for i = 1:D
%         sum_term = sum_term + (z(i)^2 - 10*cos(2*pi*z(i)) + 10);
%     end
%     y = sum_term + f_bias;
% end

function f = test_func(X)
    
    %% 提取变量
    x = X(1:4);    % x1~x4
    y = X(5:8);    % y1~y4（角度，单位：弧度）
    
    C1 = 18048;
    C2 = 13657.36315298172;
    C3 = 5497.052905295088;
    C4 = 14508.29635946082;
    C5 = 3157.294805334107;
    C6 = 105543794;
    C7 = 91598751.25016867;
    C8 = 33470578.99613227;
    
    c1 = sum(x) - C1;  
    c2 = sum(x .* cos(y)) - C2;  
    c3 = sum(x .* sin(y)) - C3;  
    c4 = sum(x .* cos(y).^2) - C4;  
    c5 = sum(x .* cos(y) .* sin(y)) - C5;  
    c6 = sum(x.^2) - C6;  
    c7 = sum(x.^2 .* cos(y)) - C7;  
    c8 = sum(x.^2 .* sin(y)) - C8;  
    
    f = c1^2 + c2^2 + c3^2 + c4^2 + c5^2 + c6^2 + c7^2 + c8^2;
end
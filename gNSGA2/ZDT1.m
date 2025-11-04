function f = ZDT1(x)
    m = length(x);  
    
    f1 = x(1);
    
    sum_x = sum(x(2:m));  
    g = 1 + 9 * (sum_x / (m - 1));  
    
    f2 = g * (1 - sqrt(f1 / g));
    
    f = [f1; f2];
end
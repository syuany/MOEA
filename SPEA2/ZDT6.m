function f = ZDT6(x)
    m = length(x);  

    f1 = 1 - exp(-4 * x(1)) * (sin(6 * pi * x(1)))^6;
    
    sum_x = sum(x(2:m));  
    g = 1 + (m - 1) * ((sum_x) / (m - 1))^0.25;  
    
    f2 = g * (1 - (f1 / g)^2);

    f = [f1; f2];
end
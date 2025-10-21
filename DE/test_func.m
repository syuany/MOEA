function y = test_func(xx, o, M)
    D = length(xx);
    f_bias=0;
    z = (xx - o) * M; 
    sum_term = 0;
    for i = 1:D
        sum_term = sum_term + (z(i)^2 - 10*cos(2*pi*z(i)) + 10);
    end
    y = sum_term + f_bias;
end
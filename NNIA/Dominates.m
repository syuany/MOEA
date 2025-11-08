function b = Dominates(x, y)
    if isstruct(x)
        x = x.Value;
    end
    if isstruct(y)
        y = y.Value;
    end
    b = all(x <= y) && any(x < y);
end
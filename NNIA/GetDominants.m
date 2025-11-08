function pop = GetDominants(pop)
    nPop = numel(pop);
    for i = 1:nPop
        for j = i + 1:nPop
            p = pop(i);
            q = pop(j);
            if Dominates(p, q)
                q.IsDominated = true;
            end
            if Dominates(q, p)
                p.IsDominated = true;
            end
            pop(i) = p;
            pop(j) = q;
        end
    end
    pop = pop(~[pop.IsDominated]);
end
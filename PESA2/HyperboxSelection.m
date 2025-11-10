function pool = HyperboxSelection(pop)
    nPop = numel(pop);
    pool = repmat(pop(1), nPop, 1);  

    for i = 1:nPop
        pool(i) = HyperboxSelect(pop);
    end
end

function selected = HyperboxSelect(pop)
    mx=max([pop.Fitness]);
    candidates=find([pop.Fitness]==mx);
    selected=pop(candidates(randi(numel(candidates))));
end
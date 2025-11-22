function varargout = NDSort(pop)
    Costs = [pop.Cost]'; 
    [N, M] = size(Costs);
    
    DominatedCount = zeros(N, 1);
    DominationSet = cell(N, 1);
    F = {[]};
    
    for i = 1:N
        for j = i+1:N
            if all(Costs(i,:) <= Costs(j,:)) && any(Costs(i,:) < Costs(j,:))
                DominationSet{i}(end+1) = j;
                DominatedCount(j) = DominatedCount(j) + 1;
            elseif all(Costs(j,:) <= Costs(i,:)) && any(Costs(j,:) < Costs(i,:))
                DominationSet{j}(end+1) = i;
                DominatedCount(i) = DominatedCount(i) + 1;
            end
        end
        if DominatedCount(i) == 0
            F{1} = [F{1}, i];
        end
    end
    
    k = 1;
    Rank = zeros(N, 1);
    Rank(F{1}) = 1;
    
    while true
        Q = [];
        for i = F{k}
            for j = DominationSet{i}
                DominatedCount(j) = DominatedCount(j) - 1;
                if DominatedCount(j) == 0
                    Q = [Q, j];
                    Rank(j) = k + 1;
                end
            end
        end
        if isempty(Q), break; end
        F{k+1} = Q;
        k = k + 1;
    end
    
    if nargout > 1
        for i = 1:N
            pop(i).DominationSet = DominationSet{i};
            pop(i).DominatedCount = 0; % 算法结束归零
            pop(i).Rank = Rank(i);
        end
        varargout{1} = pop;
        varargout{2} = F;
    else
        varargout{1} = F;
    end
end
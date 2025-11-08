function pool = TournamentSelection(pop, tournament_size)
    nPop = numel(pop);
    pool = repmat(pop(1), nPop, 1);  

    for i = 1:nPop
        % 锦标赛选择两个父代
        pool(i) = TournamentSelect(pop, tournament_size);
    end
end

function selected = TournamentSelect(pop, tournament_size)
    tournament_indices = randperm(numel(pop), tournament_size);
    tournament_members = pop(tournament_indices);
    
    best_index = 1;
    best_fitness = tournament_members(1).Fitness;
    
    for i = 2:tournament_size
        current_fitness = tournament_members(i).Fitness;
        
        if current_fitness > best_fitness
            best_index = i;
            best_fitness = current_fitness;
        end
    end
    
    selected = tournament_members(best_index);
end
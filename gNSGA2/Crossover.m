function popc = TournamentSelection(pop, crossover_params)
    nPop = numel(pop);
    popc = repmat(pop(1), nPop, 1);  

    
    for i = 1:nPop
        % 锦标赛选择两个父代
        p1 = TournamentSelect(pop, crossover_params.TournamentSize);
        p2 = TournamentSelect(pop, crossover_params.TournamentSize);
        
        % 执行模拟二进制交叉
        if rand < crossover_params.CrossoverRate
            [child1_pos, child2_pos] = SBXCrossover(p1.Position, p2.Position, crossover_params.Eta);
            popc(i).Position = child1_pos;
        else
            popc(i).Position = p1.Position;
        end
    end
end

function selected = TournamentSelect(pop, tournament_size)
    tournament_indices = randperm(numel(pop), tournament_size);
    tournament_members = pop(tournament_indices);
    
    best_index = 1;
    best_rank = tournament_members(1).Rank;
    best_crowding_distance = tournament_members(1).CrowdingDistance;
    
    for i = 2:tournament_size
        current_rank = tournament_members(i).Rank;
        current_crowding_distance = tournament_members(i).CrowdingDistance;
        
        % 如果当前个体支配等级更好
        if current_rank < best_rank
            best_index = i;
            best_rank = current_rank;
            best_crowding_distance = current_crowding_distance;
        % 如果支配等级相同，但拥挤距离更大
        elseif current_rank == best_rank && current_crowding_distance > best_crowding_distance
            best_index = i;
            best_crowding_distance = current_crowding_distance;
        end
    end
    
    selected = tournament_members(best_index);
end

function [child1, child2] = SBXCrossover(parent1, parent2, eta)
    child1 = parent1;
    child2 = parent2;
    var_min = 0;
    var_max = 1;
    
    if rand < 0.5
        % 执行SBX交叉
        for i = 1:numel(parent1)
            if rand < 0.5
                if abs(parent1(i) - parent2(i)) > 1e-6
                    y1 = min(parent1(i), parent2(i));
                    y2 = max(parent1(i), parent2(i));
                    
                    r = rand;
                    
                    if r <= 0.5
                        beta = (2*r)^(1/(eta+1));
                    else
                        beta = (1/(2*(1-r)))^(1/(eta+1));
                    end
                    
                    child1(i) = 0.5*((y1+y2) - beta*(y2-y1));
                    child2(i) = 0.5*((y1+y2) + beta*(y2-y1));

                    child1(i) = max(var_min, min(var_max, child1(i)));
                    child2(i) = max(var_min, min(var_max, child2(i)));
                else
                    child1(i) = parent1(i);
                    child2(i) = parent2(i);
                end
            else
                child1(i) = parent1(i);
                child2(i) = parent2(i);
            end
        end
    end
end
function popn = Variation(pool, crossover_params, mutate_params)
    nPop = numel(pool);
    popn = repmat(pool(1), nPop, 1);
    
    for i = 1:2:nPop-1
        p1 = pool(i);
        p2 = pool(i+1);
        
        if rand < crossover_params.CrossoverRate
            [child1_pos, child2_pos] = UniformCrossover(p1.Position, p2.Position);
            popn(i).Position = child1_pos;
            popn(i+1).Position = child2_pos;
        else
            popn(i).Position = p1.Position;
            popn(i+1).Position = p2.Position;
        end

        % Bit-flip mutation
        for j = 1:numel(p1.Position)
            if rand < mutate_params.MutationRate
                popn(i).Position(j) = BitFlipMutation(popn(i).Position(j));  
            end
            if rand < mutate_params.MutationRate
                popn(i+1).Position(j) = BitFlipMutation(popn(i+1).Position(j));  
            end
        end
    end
end

function [off1, off2] = UniformCrossover(parent1, parent2)
    off1 = parent1;
    off2 = parent2;
    
    for i = 1:numel(parent1)
        if rand < 0.5
            off1(i) = parent2(i);
            off2(i) = parent1(i);
        end
    end
end

function mutated_value = BitFlipMutation(x)
    if x == 0
        mutated_value = 1;
    else
        mutated_value = 0;
    end
end
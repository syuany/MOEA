function popn = Variation(pool, crossover_params, mutate_params)
    nPop = numel(pool);
    popn = repmat(pool(1), nPop, 1);
    
    for i = 1:2:nPop-1
        p1 = pool(i);
        p2 = pool(i+1);
        
        if rand < crossover_params.CrossoverRate
            [child1_pos, child2_pos] = SBXCrossover(p1.Position, p2.Position, crossover_params.Eta, crossover_params.lb, crossover_params.ub);
            popn(i).Position = child1_pos;
            popn(i+1).Position = child2_pos;
        else
            popn(i).Position = p1.Position;
            popn(i+1).Position = p2.Position;
        end

        for j = 1:numel(p1.Position)
            if rand < mutate_params.MutationRate
                popn(i).Position(j) = PolynomialMutation(popn(i).Position(j), mutate_params.lb, mutate_params.ub, mutate_params.Eta);  
                popn(i+1).Position(j) = PolynomialMutation(popn(i+1).Position(j), mutate_params.lb, mutate_params.ub, mutate_params.Eta);  
            end
        end
    end
end

function [off1, off2] = SBXCrossover(parent1, parent2, eta, var_min, var_max)
    off1 = parent1;
    off2 = parent2;
    
    for i = 1:numel(parent1)
        if rand < 0.5
            beta = (2 * rand)^(1 / (eta + 1));
        else
            beta = (2 - 2 * rand)^(-1 / (eta + 1));
        end
        off1(i) = 0.5 * ((1 + beta) * parent1(i) + (1 - beta) * parent2(i));
        off2(i) = 0.5 * ((1 - beta) * parent1(i) + (1 + beta) * parent2(i));

        off1(i) = max(var_min, min(var_max, off1(i)));
        off2(i) = max(var_min, min(var_max, off2(i)));
    end
end

function mutated_value = PolynomialMutation(x, xl, xu, eta)
    % Polynomial Mutation
    delta1 = (x - xl) / (xu - xl);
    delta2 = (xu - x) / (xu - xl);
    
    rand_val = rand;
    mut_pow = 1.0 / (eta + 1.0);
    
    if rand_val <= 0.5
        xy = 1.0 - delta1;
        val = 2.0 * rand_val + (1.0 - 2.0 * rand_val) * (xy^(eta + 1.0));
        deltaq = val^mut_pow - 1.0;
    else
        xy = 1.0 - delta2;
        val = 2.0 * (1.0 - rand_val) + 2.0 * (rand_val - 0.5) * (xy^(eta + 1.0));
        deltaq = 1.0 - val^mut_pow;
    end
    
    mutated_value = x + deltaq * (xu - xl);
    
    mutated_value = max(xl, min(xu, mutated_value));
end
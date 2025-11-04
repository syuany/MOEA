function popm = Mutation(pop, mutate_params)
    % 多项式变异操作
    
    nPop = numel(pop);
    popm = repmat(pop(1), nPop, 1);  % 初始化变异后代
    CostFunction=mutate_params.CostFunction;
        
    for i = 1:nPop
        % 复制个体
        popm(i) = pop(i);
        
        % 对每个维度执行多项式变异
        for j = 1:numel(pop(i).Position)
            if rand < mutate_params.MutationRate
                popm(i).Position(j) = PolynomialMutation(pop(i).Position(j), mutate_params.lb, mutate_params.ub, mutate_params.Eta);  
            end
        end
        
        % popm(i).Cost = CostFunction(popm(i).Position);
    end
end

function mutated_value = PolynomialMutation(x, xl, xu, eta)
    % 多项式变异 (Polynomial Mutation)
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
    
    % 确保变异后的值在边界内
    mutated_value = max(xl, min(xu, mutated_value));
end
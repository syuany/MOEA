clc;
clear;
close all;

%% Problem Definition

CostFunction = @(x) ZDT1(x);      % Cost Function

nVar = 30;             % Number of Decision Variables

VarSize = [1 nVar];   % Size of Decision Variables Matrix

VarMin = 0;          % Lower Bound of Variables
VarMax = 1;          % Upper Bound of Variables

% Number of Objective Functions
nObj = numel(CostFunction(unifrnd(VarMin, VarMax, VarSize)));


%% NSGA-II Parameters

MaxIt = 500;      % Maximum Number of Iterations

nPop = 50;        % Population Size

pCrossover = 0.7;                         % Crossover Percentage
nCrossover = 2*round(pCrossover*nPop/2);  % Number of Parnets (Offsprings)

pMutation = 0.4;                          % Mutation Percentage
nMutation = round(pMutation*nPop);        % Number of Mutants

mu = 0.02;                    % Mutation Rate

sigma = 0.3*(VarMax-VarMin);  % Mutation Step Size

% g-Dominance Parameters

% g=[0.4;0.2];
g=[0.6;0.5];

%% Initialization

empty_individual.Position = [];
empty_individual.Cost = [];
empty_individual.Rank = [];
empty_individual.DominationSet = [];
empty_individual.DominatedCount = [];
empty_individual.CrowdingDistance = [];

pop = repmat(empty_individual, nPop, 1);

for i = 1:nPop
    
    pop(i).Position = unifrnd(VarMin, VarMax, VarSize);
    
    pop(i).Cost = CostFunction(pop(i).Position);
    
end

% Non-Dominated Sorting
[pop, F] = NonDominatedSorting(pop, g);

% Calculate Crowding Distance
pop = CalcCrowdingDistance(pop, F);

% Sort Population
[pop, F] = SortPopulation(pop);


%% NSGA-II Main Loop

for it = 1:MaxIt
    
    % Crossover
    popc = repmat(empty_individual, nCrossover/2, 2);
    for k = 1:nCrossover/2
        
        i1 = randi([1 nPop]);
        p1 = pop(i1);
        
        i2 = randi([1 nPop]);
        p2 = pop(i2);
        
        [popc(k, 1).Position, popc(k, 2).Position] = Crossover(p1.Position, p2.Position);
        
        popc(k, 1).Cost = CostFunction(popc(k, 1).Position);
        popc(k, 2).Cost = CostFunction(popc(k, 2).Position);
        
    end
    popc = popc(:);
    
    % Mutation
    popm = repmat(empty_individual, nMutation, 1);
    for k = 1:nMutation
        
        i = randi([1 nPop]);
        p = pop(i);
        
        popm(k).Position = Mutate(p.Position, mu, sigma);
        
        popm(k).Cost = CostFunction(popm(k).Position);
        
    end
    
    % Merge
    pop = [pop
         popc
         popm]; %#ok
     
    % Non-Dominated Sorting
    [pop, F] = NonDominatedSorting(pop, g);

    % Calculate Crowding Distance
    pop = CalcCrowdingDistance(pop, F);

    % Sort Population
    pop = SortPopulation(pop);
    
    % Truncate
    pop = pop(1:nPop);
    
    % Non-Dominated Sorting
    [pop, F] = NonDominatedSorting(pop, g);

    % Calculate Crowding Distance
    pop = CalcCrowdingDistance(pop, F);

    % Sort Population
    [pop, F] = SortPopulation(pop);
    
    % Store F1
    F1 = pop(F{1});
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Number of F1 Members = ' num2str(numel(F1))]);
    
    % Plot F1 Costs
    figure(1);
    PlotCosts(F1, g);
    pause(0.01);
    
end


function b = Dominates(x, y, g)
    if nargin < 3
        g = [];  % 无g参数时使用标准Pareto支配
    end
    
    if isstruct(x)
        x = x.Cost;
    end
    
    if isstruct(y)
        y = y.Cost;
    end
    
    if isempty(g)
        % 标准Pareto支配
        b = all(x <= y) && any(x < y);
    else
        % 检查x和y是否在g点内部
        xInside = all(x <= g);
        yInside = all(y <= g);
        
        if xInside && yInside
            % 都在内部，使用标准Pareto支配
            b = all(x <= y) && any(x < y);
        elseif ~xInside && ~yInside
            % 都在外部，使用g-dominance
            dx = max(x - g);
            dy = max(y - g);
            b = (dx < dy) || (dx == dy && all(x <= y) && any(x < y));
        else
            % 一个在内部一个在外部，内部的支配外部的
            b = xInside && ~yInside;
        end
    end
end

function [pop, F] = NonDominatedSorting(pop, g)

    nPop = numel(pop);

    for i = 1:nPop
        pop(i).DominationSet = [];
        pop(i).DominatedCount = 0;
    end
    
    F{1} = [];
    
    for i = 1:nPop
        for j = i+1:nPop
            p = pop(i);
            q = pop(j);
            
            if Dominates(p, q, g)
                p.DominationSet = [p.DominationSet j];
                q.DominatedCount = q.DominatedCount+1;
            end
            
            if Dominates(q.Cost, p.Cost, g)
                q.DominationSet = [q.DominationSet i];
                p.DominatedCount = p.DominatedCount+1;
            end
            
            pop(i) = p;
            pop(j) = q;
        end
        
        if pop(i).DominatedCount == 0
            F{1} = [F{1} i];
            pop(i).Rank = 1;
        end
    end
    
    k = 1;
    
    while true
        
        Q = [];
        
        for i = F{k}
            p = pop(i);
            
            for j = p.DominationSet
                q = pop(j);
                
                q.DominatedCount = q.DominatedCount-1;
                
                if q.DominatedCount == 0
                    Q = [Q j]; %#ok
                    q.Rank = k+1;
                end
                
                pop(j) = q;
            end
        end
        
        if isempty(Q)
            break;
        end
        
        F{k+1} = Q; %#ok
        
        k = k+1;
        
    end
    

end

function pop = CalcCrowdingDistance(pop, F)

    nF = numel(F);
    
    for k = 1:nF
        
        Costs = [pop(F{k}).Cost];
        
        nObj = size(Costs, 1);
        
        n = numel(F{k});
        
        d = zeros(n, nObj);
        
        for j = 1:nObj
            
            [cj, so] = sort(Costs(j, :));
            
            d(so(1), j) = inf;
            
            for i = 2:n-1
                
                d(so(i), j) = abs(cj(i+1)-cj(i-1))/abs(cj(1)-cj(end));
                
            end
            
            d(so(end), j) = inf;
            
        end
        
        
        for i = 1:n
            
            pop(F{k}(i)).CrowdingDistance = sum(d(i, :));
            
        end
        
    end


end

function [pop, F] = SortPopulation(pop)

    % Sort Based on Crowding Distance
    [~, CDSO] = sort([pop.CrowdingDistance], 'descend');
    pop = pop(CDSO);
    
    % Sort Based on Rank
    [~, RSO] = sort([pop.Rank]);
    pop = pop(RSO);
    
    % Update Fronts
    Ranks = [pop.Rank];
    MaxRank = max(Ranks);
    F = cell(MaxRank, 1);
    for r = 1:MaxRank
        F{r} = find(Ranks == r);
    end

end

function [y1, y2] = Crossover(x1, x2)
    VarMin = 0;  
    VarMax = 1;  
    
    alpha = rand(size(x1));
    
    y1 = alpha.*x1 + (1-alpha).*x2;
    y2 = alpha.*x2 + (1-alpha).*x1;
    
    % 边界处理
    y1 = max(min(y1, VarMax), VarMin);
    y2 = max(min(y2, VarMax), VarMin);
end

function y = Mutate(x, mu, sigma)
    nVar = numel(x);
    VarMin = 0;  
    VarMax = 1;  
    
    nMu = ceil(mu*nVar);

    j = randsample(nVar, nMu);
    if numel(sigma)>1
        sigma = sigma(j);
    end
    
    y = x;
    
    y(j) = x(j) + sigma.*randn(size(j));
    
    % 边界处理
    y(j) = max(min(y(j), VarMax), VarMin);
end


function PlotCosts(pop, g)
    figure(1);  
    set(gcf, 'Position', [100 100 800 600]);
    
    f1_exact = linspace(0, 1, 200);
    f2_exact = 1 - sqrt(f1_exact);
    plot(f1_exact, f2_exact, 'k-', 'LineWidth', 2, 'DisplayName', 'Exact');
    
    hold on;
    
    Costs = [pop.Cost];
    
    validIdx = all(Costs < 1000, 1);
    if any(validIdx)
        plot(Costs(1, validIdx), Costs(2, validIdx), 'rx', ...
                'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', 'NSGA-II');
    end
    
    if length(g) >= 2
        plot(g(1), g(2), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'DisplayName', 'g');
    end
    
    xlim([0 1.2]);
    ylim([0 1.5]);
    
    xlabel('1^{st} Objective', 'FontSize', 12);
    ylabel('2^{nd} Objective', 'FontSize', 12);
    title('Pareto Front Comparison', 'FontSize', 14);
    
    legend('Location', 'northwest', 'FontSize', 10);
    
    grid on;
    box on;
    hold off;
end
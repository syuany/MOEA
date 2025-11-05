clc;
clear;
close all;

%% Problem Definition

CostFunction = @(x) ZDT1(x); % Cost Function

nVar = 30; % Number of Decision Variables
VarSize = [1, nVar]; % Size of Decision Variables Matrix
VarMin = 0; % Lower Bound of Variables
VarMax = 1; % Upper Bound of Variables

% Number of Objective Functions
nObj = numel(CostFunction(unifrnd(VarMin, VarMax, VarSize)));

%% Parameters

MaxFES = 10000;
% MaxIt = 100;      % Maximum Number of Iterations
nPop = 50; % Population Size

crossover_params.TournamentSize = 2;
crossover_params.CrossoverRate = 0.8;
crossover_params.Eta = 20;
crossover_params.CostFunction = CostFunction;

mutate_params.MutationRate = 1 / nVar;
mutate_params.Eta = 20;
mutate_params.CostFunction = CostFunction;
mutate_params.lb = VarMin;
mutate_params.ub = VarMax;

% g-Dominance Parameters

% g=[0.4;0.2];
g = [0.6; 0.5];

M = 1e10;

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

% g-Dominance
pop = gDominance(pop, g, M);

[pop, F] = NonDominatedSorting(pop);
pop = CalcCrowdingDistance(pop, F);
[pop, F] = SortPopulation(pop);

FES = 0;
it = 1;

%% NSGA-II Main Loop

while FES < MaxFES

    % Crossover
    popc = TournamentSelection(pop, crossover_params);

    % Mutation
    popm = Mutation(pop, mutate_params);

    % Evaluate offspring
    for i = 1:numel(popc)
        popc(i).Cost = CostFunction(popc(i).Position);
        FES = FES + 1;
    end

    for i = 1:numel(popm)
        popm(i).Cost = CostFunction(popm(i).Position);
        FES = FES + 1;
    end

    % Merge
    pop = [pop
           popc
           popm]; %#ok

    % g-Dominance
    pop = gDominance(pop, g, M);

    [pop, F] = NonDominatedSorting(pop);
    pop = CalcCrowdingDistance(pop, F);
    pop = SortPopulation(pop);

    % Truncate
    pop = pop(1:nPop);

    [pop, F] = NonDominatedSorting(pop);
    pop = CalcCrowdingDistance(pop, F);
    [pop, F] = SortPopulation(pop);

    F1 = pop(F{1});

    % Show Iteration Information
    disp(['Iteration ', num2str(it), ': Number of F1 Members = ', num2str(numel(F1))]);
    it = it + 1;

    % Plot F1 Costs
    PlotCosts(F1, g);
    pause(0.01);

end

%% Function Defination

function pop = gDominance(pop, g, M)
nPop = numel(pop);
for i = 1:nPop
    val = pop(i).Cost;
    isFit = all(val <= g) || all(val >= g);

    % penalize
    if ~isFit
        pop(i).Cost = pop(i).Cost + M;
    end
end
end

function b = Dominates(x, y)

if isstruct(x)
    x = x.Cost;
end
if isstruct(y)
    y = y.Cost;
end

b = all(x <= y) && any(x < y);

end

function [pop, F] = NonDominatedSorting(pop)

nPop = numel(pop);
for i = 1:nPop
    pop(i).DominationSet = [];
    pop(i).DominatedCount = 0;
end

F{1} = [];
for i = 1:nPop
    for j = i + 1:nPop
        p = pop(i);
        q = pop(j);

        if Dominates(p, q)
            p.DominationSet = [p.DominationSet, j];
            q.DominatedCount = q.DominatedCount + 1;
        end

        if Dominates(q.Cost, p.Cost)
            q.DominationSet = [q.DominationSet, i];
            p.DominatedCount = p.DominatedCount + 1;
        end

        pop(i) = p;
        pop(j) = q;
    end
    if pop(i).DominatedCount == 0
        F{1} = [F{1}, i];
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
            q.DominatedCount = q.DominatedCount - 1;

            if q.DominatedCount == 0
                Q = [Q, j]; %#ok
                q.Rank = k + 1;
            end

            pop(j) = q;
        end
    end

    if isempty(Q)
        break;
    end

    F{k+1} = Q; %#ok
    k = k + 1;

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
        for i = 2:n - 1
            d(so(i), j) = abs(cj(i+1)-cj(i-1)) / abs(cj(1)-cj(end));
        end

        d(so(1), j) = inf;
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

clc;
clear;
close all;

%% Settings
CostFunction = @(x) ZDT6(x);
nVar = 10;
VarSize = [1, nVar];
VarMin = 0;
VarMax = 1;
nObj = numel(CostFunction(unifrnd(VarMin, VarMax, VarSize)));

MaxFES = 15000;
nPop = 50;

K = 0.05;

crossover_params.CrossoverRate = 1;
crossover_params.Eta = 15;
crossover_params.lb = VarMin;
crossover_params.ub = VarMax;

mutate_params.MutationRate = 1 / nVar;
mutate_params.Eta = 20;
mutate_params.lb = VarMin;
mutate_params.ub = VarMax;

%% Initialization
empty_individual.Position = [];
empty_individual.Cost = [];
empty_individual.Fitness = [];
empty_individual.Distance = [];
empty_individual.IsDominated = [];

pop = repmat(empty_individual, nPop, 1);
for i = 1:nPop
    pop(i).Position = unifrnd(VarMin, VarMax, VarSize);
    pop(i).Cost = CostFunction(pop(i).Position);
    pop(i).Distance = zeros(1, nPop);
    pop(i).IsDominated = false;
end

[pop, c] = CalculateFitness(pop, K);
pop = EnvironmentalSelection(pop, K, c, nPop);
[pop, c] = CalculateFitness(pop, K);

FES = 0; % Function Evaluations

%% Main Loop
while FES < MaxFES
    % Binary TournamentSelection
    pool = TournamentSelection(pop, 2);

    % Recombination and Mutation
    popn = Variation(pool, crossover_params, mutate_params);

    for i = 1:numel(popn)
        popn(i).Cost = CostFunction(popn(i).Position);
        FES = FES + 1;
    end

    % Merge
    pop = [pop; popn];

    [pop, c] = CalculateFitness(pop, K);
    pop = EnvironmentalSelection(pop, K, c, nPop);
    [pop, c] = CalculateFitness(pop, K);

    pop = GetDominants(pop);

    % display
    PlotCosts(pop);
    pause(0.01);
end

%% Reults
disp(' ');
Costs = [pop.Cost];
for j = 1:nObj

    disp(['Objective #', num2str(j), ':']);
    disp(['      Min = ', num2str(min(Costs(j, :)))]);
    disp(['      Max = ', num2str(max(Costs(j, :)))]);
    disp(['    Range = ', num2str(max(Costs(j, :))-min(Costs(j, :)))]);
    disp(['    St.D. = ', num2str(std(Costs(j, :)))]);
    disp(['     Mean = ', num2str(mean(Costs(j, :)))]);
    disp(' ');

end
%% Function Defination
function [pop, c] = CalculateFitness(pop, K)
nPop = numel(pop);
Costs = [pop.Cost];
nObj = size(Costs, 1);

% Normalization
lb = zeros(1, nObj);
ub = zeros(1, nObj);
for i = 1:nObj
    lb(i) = min(Costs(i, :));
    ub(i) = max(Costs(i, :));
end
normalized_Costs = Costs;
for i = 1:nObj
    normalized_Costs(i, :) = (normalized_Costs(i, :) - lb(i)) / (ub(i) - lb(i));
end

% Calculate Indicator Costs
for i = 1:nPop
    for j = 1:nPop
        if i ~= j
            pop(i).Distance(j) = max(normalized_Costs(:, i)-normalized_Costs(:, j));
        else
            pop(i).Distance(j) = 0;
        end
    end
end
c = max([pop.Distance]);

% Calculate Fitness
for i = 1:nPop
    pop(i).Fitness = 0;
    for j = 1:nPop
        if i ~= j
            pop(i).Fitness = pop(i).Fitness - exp(-pop(j).Distance(i)/(c * K));
        end
    end
end
end

function pop = EnvironmentalSelection(pop, c, K, N)
nPop = numel(pop);
while nPop > N
    Fits = [pop.Fitness];
    [~, worst] = min(Fits);

    for i = 1:nPop
        if i ~= worst
            pop(i).Fitness = pop(i).Fitness + exp(-pop(worst).Distance(i)/(c * K));
        end
    end

    pop(worst) = [];
    nPop = nPop - 1;
end
end

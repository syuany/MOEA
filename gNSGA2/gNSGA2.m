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
crossover_params.lb = VarMin;
crossover_params.ub = VarMax;

mutate_params.MutationRate = 1 / nVar;
mutate_params.Eta = 20;
mutate_params.lb = VarMin;
mutate_params.ub = VarMax;

% g-Dominance Parameters

% g=[0.4;0.2];
g = [0.6; 0.5];
% g = [0.6; 0.22];

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
pop = Evaluate(pop, g, M);

[pop, F] = NDSort(pop);
pop = CalcCrowdingDistance(pop, F);
[pop, F] = SortPopulation(pop);

FES = 0;
it = 1;
%% NSGA-II Main Loop

while FES < MaxFES

    % Crossover
    popc = Crossover(pop, crossover_params);

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
    pop = [pop; ...
        popc; ...
        popm]; %#ok

    % g-Dominance
    pop = Evaluate(pop, g, M);

    [pop, F] = NDSort(pop);
    pop = CalcCrowdingDistance(pop, F);
    pop = SortPopulation(pop);

    % Truncate
    pop = pop(1:nPop);

    [pop, F] = NDSort(pop);
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

function pop = Evaluate(pop, g, M)
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

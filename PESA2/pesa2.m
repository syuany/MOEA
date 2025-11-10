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

MaxFES = 50000;

% PESA2 Parameters
params.f1CostMin=0;
params.f1CostMax=1;
params.f2CostMin=0;
params.f2CostMax=nVar;
params.HyperGridSize=[32 32];
nArchive=100;

nPop = 10;
nBits=30;
totalBits=nBits*nVar;

crossover_params.CrossoverRate=0.7;

mutate_params.MutationRate=1/totalBits;


%% Initialization

empty_individual.Position = [];
empty_individual.Cost = [];
empty_individual.Fitness = [];
empty_individual.IsDominated = [];
empty_individual.BoxID = [];

pop = repmat(empty_individual, nPop, 1);
for i = 1:nPop
    pop(i).Position = randi([0,1], 1, totalBits);
    x=Decode(pop(i).Position, VarMin, VarMax, nBits, nVar);
    pop(i).Cost = CostFunction(x);
    pop(i).IsDominated = false;
end

EP=[];

pop=Evaluate(pop, params);
EPT=GetDominants(pop);
EP=[EP;EPT];

FES = 0; % Function Evaluations

%% Main Loop

while FES < MaxFES
    % Binary TournamentSelection based on Hyperbox Fitness
    pool = HyperboxSelection(EP);

    % Uniform Crossover and Bit-Flip Mutation
    pop = Variation(pool, crossover_params, mutate_params);

    for i = 1:numel(pop)
        x=Decode(pop(i).Position, VarMin, VarMax, nBits, nVar);
        pop(i).Cost = CostFunction(x);
        FES = FES + 1;
    end

    pop=Evaluate(pop, params);
    EPT = GetDominants(pop);
    EP=[EP;EPT];

    EP=GetDominants(EP);
    if(numel(EP) > nArchive)
        extra = numel(EP) - nArchive;
        box_count = zeros(numel(EP), 1);
        
        for i = 1:numel(EP)
            box_id = EP(i).BoxID;
            count = 0;
            for j = 1:numel(EP)
                if EP(j).BoxID == box_id
                    count = count + 1;
                end
            end
            box_count(i) = count;
        end
        
        [~, sorted_indices] = sort(box_count, 'descend');
        to_delete = sorted_indices(1:extra);
        EP(to_delete) = [];
    end

    % display
    PlotCosts(EP);
    pause(0.01);

end


%% Function Defination

function pop = Evaluate(pop, params)
    nPop=numel(pop);
    N1=params.HyperGridSize(1);
    N2=params.HyperGridSize(2);
    boxCount=zeros(1, N1*N2);
    for k=1:nPop
        % f2: row
        i = min(floor((pop(k).Cost(2) - params.f2CostMin) * N2 / ...
                (params.f2CostMax - params.f2CostMin)) + 1, N2);
        % f1: column
        j = min(floor((pop(k).Cost(1) - params.f1CostMin) * N1 / ...
                (params.f1CostMax - params.f1CostMin)) + 1, N1);
        pop(k).BoxID=sub2ind(params.HyperGridSize, i, j);
        boxCount(pop(k).BoxID)=boxCount(pop(k).BoxID)+1;
    end

    for i=1:nPop
        pop(i).Fitness=1.0/boxCount(pop(i).BoxID);
    end

end

function x = Decode(binary_position, VarMin, VarMax, nBits, nVar)
    % Decode binary string to real values

    x = zeros(1, nVar);
    
    for i = 1:nVar
        % Extract bits for variable i
        start_idx = (i-1) * nBits + 1;
        end_idx = i * nBits;
        var_bits = binary_position(start_idx:end_idx);
        
        % Convert binary bits to decimal
        decimal_value = 0;
        for j = 1:nBits
            decimal_value = decimal_value + var_bits(j) * 2^(nBits-j);
        end
        
        % Map to real value range
        max_decimal = 2^nBits - 1;
        x(i) = VarMin + (decimal_value / max_decimal) * (VarMax - VarMin);
    end
end
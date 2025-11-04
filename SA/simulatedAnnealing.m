clc;
clear;
close all;

% Data
all_data = importdata('data200.txt');
N = all_data(1);
BC = all_data(2);
BE = all_data(3);

PF = all_data(4:4+N-1);
CO = all_data(4+N:4+2*N-1);
EC = all_data(4+2*N:4+3*N-1);

% Settings
MaxIt = 100;

nVar = N;
VarSize = [1, nVar];

T = 100;
t = T;
MinT = 1e-5;
alpha = 0.95;

% Initialization
empty_individual.Position = [];
empty_individual.Perfomance = [];
empty_individual.Cost = [];
empty_individual.EnergyConsumption = [];

solution = empty_individual;
while true
    solution.Position = randi([0, 1], VarSize) == 1;
    solution = CalcFit(solution, PF, CO, EC);
    if solution.Cost <= BC && solution.EnergyConsumption <= BE
        break;
    end
end
best = solution;

% Main Loop
while t > MinT
    for i = 1:MaxIt
        newSolution = GenerateNeighbor(solution, BC, BE, PF, CO, EC);
        newSolution = CalcFit(newSolution, PF, CO, EC);
        delta = newSolution.Perfomance - solution.Perfomance;

        if delta > 0 || rand < exp(delta/t)
            solution = newSolution;
            if solution.Perfomance > best.Perfomance
                best = solution;
            end
        end
    end
    t = t * alpha;
end

fprintf('Best Performance: \t%.2f\n', best.Perfomance);
fprintf('Cost: \t\t\t\t%.2f \t(BC: %.2f)\n', best.Cost, BC);
fprintf('Energy Consumption: %.2f \t(BE: %.2f)\n', best.EnergyConsumption, BE);
fprintf('Number of selected items: %d\n', sum(best.Position));

function newSol = GenerateNeighbor(solution, BC, BE, PF, CO, EC)
while true
    newSol = solution;
    flip = randi(length(solution.Position));
    newSol.Position(flip) = ~newSol.Position(flip);

    if newSol.Position(flip)
        newSol.Perfomance = newSol.Perfomance + PF(flip);
        newSol.Cost = newSol.Cost + CO(flip);
        newSol.EnergyConsumption = newSol.EnergyConsumption + EC(flip);
    else
        newSol.Perfomance = newSol.Perfomance - PF(flip);
        newSol.Cost = newSol.Cost - CO(flip);
        newSol.EnergyConsumption = newSol.EnergyConsumption - EC(flip);
    end

    if newSol.Cost <= BC && newSol.EnergyConsumption <= BE
        break;
    end
end
end

function solution = CalcFit(solution, PF, CO, EC)
solution.Perfomance = sum(PF(solution.Position));
solution.Cost = sum(CO(solution.Position));
solution.EnergyConsumption = sum(EC(solution.Position));
end

clc;
clear;
close all;

% Settings
TestFunction = @(x) ZDT6(x);
nVar = 10;
VarSize = [1, nVar];
VarMin = zeros(VarSize);
VarMax = ones(VarSize);
nObj = numel(TestFunction(unifrnd(VarMin, VarMax, VarSize)));

MaxIt = 200;
nDom = 100;
nActive = 20;
nClone = 100;

crossover_params.pc = 1;
crossover_params.eta = 15;
crossover_params.lb = VarMin;
crossover_params.ub = VarMax;
crossover_params.TestFunction = TestFunction;

mutate_params.pm = 1 / nVar;
mutate_params.sm = nVar;
mutate_params.eta = 20;
mutate_params.lb = VarMin;
mutate_params.ub = VarMax;
mutate_params.TestFunction = TestFunction;

% Initalization
empty_individual.Position = [];
empty_individual.Objective = [];
empty_individual.IsDominated = [];
empty_individual.CrowdingDistance = [];

B = repmat(empty_individual, nDom, 1);
for i = 1:nDom
    B(i).Position = unifrnd(VarMin, VarMax, VarSize);
    B(i).Objective = TestFunction(B(i).Position);
    B(i).IsDominated = false;
end

D = GetDominants(B);
D = CalcCrowdingDistance(D);
D = SortPopulation(D);

% Main Loop
for it = 1:MaxIt
    % Nondominated Neighbor-Based Selection
    if numel(D) > nActive
        D = CalcCrowdingDistance(D);
        D = SortPopulation(D);
        A = D(1:nActive);
    else
        A = D;
    end

    % Proporthonal Cloning
    C = ProporthonalCloning(A, nClone);

    % Recombination and Hypermutation
    R = Recombination(C, A, crossover_params);
    CM = StaticHypermutation(R, mutate_params);

    % Merge
    B = [CM; D];

    % update dominant population
    DT = GetDominants(B);
    if numel(DT) > nDom
        DT = CalcCrowdingDistance(DT);
        DT = SortPopulation(DT);
        D = DT(1:nDom);
    else
        D = DT;
    end

    % display
    figure(1);
    PlotObjs(D);
    pause(0.01);
    disp(['Iteration ', num2str(it), ': Number of D = ', num2str(numel(D))]);
end

% Reults
disp(' ');

DO = [D.Objective];
for j = 1:nObj

    disp(['Objective #', num2str(j), ':']);
    disp(['      Min = ', num2str(min(DO(j, :)))]);
    disp(['      Max = ', num2str(max(DO(j, :)))]);
    disp(['    Range = ', num2str(max(DO(j, :))-min(DO(j, :)))]);
    disp(['    St.D. = ', num2str(std(DO(j, :)))]);
    disp(['     Mean = ', num2str(mean(DO(j, :)))]);
    disp(' ');

end


function b = Dominates(x, y)
if isstruct(x)
    x = x.Objective;
end
if isstruct(y)
    y = y.Objective;
end
b = all(x <= y) && any(x < y);
end

function pop = GetDominants(pop)
nPop = numel(pop);
for i = 1:nPop
    for j = i + 1:nPop
        p = pop(i);
        q = pop(j);
        if Dominates(p, q)
            q.IsDominated = true;
        end
        if Dominates(q, p)
            p.IsDominated = true;
        end
        pop(i) = p;
        pop(j) = q;
    end
end
pop = pop(~[pop.IsDominated]);
end

function pop = CalcCrowdingDistance(pop)
nPop = numel(pop);
Objs = [pop.Objective];
nObj = size(Objs, 1);
d = zeros(nPop, nObj);
for j = 1:nObj
    [oj, so] = sort(Objs(j, :));
    for i = 2:nPop - 1
        d(so(i), j) = abs(oj(i-1)-oj(i+1)) / abs(oj(1)-oj(end));
    end
    d(so(1), j) = inf;
    d(so(end), j) = inf;
end

for i = 1:nPop
    pop(i).CrowdingDistance = sum(d(i, :));
end
end

function pop = SortPopulation(pop)
CD = [pop.CrowdingDistance];

a = isinf(CD);
[pop(a).CrowdingDistance] = deal(max(CD(~a)).*2);

[~, so] = sort(CD, 'descend');
pop = pop(so);
end

function C = ProporthonalCloning(pop, nClone)
nPop = numel(pop);
sd = sum([pop.CrowdingDistance]);
C = [];
for i = 1:nPop
    q = ceil(nClone.*pop(i).CrowdingDistance/sd);
    CT = repmat(pop(i), q, 1);
    C = [C; CT];
end
end

function R = Recombination(C, A, crossover_params)
% SBX crossover
nClone = numel(C);
nActive = numel(A);
nVar = length(C(1).Position);
eta = crossover_params.eta;
pc = crossover_params.pc;
lb = crossover_params.lb;
ub = crossover_params.ub;
TestFunction = crossover_params.TestFunction;
R = C;

for i = 1:nClone
    if rand < pc
        idx = randperm(nActive, 1);
        p = C(i).Position;
        q = A(idx).Position;
        off1 = p;
        off2 = q;
        for j = 1:nVar
            if rand < 0.5
                beta = (2 * rand)^(1 / (eta + 1));
            else
                beta = (2 - 2 * rand)^(-1 / (eta + 1));
            end
            off1(j) = 0.5 * ((1 + beta) * p(j) + (1 - beta) * q(j));
            off2(j) = 0.5 * ((1 - beta) * p(j) + (1 + beta) * q(j));

            off1(j) = max(min(off1(j), ub(j)), lb(j));
            off2(j) = max(min(off2(j), ub(j)), lb(j));
        end
        if rand < 0.5
            R(i).Position = off1;
        else
            R(i).Position = off2;
        end
        R(i).Objective = TestFunction(R(i).Position);
    end
end
end

function pop = StaticHypermutation(pop, mutate_params)

% First Constructive Mutation (FCM)
nPop = numel(pop);
pm = mutate_params.pm;
sm = mutate_params.sm;
eta = mutate_params.eta;
lb = mutate_params.lb;
ub = mutate_params.ub;
TestFunction = mutate_params.TestFunction;

for i = 1:nPop
    mu = pop(i).Position;
    ori = pop(i).Objective;
    rnd = randperm(sm);
    mutated = false;

    for k = 1:sm
        j = rnd(k);

        % polynomial mutation
        if rand < pm
            if rand < 0.5
                delta = (2 * rand)^(1 / (eta + 1)) - 1;
            else
                delta = 1 - (2 * (1 - rand))^(1 / (eta + 1));
            end

            mu(j) = mu(j) + delta * (ub(j) - lb(j));
            mu(j) = max(min(mu(j), ub(j)), lb(j));
            mut = TestFunction(mu);

            if Dominates(mut, ori)
                pop(i).Position = mu;
                pop(i).Objective = TestFunction(mu);
                mutated = true;
                break;
            end
        end

        if mutated
            break;
        end
    end

    pop(i).Position = mu;
    pop(i).Objective = TestFunction(mu);
end
end

function PlotObjs(D)
DO = [D.Objective];
plot(DO(1, :), DO(2, :), 'x');
xlabel('1^{st} Objective');
ylabel('2^{st} Objective');
grid on;
end

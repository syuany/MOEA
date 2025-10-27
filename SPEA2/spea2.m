clc;
clear;
close all;

% settings
TestFunction = @(x) ZDT6(x);
nVar = 10;
VarSize = [1, nVar];
VarMin = 0;
VarMax = 1;
nObj = numel(TestFunction(unifrnd(VarMin, VarMax, VarSize)));

MaxIt = 100;
nPop = 50;
nArchive = 50;

de_params.F = 0.5;
de_params.CR = 0.9;
de_params.xmax = VarMax;
de_params.xmin = VarMin;
de_params.nVar = nVar;
de_params.TestFunction = TestFunction;

% initalization
empty_individual.Position = [];
empty_individual.Objective = [];
empty_individual.Fitness = [];
empty_individual.Distance = [];
empty_individual.DominatedSet = [];
empty_individual.DominationCount = [];

pop = repmat(empty_individual, nPop, 1);
for i = 1:nPop
    pop(i).Position = unifrnd(VarMin, VarMax, VarSize);
    pop(i).Objective = TestFunction(pop(i).Position);
end

EP = [];

% Main Loop
for it = 1:MaxIt
    pq = [pop; EP];
    k = round(sqrt(numel(pq)));
    pq = DetermineDomination(pq);
    pq = CalculateDistance(pq);
    pq = CalculateFitness(pq, k);
    pq = RemoveDuplicates(pq);
    EP = pq([pq.Fitness] < 1);
    EP = RemoveDuplicates(EP);

    n = numel(EP);
    if n > nArchive
        EP = EnvironmentSelection(EP, nArchive, nPop);
    elseif n < nArchive
        so = SortByFitness(pq([pq.Fitness] > 1));
        EP = [EP; so(1:nArchive-n)];
        EP = DetermineDomination(EP);
        EP = CalculateDistance(EP);
        EP = CalculateFitness(EP, nArchive+nPop);
        EP = RemoveDuplicates(EP);
    end

    % binary tournament selection
    Pool = TournamentSelection(EP);

    % recombination and mutation
    pop = DE(empty_individual, Pool, de_params);

    figure(1);
    PlotFits(EP);
    pause(0.01);

    disp(['Iteration ', num2str(it), ': Number of EP = ', num2str(numel(EP))]);
end
%% Reults

disp(' ');

EPO = [EP.Objective];
for j = 1:nObj

    disp(['Objective #', num2str(j), ':']);
    disp(['      Min = ', num2str(min(EPO(j, :)))]);
    disp(['      Max = ', num2str(max(EPO(j, :)))]);
    disp(['    Range = ', num2str(max(EPO(j, :))-min(EPO(j, :)))]);
    disp(['    St.D. = ', num2str(std(EPO(j, :)))]);
    disp(['     Mean = ', num2str(mean(EPO(j, :)))]);
    disp(' ');

end


function b = Dominates(p, q)
if isstruct(p)
    p = p.Objective;
end
if isstruct(q)
    q = q.Objective;
end
b = all(p <= q) && any(p < q);
end

function pop = DetermineDomination(pop)
nPop = numel(pop);
for i = 1:nPop
    pop(i).DominatedSet = [];
    pop(i).DominationCount = 0;
end

for i = 1:nPop
    for j = i + 1:nPop
        p = pop(i);
        q = pop(j);
        if Dominates(p, q)
            q.DominatedSet = [q.DominatedSet, i];
            p.DominationCount = p.DominationCount + 1;
        end
        if Dominates(q, p)
            p.DominatedSet = [p.DominatedSet, j];
            q.DominationCount = q.DominationCount + 1;
        end

        pop(i) = p;
        pop(j) = q;
    end
end
end

function pop = CalculateDistance(pop)
nPop = numel(pop);
for i = 1:nPop
    pop(i).Distance = zeros(1, nPop);
    for j = 1:nPop
        if i ~= j
            pop(i).Distance(j) = norm(pop(i).Objective-pop(j).Objective);
        else
            pop(i).Distance(j) = Inf;
        end
    end
end
end

function pop = CalculateFitness(pop, k)
nPop = numel(pop);
for i = 1:nPop
    % calculate S(i)
    s = 0;
    for j = pop(i).DominatedSet
        s = s + pop(j).DominationCount;
    end

    % calculate D(i)
    dis = sort(pop(i).Distance);
    k = min(k, length(dis));
    d = 1 / (dis(k) + 2);

    pop(i).Fitness = s + d;
end
end

function b = comp(p, q)
pDist = sort(p.Distance);
qDist = sort(q.Distance);
nPop = length(pDist);
for i = 1:nPop - 1
    pSigma = pDist(i+1);
    qSigma = qDist(i+1);
    if pSigma < qSigma
        b = false;
        return;
    elseif pSigma > qSigma
        b = true;
        return;
    end
end
b = false;
end

function pop = EnvironmentSelection(pop, nArchive, nPop)
n = numel(pop);
indices = 1:n;

for i = 1:n - 1
    for j = 1:n - i
        if comp(pop(indices(j)), pop(indices(j+1)))
            indices([j, j + 1]) = indices([j + 1, j]);
        end
    end
end

pop = pop(indices);
pop = pop(1:nArchive);

pop = DetermineDomination(pop);
pop = CalculateDistance(pop);
pop = CalculateFitness(pop, nArchive+nPop);
pop = RemoveDuplicates(pop);
end

function pop = SortByFitness(pop)
Fits = [pop.Fitness];
[~, so] = sort(Fits);
pop = pop(so);
end

function pool = TournamentSelection(pop)
nPop = numel(pop);
pool = [];
for i = 1:nPop
    candidates = randperm(numel(pop), 2);
    candidates1 = pop(candidates(1));
    candidates2 = pop(candidates(2));
    if candidates1.Fitness < candidates2.Fitness
        pool = [pool; candidates1];
    else
        pool = [pool; candidates2];
    end
end
end

function pop = DE(empty_individual, pool, params)
F = params.F;
CR = params.CR;
xmax = params.xmax;
xmin = params.xmin;
D = params.nVar;
TestFunction = params.TestFunction;

nPop = numel(pool);
pop = repmat(empty_individual, nPop, 1);
for i = 1:nPop
    candidates = setdiff(1:nPop, i);
    idx = randperm(nPop-1, 3);
    x1 = pool(candidates(idx(1))).Position;
    x2 = pool(candidates(idx(2))).Position;
    x3 = pool(candidates(idx(3))).Position;

    % DE/rand/1
    v = x1 + F .* (x2 - x3);
    v = max(xmin, min(xmax, v));

    % crossover
    a = rand(1, D) < CR;
    jrand = randi(D);
    a(1, jrand) = 1;
    u = a .* v + (1 - a) .* pool(i).Position;

    pop(i).Position = u;
    pop(i).Objective = TestFunction(u);
end
end

function pop = RemoveDuplicates(pop, epsilon)
if nargin < 2
    epsilon = 1e-6;
end
unique_pop = [];
n = numel(pop);

for i = 1:n
    current = pop(i);
    is_duplicate = false;

    % 比较目标函数值的欧氏距离
    for j = 1:numel(unique_pop)
        dist = norm(current.Objective-unique_pop(j).Objective);
        if dist < epsilon
            is_duplicate = true;
            break;
        end
    end

    if ~is_duplicate
        unique_pop = [unique_pop; current];
    end
end
pop = unique_pop;
end

function PlotFits(EP)

EPO = [EP.Objective];
plot(EPO(1, :), EPO(2, :), 'x');
xlabel('1^{st} Objective');
ylabel('2^{nd} Objective');
grid on;

end

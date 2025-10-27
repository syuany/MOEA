clear;clc;

% init
test_func=@(x) MOP4(x);
XMax=5;
XMin=-5;
D=3;

ItMax=100;
NP=50;
nArch=50;
T=max(ceil(0.1*NP), 2);
T=min(max(T,2), 15);
crossover_params.gamma = 0.5;
crossover_params.XMin = XMin;
crossover_params.XMax = XMax;

X=unifrnd(XMin, XMax, NP, D);
nObj=numel(test_func(X(1,:)));
Fit=cell(NP,1);
z=zeros(1, nObj);
for i=1:NP
    Fit{i}=test_func(X(i,:));
    z=min(z, Fit{i});
end

lambda=rand(NP, nObj);
e=pdist2(lambda, lambda);
[~,so]=sort(e);
neighbors=so(:, 1:T);  

% weight vector
for i = 1:NP
    lambda(i,:) = lambda(i,:)/norm(lambda(i,:));
end

% Tchebycheff decomposition
g = zeros(NP, 1);
for i = 1:NP
    g(i) = max(lambda(i,:) .* abs(Fit{i} - z));
end

isDominated=DetermineDomination(Fit);
EP=find(~isDominated);

for it=1:ItMax
    for i=1:NP
        k=randsample(T,2);
        x1=X(neighbors(k(1)), :);
        x2=X(neighbors(k(2)), :);
        u=Crossover(x1, x2, crossover_params);
        f=test_func(u);
        z=min(z, f);
        
        for j=neighbors(i)
            g_u=max(lambda(j, :).*abs(f-z));
            if(g_u<g(j))
                g(j)=g_u;
                X(j,:)=u;
                Fit{j}=f;
            end
        end
    end
    
    cur_isDominated=DetermineDomination(Fit);
    cep=find(~cur_isDominated);
    EP=[EP, cep];

    if ~isempty(EP)
        EPFit=Fit(EP);
        EP_isDominated=DetermineDomination(EPFit);
        EP=EP(~EP_isDominated);
    end

    if numel(EP)>nArch
        extra = numel(EP)-nArch;
        ToBeDeleted = randsample(numel(EP), extra);
        EP(ToBeDeleted) = [];
    end

    Fitness=cell2mat(Fit(EP));

    figure(1);
    PlotCosts(EP, Fitness);
    pause(0.01);

    disp(['Iteration ' num2str(it) ': Number of Pareto Solutions = ' num2str(numel(EP))]);

end

function b = Dominates(x, y)
    b = all(x <= y) && any(x<y);
end

% Extended Line Crossover
function y = Crossover(x1, x2, params)
    gamma = params.gamma;
    XMin = params.XMin;
    XMax = params.XMax;

    alpha = unifrnd(-gamma, 1+gamma, size(x1));
    y = alpha.*x1+(1-alpha).*x2;
    y = min(max(y, XMin), XMax);
end

function isDominated = DetermineDomination(Fit)

    NP = numel(Fit);
    isDominated=false(1, NP);
    
    for i = 1:NP
        for j = i+1:NP
            if Dominates(Fit{i}, Fit{j})
                isDominated(j) = true;
            elseif Dominates(Fit{j}, Fit{i})
                isDominated(i) = true;
            end
        end
    end
end

function PlotCosts(EP, Fitness)

    plot(Fitness(:, 1), Fitness(:, 2), 'x');
    xlabel('1^{st} Objective');
    ylabel('2^{nd} Objective');
    grid on;

end
clc;clear;

% init
test_func = @(x) MOP4(x);
XMax=5;
XMin=-5;
D=3;
ItMax=100;
NP=50;

PCross=0.7;
NCross=2*round(PCross*NP/2);
PMutation=0.4;
NMutation=round(PMutation*NP);

mu=0.02;
sigma=0.1*(XMax-XMin);

X=unifrnd(XMin, XMax, NP, D);

Fit=cell(NP,1);
for i=1:NP
    Fit{i}=test_func(X(i,:));
end

[F, Rank] = NDSort(Fit);
CD = CalcCD(Fit, F);
[X, Fit, F, CD, Rank] = SortPop(X, Fit, CD, Rank);

for it=1:ItMax
    % crossover
    xc=zeros(NCross,D);
    fc=cell(NCross,1);
    
    for k=1:NCross/2
        i1=randi(NP);
        x1=X(i1,:);
        i2=randi(NP);
        x2=X(i2,:);
        [off1,off2]=Crossover(x1, x2);
        xc(2*k-1,:)=off1;
        xc(2*k,:)=off2;
        fc{2*k-1}=test_func(off1);
        fc{2*k}=test_func(off2);
    end

    % mutation
    xm=zeros(NMutation, D);
    fm=cell(NMutation,1);
    for k=1:NMutation
        i=randi(NP);
        p=X(i,:);
        xm(k,:)=Mutate(p, mu, sigma);
        fm{k}=test_func(xm(k,:));
    end

    % merge
    X=[X;xc;xm];
    Fit=[Fit;fc;fm];

    [F, Rank] = NDSort(Fit);
    CD = CalcCD(Fit, F);
    [X, Fit, F, CD, Rank] = SortPop(X, Fit, CD, Rank);

    X=X(1:NP, :);
    Fit=Fit(1:NP, :);
    Rank=Rank(1:NP);
    CD=CD(1:NP);

    [F, Rank] = NDSort(Fit);
    CD = CalcCD(Fit, F);
    [X, Fit, F, CD, Rank] = SortPop(X, Fit, CD, Rank);

    F1=X(F{1}, :);

    Fitness=cell2mat(Fit);

    disp(['Iteration ' num2str(it) ': Number of F1 Members = ' num2str(size(F1,1))]);
    
    % Plot F1 Costs
    figure(1);
    PlotCosts(Fitness);
    pause(0.01);
end

% disp(X);

function [off1, off2] = Crossover(x1, x2)
    alpha = rand(size(x1));
    off1 = alpha.*x1+(1-alpha).*x2;
    off2 = alpha.*x2+(1-alpha).*x1;
end

function y = Mutate(x, mu, sigma)

    nVar = numel(x);
    
    nMu = ceil(mu*nVar);

    j = randsample(nVar, nMu);
    if numel(sigma)>1
        sigma = sigma(j);
    end
    
    y = x;
    
    y(j) = x(j)+sigma.*randn(size(j));

end

function a = Dominates(x, y)
    a=all(x<=y) && any(x<y);
end

function [F, Rank] = NDSort(Fit)
    NP=size(Fit,1);
    n=zeros(NP,1);
    s=cell(NP,1);
    F{1}=[];
    Rank=zeros(NP,1);
    for i=1:NP
        for j=i+1:NP
            p=Fit{i};
            q=Fit{j};
            if Dominates(p, q)
                s{i}=[s{i} j]; 
                n(j)=n(j)+1;
            end
            if Dominates(q, p)
                s{j}=[s{j} i]; 
                n(i)=n(i)+1;
            end
        end
        if n(i)==0
            F{1}=[F{1} i];
            Rank(i)=1;
        end
    end

    k=1;
    while true
        q=[];
        for i=F{k}
            for j=s{i}
                n(j)=n(j)-1;
                if n(j)==0
                    q=[q j];
                    Rank(j)=k+1;
                end
            end
        end

        if isempty(q)
            break;
        end
        
        k=k+1;
        F{k}=q;
    end
end

function CD = CalcCD(Fit, F)
    NP=size(Fit,1);
    NF=numel(F);
    CD=zeros(NP,1);

    for k=1:NF
        Fits=cell2mat(Fit(F{k}));
        NOBJ=size(Fits, 2);
        n=numel(F{k});
        d=zeros(n, NOBJ);

        for j=1:NOBJ
            [cj, so]=sort(Fits(:,j));
            d(so(1), j)=inf;
            d(so(end), j)=inf;
            range = cj(end) - cj(1);
            if range == 0
                range = eps;  
            end
            for i=2:n-1
                d(so(i),j)=(cj(i+1)-cj(i-1))/range;
            end
        end

        for i=1:n
            CD(F{k}(i))=sum(d(i,:));
        end
    end
end

function [X, Fit, F, CD, Rank] = SortPop(X, Fit, CD, Rank)
    [~, idx] = sortrows([Rank, -CD]);
    
    X = X(idx, :);
    Fit = Fit(idx);
    CD = CD(idx);
    Rank = Rank(idx);
    
    MaxRank = max(Rank);
    F = cell(MaxRank, 1);
    for i = 1:MaxRank
        F{i} = find(Rank == i);
    end
end

function PlotCosts(Fit)
    
    plot(Fit(:, 1), Fit(:, 2), 'r*', 'MarkerSize', 8);
    xlabel('1^{st} Objective');
    ylabel('2^{nd} Objective');
    title('Non-dominated Solutions (F_{1})');
    grid on;

end
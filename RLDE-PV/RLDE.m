function [solution, value, NFE] = RLDE(ObjectFunction, type, Max_NFE)
%% Settings
NP=50;
CR=0.9;
NFE=0;
[LB, UB, nVar] = Parameter(type);

X=LB+(UB-LB).*rand(NP,nVar);

Fits=zeros(1, NP);
for i=1:NP
    Fits(i)=ObjectFunction(type, X(i,:));
    NFE=NFE+1;
end
[so, idx]=sort(Fits);
bestY=so(1);
bestX=X(idx(1),:);

goodArchive=[];
F=normrnd(0.5, 0.3, NP, 1);
ft=[-0.1,0,0.1];
f=ft(randi(length(ft)));

% RL Settings
alpha=0.1;
gamma=0.9;

T=1;

QAgent=zeros(2.*NP, 3);
States=2.*ones(NP, 1);
Actions=zeros(size(States));
TrialStates=States;
Reward=zeros(NP, 1);

it=1;

%% Main Loop
while NFE<Max_NFE
    for i=1:NP
        % initial stage: using standard DE strategy
        if it<2
            r1 = randi(NP);
            while r1 == i
                r1 = randi(NP);
            end

            r2 = randi(NP);
            while r2 == i || r2 == r1
                r2 = randi(NP);
            end

            r3 = randi(NP);
            while r3 == i || r3 == r1 || r3 == r2
                r3 = randi(NP);
            end

            delta = (X(r2, :) - X(r3, :));
        % subsequent stage: 50% probability of using historical successful mutations
        else
            if rand<0.5
                delta=goodArchive(randi(length(goodArchive)));
            else
                % pool=setdiff(1:NP, i);    % slow!
                % idx=randperm(NP-1, 2);    % slow!
                % r1=pool(idx(1));
                % r2=pool(idx(2));
                % delta=X(r1,:)-X(r2,:);
                r1 = randi(NP);
                while r1 == i
                    r1 = randi(NP);
                end

                r2 = randi(NP);
                while r2 == i || r2 == r1
                    r2 = randi(NP);
                end

                r3 = randi(NP);
                while r3 == i || r3 == r1 || r3 == r2
                    r3 = randi(NP);
                end

                delta = (X(r2, :) - X(r3, :));
            end
        end

        F(i)=F(i)+f;
        if F(i)<=0 | F(i)>=1
            F(i)=normrnd(0.5, 0.3);
            % F(i)=max(0.1, min(0.9, F(i)));
        end

        % DE/best/1
        [so, idx]=sort(Fits);
        v(i, :)=X(idx(1), :)+F(i).*delta;
        % a=v(i,:)<LB | v(i,:)>UB;
        % v(i,:)=a.*(LB+(UB-LB).*rand)+(1-a).*v(i,:);
        for j = 1:nVar
            if v(i, j) < LB(j)
                v(i, j) = LB(j) + rand * (UB(j) - LB(j));
            end

            if v(i, j) > UB(j)
                v(i, j) = LB(j) + rand * (UB(j) - LB(j));
            end
        end

        % binary crossover
        % a=rand(1,nVar);
        jrand=randi(nVar);
        % a(1, jrand)=1;
        % u(i,:)=a.*v(i,:)+(1-a).*X(i,:);
        for j = 1:nVar
            if rand <= CR || j == jrand
                u(i, j) = v(i, j);
            else
                u(i, j) = X(i, j);
            end
        end

        % choose and update action
        Qtemp=QAgent(2*i-1:2*i, :);
        curState=States(i);
        Qval=Qtemp(curState,:);
        temp=exp(Qval/T);
        ratio=cumsum(temp)/sum(temp);
        jtemp=find(rand<ratio);
        adjustment=jtemp(1);
        Actions(i)=adjustment;

        switch adjustment
        case 1
            f=-0.1;
        case 2
            f=0;
        case 3
            f=0.1;
        end

        % calculate the fitness of offspring
        newfit(i)=ObjectFunction(type, u(i,:));
        NFE=NFE+1;
        if newfit(i)<Fits(i)
            goodArchive(end+1,:)=u(i,:)-X(i,:);
            X(i,:)=u(i,:);
            Fits(i)=newfit(i);
            Reward(i)=1;
            TrialStates(i,:)=1;
        else
            Reward(i)=0;
            TrialStates(i,:)=2;
        end
    end
    %% Update Q table
    for k=1:NP
        Qtemp=QAgent(2*k-1:2*k, :);
        curState=States(k);
        action=Actions(k);
        nextState=TrialStates(k);
        temp=max(Qtemp(nextState, :));
        Qtemp(curState, action)=(1-alpha)*Qtemp(curState, action)+alpha*(Reward(k)+gamma*temp);
        QAgent(2*k-1:2*k,:)=Qtemp;
    end

    % state transition
    States=TrialStates;

    % random delete
    if length(goodArchive)>NP
        aa=randperm(length(goodArchive));
        bb=length(goodArchive)-NP;
        goodArchive(aa(1:bb),:)=[];
    end

    [so, idx]=sort(Fits);
    bestX=X(idx(1),:);
    bestY=so(1);
    it=it+1;
end
solution=bestX;
value=bestY;
end
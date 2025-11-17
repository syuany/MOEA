function [solution, value, nfe]=CCGDE3(Problem, nVar)
    %% Settings
    LB = 0; 
    UB = 1; 
    NumEsp=2;
    varsPerGroup=floor(nVar/NumEsp);
    Gmax=1;
    nPop=50;
    numSPop=40;
    F=0.5;
    CR=0.5;
    nfe=0;

    %% Initialization
    Index=[];
    for i=1:NumEsp;
        if i<NumEsp
            Index=[Index, ones(1, varsPerGroup)*i];
        else
            Index=[Index, ones(1, nVar-varsPerGroup*(NumEsp-1))*NumEsp];
        end
    end
    Index=Index(randperm(length(Index)));

    subDec=cell(1, NumEsp);
    Population=cell(1, NumEsp);

    Dec=unifrnd(repmat(LB, numSPop, 1), repmat(UB, numSPop, 1));
    for i=1:NumEsp
        subDec{i}=Dec(:, Index==i);
    end
    
    otherVars=[];
    for i=1:NumEsp
        for j=1:NumEsp
            % random choose numSPop individuals from subDec{j}
            otherVars=[otherVars, subDec{j}(randi(size(subDec{j}, 1), numSPop, 1))];
        end
        [Population{i}, nfe]=GetInd(subDec, otherVars, Index, j, numSPop, nVar, nfe);
    end

    pop=EnvironmentalSelection(vertcat(Population{:}), nPop);
    
    %% Main Loop
    while NotTerminated(pop, nfe)
        NDsubDec=cell(1, NumEsp);
        for i=1:NumEsp
            Front=NDSort(Population{i}.Cost, 1);
            NDsubDec{i}=subDec{i}(Front==1, :);
        end
        
        for j=1:NumEsp
            for k=1:Gmax

                % Mutation
                OffDec=CCDE(subDec{j}, subDec{j}(randi(numSPop, 1, numSPop)), subDec{j}(randi(numSPop, 1, numSPop)), LB, UB);
                otherVars=[];
                for i=1:NumEsp
                    if i~=j
                        otherVars=[otherVars, NDsubDec{j}(randi(size(NDsubDec{j}, 1), size(OffDec, 1), 1))];
                    end
                end
                [Offspring, nfe]=GetInd(OffDec, otherVars, Index, j, numSPop, nVar, nfe);

                % Selection
                Population{j}=GDE3_EnvironmentalSelection(Population{j}, Offspring, numSPop);
                
                % Update
                subDec{j}=Population{j}.decs(:, Index==j);
            end
        end
        pop=EnvironmentalSelection(vertcat(Population{:}), nPop);
    end
end

function [Population, nfe] = GetInd(subDec, otherVars, Index, j, numSPop, nVar, nfe, Problem)
    persistent empty_individual
    if isempty(empty_individual)
        empty_individual.Position = [];
        empty_individual.Cost = [];
        empty_individual.Rank = [];
        empty_individual.DominationSet = [];
        empty_individual.DominatedCount = [];
        empty_individual.CrowdingDistance = [];
    end
    
    Dec = zeros(numSPop, nVar);
    Dec(:, Index==j) = subDec;
    Dec(:, Index~=j) = otherVars;
    
    Population = repmat(empty_individual, numSPop, 1);
    
    for i = 1:numSPop
        Population(i).Position = Dec(i, :);
        Population(i).Cost = Problem(Population(i).Position);
        nfe = nfe + 1;
    end
end
function nfe=CCGDE3(nVar)
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

    TestFunction=@(x) ZDT1(x);

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

    Dec=unifrnd(repmat(LB, numSPop, nVar), repmat(UB, numSPop, nVar));
    for i=1:NumEsp
        subDec{i}=Dec(:, Index==i);
    end
    
    for i=1:NumEsp
        otherVars=[];
        for j=1:NumEsp
            if j~=i
                % random choose numSPop individuals from subDec{j}
                otherVars=[otherVars, subDec{j}(randi(size(subDec{j}, 1), numSPop, 1), :)];
            end
        end
        [Population{i}, nfe]=GetInd(subDec{i}, otherVars, Index, i, nVar, nfe, TestFunction);
    end

    pop=EnvironmentalSelection(vertcat(Population{:}), nPop);
    
    %% Main Loop
    while NotTerminated(pop, nfe)
        NDsubDec=cell(1, NumEsp);
        for j=1:NumEsp
            Front=NDSort(Population{j});
            NDsubDec{j}=subDec{j}(Front{1}, :);
        end
        
        for j=1:NumEsp
            for k=1:Gmax
                % Mutation
                OffDec=GDE3(subDec{j}, subDec{j}(randi(numSPop, 1, numSPop), :), subDec{j}(randi(numSPop, 1, numSPop), :), LB, UB);

                tempFull = zeros(size(OffDec, 1), nVar);
                for i=1:NumEsp
                    if i~=j
                         ridx = randi(size(NDsubDec{i}, 1), numSPop, 1);
                         tempFull(:, Index==i) = NDsubDec{i}(ridx, :);
                    end
                end
                otherVars = tempFull(:, Index~=j);

                [Offspring, nfe] = GetInd(OffDec, otherVars, Index, j, nVar, nfe, TestFunction);

                % GDE3 Selection
                Population{j}=GDE3_EnvironmentalSelection(Population{j}, Offspring, numSPop);
                
                % Update
                Decs = vertcat(Population{j}.Position);
                subDec{j}=Decs(:, Index==j);
                % subDec{j}=[Population{j}.Position](:, Index==j);
            end
        end
        pop=EnvironmentalSelection(vertcat(Population{:}), nPop);
        
        % 绘图优化：每 2000 nfe 画一次，或者是每代画一次但用 drawnow limitrate
        % if mod(nfe, 2000) < 100 % 简单的控制
        %      PlotCosts(pop);
        %      drawnow limitrate; % 限制刷新率，防止卡死
        % end
    end
end

function [Population, nfe] = GetInd(subDec, otherVars, Index, j, nVar, nfe, TestFunction)
    % subDec: [N x subD]
    % otherVars: [N x otherD]
    
    numSPop = size(subDec, 1);
    
    Dec = zeros(numSPop, nVar);
    Dec(:, Index==j) = subDec;
    Dec(:, Index~=j) = otherVars;
    
    Costs = TestFunction(Dec); % 输出应为 [2 x numSPop]
    
    nfe = nfe + numSPop;
    
    % 批量构建 struct （最高效）
    % 将矩阵转为 cell array 以便分配给 struct
    DecCell = mat2cell(Dec, ones(1, numSPop), nVar);
    CostCell = mat2cell(Costs, 2, ones(1, numSPop)); % 注意 Cost 是 2xN
    
    Population(numSPop, 1).Position = []; 
    
    % 批量赋值
    [Population.Position] = DecCell{:};
    [Population.Cost] = CostCell{:};
    
    % 初始化其他字段 
    [Population.Rank] = deal([]);
    [Population.DominationSet] = deal([]);
    [Population.DominatedCount] = deal([]);
    [Population.CrowdingDistance] = deal([]);
end
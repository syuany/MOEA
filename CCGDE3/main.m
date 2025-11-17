clear all
close all

rehash

format long g

warning ('off')

nObj = numel(CostFunction(unifrnd(VarMin, VarMax, VarSize)));

for type=1:7
    tic()
    if type==1
        nVar=200;
    elseif type==2
        nVar=500;
    elseif type==3
        nVar=1000;
    elseif type==4
        nVar=2000;
    elseif type==5
        nVar=3000;
    elseif type==6
        nVar=4000;
    elseif type==7
        nVar=5000;
    end
    
    [solution, value]=CCGDE3(nVar);
    fprintf('finished with the number of decision variables: %d', nVar);
    toc()
end

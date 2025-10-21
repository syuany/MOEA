
%=============================================================================================
%
% This sourse code is programmed by Zhaolu Guo
% E-mail:gzl@whu.edu.cn
% School of Science, JiangXi University of Science and Technology, Ganzhou 341000, China
%
%============================================================================================


clear all;
clc


Xmax = 100;
Xmin = -100;
D = 30;


NP = 50;    %种群规模
F  = 0.5;   %缩放因子
CR = 0.9;   %交叉概率
MAX_FES = 200000;

rand('seed', sum(100 * clock));

t = 1;
%产生初始种群
X = Xmin + rand(NP,D) .* (Xmax - Xmin);
Fits = test_func(X);

[BestFit, BestI] = min(Fits);
BestX = X(BestI, :);

FES = NP;
while FES < MAX_FES
    
    for i = 1 : NP
        %产生三个不同的数,并且也不同于i
        dx = randperm(NP);
        index = find( i == dx(1, :) );
        dx(index) = [];
        r1 = dx(1);
        r2 = dx(2);
        r3 = dx(3);
        
        %变异算子
        V = X(r1, :) + F .* ( X(r2, :) - X(r3, :) );
        
        %越界处理
        a = V > Xmax;
        V = ( Xmin + (Xmax - Xmin) .* rand(1,D) ) .* a + V .* (1 - a);
        
        a = V < Xmin;
        V = ( Xmin + (Xmax - Xmin) .* rand(1,D) ) .* a + V .* (1 - a);
        
        %交叉算子
        XI = X(i, :);
        a = rand(1, D) < CR;
        jrand = 1 + floor(rand() * D);
        a(1,jrand) = 1;
        U = V .* a + (1 - a) .* XI;
        
        %选择算子
        fitXI = Fits(i);
        fitU  = test_func(U);
        FES = FES + 1;
        if (fitU < fitXI)
            X(i, :) = U;
            Fits(i) = fitU;
            
            %保存最优个体
            if (fitU < BestFit)
                BestFit = fitU;
                BestX = X(i, :);
            end
        end
        
    end
    
    t = t + 1;    
end

BestFit
BestX


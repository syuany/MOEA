
%=============================================================================================
%
% This sourse code is programmed by Zhaolu Guo
% E-mail:gzl@whu.edu.cn
% School of Science, JiangXi University of Science and Technology, Ganzhou 341000, China
%
%============================================================================================


clear all;
clc

rand('seed', sum(100*clock));


D = 30;
Xmax = 100;
Xmin = -100;

Vmax = (Xmax - Xmin)*0.5;
Vmin = -Vmax;

NP = 30; 
C1 = 1.49618;
C2 = C1; 
W  = 0.72984;
MAX_FES = 200000;


FES = 0;
t = 0;

%产生初始种群
X = Xmin + rand(NP,D) .* (Xmax - Xmin);
V = Vmin + rand(NP,D) .* (Vmax - Vmin);
Fits = test_func(X);

[GBestFit, BestI] = min(Fits);
GBestX = X(BestI, :);

PBestX = X;
PBestFit = Fits;

while FES < MAX_FES
    
    for i = 1 : NP
        
        %更新速度
        TemV = W .* V(i, :) + C1 .* rand(1, D) .* ( GBestX - X(i, :) ) + C2 .* rand(1, D) .* ( PBestX(i, :) - X(i, :) );
        
        %速度的越界处理
        a = TemV > Vmax;
        TemV = (Vmin +(Vmax - Vmin) .* rand(1, D)) .* a + TemV .* (1 - a);
        
        a = TemV < Vmin;
        TemV = (Vmin + (Vmax - Vmin) .* rand(1, D)) .* a + TemV .* (1 - a);
        V(i, :) = TemV;
        
        %更新位置
        TemX = X(i, :) + V(i, :);
        
        %个体的越界处理
        a = TemX > Xmax;
        TemX = (Xmin + (Xmax - Xmin) .* rand(1, D)) .* a + TemX .* (1 - a);
        
        a = TemX < Xmin;
        TemX = (Xmin + (Xmax - Xmin) .* rand(1, D)) .* a + TemX .* (1 - a);
        
        TemFit  = test_func(TemX);
        FES = FES + 1;
        X(i, :) = TemX;
        
        
        %保存个体的最优位置
        if (TemFit < PBestFit(i))
            PBestX(i, :) = TemX;
            PBestFit(i) = TemFit;
            
            %保存整个种群的最优位置
            if (TemFit < GBestFit)
                GBestX = TemX;
                GBestFit = TemFit;
            end
        end
    end
    
    t = t + 1;
end
GBestX
GBestFit



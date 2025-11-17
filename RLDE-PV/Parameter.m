function [ LB,UB,Dim ] = Parameter( type )
%PARAMETER 此处显示有关此函数的摘要
%   此处显示详细说明
    if type==1
        Dim=5;
        UB=zeros(1,Dim);  %上边界
        UB(1) = 1;
        UB(2) = 1e-6;
        UB(3) = 0.5;
        UB(4) = 100;
        UB(5) = 2;

        LB=zeros(1,Dim); %下边界
        LB(5) = 1;
    elseif type==2
        Dim=7;
        UB=zeros(1,Dim);  %上边界
        UB(1) = 1;
        UB(2) = 1e-6;
        UB(3) = 0.5;
        UB(4) = 100;
        UB(5) = 2;
        UB(6) = 1e-6;
        UB(7) = 2;

        LB=zeros(1,Dim); %下边界
        LB(5) = 1;
        LB(7) = 1;
    elseif type==3
        Dim=5;
        UB=zeros(1,Dim);  %上边界
        UB(1) = 2;
        UB(2) = 5e-5;
        UB(3) = 2;
        UB(4) = 2000;
        UB(5) = 50;

        LB=zeros(1,Dim); %下边界
        LB(5) = 1;
        
    elseif type==4
        Dim=5;
        UB=zeros(1,Dim);  %上边界
        UB(1) = 2;
        UB(2) = 5e-5;
        UB(3) = 0.36;
        UB(4) = 1000;
        UB(5) = 60;

        LB=zeros(1,Dim); %下边界
        LB(5) = 1;
    elseif type==5
        Dim=5;
        UB=zeros(1,Dim);  %上边界
        UB(1) = 8;
        UB(2) = 5e-5;
        UB(3) = 0.36;
        UB(4) = 1500;
        UB(5) = 50;

        LB=zeros(1,Dim); %下边界
        LB(5) = 1;
    elseif type==6
        Dim=3;
        UB=zeros(1,Dim);  %上边界
        UB(1) = 0.5;
        UB(2) = 1000;
        UB(3) = 2;

        LB=zeros(1,Dim); %下边界
        LB(1) = 0.01;
        LB(2) = 100;
        LB(3) = 1;
    end
    

end


clc;
clear;
close all;

% Settings
nPop=100;
MaxFES=10000;
TestFunction=@(x) sphere_(x);

nVar=30;
VarMax=5;
VarMin=-5;
VarSize=[1 nVar];

empty_individual.Position=[];
empty_individual.Value=[];

pop=repmat(empty_individual, nPop, 1);
for i=1:nPop
    pop(i).Position=unifrnd(VarMin, VarMax, VarSize);
    pop(i).Value=TestFunction(pop(i).Position);
end

FES=0;      % Function Evaluations
while FES<MaxFES
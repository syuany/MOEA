clear;clc;

load('rastrigin_func_data.mat')
load('rastrigin_M_D30.mat')

D=30;
NP=10.*D;
MAX_ITE=10000;
ite=0;
w=0.6;
c1=1.8;
c2=2.2;
XMAX = 5;
XMIN = -XMAX;
k=0.3;
VMAX=k.*(XMAX-XMIN);
VMIN=-VMAX;

if length(o) >= D
    o = o(1 : D);
else
    o = XMIN + (XMAX-XMIN) * rand(1, D);
end

P=XMIN+(XMAX-XMIN).*rand(NP,D);
V=VMIN+(VMAX-VMIN).*rand(NP,D);

pbest=P;
pbest_val=zeros(NP,1);
for i=1:NP
    pbest_val(i)=test_func(P(i,:), o, M);
end
[gbest_val,index]=min(pbest_val);
gbest=P(index);

while ite<MAX_ITE
    V=w.*V+c1.*rand(NP,D).*(pbest-P)+c2.*rand(NP,D).*(gbest-P);
    P=P+V;
    
    a=P>XMAX | P<XMIN;
    b=V>VMAX | V<VMIN;
    P=a.*(XMIN+(XMAX-XMIN).*rand(NP,D))+(1-a).*P;
    V=b.*(VMIN+(VMAX-VMIN).*rand(NP,D))+(1-b).*V;
    
    for i=1:NP
        x=test_func(P(i,:), o, M);
        if x<pbest_val(i)
            pbest(i,:)=P(i,:);
            pbest_val(i)=x;
        end
        if x<gbest_val
            gbest_val=x;
            gbest=P(i,:);
        end
    end
    
    ite=ite+1;
end

disp(gbest_val);
disp(gbest);

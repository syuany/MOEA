clear;clc;

load('rastrigin_func_data.mat')
load('rastrigin_M_D30.mat')

Lb=-5;
Ub=5;
D=30;
NP=10.*D;
F=0.5;
CR=0.1;
FES_MAX=10000000;

if length(o) >= D
    o = o(1 : D);
else
    o = Lb + (Ub-Lb) * rand(1, D);
end

X=Lb+(Ub-Lb).*rand(NP, D);
Fits=zeros(NP, 1);
for i=1:NP
    Fits(i)=test_func(X(i,:), o, M);
end
FES=NP;
[best_val, best]=min(Fits);
best_x=X(best,:);

while FES<FES_MAX
    for i=1:NP
        pool=setdiff(1:NP, i);
        idx=randperm(NP-1, 3);
        r1=pool(idx(1));
        r2=pool(idx(2));
        r3=pool(idx(3));

        % DE/rand/1
        v=X(r1,:)+F.*(X(r2,:)-X(r3,:));
        a=v<Lb | v>Ub;
        v=a.*(Lb+(Ub-Lb).*rand(1,D))+(1-a).*v;

        % 二项式交叉
        a=rand(1,D)<CR;
        jrand=randi(D);
        a(1,jrand)=1;
        u=a.*v+(1-a).*X(i,:);

        f=test_func(u, o, M);
        if f<Fits(i)
            Fits(i)=f;
            X(i,:)=u;
        end
        if f<best_val
            best_val=f;
            best_x=u;
        end
        FES=FES+1;
        if FES>=FES_MAX
            break;
        end
    end
end

best_val
best_x
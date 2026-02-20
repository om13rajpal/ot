# BFS

```
clc
clear

c = [2 3 0 0];

A = [1 1 1 0;
     2 1 0 1];

b = [4;5];

z = @(x) c*x;

[m,n] = size(A);

basicSol = [];
bfSol = [];
degSol = [];

pair = nchoosek(1:n,m);

for i=1:size(pair,1)

    y=zeros(n,1);
    bv_index=pair(i,:);
    B=A(:,bv_index);

    if rank(B)<m
        continue
    end

    X=B\b;
    y(bv_index)=X;

    basicSol=[basicSol y];

    if all(X>=0)
        bfSol=[bfSol y];

        if any(X==0)
            degSol=[degSol y];
        end
    end
end

disp('Basic Solutions:')
disp(basicSol)

disp('Basic Feasible Solutions:')
disp(bfSol)

disp('Degenerate BFS:')
disp(degSol)

cost=z(bfSol);
disp(cost)

[opt_val,index]=max(cost);
opt_sol=bfSol(:,index);

disp('Optimal Solution:')
disp(opt_sol)

disp('Optimal Value:')
disp(opt_val)
```

# GRAPH

```
clc;
clear;
close all;

% Phase 1 - Input
c = [3 -5];
A = [1 1; 2 -1];
b = [6; 9];

z  = @(x1,x2) 3*x1 - 5*x2;
c1 = @(x1,x2) x1 + x2 - 6;
c2 = @(x1,x2) 2*x1 - x2 - 9;

[m, n] = size(A);

% Phase 2 - Plot Constraints
x1 = 0:0.1:10;

figure;
for i = 1:m
    if A(i,2) ~= 0
        x2 = (b(i) - A(i,1)*x1) / A(i,2);
        plot(x1,x2,'LineWidth',2);
        hold on;
    end
end

xlabel('x1');
ylabel('x2');
grid on;

% Phase 3 - Intersection and Corner Points
pt = [];
A1 = [A; -eye(2)];
b1 = [b; 0; 0];
m1 = size(A1,1);

for i = 1:m1
    for j = i+1:m1
        aa = [A1(i,:); A1(j,:)];
        bb = [b1(i); b1(j)];
        if det(aa) ~= 0
            X = aa\bb;
            if all(X >= 0)
                pt = [pt X];
            end
        end
    end
end

disp('Corner Points:');
disp(pt);

% Phase 4 - Feasible Points
FP = [];
Z  = [];

for i = 1:size(pt,2)
    PT1 = pt(1,i);
    PT2 = pt(2,i);
    if c1(PT1,PT2) <= 0 && c2(PT1,PT2) <= 0
        FP = [FP pt(:,i)];
        plot(PT1,PT2,'*r','MarkerSize',10)
        Z = [Z z(PT1,PT2)];
    end
end

disp('Feasible Points:');
disp(FP)
disp('Objective Values:');
disp(Z)

% Phase 5 - Optimal Solution (Minimization)
[optimal_value, index] = min(Z);
optimal_sol = FP(:,index);

disp('Optimal Solution (x1, x2):');
disp(optimal_sol)
disp('Optimal Value:');
disp(optimal_value)

hold off;
```

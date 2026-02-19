# DFS

```
clc
clear all

% Phase 1 input parameter
c=[6 5 0 0];
a=[1 1 1 0;
   3 2 0 1];
b=[5;12];
z  = @(x)c*x;

m=size(a,1); % no of constraints
n=size(a,2); % no of variables

% Phase 2
basicsol=[];
bfsol=[];
ncm=nchoosek(n,m);
pair=nchoosek(1:n,m);

for i = 1:ncm
    basic_index=pair(i,:);
    y=zeros(n,1);
    x=a(:,basic_index)\b;
    y(basic_index)=x;
    basicsol=[basicsol y];
    if(x>=0)
        bfsol = [bfsol y];
    end
end

disp(basicsol)
disp(bfsol)

% Phase 3 optimal solution and optimal value
cost= z(bfsol);
[optimal index]=max(cost);
optsol=bfsol(:,index);
```

# GRAPH

```
clc
clear all
%Phase 1 input 
c=[3,-5];
a=[1 1;2 -1];
b=[6;9];
z=@(x1,x2)3*x1-5*x2;
c1=@(x1,x2)x1+x2-6;
c2=@(x1,x2)2*x1-x2-9;
m=size(a,1);
n=size(a,2);
%phase 2 plotting
x1=0:max(b./a(:,1))
for i=1:m
    x2=(b(i)-a(i,1)*x1)/a(i,2)
    plot(x1,x2)
    hold on;
end
%phase 3 intersection and find corner points
a=[a;eye(2)];
b=[b;zeros(2,1)];
m=size(a,1);
pt=[];

for i=1:m
    for j=i+1:m
        aa=[a(i,:);a(j,:)];
        bb=[b(i);b(j)];
        if(det(aa)~=0)
            x=inv(aa)*bb
            if(x>=0)
             pt=[pt x]
            end
        end
    end
end
% Phase 4 
FP=[];
Z=[];
for i=1:size(pt,2)
    pt1=pt(1,i);
    pt2=pt(2,i);
    if(c1(pt1,pt2)<=0&&c2(pt1,pt2)<=0)
        FP=[FP,pt(:,i)];
        plot(pt1,pt2,'*r','MarkerSize',10);
        cost=z(pt1,pt2);
        Z=[Z cost];
    end
end
disp(FP)
disp(Z)
hold off
%Phase 5
[optimal_value index]=min(Z)
optimal_sol=FP(:,index)
```

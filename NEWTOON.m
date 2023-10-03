clc
clear all
close all
% optimization(newton)
%      min f(x)=1/2(x'Qx-b'x)
%          unconstraint problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generating problem
n=input('n = ');
%n=500;
Q=10*(rand(n)-rand(n));
lam=1e-2;
Q=Q'*Q+(lam*eye(n));
xx=5*(rand(n,1)-rand(n,1));%xx=argmin f(x)
b=Q*xx;
%******
x=5*(rand(n,1)-rand(n,1));
x0=x;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newton Method
g=Q*x-b;%gradian0
H=Q;%Hessian
d=-H\g;%d0
x5=x;
maxiter=10;
tol=1e-8;
iter=0;
alpha=1;
tic
while ((norm(g,inf)>=tol)&& maxiter>=iter)
    x=x+alpha*d;  
    g=Q*x-b;
    d=-H\g;%invH
    iter=iter+1;
end
timenew=toc;
gradian=norm(g,inf);
%%%%%%%%%%%%%%%%%
% quadprog
tic
[x1,~,~,output,~] = quadprog(Q,-b,[],[]...
    ,[],[],x5,[]);
timeqp=toc;
T=['quadprog time is :',num2str(timeqp)];
disp(T);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table part
Iters=iter;
Ngrad_inf=gradian;
Time=timenew;
Method={'Newton';};
Answer=table(n,Ngrad_inf,Time,Iters,tol,'RowNames',Method)

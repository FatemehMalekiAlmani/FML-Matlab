clear all
close all
clc
warning off
format short g
%% min .5 x'Qx-p'x
%       Ax=b
%% parameters
% n=4;
n=input('n = ');
Q=10*(rand(n)-rand(n));
lam=1e-2;
Q=Q'*Q+(lam*eye(n));
xx=5*(rand(n,1)-rand(n,1));%xx=argmin f(x)
p=Q*xx;
%******
x=5*(rand(n,1)-rand(n,1));
x0=x;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newton Method
g=Q*x-p;%gradian0
H=Q;%Hessian
d=-H\g;%d0
x5=x;
maxiter=10;
tol=1e-8;
iter=0;
alpha=1;
A=5*(rand(n)-rand(n));
b=A*xx;
lam1=ones(n,1);
lam2=zeros(n,1);
lam_total=(lam2-lam1);
c=1e-6;
tic
while (norm(lam_total,inf)>=tol && maxiter>=iter)
    while ((norm(g,inf)>=tol)&& maxiter>=iter)
        x=x+alpha*d;
        g=Q*x-p;
        d=-H\g;%invH
        iter=iter+1;
        x1=x;
    end
    h=A*x1-b;
    lam1=lam2;
    lam2=lam1+c*h;
    lam_total=(lam1-lam2);
end
timenew=toc;
gradian=norm(g,inf);
%% table part
  LAM_T=norm(lam_total,inf);
  Iters=iter;
  Ngrad_inf=gradian;
  Time=timenew;
  Method={'augmented Lagrange';};
  Answer=table(n,Ngrad_inf,LAM_T,Time,Iters,tol,'RowNames',Method)
%% quadprog part
% quadprog
tic
[x_quad,fval,~,output.iterations,~] = quadprog(Q,-p,[],[],A,b,[],[]);
timeqp=toc;
T=['quadprog time is :',num2str(timeqp)];
disp(T);
quad_detail=output.iterations;
disp('For quadprog details you can type quad_detail:')

  
  
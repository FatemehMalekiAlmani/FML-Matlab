close all
clear all
clc
warning off
%% optimazation min f(x)= .5(||(Ax-b)+||^2)

% eps=input('eps= ');
% alpha0=input('alpha0= ');
m=input('m= ');
n=input('n= ');
% m=800;n=900;
A=5*(rand(m,n)-rand(m,n));
x0=5*(rand(n,1)-rand(n,1));
b=(A*x0)+ .05*rand(m,1);
x=5*(rand(n,1)-rand(n,1));
alpha0=.5;
t1=max(0,A*x-b);
g=A'*(t1);
d=-g;
tol=1e-8;
eps=1e-4;
% mu=2;
mu=99;
k=0;
j=0;
% i=0;
maxk=n+100;
x1=x;
alpha=alpha0;
%% main part :
t1=max(0,A*x-b);
t2=max(0,A*(x1)-b);
t3=.5*(norm(t1))^2;
t4=.5*(norm(t2))^2;
%% begin mathod :
while((t4>=t3-eps*(alpha^2)*(norm(d)^2)) && maxk>=j)
    d=-g;
    if (alpha*(norm(d))<tol)
        alpha=0;
    else
        alpha=alpha*mu;
    end
    alpha1=alpha;
    j=j+1;
end
tic

while (norm(g,inf)>=tol && maxk>=k ) %&& abs(alpha-mu*alpha0)>tol )
    g=A'*t1;
    d=-g;
    x=x+alpha1*d;
    x1=x;
    t1=max(0,A*x-b);
    k=k+1;
end
time=toc;
ng=norm(g,inf);
%% quadprog
%    tic
%    N = quadplus( A,x,b,m );
%    [xx,~,~,output,~] = quadprog(N,[],[],[],[],[],x,[]);
% 800   timeqp=toc;
%    T=['quad prog time is :',num2str(timeqp)];
%    disp(T);
%% Table part
Iters=k;
Ngrad_inf=ng;
Time=time;
Method={'Armijo'};
Answer=table(m,n,Ngrad_inf,Time,Iters,alpha,tol,'RowNames',Method)
clear all
close all
clc
format short g
warning off
%% newtoon method for min c'x s.t. A*x=b & x>=0
%n=input('n = ');
%m=input('m = ');
%A=sprand(m,n);
load afiro
[m,n]=size(A)
x=rand(n,1);
c=rand(n,1);
e=ones(n,1);
lam=zeros(m,1);
s=c;
mu=x'*s/n;
%%
b=A*x;
c=A'*lam+s;
iteration=0;
gamma=1e-3;
sigmaM=.75;
sigmam=1e-2;
t=x.'*s-gamma*mu;
tic
while (x'*s>1e-10 && min(t)>=-(1e-10))
    delta_total=([zeros(n,n) A' eye(n) ;...
        A zeros(m,m) zeros(m,n);...
        diag(s) zeros(n,m) diag(x)])\...
        ([zeros(n,1);zeros(m,1);.1*mu-x.*s]);
    dx=delta_total(1:n);
    dlam=delta_total(n+1:n+m);
    ds=delta_total(m+n+1:end);%2*n+m
    t1=find(dx<0);
    t2=find(ds<0);
    alpha_1=min(-(x(t1)./dx(t1)));
    alpha_2=min(-(s(t2)./ds(t2)));
    alpha_total=.99.*min(alpha_1,alpha_2);
    alpha_total=.99*min(1,.8*alpha_total);
    mu=x'*s/n;
    x=x+alpha_total*dx;
    s=s+alpha_total*ds;
    lam=lam+alpha_total*dlam;
    t=x.'*s-gamma*mu;
    x1=x;
    lam1=lam;
    s1=s;
    iteration=iteration+1;
end
TIME=toc;
%% table part
ITERATION=iteration;
Primal=norm(A*x-b);
Dual=norm(A'*lam+s-c);
gap=x'*s;
MINS=min(s);
MINX=min(x);
Method={'LPF';};
Answer=table(n,m,TIME,ITERATION,Primal,Dual,gap ,MINX,MINS,'RowNames',Method)
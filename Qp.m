clear all
close all
clc
format short g
warning off
%% newtoon method for min .5*x'Qx+c'x s.t. A*x=b & x>=0
%n=input('n = ');
%m=input('m = ');
m=10;n=12;
A=rand(m,n);
Q=10*(rand(n)-rand(n));
lam=1e-2;
Q=Q'*Q+(lam*eye(n));
% load pilot4
% [m,n]=size(A)
x=rand(n,1);
c=rand(n,1);
e=ones(n,1);
lam=zeros(m,1);
s=c+Q*x;
mu=x'*s/n;
%%
b=A*x;
c=A'*lam+s-Q*x;
iteration=0;
eps=1e-10;
tu=(x'*s)/n;
tic
while tu>eps
    r_q=A'*lam+s-c-Q*x;
    r_p=A*x-b;
    r_c=x.*s-.1*tu*e;
    delta_total=([-Q A' eye(n) ;...
        A zeros(m,m) zeros(m,n);...
        diag(s) zeros(n,m) diag(x)])\...
        ([r_q;r_p;r_c]);
    dx=delta_total(1:n);
    dlam=delta_total(n+1:n+m);
    ds=delta_total(m+n+1:end);%2*n+m
    t1=find(dx<0);
    t2=find(ds<0);
    alpha_1=min(-(x(t1)./dx(t1)));
    alpha_2=min(-(s(t2)./ds(t2)));
    alpha_total=min(alpha_1,alpha_2);
    alpha_total=min(1,alpha_total);
    x=x-alpha_total*dx;
    s=s-alpha_total*ds;
    lam=lam-alpha_total*dlam;
    iteration=iteration+1;
    tu=(x'*s)/n;
end
TIME=toc;
%% table part
ITERATION=iteration;
ALPHA=alpha_total;
P=norm(A*x-b);
D=norm(A'*lam+s-Q*x);
gap=x'*s;
MINS=min(s);
MINX=min(x);
Method={'Qp'};
Answer=table(n,m,TIME,ITERATION,P,D,gap ,MINX,MINS,'RowNames',Method)
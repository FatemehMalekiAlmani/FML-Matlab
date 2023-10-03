clc
clear all
close all
%%      min f(x)=1/2(x'Qx-b'x)
%          unconstraint problem
%%
n=input('n = ');
m=input('m = ');
A=5*(rand(m,n)-rand(m,n));
Q=5*(rand(n)-rand(n));
%A=rand(m,n);
%Q=rand(n,n);
mu=1e-2;
Q=Q'*Q+(mu*eye(n));
xx=5*(rand(n,1)-rand(n,1));%xx=argmin f(x)
%xx=rand(n,1);
b=A*xx;
%x=5*(rand(n,1)-rand(n,1));
c=rand(n,1);
S=([Q -A'; A zeros(m,m)]);
P=([c;b]);
%lam=zeros(m,1);
c=rand(n,1);
s=c;
%% method 1
tic
z1=S\P;
x=z1(1:n);
lam=z1(n+1:m+n);
z1=[x;lam];
c=A'*lam+s;
TimB=toc;
%% method 2
tic
z2=inv(S)*P;
TimeInv=toc;
%% method 3
tic
z3=P'*inv(S)*P;
TimePIP=toc;
%% table part
Error=norm((x-xx),inf);
Method={'compare';};
Answer=table(n,m,TimB,TimeInv,TimePIP,Error,'RowNames',Method)
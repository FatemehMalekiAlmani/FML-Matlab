clear all
close all
clc
warning off
format short g
%%
% n=input('n = ');
% m=input('m = ');
m=10;n=10;
%%
A=rand(m,n);
b=rand(m,1);
c=rand(n,1);
e = ones(n,1);
lam0=(A*A')\(A*c); 
s0=c-A'*lam0; 
x0=A'*((A*A')\b);
deltax=max(-1.5*min(x0),0);
deltas=max(-1.5*min(s0),0);
deltaxbar=deltax+(0.5*(x0+deltax*e)'*(s0+deltas*e))/(sum(s0.*deltas));
deltasbar= deltas+(0.5*(x0+deltax*e)'*(s0+deltas*e))/(sum(x0.*deltax));
x=x0+deltaxbar*e;
s=s0+deltasbar*e;
lam=lam0;
%%
% A=5*(rand(m,n)-rand(m,n));
% b=rand(m,1);
% c=rand(n,1);
% x=A'*(A*A')\b;
% lam=(A*A')\A*c;
% %b=A*x;
% s=c;
% r_b=A*x-b;
% r_c=A'*lam+s-c;
% s=c-A'*lam;
%%
%gamma=1e-3;
gamma=0.5;
sigmaM=0.25;
sigmam=0.1;
sigma=(sigmaM-sigmam)/2;
mu=.1*(x'*s/n);
tol=1e-10 ;
iteration=0;
d_1=x<0;
x(d_1)=.1;
d_2=s<0;
s(d_2)=.1;
%% METHOD IP
r_b=A*x-b;
r_c=A'*lam+s-c;
tic
while x'*s>tol  %&& norm(r_c,inf)>=tol && norm(r_b,inf)>=tol 
    delta_total=([zeros(n,n) A' eye(n) ;....
        A zeros(m,m) zeros(m,n);...
        diag(s) zeros(n,m) diag(x)])\...
        ([-r_c;-r_b;.1*mu-x.*s]);
    dx=delta_total(1:n);
    dlam=delta_total(n+1:n+m);
    ds=delta_total(m+n+1:end);%2*n+m
    t1=find(dx<0);
    t2=find(ds<0);
    alpha_1=min(-(x(t1)./dx(t1)));
    alpha_2=min(-(s(t2)./ds(t2)));
    alpha_total=.99.*min(alpha_1,alpha_2);
    alpha_total=.99*min(1,.8*alpha_total);
    mu=.1*(x'*s/n);
    x=x+alpha_total*dx;
    s=s+alpha_total*ds;
    lam=lam+alpha_total.*dlam;
    iteration=iteration+1
    
    r_b=A*x-b;
    r_c=A'*lam+s-c;
    gap=x'*s
end
Time=toc;
%% Table part
ITERATION=iteration;
Primal=norm(A*x-b);
Dual=norm(A'*lam+s-c);
gap=x'*s
% MINS=min(s);
% MINX=min(x);
% NormRC=norm(r_c,inf);
% NormRB=norm(r_b,inf);
% TIME=Time;
% Method={'IP';};
% Answer=table(n,m,TIME,ITERATION,Primal,Dual,gap ,NormRC,NormRB,MINX,MINS,'RowNames',Method)


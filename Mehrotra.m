clear all
close all
clc
format short g
warning off
%% PRIMAL problem
%n=input('n = ');
%m=input('m = ');
%  m=100;
%  n=200;
load afiro
%  load lpi_cplex2
%  A=Problem.A;
 [m,n]=size(A);
%  b=Problem.b;
% A=5*(rand(m,n)-rand(m,n));
b=A*rand(n,1);
c=rand(n,1);
lam=(A*A')\(A*c);
s=c-A'*lam;
x=A'*((A*A')\b);
r_b=A*x-b;
r_c=A'*lam+s-c;
maxiter=100;
%% PARAMETERS PART
mu=((x'*s)/n);
tol=1e-10 ;
iteration=0;
%t=x.'*s-gamma*mu;
d_1=x<0;
x(d_1)=.1;
d_2=s<0;
s(d_2)=.1;
t5=(norm(r_b,inf))/(1+norm(b,inf));
t6=(norm(r_c,inf))/(1+norm(c,inf));
%% BEGIN METHOD
tic
while (x'*s>tol || t5>=tol || t6>=tol) && iteration<maxiter
    delta_total1=([zeros(n,n) A' eye(n) ;...
        A zeros(m,m) zeros(m,n);...
        diag(s) zeros(n,m) diag(x)])\...
        ([-r_c;-r_b;-x.*s]);
    dx_off=delta_total1(1:n);
    dlam_off=delta_total1(n+1:n+m);
    ds_off=delta_total1(m+n+1:end);%2*n+m 
    t1=find(dx_off<0);
    t2=find(ds_off<0);
    alpha_1=min(-(x(t1)./dx_off(t1)));
    alpha_2=min(-(s(t2)./ds_off(t2)));
            if isempty(alpha_1)==1
            alpha_1=1 ;
         else
     alpha_1=min(-(x(t1)./dx_off(t1)));
            end
            if isempty(alpha_2)==1
            alpha_2=1 ;
         else
     alpha_2=min(-(s(t2)./ds_off(t2)));
            end
%      alpha_prioff=min(1,.9999*alpha_1);
%      alpha_dualoff=min(1,.9999*alpha_2);
    mu_off=((x+alpha_1*dx_off)'*(s+alpha_2*ds_off))/n;
    sigma=(mu_off/mu)^3;
%=============================================================================================
    delta_total2=([zeros(n,n) A' eye(n) ;...
        A zeros(m,m) zeros(m,n);...
        diag(s) zeros(n,m) diag(x)])\...
        ([zeros(n,1);zeros(m,1);(sigma*mu)*ones(n,1)]);
    dx=delta_total2(1:n);
    dlam=delta_total2(n+1:n+m);
    ds=delta_total2(m+n+1:end);%2*n+m
    dx_total=dx+dx_off;
    ds_total=ds+ds_off;
    dlam_total=dlam+dlam_off;
     t3=find(dx_total<0);
     t4=find(ds_total<0);
      alpha_3=min(-(x(t3)./dx_total(t3)));
      alpha_4=min(-(s(t4)./ds_total(t4)));
            if isempty(alpha_3)==1
            alpha_3=1 ;
         else
     alpha_3=min(-(x(t3)./dx_total(t3)));
            end
         if isempty(alpha_4)==1
            alpha_4=1 ;
        else
    alpha_4=min(-(s(t4)./ds_total(t4)));
         end
    alpha_pritotal=min(1,.999*alpha_3);
    alpha_dualtotal=min(1,.999* alpha_4);
    x1=x+alpha_pritotal*dx_total;
    s1=s+alpha_dualtotal*ds_total;
    lam1=lam+alpha_dualtotal*dlam_total;
    r_b=A*x-b;
    r_c=A'*lam+s-c;
    iteration=iteration+1;
    %alpha_4=1;
         x=x1;
         s=s1;
         lam=lam1;
         mu=((x'*s)/n);
end
 d_1=x<0;
 x(d_1)=.1;
 d_2=s<0;
s(d_2)=.1;
Time=toc;
%% Table part
ITERATION=iteration;
% Primal=norm(A*x-b);
% Dual=norm(A'*lam+s-c);
gap=x'*s;
MINS=min(s);
MINX=min(x);
RC=norm(r_c,inf);
RB=norm(r_b,inf);
TIME=Time;
Method={'MPC'};
Answer=table(n,m,TIME,ITERATION,gap ,RC,RB,MINX,MINS,'RowNames',Method)
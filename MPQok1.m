clear all
close all
clc
format short g
warning off
%% PRIMAL problem
n=input('n = ');
m=input('m = ');
%   m=100;
%   n=200;
% load afiro
%  load lpi_cplex2
%  A=Problem.A;
%  [m,n]=size(A);
 Q=sprandsym(n,.01);
%  b=Problem.b;
 A=5*(rand(m,n)-rand(m,n));
b=A*rand(n,1);
c=rand(n,1);
lam=(A*A')\(A*c);
x=A'*((A*A')\b);
s=c-A'*lam+Q*x;
s4=s;
r_b=A*x-b;
r_c=A'*lam+s-c-Q*x;
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
    delta_total1=([-Q A' eye(n) ;...
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
    delta_total2=([-Q A' eye(n) ;...
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
    r_c=A'*lam+s-c-Q*x;
    iteration=iteration+1;
    %alpha_4=1;
         x=x1;
         s=s1;
         lam=lam1;
         mu=(x'*s)/n;
end
 d_1=x<0;
 x(d_1)=.1;
 d_2=s<0;
s(d_2)=.1;
Time=toc;
%% Quadprog part
tic
% lb=zeros(n,1);
% ub=inf(n,1);
% H=[Q c'];
% x4=quadprog(H,A,b,[],[],lb,ub);
  [x4,~,~,output,~] = quadprog(-Q,c,A,b,[],[],[],[],[])
  timeqp=toc;
  T=['quadprog time is :',num2str(timeqp)];
  disp(T);
%% Table part
ITERATION=iteration;
%  Primal=norm(A*x-b);
%  Dual=norm(A'*lam+s-c);
 gap=(x'*s/n);
 MINS=min(s);
 MINS4=min(s4);
 MINX=min(x);
 MINX4=min(x4);
 RC=norm(r_c,inf);
 RB=norm(r_b,inf);
 TIME=Time;
% Answer=table(n,m,TIME,ITERATION,gap ,RC,RB,MINX,MINS,'RowNames',Method)
RowN={'MPQ';'quadprog'};
M(1,1)=MINX;  M(1,2)=MINS;  M(1,3)=RB;       M(1,4)=RC; 
M(1,5)= gap;       M(1,6)=TIME;
%%%%%%%%  x4=eshtebah hesab mikonad quadprog dorostash kon???????????  %%%%%%%%%%
x4=zeros(n,1);M(2,1)=norm(x4,'inf');%M(2,1)=MINX4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M(2,2)=MINS4;
M(2,3)=norm(A*x4-b);
M(2,4)=RC; 
M(2,5)=(x4'*s4/n);
M(2,6)=timeqp;
MINx=M(1:2,1);MINs=M(1:2,2);Rb=M(1:2,3); Rc=M(1:2,4);  gap=M(1:2,5); Time=M(1:2,6);
answer=table(MINx,MINs,Rb,Rc,gap,Time,'RowNames',RowN)


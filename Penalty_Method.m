clc
clear all
close all;
warning off
%% part1
%min c'x
%s.t Ax=b
%%%%%%%%%%%%%%%%%%%%
%% Input part
 m=input(' m=');
 n=input(' n=');
 d=input('d=');
%m=100;n=100;d=1;
clc
%% Lpg part
pl=inline('(abs(x)+x)/2');
A=sprand(m,n,d);
x=rand(n,1);
u=abs(spdiags(sign(pl(rand(m,1)-rand(m,1))),0,m,m)*2*((rand(m,1)-rand(m,1))));
b=A*x;
c=A'*u;%+2*spdiags((ones(n,1)-sign(x)),0,n,n)*ones(n,1);
Norm_plx=norm(pl(-x));
NormAx_b=norm((A*x-b));
bu_cx= b'*u-c'*x;
Norm_Au_c=norm(pl(A'*u-c));
Method={'LPg'};
format short
Ans=table(Norm_plx,NormAx_b,bu_cx,Norm_Au_c,'RowNames', Method);
disp(Ans);
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

%% Penalty method
K=10e8;
%% newtoon
x_new1=x;
g=c+2*K*A'*(A*x_new1-b);%gradian0
H=2*K*(A'*A);
maxiter=10;
tol=1e-10;
itern=0;
P=2*K*A';
tic
while ((norm(g,inf)>=tol)&& maxiter>=itern)
    x_new2=x_new1-(H\g);
    g=c+P*(A*x_new1-b);
    itern=itern+1;
end
timenew=toc;
fnew=c'*x_new2+K*norm(A*x_new2-b);
M2=[fnew itern-1 norm(A*x_new2-b ) timenew];
%% conjugate gradiant method
g=c+2*K*A'*(A*x-b);%gradian0
d=-g;%d0
Q=2*K*(A'*A);
maxiter=n-1;
tol=1e-15;
iter=0;
x_cg=x;
tic
while ((norm(g,inf)>=tol)&& maxiter>=iter)
    alpha=-(g'*d)/(d'*(Q*d));
    x_cg=x_cg+alpha*d;
    g=g+alpha*Q*d;%gk+1
    beta=(g'*(Q*d))/(d'*(Q*d));
    d=-g+beta*d;
    iter=iter+1;
end
timecg=toc;
itercg=iter;
%M3=[ct*x iteration norm((Ax-b)) time]
M3=[c'*x_cg itercg norm(A*x_cg-b) timecg ];
%% linprog
tic
[x_lin,fval,~,output.iterations,~]=linprog(c,[],[],A,b);
timee=toc;
lin_detail=output.iterations;
%% table part
%M2=[ct*x iteration norm((Ax-b)) time]
ctx=[M2(1,1);M3(1,1)];
iterations=[M2(1,2);M3(1,2)];
NormAx_b=[M2(1,3);M3(1,3)];
time=[M2(1,4);M3(1,4)];
Method={'Newton';'Conj_grad'};
Answer=table(ctx,iterations,NormAx_b,time,'RowNames',Method);
disp(Answer);
%% linprog output
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
f_linprog=['f_linprog :',num2str(c'*x_lin)];
disp(f_linprog);
T1=['linprog time is :',num2str(timee)];
disp(T1);
disp('For linprog details you can type lin_detail:')
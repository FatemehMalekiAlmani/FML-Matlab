clear all
close all
clc
%==============Piecewise Linear Finite Element Approximations=================
%===============================Ac=f_h================================
%% 
% generating problem
n=input('n = ');
%n=50;
h=1/(n+1);
f=.5; 
%-----------------------------
d0=rand(1,n-1);
dn0=(1/h)+(h/3);
aii0=(2/h)+(2*h/3);
aij0=(h/6)-(1/h);
D0=[d0,dn0];
G0=diag(aii0*D0);
E0=diag(aij0*D0(1:n-1),1);
P0=diag(aij0*D0(1:n-1),-1);
A0=G0+E0+P0;%_____________________________
%C=5*(rand(n,1)-rand(n,1));
Fh=f*h*rand(1,n-1);
Fhn=(f*h/2);
FT=[Fh Fhn];
C22=FT/A0;
%FT=A0*C;
%%
%=====================practical situation==================
tic
d=ones(1,n-1);
dn=(1/h)+(h/3);
aii=(2/h)+(2*h/3);
aij=(h/6)-(1/h);
D=[d,dn];
G=diag(aii*D);
E=diag(aij*D(1:n-1),1);
P=diag(aij*D(1:n-1),-1);
A=G+E+P;
Fh=f*h*ones(1,n-1);
Fhn=(f*h/2);
FT=[Fh Fhn];
C=FT/A;
T=toc;
%%
% Table part
C1=norm(C,inf);
C2=norm(C22,inf);
gap=abs(C1-C2);
 Time=T;
 Method={'PLFEA';};
 Answer=table(n,C1,C2,gap,h,Time,'RowNames',Method)


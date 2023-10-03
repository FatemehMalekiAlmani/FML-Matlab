clear all
close all
clc
%==============Piecewise Linear Finite Element Approximations=================
%% 
%=====================practical situations==================
n=input("n=");
h=1/(n+1);
f=.5;
tic
d=ones(1,n-1);
dn0=(1/h)+(h/3);
aii0=(2/h)+(2*h/3);
aij0=(h/6)-(1/h);
D0=[d,dn0];
G0=diag(aii0*D0);
E0=diag(aij0*D0(1:n-1),1);
P0=diag(aij0*D0(1:n-1),-1);
A0=G0+E0+P0;
Fh0=f*h*ones(1,n-1);
Fhn0=(f*h/2);
FT0=[Fh0 Fhn0];
%Fh=f*(2-h)*ones(n,1);
C0=FT0/A0;
T0=toc;
%%
%=======================new situation=====================
tic
dn=(h/3)-1;
aii=(4/h)+(2*h/3)-2;
aij0=(h/6)-1;
D=[d,dn];
G=diag(aii*D);
E=diag(aij0*D(1:n-1),1);
P=diag(aij0*D(1:n-1),-1);
A=G+E+P;
Fh=f*(2-h)*ones(n,1);
C=Fh'/A;
T=toc;
%%
% Table part
C=norm(C0);
Cnew=norm(C);
gap=abs(C-Cnew);
Time0=T0;
Timenew=T;
Method={'PLFEA';};
Answer=table(n,C,Cnew,gap,h,Time0,Timenew,'RowNames',Method)


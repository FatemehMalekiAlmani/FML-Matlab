clear all
close all
clc
%%
%              (4) ========= -u''+u=x^2  ==============
%                  ========= u(0)=1      ==============
%                  ========= u(1)=2      ==============
%                  ========= 0<x<1       ==============
%%
syms u(x);
du=diff(u);
ode=diff(u,x,2)==-x^2+u;
cond1=u(0)==1;
cond2=u(1)==2;
cond=[cond1 cond2];
usol(x)=dsolve(ode,cond);
usol=simplify(usol);
u0=(-1/(exp(1)+1))*exp(x)-(exp(1)/(exp(1)+1))*exp(-x)+x^2+2;
%% 
%```````````````````````plot part````````````````
fplot(x,u0,[0,1],'--or');
hold on
fplot(x,usol,[0,1],'b','Linewidth',1.25);
grid on
clear all
close all
clc
%%
%              (2) ========= -u''+u=sin(x)+2  ==============
%                  ========= u(0)=0      ==============
%                  ========= u(1)=0      ==============
%                  ========= 0<x<1       ==============
%%
syms u(x);
du=diff(u);
ode=-diff(u,x,2)+u==sin(x)+2;
cond1=u(0)==0;
cond2=u(1)==0;
cond=[cond1 cond2];
usol(x)=dsolve(ode,cond);
usol=simplify(usol);
u0=(-2-((2*exp(1)-2-.5*sin(1))/((1/exp(1))-exp(1))))*exp(x)+((2*exp(1)-2-.5*sin(1))/((1/exp(1))-exp(1)))*exp(-x)+.5*sin(x)+2;
%% 
%```````````````````````plot part````````````````
fplot(x,u0,[0,1],'--or');
hold on
fplot(x,usol,[0,1],'b','Linewidth',1.25);
grid on
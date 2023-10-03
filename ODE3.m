clear all
close all
clc
%%
%              (3) ========= u''+u=x^2  ==============
%                  ========= u(0)=1      ==============
%                  ========= u'(1)=1      ==============
%                  ========= 0<x<1       ==============
%%
syms u(x);
du=diff(u);
ode=diff(u,x,2)==x^2-u;
cond1=u(0)==1;
cond2=du(1)==1;
cond=[cond1 cond2];
usol(x)=dsolve(ode,cond);
usol=simplify(usol);
u0=3*cos(x)+((-1+3*sin(1))/(cos(1)))*sin(x)+x^2-2;
%% 
%```````````````````````plot part````````````````
fplot(x,u0,[0,1],'--or');
hold on
fplot(x,usol,[0,1],'b','Linewidth',1.25);
grid on
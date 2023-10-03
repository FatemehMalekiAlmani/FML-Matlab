clc
clear all
close all
m=800;n=900;
A=5*(rand(m,n)-rand(m,n));
x0=5*(rand(n,1)-rand(n,1));
b=(A*x0)+ .05*rand(m,1);
x=rand(n,1);
for i=(1:1:m)
a=A(i,:);
if (a*x-b)>0
    H1=a*x-b;
elseif (a*x-b)<0
        H2=zeros(m,n);
elseif (a*x-b)==0
            H3=1e-2*ones(m,n);   
end
end
H1
H2
H3
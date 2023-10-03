function [ f1 ] = nonlin( A,x0 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
A=5*(rand(m,n)-rand(m,n));
x0=5*(rand(n,1)-rand(n,1));
b_1=A*x0;
b=b_1+ones(m,1);
f1=A*x0-b;
f=find(f1<0);
f==0;
f1=f;
end


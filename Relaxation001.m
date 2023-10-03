clc
clear
close all

n=1;
k=10;
c_beta=rand(1,n);
mu=rand(1,n);
e=ones(n,1);
Q=rand(n);

x(1)=rand(n,1);
x(2)=rand(n,1);
x(3)=rand(n,1);
x(4)=rand(n,1);

X=[x(1) x(2) x(3) x(4)]' ;

fun=@(X) c_beta*x(4)-mu'*x(1);

A=[zeros(1,n) -e'*x(2) zeros(1,n) zeros(1,1)];

b=k-n;

Aeq=[e'*x(1) zeros(1,n) zeros(1,n) zeros(1,1) ; ...
    sqrt(Q)*x(1) zeros(n,n) -x(3) zeros(n,1)];

beq=[1 zeros(n,1)]';

lb=[zeros(n,1) zeros(n,1) [] zeros(1,1) ]';
ub=[[] ones(n,1) [] []]';



X1=fmincon(fun,X,A,b, Aeq ,beq,lb,ub,@NLC)














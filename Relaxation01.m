clc
clear
close all

k=10;
n=4;

X=rand(1,4*n);

Q=rand(n);

e=ones(n,1);


mu=rand(1,n);

x=X(1:n);
y=X(n+1:2*n);
w=X(2*n+1:3*n);
v=X(3*n+1:4*n);

c_beta=rand(1,n);


fun=@(X) c_beta*v'-mu*x';
A=[zeros(1,n) -e' zeros(1,n) zeros(1,n)];
b=k-n;

Aeq=[e' zeros(1,n) zeros(1,n) zeros(1,n) ; ...
    sqrt(Q) zeros(n,n) -ones(n,n) zeros(n,n)];

beq=[1 ; zeros(n,1)];
lb=[zeros(n,1); zeros(n,1); [] ;zeros(1,1) ];
ub=[[]; ones(n,1); []; []];


X1=fmincon(fun,X,A,b, Aeq ,beq,lb,ub,@NLP)









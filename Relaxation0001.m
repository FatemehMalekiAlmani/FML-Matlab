clc
clear
close all
warning off
%% OPTIMAL VALUE FOR RELAXATION PROBLEM
k=10;
n=1;

X=rand(1,4*n); %initial point
Q=rand(n);
e=ones(n,1);
mu=rand(n,1);

x=X(1:n);
y=X(n+1:2*n);
w=X(2*n+1:3*n);
v=X(3*n+1:4*n);

alpha=1;
c_beta=rand(1,n);
t=0;

tol1=1e-6;
tol2=1e-8;
%% calculation linprog for d

f=[-mu zeros(1,n) zeros(1,n) c_beta];

A=[zeros(1,n) zeros(1,n) 2*w -2*v ;...
    zeros(1,n) -e' zeros(1,n) zeros(1,n)];

b=[0 n]';

Aeq=[y x zeros(1,n) zeros(1,n) ;...
    e' zeros(1,n) zeros(1,n) zeros(1,n);...
    sqrt(Q) zeros(n) ones(n) zeros(n)];

beq=[0 0 zeros(1,n) ]';

ub=[[] e' [] []]';
lb=[-e' [] [] zeros(1,n)]';

d=linprog(f,A,b,Aeq,beq,lb,ub,X);


iter=0;
maxiter=100;
tic
while (norm(x.*y,inf)<= tol1 || t<tol2 && maxiter>=iter )
    
    if d==zeros(4*n,1);
        disp('this x is your answer')
        break
    else
          
        Xt=X'+alpha*d;
        x=Xt(1:n);
        y=Xt(n+1:2*n);
        w=Xt(2*n+1:3*n);
        v=Xt(3*n+1:4*n);
        t=.01*t;
        
        A=[zeros(1,n) zeros(1,n) 2*w' -2*v' ;...
            zeros(1,n) -e' zeros(1,n) zeros(1,n)];
        
        Aeq=[y' x' zeros(1,n) zeros(1,n) ;...
            e' zeros(1,n) zeros(1,n) zeros(1,n);...
            sqrt(Q) zeros(n) ones(n) zeros(n)];
        
        d=linprog(f,A,b,Aeq,beq,lb,ub,Xt);
        
        iter=1+iter;
        
    end
    
end
Time=toc;
yep=norm(x.*y,inf)

f*d

















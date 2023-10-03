clc
clear
close all
warning off


X=rand(1,4*n)-rand(1,4*n);


e=ones(n,1);

x1=X(1:n);
y=X(n+1:2*n);
w=X(2*n+1:3*n);
v=X(3*n+1:4*n);

lamda=rand(1,3*n)-rand(1,3*n);
lambda1=lambda(1:n);
lambda2=lambda(n+1:2*n);
lambda3=lambdaX(2*n+1:3*n);


mu=rand(1,5*n)-rand(1,5*n);
mu1=mu(1:n);
mu2=mu(n+1:2*n);
mu3=mu(2*n+1:3*n);
mu4=mu(3*n+1:4*n);
mu5=mu(4*n+1:5*n);

iter=0;

maxiter=10;
tol=1e-8;
c=1e-4;
k=1;
while (norm(x-x1,inf)<tol && maxiter>=iter)
    h1=x1.*y;
    h2=e'*x1-1;
    h3=w-sqrt(Q)*x;
    
    H=[h1 h2 h3];
    
    gh1=[y x zeros(n,1) zeros(n,1)];
    gh2=[e' zeros(n,1) zeros(n,1) zeros(n,1)];
    gh3=[-sqrt(Q) zeros(n,n) ones(n,n) zeros(n,n) ];
    
    GH=[gh1 gh2 gh3];
    
    %       gL=[-mu+lamda1'*y
    g1=w'*w-v.*2;
    g2=-e'*y-k+n;
    g3=-y;
    g4=-v;
    g5=y;
    
    G=[g1 g2 g3 g4 g5];
    
    gg1=[zeros(1,n) zeros(1,n) 2*w -2*v] ;
    gg2=[zeros(1,n) -e' zeros(1,n) zeros(1,n)];
    gg3=[zeros(1,n) -ones(1,n) zeros(1,n) zeros(1,n)];
    gg4=[zeros(1,n) zeros(1,n) zeros(1,n) -ones(1,n)];
    gg5=[zeros(1,n) ones(1,n) zeros(1,n) zeros(1,n)];
    
    
    GG=[gg1 gg2 gg3 gg4 gg5];
    
    grad2_Lxx=[zeros(1,n) lambda1 zeros(1,n) zeros(1,n);...
    lambda1 zeros(1,n) zeros(1,n) zeros(1,n);...
    zeros(1,n) zeros(1,n) 2*mu1 zeros(1,n);...
    zeros(1,n) zeros(1,n) zeros(1,n) -2*mu1];


    grad2_LTOTAL=[grad2_Lxx GH GG;...
         GH zeros(3) zeros(3,5);...
         GG zeros(5,3) zeros(5)];

    f=[-mu0 zeros(1,n) zeros(1,n) C_beta];
    F=f*X;
    Aeq=[H G];
    beq=
    
d=fmincon(F,X,[],[],Aeq,






    iter=iter+1;
end











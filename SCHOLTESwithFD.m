clear all
close all
clc
format short g
warning off
%%
% m=20;

n=10;
beta=.9;

% n=input('n=');%113s
% beta=input('beta=');

k=10;
% X=rand(1,4*n);
X=[ones(1,n) ones(1,n) ones(1,n) ones(1,n)];%perfect

x1=X(1:n);
y=X(n+1:2*n);
w=X(2*n+1:3*n);
v=X(3*n+1:4*n);

t=1;
e=ones(n,1);
Q=rand(n);
tol=1e-6;
% beta=input('beta=');
mu=rand(n,1);
%lam=1e-10;
Q=rand(n);
% beta=input('beta=');
x=rand(n,1);
%  Q=Q'*Q+(lam*eye(n));
%% RISK MEASER FUNCTION

VaR=VaRfun( mu,Q,beta,x );

CVaR=CVaRfun( mu,Q,beta,x );

RVaR=RVAR( mu,Q,beta,x );

RCVaR=RCVARfun( mu,Q,beta,x );


%% TABLE PART
VAR=VaR;
CVAR=CVaR;
RVAR=RVaR;
RCVAR=RCVaR;
Method={'RISK(c_beta)';};
Answer=table(VAR,CVAR,RVAR,RCVAR,'RowNames',Method)

c_beta=[VAR CVAR RVAR RCVAR];

C=[-mu' zeros(1,n) zeros(1,n) c_beta zeros(1,n-4)]';
C_beta=[VAR CVAR RVAR RCVAR zeros(1,n-4)];
%% ~~~~~~~~~~~~~~Find initial feasible point(faaz1)~~~~~~~~~~~~~~~~~~~
Aeq=[sqrt(Q) zeros(n) ones(n) zeros(n);...
    e' zeros(1,n) zeros(1,n) zeros(1,n)];

A=[zeros(1,n) zeros(1,n) w -v;...
    zeros(1,n) x1 zeros(1,n) zeros(1,n);...
    zeros(1,n) -e' zeros(1,n) zeros(1,n)];

beq=[zeros(1,n) 1]';
b=[0 t k-n]';

lb=[[] 0 [] 0];
ub=[[] e [] []];

Mu=rand(1,n);
f=[-Mu zeros(1,n) zeros(1,n) C_beta];

X=linprog(f,A,b,Aeq,beq,lb,ub);
X=X/norm(X);
maxiter=30;
iteration=0;


%% ~~~~~~~~~~~~~~~~~find direction with linprog (Faaz2)~~~~~~~~~~~~~~~~~~~~~~~
FF=-Mu*x1'+C_beta*v';
tic
while f*X <= tol && iteration<maxiter
    
    beqd=zeros(n+1,1);
    bd=zeros(3,1);
    
    ubd=[[] 0 [] []];
    d=linprog(f,A,bd,Aeq,beqd,lb,ubd);
    D=d/norm(d);
    
    if norm(D)<tol
        
        disp('this X is your answer')
        break
    else
        
        if isempty(D)==1
            disp('this X is your answer')
            break
        else
            %% ~~~~~~~~~~~~~~~calculate alpha (faaz3)~~~~~~~~~~~~~~~~~~~~~~
            alpha0=1;
            F=@(X) -Mu*x1'+C_beta*v';
            lb=0;%zeros(4*n,1);
            ub=1;%ones(4*n,1);
            
            alpha=fmincon(F,alpha0,[],[],[],[],lb,ub)
            %% keep on while
            X=X+alpha*D;
            x1=X(1:n);
            y=X(n+1:2*n);
            w=X(2*n+1:3*n);
            v=X(3*n+1:4*n);
            X=X/norm(X);
            FF=-Mu*x1+C_beta*v;
            iteration=iteration+1;
        end
    end
end
Time=toc;

[norm(X)
norm(D)
alpha
Time
FF]










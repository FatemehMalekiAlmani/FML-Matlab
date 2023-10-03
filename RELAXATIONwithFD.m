clc
clear
close all
warning off
%%
n=10;
% n=input('n=');%124s
%   load EuroStoxx50
%   Q=problem.Assets;
%  [n,n]=size(Q);
%   mu=problem.Index;

mu=rand(n,1);
Q=rand(n);
lam=1e-10;
beta=.9;
% beta=input('beta=');
x=rand(n,1);
Q=Q'*Q+(lam*eye(n));
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

%%
c_beta=[VAR CVAR RVAR RCVAR];

%% OPTIMAL VALUE FOR RELAXATION PROBLEM
k=10;


% X=rand(1,4*n)-rand(1,4*n); %initial point ok
% Q=rand(n);
e=ones(n,1);
Mu=rand(1,n);
%  X=[zeros(1,n) e' zeros(1,n) zeros(1,n)]; % d=[] & no answer ~~~~~~~~~~~~~~~~~~~~~~~~~~
%  X=[zeros(1,n) e' ones(1,n) ones(1,n)]; perfect



X0=[zeros(1,n) ones(1,n) ones(1,n) ones(1,n)];%perfect


x1=X0(1:n);
y=X0(n+1:2*n);
w=X0(2*n+1:3*n);
v=X0(3*n+1:4*n);

alpha=1;
C_beta=[VAR CVAR RVAR RCVAR zeros(1,n-4)];
t=0;

tol1=1e-6;
tol2=1e-8;
%% ~~~~~~~~~~~~~calculation linprog for initial point X~~~~~~~~~~~~~~

f=[-Mu zeros(1,n) zeros(1,n) C_beta];

A1=[x1 zeros(1,n) zeros(1,n) zeros(1,n);...
    sqrt(Q) zeros(n) ones(n) zeros(n);...
    e' zeros(1,n) zeros(1,n) zeros(1,n)];

A2=[zeros(1,n) zeros(1,n) w v;...
    zeros(1,n) e' zeros(1,n) zeros(1,n);...
    zeros(1,n) -e' zeros(1,n) zeros(1,n);...
    zeros(1,n) zeros(1,n) zeros(1,n) -e';...
    zeros(1,n) e' zeros(1,n) zeros(1,n)];

b1=[zeros(n+1,1);1];
b2=[0 k-n 0 0 1]';

% A=[A1;A2];
% b=[b1;b2];

X=linprog(f,A1,b1,A2,b2,[],[],X0);

%% ~~~~~~~~~~~ calculate linprog for direction d  ~~~~~~~~~~~~~~~~~~~~~~~~~~~
bd1=zeros(n+2,1);
bd2=zeros(5,1);

ub=[e' e' e' e']';

lb=[-e' -e' -e' -e']';
d=linprog(f,A1,bd1,A2,bd2,lb,ub,X);

D=d/norm(d);
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
g=f;

X1=X/norm(X);

X2=X1+D;
%  F1=f*X1';
% F1=f*D;
F1=f*X2;
% b=A*d;
% beq=Aeq*d;
tol=1e-8;
iter=0;
maxiter=10;
% t=(f*X')./norm(f.*X,inf);
tic
while ((f*X1<=tol)&& maxiter>=iter)
    
    if f*X1==0
        disp('this X is your answer')
        break
    else
        %         alpha=.99*alpha;
        % alpha=-(f*d)/(d'*(Q*d));
        AD=A2*D;
        d_1=AD<=0;
        AD(d_1)=.1;
        bbar=b2-A2*X1;
        alpha=min(bbar./AD);
        if isempty(alpha)==1
            alpha=1 ;
        else
            alpha=min(bbar./ AD);
        end
        
        X=X1+alpha*D;
        
        X1=X/norm(X);
        
        
        x1=X1(1:n);
        y=X1(n+1:2*n);
        w=X1(2*n+1:3*n);
        v=X1(3*n+1:4*n);
        
        
        
        A1=[x1 zeros(1,n) zeros(1,n) zeros(1,n);...
            sqrt(Q) zeros(n) ones(n) zeros(n);...
            e' zeros(1,n) zeros(1,n) zeros(1,n)];
        
        A2=[zeros(1,n) zeros(1,n) w v;...
            zeros(1,n) e' zeros(1,n) zeros(1,n);...
            zeros(1,n) -e' zeros(1,n) zeros(1,n);...
            zeros(1,n) zeros(1,n) zeros(1,n) -e';...
            zeros(1,n) e' zeros(1,n) zeros(1,n)];
        
        b1=zeros(n+2,1);
        b2=zeros(5,1);
        
        
        
        d=linprog(f,A1,b1,A2,b2,[],[],X1);
        
        D=d/norm(d);
        
        
        d1=D(1:n);
        d2=D(n+1:2*n);
        d3=D(2*n+1:3*n);
        d4=D(3*n+1:4*n);
        
 %% ~~~~~~~~~~~~~~~~~~~~~~~~~F MIN CON ~~~~~~~~~~~~~~~~~~~~~~~~~
        
        fun=-Mu*x1+C_beta*v;
        
        A=[A1;A2];
        b=[b1;b2];
        
        
        
        xx= fmincon(fun,x0,A1,b1,A2,b2,lb,ub)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        iter=1+iter;
        
    end
end

Time=toc;

f*X1'


%% TABLE PART
% XOY=norm(g,inf);
% T=Time;
% Ft=abs((F1-F2)/F2);
% XoY=XOY;
% TIME=T;
% gap=Ft;
% Mx=min(x1);
% My=min(y);
% Fin=F1;
% FT=F2;
% Method={'RELAXATION';};
% Answer=table(XoY,TIME,gap,Fin,Mx,My,FT,'RowNames',Method)

%% PLOT
% line(iter,yep)
% figure()
% yep=norm(x1.*y,inf);
% Y=iter;
% plot(X,FT)



clc
clear
close all
warning off
%%

n=10;
mu=rand(n,1);
lam=1e-10;
Q=rand(n);
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


 X=rand(1,4*n)-rand(1,4*n); %initial point
% Q=rand(n);
e=ones(n,1);
Mu=rand(1,n);

%  X=[zeros(1,n) e' zeros(1,n) zeros(1,n)];

x1=X(1:n);
y=X(n+1:2*n);
w=X(2*n+1:3*n);
v=X(3*n+1:4*n);

alpha=1;
C_beta=[VAR CVAR RVAR RCVAR zeros(1,n-4)];
t=0;

tol1=1e-6;
tol2=1e-8;
%% calculation linprog for d

f=[-Mu zeros(1,n) zeros(1,n) C_beta];

A=[zeros(1,n) zeros(1,n) 2*w -2*v ;...
    zeros(1,n) -e' zeros(1,n) zeros(1,n)];

 b=[0 k]';
%  b=[v*v'-w*w' e'*y'+k-n]';
 
 
Aeq=[y x1 zeros(1,n) zeros(1,n) ;...
    e' zeros(1,n) zeros(1,n) zeros(1,n);...
    -sqrt(Q) zeros(n) ones(n) zeros(n)];

 beq=[0 1 zeros(1,n) ]';
% beq=[-x1*y' 1-e'*x1' (sqrt(Q)*x1'-w')']';


 ub=[[] e' [] []]';
% ub=[[] -y [] []]';

%  lb=[-e' [] [] zeros(1,n)]';
 lb=[[] -e' [] zeros(1,n)]';


d=linprog(f,A,b,Aeq,beq,lb,ub,X);

F1=f*X';


% b=A*d;
% beq=Aeq*d;


iter=0;
maxiter=10;
tic
while (norm(x1.*y,inf)<= tol1 || t<tol2) && maxiter>=iter
    
    if d==zeros(4*n,1);
        disp('this X is your answer')
        break
    else
        
        if isempty(d)==1
            disp('this X is your answer')
            break
        else
            alpha=.99*alpha;
            
            X=X+alpha*d';
            
            x1=X(1:n);
            y=X(n+1:2*n);
            w=X(2*n+1:3*n);
            v=X(3*n+1:4*n);
            t=.01*t;
   
            
            
            Aeq=[zeros(1,n) zeros(1,n) 2*w -2*v ;...
                zeros(1,n) -e' zeros(1,n) zeros(1,n);...
                y x1 zeros(1,n) zeros(1,n) ;...
                e' zeros(1,n) zeros(1,n) zeros(1,n);...
                -sqrt(Q) zeros(n) ones(n) zeros(n)];

             grad_f=[-Mu zeros(1,n) zeros(1,n) C_beta];

            % grad2_L=

            f=grad_f*d%+.5*d'*grad2_L*d;

            beq=[v*v'-w*w' e'*y'+k-n -x1*y' 1-e'*x1' (sqrt(Q)*x1'-w')']';
            lb=[-e' [] [] zeros(1,n)]';
            ub=[[] e' [] []]'; 

            d=fmincon(f,[],[],Aeq,beq,lb,ub,X);
            
            d1=d(1:n);
            d2=d(n+1:2*n);
            d3=d(2*n+1:3*n);
            d4=d(3*n+1:4*n);
            
%             d=d.*d1;
            
%             F2=C_beta*v'-mu'*x1'-mu'*d1+C_beta*d4;
             F2=f*X';
            iter=1+iter;
            
        end
    end
    
end
Time=toc;


%% TABLE PART
yep=norm(x1.*y,inf);
T=Time;
Ft=(F1-F2)/F2;
YEP=yep;
TIME=T;
FT=Ft;
Mx=min(x1);
My=min(y);
Fin=F1;
FT=F2;
Method={'RELAXATION';};
Answer=table(YEP,TIME,FT,Fin,Mx,My,FT,'RowNames',Method)



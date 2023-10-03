function [ Time, gap, iter ] = relaxtion_fun ( n, Q, VaR, CVaR, RVaR, RCVaR  ) 
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%% TABLE PART
VAR=VaR;
CVAR=CVaR;
RVAR=RVaR;
RCVAR=RCVaR;
Method={'RISK(c_beta)';};
Answer=table(VAR,CVAR,RVAR,RCVAR,'RowNames',Method)

%%
c_beta=[VAR CVAR RVAR RCVAR];
% `````````````````````````````````````````````for o=1:10
%% OPTIMAL VALUE FOR RELAXATION PROBLEM
k=10;
% X=rand(1,4*n)-rand(1,4*n); %initial point ok
% Q=rand(n);
e=ones(n,1);
Mu=rand(1,n);
%  X=[zeros(1,n) e' zeros(1,n) zeros(1,n)]; % d=[] & no answer ~~~~~~~~~~~~~~~~~~~~~~~~~~
%  X=[zeros(1,n) e' ones(1,n) ones(1,n)]; perfect

X=[zeros(1,n) ones(1,n) ones(1,n) ones(1,n)];%perfect


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

D=d/norm(d);

X1=X/norm(X);

X2=X1+D';
%  F1=f*X1';
% F1=f*D;
F1=f*X2';
% b=A*d;
% beq=Aeq*d;

iter=0;
maxiter=10;
tic
while (norm(x1.*y,inf)>= tol1 || t<tol2) && maxiter>=iter
    
    if norm(D)<tol1
        %         d==zeros(4*n,1);
        disp('this X is your answer')
        break
    else
        
        if isempty(D)==1
            disp('this X is your answer')
            break
        else
            alpha=.99*alpha;
            % alpha=-(f*d)/(d'*(Q*d));
            
            X=X+alpha*D';
            
            
            X1=X/norm(X);
            
            
            x1=X1(1:n);
            y=X1(n+1:2*n);
            w=X1(2*n+1:3*n);
            v=X1(3*n+1:4*n);
            t=.01*t;
            
            A=[zeros(1,n) zeros(1,n) 2*w -2*v ;...
                zeros(1,n) -e' zeros(1,n) zeros(1,n)];
            
            b=[v*v'-w*w' e'*y'+k-n]';
            
            
            Aeq=[y x1 zeros(1,n) zeros(1,n) ;...
                e' zeros(1,n) zeros(1,n) zeros(1,n);...
                -sqrt(Q) zeros(n) ones(n) zeros(n)];
            
            
            beq=[-x1*y' 1-e'*x1' (sqrt(Q)*x1'-w')']';
            
            
            d=linprog(f,A,b,Aeq,beq,lb,ub,X1);
            D=d/norm(d);
            d1=D(1:n);
            d2=D(n+1:2*n);
            d3=D(2*n+1:3*n);
            d4=D(3*n+1:4*n);
            
            %             d=d.*d1;
            
            %             F2=C_beta*v'-mu'*x1'-mu'*d1+C_beta*d4;
            F2=f*X1';
            iter=1+iter;
            
        end
    end
end

%     Time(o)=toc;
Ft=abs((F1-F2)/F2);
% FV(o)=Ft;
% [o Time Ft]
% plot(Time,FV,'o-')
% end
Time=toc;
XOY=norm(x1.*y,inf);
gap=XOY;

end


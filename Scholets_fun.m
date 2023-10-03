function [ Time, gap, iter ] = Scholets_fun( n, Q, VaR, CVaR, RVaR, RCVaR  )
%UNTITLED7 Summary of this function goes here
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

%% OPTIMAL VALUE FOR RELAXATION PROBLEM
k=10;

    %% RISK MEASER FUNCTION


%   X=rand(1,4*n)-rand(1,4*n); %initial point
e=ones(n,1);
Mu=rand(1,n);
% X=[zeros(1,n) e' zeros(1,n) zeros(1,n)];%NO answer ~~~~~~~~~~~~~~~~~~~
X=ones(1,4*n);
% X=[zeros(1,n) ones(1,n) ones(1,n) ones(1,n)];%perfect

x1=X(1:n);
y=X(n+1:2*n);
w=X(2*n+1:3*n);
v=X(3*n+1:4*n);

alpha=1;
C_beta=[VAR CVAR RVAR RCVAR zeros(1,n-4)];
t=1;

tol=1e-6;
%% calculation linprog for d
T=.1*rand(n,1);
f=[-Mu zeros(1,n) zeros(1,n) C_beta];

A=[zeros(1,n) zeros(1,n) 2*w -2*v ;...
    y x1 zeros(1,n) zeros(1,n);...
    -y -x1 zeros(1,n) zeros(1,n);...
    zeros(1,n) -e' zeros(1,n) zeros(1,n)];

b=[w*w'+v*v' -x1*y'+T'*e x1*y'+T'*e k-n+e'*y']';
Aeq=[sqrt(Q) zeros(n) -ones(n) zeros(n);...
    e' zeros(1,n) zeros(1,n) zeros(1,n)];
beq=[(-sqrt(Q)*x1'+w')' 1-e'*x1']';
ub=[[] -y+e' [] []]';
lb=[[] -y [] -v]';
d=linprog(f,A,b,Aeq,beq,lb,ub,X);
D=d/norm(d);
X1=X/norm(X);
F1=f*X1';
iter=0;
maxiter=10;
tic
%``````````````````````````  for o=1:1:n
while (norm(x1.*y,inf)>= tol || t<tol) && maxiter>=iter
%     (x1*y'>= tol || t<tol) && maxiter>=iter
    if norm(D)<tol
        %         d==zeros(4*n,1);
        disp('this X is your answer')
        break
    else
        
        if isempty(D)==1
            disp('this X is your answer')
            break
        else
            alpha=.99*alpha;
            X=X+alpha*D';
            X1=X/norm(X);
            x1=X1(1:n);
            y=X1(n+1:2*n);
            w=X1(2*n+1:3*n);
            v=X1(3*n+1:4*n);
            t=.01*t;
            A=[zeros(1,n) zeros(1,n) 2*w -2*v ;...
                y x1 zeros(1,n) zeros(1,n);...
                -y -x1 zeros(1,n) zeros(1,n);...
                zeros(1,n) -e' zeros(1,n) zeros(1,n)];
            
            b=[w*w'+v*v' -x1*y'+T'*e x1*y'+T'*e k-n+e'*y']';
            
            Aeq=[sqrt(Q) zeros(n) -ones(n) zeros(n);...
                e' zeros(1,n) zeros(1,n) zeros(1,n)];
            beq=[(-sqrt(Q)*x1'+w')' 1-e'*x1']';
            ub=[[] -y+e' [] []]';
            lb=[[] -y [] -v]';
            d=linprog(f,A,b,Aeq,beq,lb,ub,X1);
            D=d/norm(d);
            d1=D(1:n);
            d2=D(n+1:2*n);
            d3=D(2*n+1:3*n);
            d4=D(3*n+1:4*n);
            F2=f*X1';
            iter=1+iter;
        end
    end
end
Time=toc;
%```````````````````````````  Time(o)=toc;
%Ft=abs((F1-F2)/F2);
%``````````````````````````` FV(o)=Ft;
%``````````````````````````` [o Time Ft];
%``````````````````````````` o=o+1;
%``````````````````````````` plot(Time,FV,'o-')
%end
%% TABLE PART
XOY=norm(x1.*y,inf);
gap=XOY;

end


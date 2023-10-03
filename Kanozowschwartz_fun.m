function [ Time, gap, iter ] = Kanozowschwartz_fun ( n, Q, VaR, CVaR, RVaR, RCVaR  ) 
%UNTITLED5 Summary of this function goes here
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


%   X=rand(1,4*n)-rand(1,4*n); %initial point
% Q=rand(n);
e=ones(n,1);
Mu=rand(1,n);
% X=[zeros(1,n) e' zeros(1,n) zeros(1,n)];
% X=[zeros(1,n) e' ones(1,n) ones(1,n)];% perfect

X=[zeros(1,n) ones(1,n) ones(1,n) ones(1,n)];


x1=X(1:n);
y=X(n+1:2*n);
w=X(2*n+1:3*n);
v=X(3*n+1:4*n);

alpha=1;
C_beta=[VAR CVAR RVAR RCVAR zeros(1,n-4)];
t=0;

tol1=1e-6;
tol2=1e-6;



%% calculation linprog for d
T=.1*rand(n,1);
if norm(x1)<= norm(y)
    B=x1'-T;
    
    f=[-Mu zeros(1,n) zeros(1,n) C_beta];
    
    A=[zeros(1,n) zeros(1,n) 2*w -2*v ;...
        zeros(1,n) -e' zeros(1,n) zeros(1,n)];
    
    b=[v*v'-w*w' e'*y'+k-n]';
    
    
    Aeq=[-sqrt(Q) zeros(n) ones(n) zeros(n);...
        e' zeros(1,n) zeros(1,n) zeros(1,n)];
    
    beq=[(sqrt(Q)*x1'-w')' 1-e'*x1']';
    
    
    ub=[-B -y'+e [] []]';
    
    lb=[[] [] [] -v]';
    
    
    d=linprog(f,A,b,Aeq,beq,lb,ub,X);
    
    
else if norm(x1)>norm(y)
        B=y-T;
        f=[-Mu zeros(1,n) zeros(1,n) C_beta];
        
        A=[zeros(1,n) zeros(1,n) 2*w -2*v ;...
            zeros(1,n) -e' zeros(1,n) zeros(1,n)];
        
        b=[v*v'-w*w' e'*y'+k-n]';
        
        
        Aeq=[-sqrt(Q) zeros(n) ones(n) zeros(n);...
            zeros(1,n) e' zeros(1,n) zeros(1,n)];
        
        beq=[(sqrt(Q)*x1'-w')' 1-e'*x1']';
        
        
        ub=[-B -y'+e [] []]';
        
        lb=[[] [] [] -v]';
        
        
        d=linprog(f,A,b,Aeq,beq,lb,ub,X);
        
        
    end
end

D=d/norm(d);

X1=X/norm(X);

X2=X1+D';
%  F1=f*X1';
% F1=f*D;
F1=f*X1';
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
            
            if norm(x1)<= norm(y)
                B=x1'-T;
                
                f=[-Mu zeros(1,n) zeros(1,n) C_beta];
                
                A=[zeros(1,n) zeros(1,n) 2*w -2*v ;...
                    zeros(1,n) -e' zeros(1,n) zeros(1,n)];
                
                b=[v*v'-w*w' e'*y'+k-n]';
                
                
                Aeq=[-sqrt(Q) zeros(n) ones(n) zeros(n);...
                    e' zeros(1,n) zeros(1,n) zeros(1,n)];
                
                beq=[(sqrt(Q)*x1'-w')' 1-e'*x1']';
                
                
                ub=[-B -y'+e [] []]';
                
                lb=[[] [] [] -v]';
                
                
                d=linprog(f,A,b,Aeq,beq,lb,ub,X1);
                
                
            else if norm(x1)>norm(y)
                    B=y-T;
                    f=[-Mu zeros(1,n) zeros(1,n) C_beta];
                    
                    A=[zeros(1,n) zeros(1,n) 2*w -2*v ;...
                        zeros(1,n) -e' zeros(1,n) zeros(1,n)];
                    
                    b=[v*v'-w*w' e'*y'+k-n]';
                    
                    
                    Aeq=[-sqrt(Q) zeros(n) ones(n) zeros(n);...
                        zeros(1,n) e' zeros(1,n) zeros(1,n)];
                    
                    beq=[(sqrt(Q)*x1'-w')' 1-e'*x1']';
                    
                    
                    ub=[-B -y'+e [] []]';
                    
                    lb=[[] [] [] -v]';
                    
                    
                    d=linprog(f,A,b,Aeq,beq,lb,ub,X1);
                    
                    
                end
            end
            D=d/norm(d);
            d1=D(1:n);
            d2=D(n+1:2*n);
            d3=D(2*n+1:3*n);
            d4=D(3*n+1:4*n);
            
            %             d=d.*d1;
            
            %             F2=C_beta*v'-mu'*x1'-mu'*d1+C_beta*d4;
            F2=-f*X1';
            iter=1+iter;
            
        end
    end
    
end
Time=toc;

XOY=norm(x1.*y,inf);
gap=XOY;

end


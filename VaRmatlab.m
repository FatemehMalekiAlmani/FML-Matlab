clc;
clear;
close all;
n=5;
mu=rand(n,1);
%lam=1e-10;
Q=rand(n);
beta=input('beta=');
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

%%
c_beta=[VAR CVAR RVAR RCVAR]



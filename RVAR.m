function [ Rvar ] = RVAR( mu,Q,beta,x )
Rvar=((2*beta-1)/(2*(sqrt(beta*(1-beta)))))*sqrt(x'*Q*x)-mu'*x;
end
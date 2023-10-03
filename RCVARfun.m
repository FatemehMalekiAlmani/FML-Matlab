function [ Rcvar ] = RCVARfun( mu,Q,beta,x )
 Rcvar=(sqrt(beta/(1-beta)))*sqrt(x'*Q*x)-mu'*x;
end
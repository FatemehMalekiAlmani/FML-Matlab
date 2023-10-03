function [ VaR ] = VaRfun( mu,Q,beta,x )
if beta==.9
    stigma=1.282;
elseif beta==.95
    stigma=1.645;
elseif beta==.99
    stigma=2.326;
end
      

VaR=stigma*(sqrt(x'*Q*x))-mu'*x;
end


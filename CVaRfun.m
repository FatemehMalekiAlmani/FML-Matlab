function [ CVaR ] = CVaRfun( mu,Q,beta,x )

if beta==.9
    B=.1;
elseif beta==.95
    B=.05;
elseif beta==.99
    B=.01;
end


if beta==.9
    stigma=1.282;
elseif beta==.95
    stigma=1.645;
elseif beta==.99
    stigma=2.326;
end


T=exp(-(stigma^2)/2);
eta=((1/sqrt(2*pi))*T)/B ;
CVaR=eta*(sqrt(x'*Q*x))-mu'*x;

end


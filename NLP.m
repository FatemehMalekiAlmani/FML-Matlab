function [c , ceq ] = NLP(x,y,w,v)
% w=sqrt(Q)*x;
% X=rand(1,4*n);
% 
%   x=X(1:n);
%  y=X(n+1:2*n);
%  w=X(2*n+1:3*n);
%  v=X(3*n+1:4*n);

  c=w'*w-v'*v;


  ceq=x'*y;

end

%  c=[zeros(1,n) zeros(1,n) w' -v']*[x;y;w;v];
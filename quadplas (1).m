function N = quadplas( A,x,b,m )
[~,n]=size(A);
N=zeros(m,n);
H=zeros(m,1);
for i=(1:1:m)
a=A(i,:);
if (a*x-b)>0
    H=a*x-b;
elseif (a*x-b)<0
        H=zeros(m,1);
elseif (a*x-b)==0
            H=100*rand(m,1);
        
end
N(i,1:end)=H(i,1:end);
% i=i+1;
end
end


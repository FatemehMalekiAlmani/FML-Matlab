function H = quadplas( A,x,b,m )
for i=(1:1:m)
a=A(i,:);
if (a*x-b)>0
    H=a*x-b;
else if (a*x-b)<0
        H=0;
    else if (a*x-b)==0
            H=1e-2;
        end
    end
% i=i+1;
end
end


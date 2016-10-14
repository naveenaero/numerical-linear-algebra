function [ Q,R,U ] = houseqr( A )
%HOUSEQR Summary of this function goes here
%   Detailed explanation goes here

[m,n]=size(A);
U = zeros(m,n);
Q = eye(m);

    function [v, beta] = house(x)
        sigma = x(2:end)'*x(2:end);
        v = [1; x(2:end)];
        if sigma==0 && x(1)>=0
            beta = 0;
        elseif sigma==0 && x(1)<0
            beta=-2;
        else
            mu = sqrt(x(1)^2 + sigma);
            if x(1)<=0
                v(1)=x(1)-mu;
            else
                v(1)=-sigma/(x(1)+mu);
            end
            beta = 2*v(1)^2/(sigma+v(1)^2);
            v = v/v(1);
        end
    end

for j = 1:n
    [v,beta] = house(A(j:m,j));
    U(j:end,j) = v;
    temp = v'*A(j:m,j:n);
    A(j:m,j:n) = A(j:m,j:n) - beta*v*temp;
    H = eye(m);
    H(j:m,j:m) = H(j:m,j:m) - beta*(v*v');
    Q = Q*H;
end
R = A;

end


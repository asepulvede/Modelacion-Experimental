function [phi, Y, param, cova] = MCO(data,na, nb, nk)
% Function for Least Squares.
% Arguments:
% - data: matrix as [t, u, y].
% - nk: delay of the system. nk = d + 1.
% Output:
% - phi: matrix with the autorregressive model in each time.
% - Y: measured observations from t = na + nk - 1 to t = N.
% - param: vector with the estimation of the parameters.

y = data(:,3);
u = data(:,2);
N = length(data);
Q = max(na, nb + nk - 1);
size_ = N - Q;
phi = zeros(size_, na + nb);
count_aux = 0;
aux = max(na, nb);
for i = 1:size_
    count = 1;
    for j = 1:na
        phi(i,count) = -y(Q - j + 1 + count_aux);
        count = count + 1;
    end 
    for j = 1:nb
        phi(i,count) = u(nb - j + 1 + count_aux);
        count = count + 1;
    end
    count_aux = count_aux + 1;
end
Y = data((Q+1):end,3);
param = inv(phi'*phi)*(phi')*Y;
V = 0.5*(Y'*Y - Y'*phi*inv(phi'*phi)*phi'*Y);
lambda = 2*V/(N-nk-na-nb);
cova = lambda*inv(phi'*phi);
end

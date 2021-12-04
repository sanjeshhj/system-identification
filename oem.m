function [ A,B,theta ] = oem( y, u, na, nb, d )
%Output Error Method

N=length(y);
y=reshape(y,N,1);
u=reshape(u,N,1);

%Initialisation of the algorithm
theta=[zeros(na,1);zeros(nb,1)];
D=1000*eye(na+nb);
y0 = zeros(N,1);
for i=max([na,nb]):N-1
    phi=[-y0(i:-1:i-na+1);u(i:-1:i-nb+1)];
    y0(i+1) = theta' * phi;
    epr = y(i+1) - y0(i+1);
    epr = epr/(1+ phi' * D * phi);
    theta = theta + (D * phi * epr);
    D = D - (D * (phi) * phi' * D) / (1 + phi' * D * phi);
end
A = [ 1 theta(1:na)'];
B = [ 0 theta(na+1:end)'];
end
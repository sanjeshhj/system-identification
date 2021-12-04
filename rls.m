function [ A,B,theta ] = rls( y, u, na, nb, d )
% Recursive Least Square Algorithm


N=length(u);
y=reshape(y,N,1);
u=reshape(u,N,1);

%Initialisation of the algorithm
theta=[zeros(na,1);zeros(nb,1)];
D=100000*eye(na+nb);


for i=max(na,nb):N-1
    phi=[ -y(i: -1: i-na+1) ; u(i : -1 : i -nb+1)];
    yh0=theta'*phi;
    e0=y(i+1)-yh0;
    e=e0/(1+phi'*D*phi);
    theta=theta+D*phi/(1+phi'*D*phi)*e;
    D=D - (D*(phi*phi')*D/(1+phi'*D*phi));
end

A = [1 theta(1:na)'];
B = [0 theta(na+1:end)'];

end

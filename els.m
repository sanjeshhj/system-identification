function [ A,B,C,theta,e] = els( y, u, na, nb, nc, d )
%  Extended Least Square Algorithm

N=length(u);
y=reshape(y,N,1);
u=reshape(u,N,1);
e=zeros(N,1);

%Initialisation of the algorithm
theta=[zeros(na,1);zeros(nb,1);zeros(nc,1)];
D=10000000*eye(na+nb+nc);


for i=max(max(na,nb),nc):N-1
    phi=[ -y(i: -1: i-na+1) ; u(i : -1 : i-nb+1); e(i : -1 : i-nc+1)];
    yh0=theta'*phi;
    e0=y(i+1)-yh0;
    e1=e0/(1+phi'*D*phi);
    theta=theta+D*phi*e1;
    e(i+1)=e1;
    D=D - (D*(phi*phi')*D/(1+phi'*D*phi));
end

A = [1 theta(1:na)'];
B = [0 theta(na+1:na+nb)'];
C = [1 theta(na+nb+1:end)'];

end

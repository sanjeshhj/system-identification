clear all
%run section wise for better view of figures and analyse data

%first parameter set
na=3;nb=2;nc=5;

%Second parameter set (uncomment for the second run)
%na=5;nb=3;nc=10;

%data set for estimation
load('quanser/TestDoubleStep2.mat')
y=tach(1000:6000,1);
y=y-mean(y);
yh(1:na)=y(1:na);
u=V1(1000:6000);
u=u-mean(u);

%data set for validation
load('quanser/TestStep2.mat')
yv=tach(1000:6000,1);
yv=yv-mean(yv);
yh(1:na)=yv(1:na);
uv=V1(1000:6000);
uv=uv-mean(uv);
%% RLS 

[Ar,Br,thetao]=rls(y,u,na,nb,0);

%yh and error computation
for i =max(na,nb):length(uv)-1
    phi=[-yv(i:-1:i-na+1); uv(i:-1:i-nb+1)];
    yh(i+1)=thetao'*phi;  
end
err=yv-yh';
plot(err)

%using matlab resid function for residual analysis
sys3=idpoly(Ar,Br);
data3=iddata(yv,uv);
figure
resid(sys3,data3)

%error computation using inbuilt algorithm
errn=pe(sys3,data3);
figure
autocorr(errn.y)

%using the matlab arx fucntion which uses rls to estimate parameters
armodel=arx(data3,[na nb 0]);
resid(armodel,data3)
figure
present(armodel)


%% ELS
[Ae,Be,Ce,thetae]=els(y,u,na,nb,nc,0);


%yh and error computation
e=zeros(length(uv),1);
D=1000000*eye(na+nb+nc);
for i =max(max(na,nb),nc):length(uv)-1
    phi=[-yv(i:-1:i-na+1); uv(i:-1:i-nb+1);  e(i : -1 : i -nc+1)];
    yh(i+1)=thetae'*phi;
    e(i+1)=(yv(i+1)-yh(i+1))/(1+phi'*D*phi);
    D=D-(D*(phi*phi')*D/(1+phi'*D*phi));
end
err=yv-yh';
figure
plot(err)

%using matlab resid function for residual analysis
sys4=idpoly(Ae,Be,Ce);
data4=iddata(yv,uv);
figure
resid(sys4,data4)

%error computation using inbuilt algorithm
errn4=pe(sys4,data4);
figure
autocorr(errn4.y)

%using the matlab arxmax fucntion which uses els to estimate parameters
armmodel=armax(data4,[na nb nc 0]);
present(armmodel)
figure
resid(armmodel,data4)


%% GLS
[Ag,Bg,Cg,thetag]=gls(y,u,na,nb,nc,0);

%yh and error computation
e=zeros(length(u),1);
a=zeros(length(u),1);
D=1000000*eye(na+nb+nc);
for i=max(max(na,nb),nc):length(u)-1
    A=Ag;
    B=Bg;
    theta=thetag;
    a(i+1)= (A)*[yv(i+1:-1:i-na+1)]-B*[uv(i+1:-1:i-nb+1)];
    phi=[ -yv(i: -1: i-na+1) ; uv(i : -1 : i -nb+1); -a(i : -1 : i -nc+1)];
    yh(i+1)=theta'*phi;
    e0=yv(i+1)-yh(i+1);
    e1=e0/(1+phi'*D*phi);
    theta=theta+D*phi*e1;
    D=D - (D*(phi*phi')*D/(1+phi'*D*phi));
end
err=yv-yh';
figure
plot(err)

%using matlab resid function for residual analysis
sys5=idpoly(Ag,Bg,Cg);
data5=iddata(y,u);
figure
resid(sys5,data5)

%error computation using inbuilt algorithm
errn5=pe(sys5,data5);
figure
autocorr(errn5.y)

%using the matlab polyest fucntion which uses gls to estimate parameters
glsmodel = polyest(data5,[na nb 0 nc 0 0])
figure
resid(glsmodel,data5)


%% output error model

[Ao,Bo,thetao]=oem(y,u,na,nb,0);

%yh and error computation
for i =max(na,nb):length(uv)-1
    phi=[-yv(i:-1:i-na+1); uv(i:-1:i-nb+1)];
    yh(i+1)=thetao'*phi;  
end
err=yv-yh';
figure
plot(err)

%using matlab resid function for residual analysis
sys6=idpoly(Ao,Bo);
data6=iddata(yv,uv);
figure
resid(sys6,data6)

%error computation using inbuilt algorithm
errn=pe(sys6,data6);
figure
autocorr(errn.y)


%using the matlab polyest fucntion which uses gls to estimate parameters
oemodel=oe(data6,[na nb 0]);
figure
resid(oemodel,data6)
present(oemodel)





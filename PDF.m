clear;
clc;

%%% system parameters
r0=0.2;
miu0=0.25;
c0=0.02;
k0=10;
E0=0.4;
C0=1;
sigmax=0.1;
sigmay=0.1;

%%% parameters
xmin=0;
xmax=30;
ymin=0;
ymax=30;
J=100;
v=linspace(-2,2,4*J+1);
v(1)=[];
v(end)=[];
hx=v(2)-v(1);
w=linspace(-2,2,4*J+1);
w(1)=[];
w(end)=[];
hy=w(2)-w(1);
deltat=1e-3;
nT=1e4;
alpha=1.5;
% sigma=0.001;
% sigmax=(epsilong*sigma)^(1/alpha);
% sigmay=sigma^(1/alpha);
Calpha=alpha*gamma((1+alpha)/2)/(2^(1-alpha)*sqrt(pi)*gamma(1-alpha/2));
sigmaB=0;
Chx=sigmaB/2-(2*sigmax/(xmax-xmin))^alpha*Calpha*zeta(alpha-1)*hx^(2-alpha);
Chy=sigmaB/2-(2*sigmay/(ymax-ymin))^alpha*Calpha*zeta(alpha-1)*hy^(2-alpha);
C1x=Chx/hx^2;
C2x=Calpha/alpha*(2*sigmax/(xmax-xmin))^alpha;
C3x=Calpha*(2*sigmax/(xmax-xmin))^alpha*hx;
C1y=Chy/hy^2;
C2y=Calpha/alpha*(2*sigmay/(ymax-ymin))^alpha;
C3y=Calpha*(2*sigmay/(ymax-ymin))^alpha*hy;

%%% initialization
x=(xmax-xmin)/2*v(J+1:3*J-1)+(xmax+xmin)/2;
y=(ymax-ymin)/2*w(J+1:3*J-1)+(ymax+ymin)/2;
P1=normpdf(x,(xmax+xmin)/2,0.01*(xmax-xmin));
P2=normpdf(y,(ymax+ymin)/2,0.01*(ymax-ymin));
P0=P1'*P2;
P=reshape(P0,(2*J-1)*(2*J-1),1);

A1x=zeros(2*J-1,2*J-1);
for i=1:2*J-1
    A1x(i,:)=1./abs(v(2*J+1-i:4*J-1-i)).^(1+alpha);
    A1x(i,i)=0;
end
v1=0.5./abs(v(1:2*J-1)).^(1+alpha);
v2=0.5./abs(v(2*J+1:4*J-1)).^(1+alpha);
v1=flipud(v1);
v2=flipud(v2);
A0=[v1' A1x v2'];
for i=1:2*J-1
    A1x(i,i)=-sum(A0(i,:));
end
A1x=C3x*A1x;
v0=C2x./(1+v(J+1:3*J-1)).^alpha+C2x./(1-v(J+1:3*J-1)).^alpha;
A1x=A1x-2*C1x*eye(2*J-1)-diag(v0);
A1x(1:2*J-2,2:2*J-1)=A1x(1:2*J-2,2:2*J-1)+C1x*eye(2*J-2);
A1x(2:2*J-1,1:2*J-2)=A1x(2:2*J-1,1:2*J-2)+C1x*eye(2*J-2);

A1y=zeros(2*J-1,2*J-1);
for i=1:2*J-1
    A1y(i,:)=1./abs(w(2*J+1-i:4*J-1-i)).^(1+alpha);
    A1y(i,i)=0;
end
w1=0.5./abs(w(1:2*J-1)).^(1+alpha);
w2=0.5./abs(w(2*J+1:4*J-1)).^(1+alpha);
w1=flipud(w1);
w2=flipud(w2);
A0=[w1' A1y w2'];
for i=1:2*J-1
    A1y(i,i)=-sum(A0(i,:));
end
A1y=C3y*A1y;
w0=C2y./(1+w(J+1:3*J-1)).^alpha+C2y./(1-w(J+1:3*J-1)).^alpha;
A1y=A1y-2*C1y*eye(2*J-1)-diag(w0);
A1y(1:2*J-2,2:2*J-1)=A1y(1:2*J-2,2:2*J-1)+C1y*eye(2*J-2);
A1y(2:2*J-1,1:2*J-2)=A1y(2:2*J-1,1:2*J-2)+C1y*eye(2*J-2);

% L1=1:(2*J-1);
% L2=0:(2*J-1):(2*J-2)*(2*J-1);
% A1=sparse((2*J-1)*(2*J-1),(2*J-1)*(2*J-1));
% for i=1:2*J-1
%     A1((i-1)*(2*J-1)+L1,(i-1)*(2*J-1)+L1)=A1((i-1)*(2*J-1)+L1,(i-1)*(2*J-1)+L1)+A1x;
%     A1(L2+i,L2+i)=A1(L2+i,L2+i)+A1y;
% end

%%% drift
[X,Y]=meshgrid(x,y);
X=X';
Y=Y';
fXY=C0*X.^2./(X.^2+100);
fx=r0*X-c0*X.^2-Y.*fXY;
Mfx1=(fx+max(max(abs(fx))))/2;
Mfx2=(fx-max(max(abs(fx))))/2;
fy=(E0*fXY-miu0).*Y;
Mfy1=(fy+max(max(abs(fy))))/2;
Mfy2=(fy-max(max(abs(fy))))/2;

L=nT/10;
Ptotal=[];

%%% integral
for i=1:nT
    U=P;
    
    Ux=reshape(U,2*J-1,2*J-1);
    phi=Mfx1.*Ux;
    phixp=2/(xmax-xmin)*leftbiased(phi,hx);
    phixp=reshape(phixp,(2*J-1)*(2*J-1),1);
    phi=Mfx2.*Ux;
    phixn=2/(xmax-xmin)*rightbiased(phi,hx);
    phixn=reshape(phixn,(2*J-1)*(2*J-1),1);
    
    phi=Mfy1.*Ux;
    phi=phi';
    phiyp=2/(ymax-ymin)*leftbiased(phi,hy);
    phiyp=reshape(phiyp',(2*J-1)*(2*J-1),1);
    phi=Mfy2.*Ux;
    phi=phi';
    phiyn=2/(ymax-ymin)*rightbiased(phi,hy);
    phiyn=reshape(phiyn',(2*J-1)*(2*J-1),1);
    U0=A1x*Ux+(A1y*Ux')';
    U0=reshape(U0,(2*J-1)*(2*J-1),1);
    U1=U+deltat*(U0-phixp-phixn-phiyp-phiyn);
    
    
    Ux=reshape(U1,2*J-1,2*J-1);
    phi=Mfx1.*Ux;
    phixp=2/(xmax-xmin)*leftbiased(phi,hx);
    phixp=reshape(phixp,(2*J-1)*(2*J-1),1);
    phi=Mfx2.*Ux;
    phixn=2/(xmax-xmin)*rightbiased(phi,hx);
    phixn=reshape(phixn,(2*J-1)*(2*J-1),1);
    
    phi=Mfy1.*Ux;
    phi=phi';
    phiyp=2/(ymax-ymin)*leftbiased(phi,hy);
    phiyp=reshape(phiyp',(2*J-1)*(2*J-1),1);
    phi=Mfy2.*Ux;
    phi=phi';
    phiyn=2/(ymax-ymin)*rightbiased(phi,hy);
    phiyn=reshape(phiyn',(2*J-1)*(2*J-1),1);
    U0=A1x*Ux+(A1y*Ux')';
    U0=reshape(U0,(2*J-1)*(2*J-1),1);
    U2=3/4*U+1/4*U1+deltat/4*(U0-phixp-phixn-phiyp-phiyn);
    
    
    Ux=reshape(U2,2*J-1,2*J-1);
    phi=Mfx1.*Ux;
    phixp=2/(xmax-xmin)*leftbiased(phi,hx);
    phixp=reshape(phixp,(2*J-1)*(2*J-1),1);
    phi=Mfx2.*Ux;
    phixn=2/(xmax-xmin)*rightbiased(phi,hx);
    phixn=reshape(phixn,(2*J-1)*(2*J-1),1);
    
    phi=Mfy1.*Ux;
    phi=phi';
    phiyp=2/(ymax-ymin)*leftbiased(phi,hy);
    phiyp=reshape(phiyp',(2*J-1)*(2*J-1),1);
    phi=Mfy2.*Ux;
    phi=phi';
    phiyn=2/(ymax-ymin)*rightbiased(phi,hy);
    phiyn=reshape(phiyn',(2*J-1)*(2*J-1),1);
    U0=A1x*Ux+(A1y*Ux')';
    U0=reshape(U0,(2*J-1)*(2*J-1),1);
    P=1/3*U+2/3*U2+deltat*2/3*(U0-phixp-phixn-phiyp-phiyn);
    
    if rem(i,L)==0
        Ptotal=[Ptotal P];
    end
end

% for i=1:nT
%     U=P;
%     
%     Ux=reshape(U,2*J-1,2*J-1);
%     phi=Mfx1.*Ux;
%     phixp=2/(xmax-xmin)*leftbiased(phi,hx);
%     phixp=reshape(phixp,(2*J-1)*(2*J-1),1);
%     phi=Mfx2.*Ux;
%     phixn=2/(xmax-xmin)*rightbiased(phi,hx);
%     phixn=reshape(phixn,(2*J-1)*(2*J-1),1);
%     
%     phi=Mfy1.*Ux;
%     phi=phi';
%     phiyp=2/(ymax-ymin)*leftbiased(phi,hy);
%     phiyp=reshape(phiyp',(2*J-1)*(2*J-1),1);
%     phi=Mfy2.*Ux;
%     phi=phi';
%     phiyn=2/(ymax-ymin)*rightbiased(phi,hy);
%     phiyn=reshape(phiyn',(2*J-1)*(2*J-1),1);
%     U1=U+deltat*(A1*U-phixp-phixn-phiyp-phiyn);
%     
%     
%     Ux=reshape(U1,2*J-1,2*J-1);
%     phi=Mfx1.*Ux;
%     phixp=2/(xmax-xmin)*leftbiased(phi,hx);
%     phixp=reshape(phixp,(2*J-1)*(2*J-1),1);
%     phi=Mfx2.*Ux;
%     phixn=2/(xmax-xmin)*rightbiased(phi,hx);
%     phixn=reshape(phixn,(2*J-1)*(2*J-1),1);
%     
%     phi=Mfy1.*Ux;
%     phi=phi';
%     phiyp=2/(ymax-ymin)*leftbiased(phi,hy);
%     phiyp=reshape(phiyp',(2*J-1)*(2*J-1),1);
%     phi=Mfy2.*Ux;
%     phi=phi';
%     phiyn=2/(ymax-ymin)*rightbiased(phi,hy);
%     phiyn=reshape(phiyn',(2*J-1)*(2*J-1),1);
%     U2=3/4*U+1/4*U1+deltat/4*(A1*U1-phixp-phixn-phiyp-phiyn);
%     
%     
%     Ux=reshape(U2,2*J-1,2*J-1);
%     phi=Mfx1.*Ux;
%     phixp=2/(xmax-xmin)*leftbiased(phi,hx);
%     phixp=reshape(phixp,(2*J-1)*(2*J-1),1);
%     phi=Mfx2.*Ux;
%     phixn=2/(xmax-xmin)*rightbiased(phi,hx);
%     phixn=reshape(phixn,(2*J-1)*(2*J-1),1);
%     
%     phi=Mfy1.*Ux;
%     phi=phi';
%     phiyp=2/(ymax-ymin)*leftbiased(phi,hy);
%     phiyp=reshape(phiyp',(2*J-1)*(2*J-1),1);
%     phi=Mfy2.*Ux;
%     phi=phi';
%     phiyn=2/(ymax-ymin)*rightbiased(phi,hy);
%     phiyn=reshape(phiyn',(2*J-1)*(2*J-1),1);
%     P=1/3*U+2/3*U2+deltat*2/3*(A1*U2-phixp-phixn-phiyp-phiyn);
% end

P1=reshape(P,2*J-1,2*J-1);
figure;
mesh(X,Y,P1);

P2=reshape(P0,2*J-1,2*J-1);
figure;
mesh(X,Y,P2);


% for i=1:10
%     P3=reshape(Ptotal(:,i),2*J-1,2*J-1);
%     figure;
%     mesh(X,Y,P3);
% end


function phixp=leftbiased(phi,h)

% delta=1e-6;
% dpphij=[phi(2:end,:);phi(end,:)]-phi;
% dpphij_1=phi-[phi(1,:);phi(1:end-1,:)];
% dpphij_2=[phi(1,:);phi(1:end-1,:)]-[phi(1,:);phi(1,:);phi(1:end-2,:)];
% dndpphij=dpphij-[dpphij(1,:);dpphij(1:end-1,:)];
% dndpphij_1=dpphij_1-[dpphij_1(1,:);dpphij_1(1:end-1,:)];
% 
% r=(delta+dndpphij_1.^2)./(delta+dndpphij.^2);
% ommega=1./(1+2*r.^2);
% phixp=(dpphij_1+dpphij)/2/h-(dpphij_2-2*dpphij_1+dpphij).*ommega/2/h;


delta=1e-6;
L=zeros(1,length(phi(1,:)));
dpphij=[phi(2:end,:);L]-phi;
dpphij_1=phi-[L;phi(1:end-1,:)];
dpphij_2=[L;phi(1:end-1,:)]-[L;L;phi(1:end-2,:)];
dndpphij=dpphij-[L;dpphij(1:end-1,:)];
dndpphij_1=dpphij_1-[L;dpphij_1(1:end-1,:)];

r=(delta+dndpphij_1.^2)./(delta+dndpphij.^2);
ommega=1./(1+2*r.^2);
phixp=(dpphij_1+dpphij)/2/h-(dpphij_2-2*dpphij_1+dpphij).*ommega/2/h;

function phixn=rightbiased(phi,h)

% delta=1e-6;
% dpphij=[phi(2:end,:);phi(end,:)]-phi;
% dpphij_1=phi-[phi(1,:);phi(1:end-1,:)];
% dpphij_p1=[phi(3:end,:);phi(end,:);phi(end,:)]-[phi(2:end,:);phi(end,:)];
% dndpphij=dpphij-[dpphij(1,:);dpphij(1:end-1,:)];
% dndpphij_p1=dpphij_p1-[dpphij_p1(1,:);dpphij_p1(1:end-1,:)];
% 
% r=(delta+dndpphij_p1.^2)./(delta+dndpphij.^2);
% ommega=1./(1+2*r.^2);
% phixn=(dpphij_1+dpphij)/2/h-(dpphij_p1-2*dpphij+dpphij_1).*ommega/2/h;


delta=1e-6;
L=zeros(1,length(phi(1,:)));
dpphij=[phi(2:end,:);L]-phi;
dpphij_1=phi-[L;phi(1:end-1,:)];
dpphij_p1=[phi(3:end,:);L;L]-[phi(2:end,:);L];
dndpphij=dpphij-[L;dpphij(1:end-1,:)];
dndpphij_p1=dpphij_p1-[L;dpphij_p1(1:end-1,:)];

r=(delta+dndpphij_p1.^2)./(delta+dndpphij.^2);
ommega=1./(1+2*r.^2);
phixn=(dpphij_1+dpphij)/2/h-(dpphij_p1-2*dpphij+dpphij_1).*ommega/2/h;

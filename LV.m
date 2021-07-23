%  Competing species model.
clear
hold on
sys = @(t,x) [(2./3).*x(1)-(4./3).*x(1).*x(2);
x(1).*x(2)-x(2)];
options = odeset('RelTol',1e-4,'AbsTol',1e-4);
[t,xa]=ode45(sys,[0 100],[7.1 .1],options);
plot(xa(:,1),xa(:,2))
axis([0 8 0 8])
fsize=15;
set(gca,'XTick',0:1:8,'FontSize',fsize)
set(gca,'YTick',0:1:8,'FontSize',fsize)
xlabel('x(t)','FontSize',fsize)
ylabel('y(t)','FontSize',fsize)
hold off
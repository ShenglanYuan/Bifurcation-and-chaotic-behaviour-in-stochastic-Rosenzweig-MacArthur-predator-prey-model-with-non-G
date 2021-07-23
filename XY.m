% Predator-prey model.
% Time series of a predator/prey model.
clear
hold on
sys = @(t,x) [2.5.*x(1)-0.02.*(x(1).^2)-((x(1).^2).*x(2))./(x(1).^2+100);
((0.4.*(x(1).^2))./(x(1).^2+100)-0.2).*x(2)];
options = odeset('RelTol',1e-4,'AbsTol',1e-4);
[t,xa]=ode45(sys,[0 200],[7.1 .1],options);
plot(t,xa(:,1),'r')
plot(t,xa(:,2),'b')
legend('prey','predator')
axis([0 200 0 150])
fsize=15;
set(gca,'XTick',0:50:200,'FontSize',fsize)
set(gca,'YTick',0:30:150,'FontSize',fsize)
xlabel('time','FontSize',fsize)
ylabel('population','FontSize',fsize)
hold off
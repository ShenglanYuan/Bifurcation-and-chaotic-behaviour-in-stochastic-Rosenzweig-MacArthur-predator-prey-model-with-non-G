% Program: Rosenzweig-MacArthur.
clear
sys = @(t,x) [0.2.*x(1)-0.02.*(x(1)).^2-x(2).*((x(1)).^2./((x(1))^2+100));(0.4.*((x(1)).^2./((x(1)).^2+100))-0.25).*x(2)];
vectorfield(sys,0:1:80,0:1:60);
hold on
sep=1;
for x0=0:sep:20
for y0=0:sep:20
[ts,xs] = ode45(sys,[0 40],[x0 y0]);
plot(xs(:,1),xs(:,2))
end
end
for x0=0:sep:20
for y0=0:sep:20
[ts,xs] = ode45(sys,[0 -40],[x0 y0]);
plot(xs(:,1),xs(:,2))
end
end
hold off
axis([0 20 0 20])
fsize=15;
set(gca,'XTick',0:4:20,'FontSize',fsize)
set(gca,'YTick',0:4:20,'FontSize',fsize)
xlabel('X(t)','FontSize',fsize)
ylabel('Y(t)','FontSize',fsize)
hold off
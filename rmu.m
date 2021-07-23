clear all
fsize=15;
x=0:.001:4;
figure;
plot(x,(0.4.*(x.^2))./(x.^2+0.04))
text(1,0.37,'\leftarrow curve T','Fontsize',fsize)
axis([0 4 0 0.5])
set(gca,'XTick',0:1:4,'FontSize',fsize)
set(gca,'YTick',0:0.1:0.5,'FontSize',fsize)
xlabel('r','Fontsize',fsize);
ylabel('\mu','Fontsize',fsize);
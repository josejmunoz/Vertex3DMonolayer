function PlotWoundEvolution(Set,Ablated)
nIncr=length(Set.t);
Height=Ablated.Height(1:nIncr);
Volume=Ablated.Volume(1:nIncr);
AreaTop=Ablated.AreaTop(1:nIncr);
AreaBottom=Ablated.AreaBottom(1:nIncr);
AreaLateral=Ablated.AreaLateral(1:nIncr);

figure(1)
plot(Set.t,Height/Height(1),'-','LineWidth',2,'DisplayName','Height'); 
grid on;
xlabel('Time','FontSize',20);
legend('show')

figure(2)
plot(Set.t,Volume/Volume(1),'>-','LineWidth',2,'DisplayName','Volume'); 
grid on;
xlabel('Time','FontSize',20);
legend('show')

figure(3)
plot(Set.t,(AreaTop/AreaTop(1))*100,'O-','LineWidth',2,'DisplayName','AreaTop'); 
grid on;
xlabel('Time','FontSize',20);
legend('show')

figure(4)
plot(Set.t,AreaBottom/AreaBottom(1),'^-','LineWidth',2,'DisplayName','AreaBottom'); 
grid on;
xlabel('Time','FontSize',20);
legend('show')

figure(5)
plot(Set.t,AreaLateral,'*-','LineWidth',2,'DisplayName','AreaLateral'); 
grid on;
xlabel('Time','FontSize',20);
legend('show')

end 
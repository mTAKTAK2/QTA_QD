function onlineDisplay(x1,x2,x3,x4,PHU1,PHL1,PHU2,PHL2,PHU3,PHL3,PHU4,PHL4,qw)    
figure(1)
subplot(221)
x=plot(x1); hold on,
y=plot(PHU1,':');hold on,z=plot(PHL1,'--r');
legend([x y z],'Res1','PHPaaU','PHPaaL','location','best')
xlabel(['slice position ' num2str(qw)])
subplot(222)
x=plot(x2); hold on,
y=plot(PHU2,':');hold on,z=plot(PHL2,'--r');
legend([x y z],'Res2','PHPaaU','PHPaaL','location','best')
xlabel(['slice position ' num2str(qw)])
subplot(223)
x=plot(x3);hold on,
y=plot(PHU3,':');hold on,z=plot(PHL3,'--r');
legend([x y z],'Res3','PHPaaU','PHPaaL','location','best')
xlabel(['slice position ' num2str(qw)])
subplot(224)
x=plot(x4);hold on,
y=plot(PHU4,':');hold on,z=plot(PHL4,'--r');
legend([x y z],'Res4','PHPaaU','PHPaaL','location','best')
xlabel(['slice position ' num2str(qw)])
%pause 
end

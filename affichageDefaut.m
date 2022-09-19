clear all
%load sans_fault_periode04_.mat
%load fault_C1_moins_periode04_td035.mat    
%load fault_C1_moins_periode04_td045.mat    
%load fault_C1_plus_periode04_td035.mat     
%load fault_C1_plus_periode04_td045.mat     
%load fault_C2_moins_periode04_td035.mat    
%load fault_C2_moins_periode04_td045.mat    
%load fault_C2_plus_periode04_td035.mat     
%load fault_C2_plus_periode04_td045.mat     
%load fault_IMs_moins_periode04_td035.mat   
%load fault_IMs_moins_periode04_td045.mat   
%load fault_IMs_plus_periode04_td035.mat    
%load fault_IMs_plus_periode04_td045.mat    
%load fault_IMu_moins_periode04_td035.mat   
%load fault_IMu_moins_periode04_td045.mat   
%load fault_IMu_plus_periode04_td035.mat    
load fault_IMu_plus_periode04_td045.mat    
%load fault_MR1_moins_periode04_td0p35.mat  
%load fault_MR1_moins_periode04_td0p45.mat  
%load fault_MR1_plus_periode04_td0p35.mat   
%load fault_MR1_plus_periode04_td0p45.mat   
%load fault_MR2_moins_periode04_td0p35.mat  
%load fault_MR2_moins_periode04_td0p45.mat  
%load fault_MR2_plus_periode04_td0p35.mat   
%load fault_MR2_plus_periode04_td0p45.mat   
%load fault_MR3_moins_periode04_td0p35.mat  
%load fault_MR3_moins_periode04_td0p45.mat  
%load fault_MR3_plus_periode04_td0p35.mat   
%load fault_MR3_plus_periode04_td0p45.mat   
%load fault_MR4_moins_periode04_td0p35.mat  
%load fault_MR4_moins_periode04_td0p45.mat  
%load fault_MR4_plus_periode04_td0p35.mat   
%load fault_MR4_plus_periode04_td0p45.mat   
nom='fault IMu plus periode04 td045';
[mesure,nSamp]=size(res); % 5x805
res1=res(3,:);%P1
res2=res(2,:);%P2
res3=res(5,:);%Vs
res4=res(4,:);%Vu
temp=res(1,:);%time
tsamp=1:nSamp;%sample number
figure(1)
%Title(nom)
%set(gcf,'Title',nom)

subplot(221),plot(res1),xlabel('samples'),ylabel('res1'),Title(nom),set(gca,'XGrid','on','XTick',[0 200 400 600 800 1000])
subplot(222),plot(res2),xlabel('samples'),ylabel('res2'),set(gca,'XGrid','on','XTick',[0 200 400 600 800 1000])
subplot(223),plot(res3),xlabel('samples'),ylabel('res3'),set(gca,'XGrid','on','XTick',[0 200 400 600 800 1000])
subplot(224),plot(res4),xlabel('samples'),ylabel('res4'),set(gca,'XGrid','on','XTick',[0 200 400 600 800 1000])


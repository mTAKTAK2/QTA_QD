
clear all

signature={'-','0','+'};
message1=['No Abrupt change'];
%%%%%%%%%%%%%%%%%%%%%%%%%
% .: Data generation :. %
%%%%%%%%%%%%%%%%%%%%%%%%%
load fault_C1_moins_periode04_td035.mat    
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
%load fault_IMu_plus_periode04_td045.mat    
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


fid=fopen('essai.xls','w');
fprintf(fid,'%s\n','essai');
fprintf(fid,'%s\t %s\t %s\t %s\t %s\t %s\n','','slice','res1','res2','res3','res4');

[mesure,nSamp]=size(res); % 5x805
res1=res(3,:);%P1
res2=res(2,:);%P2
res3=res(5,:);%Vs
res4=res(4,:);%Vu
temp=res(1,:);%time
tsamp=1:nSamp;%sample number

zx=[];

N1=100;
N2=20;
alpha=2;
qw=N1;
c=1;

while qw < nSamp
    
    indiceU1(1)=1;indiceU2(1)=1;indiceU3(1)=1;indiceU4(1)=1;
    indiceL1(1)=1;indiceL2(1)=1;indiceL3(1)=1;indiceL4(1)=1;
    r=[];tend=[];st=[];  
    x1=res1(qw-N1+1:qw);
    x2=res2(qw-N1+1:qw);
    x3=res3(qw-N1+1:qw);
    x4=res4(qw-N1+1:qw);
    t_sli_wind=qw-N1+1:qw;
    %%%%%%%%%%%%%%%%%%%%%%%% Z NORMALISATION
    x1 = (x1 - mean(x1))/std(x1);
    x2 = (x2 - mean(x2))/std(x2);
    x3 = (x3 - mean(x3))/std(x3);
    x4 = (x4 - mean(x4))/std(x4);
    x=[x1 x2 x3 x4];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIN ZNORM
    %%%%%%%%%%%%%%%%%%%%%%%%% PAA + PHT
    [lin_paa1,PHU1,PHL1,deltaPHU1,deltaPHL1]=PiAA(x1,N1,25,1,0);
    [lin_paa2,PHU2,PHL2,deltaPHU2,deltaPHL2]=PiAA(x2,N1,25,1,0);
    [lin_paa3,PHU3,PHL3,deltaPHU3,deltaPHL3]=PiAA(x3,N1,25,1,0);
    [lin_paa4,PHU4,PHL4,deltaPHU4,deltaPHL4]=PiAA(x4,N1,25,1,0);
    %%%%%%%%%%%%%%%%%%%%%%%%% Fin PAA + PHT
    maxPHU1=max(find(PHU1==max(PHU1)));%sample
    maxPHL1=max(find(PHL1==max(PHL1)));
    maxPHU2=max(find(PHU2==max(PHU2)));
    maxPHL2=max(find(PHL2==max(PHL2)));
    maxPHU3=max(find(PHU3==max(PHU3)));
    maxPHL3=max(find(PHL3==max(PHL3)));
    maxPHU4=max(find(PHU4==max(PHU4)));
    maxPHL4=max(find(PHL4==max(PHL4)));%sample
    SampMaxPHL=[maxPHL1,maxPHL2,maxPHL3,maxPHL4];
    SampMaxPHU=[maxPHU1,maxPHU2,maxPHU3,maxPHU4];
    PHL=[PHL1(maxPHL1),PHL2(maxPHL2),PHL3(maxPHL3),PHL4(maxPHL4)];
    PHU=[PHU1(maxPHU1),PHU2(maxPHU2),PHU3(maxPHU3),PHU4(maxPHU4)];
    lamda=7;
    ab_PHU=find(PHU>lamda);
    ab_PHL=find(PHL>lamda);
    if ~isempty(ab_PHU)
        if ~isempty(ab_PHL)
            min_ab_PHL=min(PHL(ab_PHL));
            min_ab_PHU=min(PHU(ab_PHU));
            min_ab_samp_PHL=min(SampMaxPHL(ab_PHL));
            min_ab_samp_PHU=min(SampMaxPHU(ab_PHU));
            if min_ab_samp_PHU < min_ab_samp_PHL
                message1=['(-) change in X_' num2str(find(SampMaxPHU==min_ab_samp_PHU(1)))];
                onlineDisplay(x1,x2,x3,x4,PHU1,PHL1,PHU2,PHL2,PHU3,PHL3,PHU4,PHL4,qw)
                indiceU1=find(deltaPHU1<0);indiceU2=find(deltaPHU2<0);indiceU3=find(deltaPHU3<0);indiceU4=find(deltaPHU4<0);
                indiceL1=find(deltaPHL1<0);indiceL2=find(deltaPHL2<0);indiceL3=find(deltaPHL3<0);indiceL4=find(deltaPHL4<0);
            else
                message1=['(+) change in X_' num2str(find(SampMaxPHL==min_ab_samp_PHL(1)))];
                onlineDisplay(x1,x2,x3,x4,PHU1,PHL1,PHU2,PHL2,PHU3,PHL3,PHU4,PHL4,qw)
                indiceU1=find(deltaPHU1<0);indiceU2=find(deltaPHU2<0);indiceU3=find(deltaPHU3<0);indiceU4=find(deltaPHU4<0);
                indiceL1=find(deltaPHL1<0);indiceL2=find(deltaPHL2<0);indiceL3=find(deltaPHL3<0);indiceL4=find(deltaPHL4<0);
            end
        else
            min_ab_PHU=min(PHU(ab_PHU));
            min_ab_samp_PHU=min(SampMaxPHU(ab_PHU));
            message1=['(-) change in X_' num2str(find(SampMaxPHU==min_ab_samp_PHU(1)))];
            onlineDisplay(x1,x2,x3,x4,PHU1,PHL1,PHU2,PHL2,PHU3,PHL3,PHU4,PHL4,qw)
            indiceU1=find(deltaPHU1<0);indiceU2=find(deltaPHU2<0);indiceU3=find(deltaPHU3<0);indiceU4=find(deltaPHU4<0);
            indiceL1=find(deltaPHL1<0);indiceL2=find(deltaPHL2<0);indiceL3=find(deltaPHL3<0);indiceL4=find(deltaPHL4<0);
        end
    elseif ~isempty(ab_PHL)
        if ~isempty(ab_PHU)
            min_ab_PHL=min(PHL(ab_PHL));
            min_ab_PHU=min(PHU(ab_PHU));
            min_ab_samp_PHL=min(SampMaxPHL(ab_PHL));
            min_ab_samp_PHU=min(SampMaxPHU(ab_PHU));
            if min_ab_samp_PHU < min_ab_samp_PHL
                message1=['(-) change in X_' num2str(find(SampMaxPHU==min_ab_samp_PHU(1)))];
                onlineDisplay(x1,x2,x3,x4,PHU1,PHL1,PHU2,PHL2,PHU3,PHL3,PHU4,PHL4,qw)
                indiceU1=find(deltaPHU1<0);indiceU2=find(deltaPHU2<0);indiceU3=find(deltaPHU3<0);indiceU4=find(deltaPHU4<0);
                indiceL1=find(deltaPHL1<0);indiceL2=find(deltaPHL2<0);indiceL3=find(deltaPHL3<0);indiceL4=find(deltaPHL4<0);
            else
                message1=['(+) change in X_' num2str(find(SampMaxPHL==min_ab_samp_PHL(1)))];
                onlineDisplay(x1,x2,x3,x4,PHU1,PHL1,PHU2,PHL2,PHU3,PHL3,PHU4,PHL4,qw)
                indiceU1=find(deltaPHU1<0);indiceU2=find(deltaPHU2<0);indiceU3=find(deltaPHU3<0);indiceU4=find(deltaPHU4<0);
                indiceL1=find(deltaPHL1<0);indiceL2=find(deltaPHL2<0);indiceL3=find(deltaPHL3<0);indiceL4=find(deltaPHL4<0);
            end
        else
            min_ab_PHL=min(PHL(ab_PHL));
            min_ab_samp_PHL=min(SampMaxPHL(ab_PHL));
            message1=['(+) change in X_' num2str(find(SampMaxPHL==min_ab_samp_PHL(1)))];
            onlineDisplay(x1,x2,x3,x4,PHU1,PHL1,PHU2,PHL2,PHU3,PHL3,PHU4,PHL4,qw)
            indiceU1=find(deltaPHU1<0);indiceU2=find(deltaPHU2<0);indiceU3=find(deltaPHU3<0);indiceU4=find(deltaPHU4<0);
            indiceL1=find(deltaPHL1<0);indiceL2=find(deltaPHL2<0);indiceL3=find(deltaPHL3<0);indiceL4=find(deltaPHL4<0);
        end
    end
    
    
    
    %if ~isempty(ab_PHL) | ~isempty(ab_PHU)
    %    min_ab_PHL=min(PHL(ab_PHL));
    %    min_ab_samp_PHL=SampMaxPHL(find(min(PHL(ab_PHL))));
    %    min_ab_PHU=min(PHU(ab_PHU));
    %    min_ab_samp_PHU=SampMaxPHU(find(min(PHU(ab_PHU))));
    %    if min_ab_samp_PHU < min_ab_samp_PHL
    %        message1=['(-) abrupt_change in X_' num2str(find(SampMaxPHU==min_ab_samp_PHU))];
    %    else
    %        message1=['(+) abrupt_change in X_' num2str(find(SampMaxPHL==min_ab_samp_PHL))];
    %    end
    %end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(fid,'%s\t %1.4f\t %1.4f\t %1.4f\t %1.4f\t %1.4f\n','PHU',qw,PHU1(maxPHU1),PHU2(maxPHU2),PHU3(maxPHU3),PHU4(maxPHU4));
    %fprintf(fid,'%s\t %1.4f\t %1.4f\t %1.4f\t %1.4f\t %1.4f\n','DeltaPHU',[],deltaPHU1(maxPHU1),deltaPHU2(maxPHU2),deltaPHU3(maxPHU3),deltaPHU4(maxPHU4));
    fprintf(fid,'%s\t %1.4f\t %1.4f\t %1.4f\t %1.4f\t %1.4f\n','DeltaPHU',[],deltaPHU1(indiceU1(1)),deltaPHU2(indiceU2(1)),deltaPHU3(indiceU3(1)),deltaPHU4(indiceU4(1)));
    fprintf(fid,'%s\t %1.4f\t %1.4f\t %1.4f\t %1.4f\t %1.4f\n','time in slice PHU','',maxPHU1,maxPHU2,maxPHU3,maxPHU4);
    fprintf(fid,'%s\t %1.4f\t %1.4f\t %1.4f\t %1.4f\t %1.4f\n','PHL',qw,PHL1(maxPHL1),PHL2(maxPHL2),PHL3(maxPHL3),PHL4(maxPHL4));
    %fprintf(fid,'%s\t %1.4f\t %1.4f\t %1.4f\t %1.4f\t %1.4f\n','DeltaPHL',[],deltaPHL1(maxPHL1),deltaPHL2(maxPHL2),deltaPHL3(maxPHL3),deltaPHL4(maxPHL4));
    fprintf(fid,'%s\t %1.4f\t %1.4f\t %1.4f\t %1.4f\t %1.4f\n','DeltaPHL',[],deltaPHL1(indiceL1(1)),deltaPHL2(indiceL2(1)),deltaPHL3(indiceL3(1)),deltaPHL4(indiceL4(1)));
    fprintf(fid,'%s\t %1.4f\t %1.4f\t %1.4f\t %1.4f\t %1.4f\n','time in slice PHL','',maxPHL1,maxPHL2,maxPHL3,maxPHL4);
    fprintf(fid,'%s\n',message1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pause
    close all
    clear x1 x2 x3 x4 t_sli_wind PHU1 PHL1 PHU2 PHL2 PHU3 PHL3 PHU4 PHL4 lin_paa1 lin_paa2 lin_paa3 lin_paa4 
    clear ab_PHL ab_PHU
    clear indiceU1 indiceU2 indiceU3 indiceU4 indiceL1 indiceL2 indiceL3 indiceL4
    qw=qw+N2;
end%while
ST=fclose(fid);


%figure(1)
%for i=2:5
%    j=i-1;
%   subplot(4,1,j)
%   plot(res(i,:))
%end
%[lin_paa1,PHU1,PHL1]=PiAA(res(2,1:800),800,10,1,0);
%[lin_paa2,PHU2,PHL2]=PiAA(res(3,1:800),800,10,1,0);
%[lin_paa3,PHU3,PHL3]=PiAA(res(4,1:800),800,10,1,0);
%[lin_paa4,PHU4,PHL4]=PiAA(res(5,1:800),800,10,1,0);
%figure(2)
%subplot(411)
%plot(PHU1,':'),hold on,plot(PHL1,'--r')
%subplot(412)
%plot(PHU2,':'),hold on,plot(PHL2,'--r')
%subplot(413)
%plot(PHU3,':'),hold on,plot(PHL3,'--r')
%subplot(414)
%plot(PHU4,':'),hold on,plot(PHL4,'--r')


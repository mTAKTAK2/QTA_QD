
function PHT
clear all
%N=300;  % Number of samples
%nl=2;    % Number of latent variables
%nx=6;    % Number of observable variables
%[tt,x,latgen]=datagenPH(N,nl,nx);
%r1=latgen(:,2);

%load x
%r1=x';
%w=1:1100;
%r1=[0.1*ones(1,400) 0.3*ones(1,300) 0.1*ones(1,300)]';
r1=[0.1*ones(1,450) 0.4*ones(1,250) 0.1*ones(1,300)]';
r=r1;
%r1=[0.1*ones(1,200) 5 0.1*ones(1,199) 4*ones(1,300) 0.1*ones(1,150) 5 0.1*ones(1,149)]';%w
%r1=[0.1*ones(1,100) 0.11*ones(1,200) 0.1*ones(1,400).*w(1:400) 0.3*ones(1,300) -0.1*ones(1,100)]';
r1n=sqrt(2.25)*rand(1000,1);%sqrt(variance)
%SNR=(sqrt(mean(r1.^2))/sqrt(mean(r1n.^2)))^2;%signal to noise ratio
r2n=sqrt(3.25)*rand(1000,1);
SNR=mean(r1.^2)/mean(r1n.^2);
SNR2=mean(r1.^2)/mean(r2n.^2);
DB1=10*log10(var(r1n));
DB2=10*log10(var(r1));
DB=DB2-DB1
SNR_bd=10*log10(SNR)
SNR_bd2=10*log10(SNR2)
r1=r+r1n; %signal with noise
r2=r+r2n;
%r2=r2+r2n;
r1=zscore(r1);
r2=zscore(r2);
t=1:length(r1);
data=[t' r1];
data2=[t' r2];
[lin_paa,PHU,PHL,deltaPHU,deltaPHL]=PiAA(r1,1000,20,1,0);
[lin_paa2,PHU2,PHL2,PPHU2,PPHL2]=PiAA(r2,1000,20,1,0);
%disp('max PHU')
%max(find(PHU==max(PHU(100:length(PHU)-100))))
%disp('amp max PHU')
%PHU(max(find(PHU==max(PHU))))
%disp('max PHL')
%max(find(PHL==max(PHL(100:length(PHL)-100))))
%disp('amp max PHL')
%PHL(max(find(PHL==max(PHL))))
%detect_indice2=p_h(data2,0.999,0.375,0.62);%with noise
detect_indice2=p_h(data2,0.999,0.237,0.645);%with noise
detect_indice=p_h(data,0.999,0.237,0.645);%with noise
%detect_indice=p_h(data,0.999,0.005,0.02);%without noise

subplot(221)
plot(data(:,1),data(:,2),'-.b'),
hold on
plot(lin_paa,'k')
line([data(1,1) data(end,1)],[mean(data(:,2)) mean(data(:,2))])
ylabel('data (variance 2.25)')
subplot(222)
plot(data2(:,1),data2(:,2),'-.b'),
hold on
plot(lin_paa2,'k')
line([data2(1,1) data2(end,1)],[mean(data2(:,2)) mean(data2(:,2))])
ylabel('data (variance 3.25)')
%subplot(412)
%plot(lin_paa,'r'),
%ylabel('PAA')
%subplot(413)
%x=plot(PHU,'--b');
%ylabel('--PHU ..PHL')
%hold on
%z=plot(PPHU,'+g'),
%y=plot(PHL,':k');
%w=plot(PPHL,'oy'),
%legend([x y z w],'PHU','PHL','PPHU','PPHL','location','best')
%legend([x y],'PHU','PHL','location','best')
subplot(223)
plot(PHL,':k')
hold on
plot(PHU,'--k')
ylabel('PHT+PAA')
xlabel('samples')
%plot(detect_indice,'b*')
%ylabel('drift indice')
%set(gca,'ytick',[-1 0 1])
%set(gca,'yticklabel',{'Down';'None';'Up'})
subplot(224)
plot(PHL2,':k')
hold on
plot(PHU2,'--k')
ylabel('PHT+PAA')
xlabel('sample')
%plot(detect_indice2,'b*')
%set(gca,'ytick',[-1 0 1])
%set(gca,'yticklabel',{'Down';'None';'Up'})

end



function [detect_indice]=p_h(res,alpha,gamma,delta)

detect_indice=[];
minsize=10;
raw=res(:,2);
%raw=(raw-mean(raw))/std(raw);
time=res(:,1);
dtime=diff(time);
T=dtime(1);
X=[];
t=1;
U=0;L=0;x_m=0;MinU=0;MaxL=0;
%MinU1=0;MinU2=0;MaxL1=0;MaxL2=0;
while t<length(time)
    if t < minsize
        X=[X;raw(t)];
        detect_indice(t)=0;
    end
    
    if t==minsize 
        x_m=mean(X);
        U=0;
        L=0;
        MinU=0;
        MaxL=0;
        detect_indice(t)=0;
    end
    
    if t>minsize
        U=alpha*U+(raw(t)-x_m-gamma);
        %MinU=MinU+abs(MinU*(1-alpha));
        if U<MinU
            MinU=U;
        end
        %if MinU<MinU2
        %    MinU2=MinU;
        %end
        %PHU=U-min(MinU1,MinU2);
        PHU=U-MinU;
        L=alpha*L+(raw(t)-x_m+gamma);
        %MaxL=MaxL-abs(MaxL*(1-alpha));
        if L>MaxL
            MaxL=L;
        end
        %if MaxL>MaxL2
        %    MaxL2=MaxL;
        %end
        PHL=MaxL-L;
        if PHL>delta 
            detect_indice(t)=-1;%down
            L=0;
            MaxL=0;
        elseif PHU>delta
            detect_indice(t)=1;%up
            U=0;
            MinU=0;
        else
            x_m=((T-1)*x_m+raw(t))/T;
            detect_indice(t)=0;
        end
    end
    t=t+1;
    end
end

function [p,PHU,PHL,deltaPHU,deltaPHL] =  PiAA(data, N, n,alpha,gamma)
% Input:
%   data              is the raw time series. 
%   N                 is the length of sliding window (use the length of the raw time series
%                     instead if you don't want to have sliding windows)
%   n                 is the number of symbols in the low dimensional approximation of the sub sequence.
% The variable "win_size" is assigned to N/n, this is the number of data points on the raw 
% time series that will be mapped to a single symbol, and can be imagined as the 
% "compression rate".
% win_size is the number of data points on the raw time series that will be
% mapped to a single symbol
MinUpaa=0;MaxLpaa=0;
win_size = floor(N/n);                         
% Scan accross the time series extract sub sequences, and converting them to strings.
for i = 1 : length(data) - (N -1)                                       
    
    % Remove the current subsection.
    sub_section = data(i:i + N -1); 
    GM=mean(sub_section);% Global mean slim
    
    % Z normalize it.
    %sub_section = (sub_section - mean(sub_section))/std(sub_section);     
    
    % take care of the special case where there is no dimensionality reduction
    if N == n
        PAA = sub_section;
        Upaa=cumsum(PAA-GM);%slim
        PHU=Upaa-min(Upaa);
        PHL=max(Upaa)-Upaa;
        %detec_ind=0;% no change (slim)
    % Convert to PAA.    
    else
        % N is not dividable by n
        if (N/n - floor(N/n))                               
            temp = zeros(n, N);
            for j = 1 : n
                temp(j, :) = sub_section;
            end
            expanded_sub_section = reshape(temp, 1, N*n);
            PAA = [mean(reshape(expanded_sub_section, N, n))];
            UUpaa=cumsum(PAA-GM);%up slim
            Upaa(1)=PAA(1)-GM-gamma;
            Lpaa(1)=PAA(1)-GM+gamma;
            PHU(1)=Upaa(1);PHL(1)=-Lpaa(1);
            deltaPHL(1)=PHU(1);%%%%%
            deltaPHU(1)=PHL(1);%%%%%%
            for h=2:length(PAA)
                Upaa(h)=alpha*Upaa(h-1)+(PAA(h)-GM-gamma);
                if Upaa(h)<MinUpaa
                    MinUpaa=Upaa(h);
                end
                Lpaa(h)=alpha*Lpaa(h-1)+(PAA(h)-GM+gamma);
                if Lpaa(h)>MaxLpaa
                    MaxLpaa=Lpaa(h);
                end
                PHU(h)=Upaa(h)-MinUpaa;
                PHL(h)=MaxLpaa-Lpaa(h);
                if PHU(h)<PHU(h-1)
                    MinUpaa=0;
                    Upaa(h)=0;
                    %MaxLpaa=0;
                    %Lpaa(h)=0;
                end
                if PHL(h)<PHL(h-1)
                    MaxLpaa=0;
                    Lpaa(h)=0;
                    %MinUpaa=0;
                    %Upaa(h)=0;
                end
            end
            PPHU=Upaa-min(Upaa);%slim
            PPHL=max(Upaa)-Upaa;%slim
        % N is dividable by n
        else                                                
            PAA = [mean(reshape(sub_section,win_size,n))];
            UUpaa=cumsum(PAA-GM);% slim
            Upaa(1)=PAA(1)-GM-gamma;
            Lpaa(1)=PAA(1)-GM+gamma;
            PHU(1)=Upaa(1);PHL(1)=-Lpaa(1);
            deltaPHL(1)=PHU(1);%%%%%
            deltaPHU(1)=PHL(1);%%%%%%
            for h=2:length(PAA)
                Upaa(h)=alpha*Upaa(h-1)+(PAA(h)-GM-gamma);
                if Upaa(h)<MinUpaa
                    MinUpaa=Upaa(h);
                end
                Lpaa(h)=alpha*Lpaa(h-1)+(PAA(h)-GM+gamma);
                if Lpaa(h)>MaxLpaa
                    MaxLpaa=Lpaa(h);
                end
                PHU(h)=Upaa(h)-MinUpaa;
                deltaPHU(h)=PHU(h)-PHU(h-1);%%%%%%%%%%%
                PHL(h)=MaxLpaa-Lpaa(h);
                deltaPHL(h)=PHL(h)-PHL(h-1);%%%%%%%%%%%%
                %if PHU(h)<PHU(h-1)
                if deltaPHU(h)<0
                    MinUpaa=0;
                    Upaa(h)=0;
                end
                %if PHL(h)<PHL(h-1)
                if deltaPHL(h)<0
                    MaxLpaa=0;
                    Lpaa(h)=0;
                end

            end

            
            PPHU=Upaa-min(Upaa);% slim
            PPHL=max(Upaa)-Upaa;%slim
        end
    end
    PAA = repmat(PAA, win_size, 1);
    PHU = repmat(PHU, win_size, 1);
    PHL = repmat(PHL, win_size, 1);
    %PPHU = repmat(PPHU, win_size, 1);
    %PPHL = repmat(PPHL, win_size, 1);
    deltaPHU = repmat(deltaPHU, win_size, 1);
    deltaPHL = repmat(deltaPHL, win_size, 1);

    p = reshape(PAA, N, 1);
    PHU = reshape(PHU, N, 1);
    PHL = reshape(PHL, N, 1);
    %PPHU = reshape(PPHU, N, 1);
    %PPHL = reshape(PPHL, N, 1);
    deltaPHU = reshape(deltaPHU, N, 1);
    deltaPHL = reshape(deltaPHL, N, 1);

end
   %plot(p),hold on,plot(PHU,':b'),plot(PHL,':r')
   %pause
end



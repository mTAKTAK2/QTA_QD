
function testPHTPAA
clear all
Abr_amp=[2.5 2.75 3 3.25 3.5];
vari=[2.25 2.5 2.75 3 3.25];
TP=0;FP=0;FN=0;TN=0;
fid=fopen('testPHPAA.xls','w');
fprintf(fid,'%s\n','PHT vs PH-PAA');
fprintf(fid,'%s\t %s\t %s\t %s\t %s\t %s\t %s\n','change amp.','variance','TP','FN','FP','maxPHPaa_U','maxPHPaa_L')

for i=1:length(Abr_amp)
    for j=1:length(vari)
        r1=[0.1*ones(1,400) Abr_amp(i)*ones(1,300) 0.1*ones(1,300)]';
        r1n=sqrt(vari(j))*rand(1000,1);%sqrt(variance)
        %SNR=mean(r1.^2)/mean(r1n.^2);
        %DB1=10*log10(var(r1n));
        %DB2=10*log10(var(r1));
        %DB=DB2-DB1;
        %SNR_bd=10*log10(SNR);
        r1=r1+r1n; %signal with noise
        r1=zscore(r1);
        t=1:length(r1);
        data=[t' r1];
        [lin_paa,PHU,PHL,PPHU,PPHL]=PiAA(r1,1000,20,1,0);
        maxPHU=max(find(PHU==max(PHU)));
        maxPHL=max(find(PHL==max(PHL)));
        detect_indice=p_h(data,0.999,0.375,0.62);%with noise
        FP1=find(detect_indice(1:399)~=0);
        FP2=find(detect_indice(701:999)~=0);
        FP=[FP1 FP2];
        TP=find(detect_indice(400:700)==1);
        if isempty(FP)
            FP=0;
        else
            FP=length(FP);
        end
        if isempty(TP)
            TP=0;
        else
            TP=length(TP);
            FN=(300-TP);
        end
        fprintf(fid,'%1.4f\t %1.4f\t %1.4f\t %1.4f\t %1.4f\t %1.4f\t %1.4f\n',Abr_amp(i),vari(j),TP,FN,FP,maxPHU,maxPHL);
        %detect_indice=p_h(data,0.999,0.005,0.02);%without noise
    end
end
ST=fclose(fid);
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

function [p,PHU,PHL,PPHU,PPHL] =  PiAA(data, N, n,alpha,gamma)
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
                end
                if PHL(h)<PHL(h-1)
                    MaxLpaa=0;
                    Lpaa(h)=0;
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
                end
                if PHL(h)<PHL(h-1)
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
    PPHU = repmat(PPHU, win_size, 1);
    PPHL = repmat(PPHL, win_size, 1);

    p = reshape(PAA, N, 1);
    PHU = reshape(PHU, N, 1);
    PHL = reshape(PHL, N, 1);
    PPHU = reshape(PPHU, N, 1);
    PPHL = reshape(PPHL, N, 1);

end
   %plot(p),hold on,plot(PHU,':b'),plot(PHL,':r')
   %pause
end



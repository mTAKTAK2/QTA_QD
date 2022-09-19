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
        deltaPHU=0;
        deltaPHL=0;
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
            deltaPHU(1)=PHU(1);%%%%%
            deltaPHL(1)=PHL(1);%%%%%%
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
                deltaPHU=PHU(h)-PHU(h-1);%%%%%%%%%%%%%%
                PHL(h)=MaxLpaa-Lpaa(h);
                deltaPHL(h)=PHL(h)-PHL(h-1);%%%%%%%%%%%%%
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
            PPHU=Upaa-min(Upaa);%slim
            PPHL=max(Upaa)-Upaa;%slim
        % N is dividable by n
        else                                                
            PAA = [mean(reshape(sub_section,win_size,n))];
            UUpaa=cumsum(PAA-GM);% slim
            Upaa(1)=PAA(1)-GM-gamma;
            Lpaa(1)=PAA(1)-GM+gamma;
            PHU(1)=Upaa(1);PHL(1)=-Lpaa(1);
            deltaPHU(1)=PHU(1);%%%%%
            deltaPHL(1)=PHL(1);%%%%%%
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
                deltaPHU(h)=PHU(h)-PHU(h-1);%%%%%%%%%%%%%%
                PHL(h)=MaxLpaa-Lpaa(h);
                deltaPHL(h)=PHL(h)-PHL(h-1);%%%%%%%%%%%%%
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

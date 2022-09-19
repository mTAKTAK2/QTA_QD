

function PAA
N=300;  % Number of samples
nl=2;    % Number of latent variables
nx=6;    % Number of observable variables

[tt,x,latgen]=datagen(N,nl,nx);


%load x
%r1=x;



r1=x(:,nx)';
r1=zscore(r1);

lin_paa=PiAA(r1,N,10);
%hold on,
plot(lin_paa,'b'),hold on,plot(r1,':r')
end







function [p] =  PiAA(data, N, n)
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
win_size = floor(N/n);                         


% Scan accross the time series extract sub sequences, and converting them to strings.
for i = 1 : length(data) - (N -1)                                       
    
    % Remove the current subsection.
    sub_section = data(i:i + N -1); 
    
    % Z normalize it.
    sub_section = (sub_section - mean(sub_section))/std(sub_section);     
    
    % take care of the special case where there is no dimensionality reduction
    if N == n
        PAA = sub_section;
        
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
        % N is dividable by n
        else                                                
            PAA = [mean(reshape(sub_section,win_size,n))];
        end
    end
    PAA = repmat(PAA, win_size, 1);
    p = reshape(PAA, N, 1);
   end
end


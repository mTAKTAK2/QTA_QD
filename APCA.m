function apca
N=300;  % Number of samples
nl=2;    % Number of latent variables
nx=6;    % Number of observable variables

[tt,x,latgen]=datagen(N,nl,nx);


%load x
%r1=x;



r1=x(:,nx)';
r1=zscore(r1);
seg_num=9;
lin_apca=ts2apca(r1,seg_num);
%hold on,plot(lin_apca,'b')
end


function TS = ts2apca(data,num_segment)
% Convert segment from ACPA to discrete time series
% ======================================
% == Input ==  a time series in a row vector; 
%              number of segments;
%===============================================
% segment : an output from apca()
%       - segment(i).lx  -> left position on x axis
%       - segment(i).rx  -> right position on x axis
%       - segment(i).mc  -> error (don't sure) >>>>> MERGE COST
%       - segment(i).y   -> y position of this interval
% == Output ==
% TS     : a discrete ts which has the same length as input ts (in high dimensional space) before using apca()
% ======================================
% last updated: Feb 26  2011
   segment=apca1(data,num_segment,1);
    TS(segment(end).rx)=0;
    for i=1:length(segment)
       TS(segment(i).lx:segment(i).rx)=round(segment(i).y);
       %TS(segment(i).lx:segment(i).rx)=max(segment(i).y);
    end

end




function segment = apca1(data,num_segments,dummy)

data=data';
if num_segments==1
    segment(1).lx=1;
    segment(1).rx=length(data);
    segment(1).y=0;
    segment(1).mc=0;
    segment(1).mean=mean(data);%by slim
    
else
    
left_x = [1 : 2 : size(data,1)-1];
right_x      = left_x + 1;
right_x(end) = size(data,1);
number_of_segments = length(left_x );

for i = 1 : number_of_segments 
   
   segment(i).lx = left_x(i);
   segment(i).rx = right_x(i);
   segment(i).mc = inf;
   segment(i).y = inf;
   segment(i).mean=inf;%by slim
end;

for i = 1 : number_of_segments - 1
   coef = polyfit([segment(i).lx :segment(i+1).rx]',data(segment(i).lx :segment(i+1).rx),0);
   best = (0*( [segment(i).lx :segment(i+1).rx]' ))+coef(1);
   segment(i).mc = sum((data([segment(i).lx :segment(i+1).rx]')-best).^2);
end;

while  length(segment) > num_segments  
      
   [value, i ] = min([segment(:).mc]);
      
     if i > 1 && i < length(segment) -1								
      
    	coef = polyfit([segment(i).lx :segment(i+2).rx]',data(segment(i).lx :segment(i+2).rx),0);        
    	best = (0 *( [segment(i).lx :segment(i+2).rx]' ))+coef(1);       
        segment(i).mc = sum((data([segment(i).lx :segment(i+2).rx]')-best).^2);
	 	segment(i).rx = segment(i+1).rx;
    	segment(i+1) = [];
        segment(i).mean=coef(1);%by slim
       	i = i - 1;    
        coef = polyfit([segment(i).lx :segment(i+1).rx]',data(segment(i).lx :segment(i+1).rx),0);
      	best = (0*( [segment(i).lx :segment(i+1).rx]' ))+coef(1);    
        segment(i).mc = sum((data([segment(i).lx :segment(i+1).rx]')-best).^2);
        segment(i).mean=coef(1);%by slim
       
   elseif i == 1
       
	   coef = polyfit([segment(i).lx :segment(i+2).rx]',data(segment(i).lx :segment(i+2).rx),0);
       best = (0*( [segment(i).lx :segment(i+2).rx]' ))+coef(1);
       segment(i).mc = sum((data([segment(i).lx :segment(i+2).rx]')-best).^2);
	   segment(i).rx = segment(i+1).rx;
      segment(i+1) = [];
      segment(i).mean=coef(1);%by slim
              
   else
     
      segment(i).rx = segment(i+1).rx;
      segment(i).mc = inf;
      segment(i+1) = [];    
      i = i - 1;       
      coef = polyfit([segment(i).lx :segment(i+1).rx]',data(segment(i).lx :segment(i+1).rx),0);
      best = (0*( [segment(i).lx :segment(i+1).rx]' ))+coef(1);    
     segment(i).mc = sum((data([segment(i).lx :segment(i+1).rx]')-best).^2);
     segment(i).mean=coef(1);%by slim
  end; 
          
end;


for i = 1 : length(segment) 
   
      coef = polyfit([segment(i).lx :segment(i).rx]',data(segment(i).lx :segment(i).rx),0);
      best = (0*( [segment(i).lx :segment(i).rx]' ))+coef(1);
      segment(i).y = best(1);
      segment(i).mean=coef(1);% slim
end;

residuals=[];
for i = 1 : length(segment)    
      coef = polyfit([segment(i).lx :segment(i).rx]',data(segment(i).lx :segment(i).rx),0);     
      best = (0*( [segment(i).lx :segment(i).rx]' ))+coef(1);
      residuals = [    residuals ; sum((data([segment(i).lx :segment(i).rx]')-best).^2)];
end;

if nargin == 3

hold on;
plot(data+0,'r');
temp = [];
for i = 1 : length(segment) 
   
      coef = polyfit([segment(i).lx :segment(i).rx]',data(segment(i).lx :segment(i).rx),0);
      best = (0*( [segment(i).lx :segment(i).rx]' ))+coef(1);
      segment(i).mean=coef(1);% slim
      plot([segment(i).lx :segment(i).rx]', best,'b');
      temp = [temp; [best(1) best(end)]]; 
end;

[dMoyenne, indice ] = max(diff([segment(:).mean]))%slim
segment(indice).rx

for i = 1 : length(segment)  - 1 
        plot([segment(i).rx :segment(i+1).lx]', [ temp(i,2) temp(i+1,1)  ],'g');
end;
end

end

end













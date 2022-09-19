% ce programme permet de segmenter une série de points de mesures
% suivant la méthode de la fenêtre glissante (regression), et représente par la suite la
% dynamique de l'évolution en appliquant une classification floue des segments.

%%%%%%%%%%%%%%%%%%%%%%%%%
% .: Data generation :. %
%%%%%%%%%%%%%%%%%%%%%%%%%
N=1000;  % Number of samples
nl=2;    % Number of latent variables
nx=6;    % Number of observable variables

[tt,x,latgen]=datagen(N,nl,nx);

r1=x(:,nx);

%%%%%%%%%%%%%%%%%%%%%%%%
%load x                 %
%t=1:length(x);         %
%tt=t';                 %
%r1=x(1,:)';            %
%%%%%%%%%%%%%%%%%%%%%%%%

zx=[];

N1=80;
N2=10;
alpha=2;
qw=N1;
c=1;

while qw < length(r1)
    
    r=[];tend=[];st=[];

    
    x1=r1(qw-N1+1:qw);
    xx1=x1;
    %x1=r1(1:qw);
    t=qw-N1+1:qw;
    %t=1:qw;
    %%%%%%%%%%%%%%%%%%%%%%%%%%% normalization MaxMin
    Max=max(x1);Maxi=Max;Min=min(x1);Mini=Min;
    cst1=(Maxi+Mini)/2;cst2=(Maxi-Mini)/2;
    for i=1:length(x1)
       x1(i)=(x1(i)-cst1)/cst2;
    end
    x=[t' x1];
    %%%%%%%%%%%%%%%%%%%%%%FIN NORMALISATION MaxMin
    
    %%%%%%%%%%%%%%%%%%%%%%%% Z NORMALISATION
    x1 = (x1 - mean(x1))/std(x1);
    zx=[t' x1];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIN ZNORM
    
    
    %%%%%%%%%%%%% seg MaxMin
    [xseg,xer]=segmentation(x,1);
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%% seg Znorm
    [zxseg,zxer]=segmentation(zx,2);
    %%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calvelo
    seuil_std=80;
    st=[st,std(xx1)];
    z=corrcoef(xx1,t);
    r=[r,z(1,2)];
    for i=1:length(r)
        if st(i)>seuil_std,
        tend=[tend;i,r(i) st(i) 3];
        else
            if (-0.8<r(i))&(r(i)<0.8),
                tend=[tend;i r(i) st(i) 0];

            elseif r(i)<=-0.8

                tend=[tend;i r(i) st(i) -1];
            elseif r(i)>=0.8,

                tend=[tend; i r(i) st(i) 1];
            end
        end
    end
    %if size(tend,1)>1
    %    for w=1:length(tend)-1,
    %        if tend(w+1,4)~=tend(w,4)
    %            rupture_derivee=[rupture_derivee;tend(w+1,1)-floor(N1/2)];
    %        end
    %    end
    %end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tend
    
    subplot(311)
    plot(tt,r1)
    ylabel('res1')
    xlabel('samples')

    subplot(312)
    plot(x(:,1),x(:,2),':b','linewidth',1.25);
    ylabel('MaxMin Normal')
    hold on
    for k=1:2:length(xseg)
        line([xseg(k,1) xseg(k+1,1)],[xseg(k,2) xseg(k+1,2)],'Color','r','LineWidth',1.6)
    end

    subplot(313)
    plot(zx(:,1),zx(:,2),':r','linewidth',1.25);
    ylabel('Z normalized data')
    xlabel('Samples in Sliding Window')
    hold on
    for k=1:2:length(zxseg)
        line([zxseg(k,1) zxseg(k+1,1)],[zxseg(k,2) zxseg(k+1,2)],'Color','r','LineWidth',1.6)
    end

    
    %subplot(414)
    %plot(t,xx1,':b','linewidth',1.25)
    %hold on
    %plot(t(rupture_derivee,1),xx1(rupture_derivee,2),'r.','markersize',20)
    %plot(t,tend(:,4)*ones(1,N1),'r','linewidth',2)
    %ylabel('Colveol')

    
    
    
    
    %d=mat_info(j_fin);
    %am_du=readfis('ampl_dure');
    %cin=size(d,2);s=[];tt=[];
    %figure
    %for e=1:cin
    %   sortie=evalfis([d(5,e) d(6,e)],am_du);
    %   s=[s sortie];
    %   tt=(d(2,e)+d(1,e))/2;
    %   plot([d(1,e) d(2,e)],[sortie sortie],'m','linewidth',2.5);hold on;
    %   plot(tt,sortie,'or');
    %end
    %plot([v(1) v(2)],[0.12 0.12],':b');
    %plot([v(1) v(2)],[0.38 0.38],':b');
    %plot([v(1) v(2)],[0.62 0.62],':b');
    %plot([v(1) v(2)],[0.12 0.12],':b');
    %plot([v(1) v(2)],[0.88 0.88],':b');

    pause
    close all
    clear x1 t j jj j_fin d r st tend
    qw=qw+N2;
end%while







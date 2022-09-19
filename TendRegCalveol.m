d=80; %longueur de la fenêtre d''analyse
seuil_std=2; %seuil sur l''ecart type

N=1000;  % Number of samples
nl=2;    % Number of latent variables
nx=6;    % Number of observable variables

[tt,x,latgen]=datagen(N,nl,nx);


X=[tt x(:,nx)];
    for i=2,
    X(:,i)=X(:,i)-mean(X(1:10,i));
    end

    clear s r
temps=1:d+1;
for i=d+1:length(X(:,2)),
    st(i)=std(X(i-d:i,2));
    z=corrcoef(X(i-d:i,2),temps');
    r(i)=z(1,2);
end
    
tend=[];
for i=1:length(r),
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

rupture_derivee=[];
for i=1:length(tend)-1,
    if tend(i+1,4)~=tend(i,4)
        rupture_derivee=[rupture_derivee;tend(i+1,1)-floor(d/2)];
    end
end

figure(3)
clf
subplot(2,1,1)
plot(X(:,1),X(:,2),'linewidth',2)
axis tight
hold
plot(X(rupture_derivee,1),X(rupture_derivee,2),'r.','markersize',20)
grid on
axis tight
subplot(2,1,2)
plot(tend(:,1)-floor(d/2),tend(:,4),'r','linewidth',2)
grid on
axis tight
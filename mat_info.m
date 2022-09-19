function d=mat_info(x)

[l,c]=size(x);
g=1;
for i=1:2:l
   d(1,g)=x(i,1);%temps du 1 point
   d(2,g)=x(i+1,1);%temps du 2 point
   d(3,g)=x(i,2);%mesure du 1 point
   d(4,g)=x(i+1,2);%mesure du 2 point
   d(5,g)=abs(x(i+1,2)-x(i,2));%l'amplitude
   d(6,g)=x(i+1,1)-x(i,1);%la duree
   g=g+1;
end

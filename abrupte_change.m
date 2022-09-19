function [ecartt,instant,signature,dd]=abrupte_change(d)

[l,c]=size(d);
signature=[0 0 0];
%max_ecart=0.39;
max_ecart=0;
j=1;dd=[];ecartt=[];
% éliminer les points inutiles
for i=1:c
   if d(1,i)~=d(2,i) & d(6,i)~=0 & d(7,i)~=0.5
      dd(1,j)=d(1,i);
      dd(2,j)=d(2,i);
      dd(3,j)=d(3,i);
      dd(4,j)=d(4,i);
      dd(5,j)=d(5,i);
      dd(6,j)=d(6,i);
      dd(7,j)=d(7,i);%coeff. de dynamique
      j=j+1;
   end
end
le=size(dd,2);
nom=char('mat_info');
fichier(dd,dd,nom);
% recherche d'un éventuel abrupt change
for k=2:le
   ecart=abs(dd(7,k)-dd(7,k-1));
   ecartt=[ecartt ecart];
   if ecart>=max_ecart
      max_ecart=ecart;
      instant=dd(2,k-1);
      signature(1)=sign(dd(4,k+1)-dd(4,k));
      signature(2)=sign(dd(4,k+1)-dd(3,k+1));
      %signature(3)=sign(((dd(4,k+2)-dd(3,k+2))/(dd(2,k+2)-dd(1,k+2)))-((dd(4,k+1)-dd(3,k+1))/(dd(2,k+1)-dd(1,k+1))));
      %signature(3)=sign((dd(4,k+2)-dd(3,k+2))-(dd(4,k+1)-dd(3,k+1)));
   end
end
% calcul de la signature pour une surveillance progressive
%ind=find(dd(1,:)==max_ecart);
%while dd(7,ind) > 0.12
%   
%end
% calcul de la signature pour chaque segment de la courbe
j=1;
for k=1:le-2
   sig_entier(j,1)=dd(2,k);%instant
   %if abs(dd(4,k+1)-dd(4,k))<0.07
   if dd(7,k+1)<=0.2
      sig_entier(j,2)=0;
   else
      sig_entier(j,2)=sign(dd(4,k+1)-dd(4,k));%amplitude
   end
   
   sig_entier(j,3)=sign(dd(4,k+1)-dd(3,k+1));%pente
   %sig_entier(j,4)=sign((dd(4,k+2)-dd(3,k+2))-(dd(4,k+1)-dd(3,k+1)));%courbure
   j=j+1;
end


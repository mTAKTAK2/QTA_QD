function d=logic_floue(dur_max2,dur_max1,dur_min2,dur_min1,dur_moy,am_max2,am_max1,am_min2,am_min1,am_moy,d,vv)

%range1=[nn mm];
%range2=[n m];

range1=[am_min2 am_max2];
range2=[dur_min2 dur_max2];

M=[am_min2 am_min2 (am_min1+am_moy)/2 am_moy];
N=[(am_min1+am_moy)/2 am_moy (am_max1+am_moy)/2];
P=[am_moy (am_max1+am_moy)/2 am_max2 am_max2];

%M=[nn nn (nn+oo+(pp-oo))/2 oo+(pp-oo)];
%N=[(nn+oo+(pp-oo))/2 oo+(pp-oo) (oo+(pp-oo)+mm)/2];
%P=[oo+(pp-oo) (oo+(pp-oo)+mm)/2 mm mm];

Mm=[dur_min2 dur_min2 (dur_min1+dur_moy)/2 dur_moy];
Nn=[(dur_min1+dur_moy)/2 dur_moy (dur_max1+dur_moy)/2];
Pp=[dur_moy (dur_max1+dur_moy)/2 dur_max2 dur_max2];

%Mm=[n n (n+o+(p-o))/2 o+(p-o)];
%Nn=[(n+o+(p-o))/2 o+(p-o) (o+(p-o)+m)/2];
%Pp=[o+(p-o) (o+(p-o)+m)/2 m m];
%pause
a=newfis('variation_amp_dur','mamdani','prod','probor','min','max','centroid');
a=addvar(a,'input','amplitude',range1);
a=addmf(a,'input',1,'small','trapmf',M);
a=addmf(a,'input',1,'mean','trimf',N);
a=addmf(a,'input',1,'large','trapmf',P);
a=addvar(a,'input','duree',range2);
a=addmf(a,'input',2,'small','trapmf',Mm);
a=addmf(a,'input',2,'mean','trimf',Nn);
a=addmf(a,'input',2,'large','trapmf',Pp);
a=addvar(a,'output','variation',[0 1]);
a=addmf(a,'output',1,'lente','trapmf',[-2 -1 0 0.25]);
a=addmf(a,'output',1,'assez_lente','trimf',[0 0.25 0.5]);
a=addmf(a,'output',1,'moyenne','trimf',[0.25 0.5 0.75]);
a=addmf(a,'output',1,'assez_rapide','trimf',[0.5 0.75 1]);
a=addmf(a,'output',1,'rapide','trapmf',[0.75 1 1.25 2]);
%rulelist=[1 1 3 1 1;1 2 2 1 1;1 3 1 1 1;2 1 4 1 1;2 2 3 1 1;2 3 2 1 1;3 1 5 1 1;3 2 4 1 1;3 3 3 1 1]; 
%a=addrule(a,rulelist);
%rule1='If (amplitude is small) and (duree is small) then (variation is moyenne)    ';    
rule1='If (amplitude is small) and (duree is small) then (variation is assez_lente)';
rule2='If (amplitude is small) and (duree is mean) then (variation is assez_lente) '; 
rule3='If (amplitude is small) and (duree is large) then (variation is lente)      ';      
rule4='If (amplitude is mean) and (duree is small) then (variation is assez_rapide)';
rule5='If (amplitude is mean) and (duree is mean) then (variation is moyenne)      ';       
rule6='If (amplitude is mean) and (duree is large) then (variation is assez_lente) '; 
rule7='If (amplitude is large) and (duree is small) then (variation is rapide)     ';      
rule8='If (amplitude is large) and (duree is mean) then (variation is assez_rapide)'; 
rule9='If (amplitude is large) and (duree is large) then (variation is moyenne)    '; 
rule=[rule1;rule2;rule3;rule4;rule5;rule6;rule7;rule8;rule9];
a=parsrule(a,rule);
showrule(a);
%pause
WRITEFIS(a,'variation_amp_dur.fis');
am_du=readfis('variation_amp_dur.fis');
[lin,cin]=size(d);s=[];ttt=[];cusum1=0;cusum2=[];mf_class=[];
for e=1:cin
   [sortie,IRR]=evalfis([d(5,e) d(6,e)],am_du);
   %sortie_before=[sortie_before ARR];
   %1er agrégation (ET)
   sor_first_agr=IRR(:,1).*IRR(:,2);
   %2eme agrégation (OU)
   class_moyenne1=sor_first_agr(1,1)+sor_first_agr(5,1)-(sor_first_agr(1,1)*sor_first_agr(5,1));
   class_moyenne2=sor_first_agr(9,1)+class_moyenne1(1,1)-(sor_first_agr(9,1)*class_moyenne1(1,1));
   class_ass_lente=sor_first_agr(2,1)+sor_first_agr(6,1)-(sor_first_agr(2,1)*sor_first_agr(6,1));
   class_ass_rapide=sor_first_agr(4,1)+sor_first_agr(8,1)-(sor_first_agr(4,1)*sor_first_agr(8,1));
   sor_final_agr=[sor_first_agr(3,1) class_ass_lente(1,1) class_moyenne2(1,1) class_ass_rapide(1,1) sor_first_agr(7,1)];
   mf_class=[mf_class;sor_final_agr];
   s=[s sortie];
   cusum1=cusum1+sortie;
   cusum2=[cusum2 cusum1];
   d(7,e)=sortie;% stocker le coeff. dans la 7eme ligne de la matrice d
   %tt=[tt (d(2,e)+d(1,e))/2];
   tt=(d(2,e)+d(1,e))/2;
   ttt=[ttt tt];
   %figure(3)
   %set(gcf,'position',[188   435   821   273])
   subplot(212)
   axis=[vv(1) vv(2) 0 1];
   %if d(6,e)~=0 & sortie~=0.5
      plot([d(1,e) d(2,e)],[sortie sortie],'m','linewidth',2.5);hold on;
      plot(tt,sortie,'or');
      %plot(tt,cusum2);
   %end 
end
v=axis;
plot([v(1) v(2)],[0.12 0.12],':b','linewidth',2.5);
plot([v(1) v(2)],[0.38 0.38],':b','linewidth',2.5);
plot([v(1) v(2)],[0.62 0.62],':b','linewidth',2.5);
plot([v(1) v(2)],[0.12 0.12],':b','linewidth',2.5);
plot([v(1) v(2)],[0.88 0.88],':b','linewidth',2.5);
detec=[ttt' s'];
%variation maxi

max_var=max(detec(:,2));
ind=find(detec(:,2)==max_var);
if size(ind,1)>1
   for i=1:size(ind,1)
      plot(detec(ind(i,1),1),detec(ind(i,1),2),'*r');
      plot([detec(ind(i,1),1) detec(ind(i,1),1)],[0 v(4)],':k');
   end   
else
   plot(detec(ind,1),detec(ind,2),'*g');
   plot([detec(ind,1) detec(ind,1)],[0 v(4)],':k');
end

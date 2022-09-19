function [j,jj,er]=pwreg_bar(x,max_erreur)
% piecewise regression method of segmentation

n=length(x);
inc=1;
j=[];jj=[];er=[];
%j(1,1)=x(1,1);
%j(1,2)=x(1,2);
comp=1;
while inc < n
   i=1;
   %[erreur,j,jj]=cal_err(x,inc,inc+i,comp);
   while (inc+i<n) & ( cal_err(x,inc,inc+i)< max_erreur*(x(inc+i,1)-x(inc,1)))
   %while (inc+i<n) & ( cal_err(x,inc,inc+i)< max_erreur)
      i=i+1;
   end % while
   %comp=comp+1;
   %if inc+i<n
   app=polyfit(x(inc:inc+i,1),x(inc:inc+i,2),1);
   val=polyval(app,x(inc:inc+i,1));
   err=max(abs(val-x(inc:inc+i,2)));%calcul de l'erreur 
   j(comp,1)=x(inc,1);
   j(comp,2)=polyval(app,x(inc,1));
   er=[er;err];
   %comp=comp+1;
   jj(comp,1)=x(inc+i,1);
   jj(comp,2)=polyval(app,x(inc+i,1));
   %j(comp,1)=x(inc+(i-1),1);
   %j(comp,2)=x(inc+(i-1),2);
   comp=comp+1;
   inc=inc+i;
   %app=polyfit(x(inc:inc+1,1),x(inc:inc+1,2),1);
   %val=polyval(app,x(inc:inc+1,1));
   %err=max(abs(val-x(inc:inc+1,2)));%calcul de l'erreur 
   %er=[er;err];
end % while
% Completer la segmentation (ceci n'est valable qu'en travaille hors ligne)
app=polyfit(x(inc:n,1),x(inc:n,2),1);
val=polyval(app,x(inc:n,1));
err=max(abs(val-x(inc:n,2)));%max de l'erreur par segment 
j(comp,1)=x(inc,1);
j(comp,2)=polyval(app,x(inc,1));
er=[er;err];
jj(comp,1)=x(n,1);
jj(comp,2)=polyval(app,x(n,1));

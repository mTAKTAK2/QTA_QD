function [j_fin,er]=segmentation(data,max_erreur)

    x=data;
    %sigma=std(x1);
    %max_erreur=1.5*std(x1);
    [j,jj,er]=pwreg_bar(x,max_erreur);
    n_j=length(j);
    comp=1;k=1;
    while comp < n_j
       j_fin(k,1)=j(comp,1);
       j_fin(k,2)=j(comp,2);
       comp=comp+1;k=k+2;
    end
    comp=1;kk=2;
    while comp < n_j
       j_fin(kk,1)=jj(comp,1);
       j_fin(kk,2)=jj(comp,2);
       comp=comp+1;kk=kk+2;
    end
    delta=min(diff(x(:,1)));

    %%%%%%%%%%%%%%%%%%%%%
    % calcul de l'écart d'approximation
    compteur=x(1,1);comp=1;
    nx=length(j);er=[];
    for i=2:nx   								
       err_suiv=0;				
       xx(comp,1)=compteur;									
       xx(comp,2)=x(comp,2);						
       pent=(j(i,2)-j(i-1,2))/(j(i,1)-j(i-1,1));		
       erreur2=0;												
       for k=j(i-1,1)+1:j(i,1)				
          compteur=compteur+1;						
          comp=comp+1;										
          xx(comp,1)=compteur;		
          xx(comp,2)=pent*(xx(comp,1)-j(i-1,1))+j(i-1,2);
          %erreur2=erreur2+(abs(xx(comp,2)-x(comp,2)));
          erreur2=erreur2+(xx(comp,2)-x(comp,2))^2;
       end
       er=[er;erreur2];			
    end
    er=sqrt(sum(er)/length(er));

end

function erreur=cal_err(x,inc,incc)

app=polyfit(x(inc:incc,1),x(inc:incc,2),1);%approximer la courbe par une regression d'ordre 1 
                                           %(lineaire)
val=polyval(app,x(inc:incc,1));
err=abs(val-x(inc:incc,2));%calcul de l'erreur sur tout l'axe de temps compris entre (inc,incc)
erreur=err*(x(incc,1)-x(inc,1));

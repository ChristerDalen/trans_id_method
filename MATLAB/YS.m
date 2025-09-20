function [Km, tm, dm]=YS(Y,R,T,Kc)

[cp1,cp2,A,cm1,tm1,tp1]=id_6_param3(Y,R,T);
cinf=(cp1*cp2-cm1^2)/(cp1+cp2-2*cm1);
v=(cinf-cm1)/(cp1-cinf);
xi=-log(v)/sqrt(pi^2+log(v)^2);

Km=cinf/(Kc*(A-cinf));
K=Km*Kc;
s1=xi*sqrt(K+1)+sqrt(xi^2*(K+1)+K-1); %the correct 
s1=real(s1);
s2=sqrt((1-xi^2)*(K+1));
dm=2*(tm1-tp1)*s2/(pi*s1);
tm=(tm1-tp1)*s1*s2/pi;
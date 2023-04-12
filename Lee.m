function [Km, tm, dm]=Lee(Y,R,T,Kc)

[cp1,cp2,A,cm1,tm1,tp1]=id_6_param3(Y,R,T);

cinf=(cp1*cp2-cm1^2)/(cp1+cp2-2*cm1);

v=(cinf-cm1)/(cp1-cinf);
xi=-log(v)/sqrt(pi^2+log(v)^2);
tau=(tm1-tp1)*sqrt(1-xi^2)/pi;
Km=cinf/(Kc*(A-cinf));
alpha=xi/tau;
beta=sqrt(1-xi^2)/tau;
v=atan(beta/alpha);

dm=(v+pi/4)/beta;
for i=1:10

dm=1/beta*(v+atan(beta*exp(-alpha*dm)/(Kc*Km*sqrt(alpha^2+beta^2)*cos(beta*dm-v))));

end

tm=1/alpha*(1+Kc*Km*exp(alpha*dm)*cos(beta*dm));


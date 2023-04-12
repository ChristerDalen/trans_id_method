function [Km, tm, dm]=Chen(Y,R,T,Kc)

global d xi tau

[cp1,cp2,A,cm1,tm1,tp1]=id_6_param3(Y,R,T);

cinf=(cp1*cp2-cm1^2)/(cp1+cp2-2*cm1);

Km=cinf/(Kc*(A-cinf));
K=cinf/A;
H=1/3*((cp1-cinf)/cinf+(cinf-cm1)/(cp1-cinf)+(cp2-cinf)/(cinf-cm1));
xi=-log(H)/sqrt(pi^2+(log(H)^2));
tau=(tm1-tp1)*sqrt(1-xi^2)/pi;
d=2*tp1-tm1;

wu=1;
wu=fsolve(@fun,wu,optimset('Display','off'));
Gcl=K/sqrt((1-tau^2*wu^2)^2+(2*xi*tau*wu)^2);
GcGp=Gcl/sqrt(1+2*Gcl+Gcl^2);
GM=1/GcGp;
Kcu=Kc*GM;
tm=1/wu*sqrt(Kcu^2*Km^2-1);
dm=1/wu*(pi-atan(tm*wu));
tm=real(tm); dm=real(dm);
function [Km, tm, dm]=JR(Y,R,T,Kc)

[cp1,cp2,A,cm1,tm1,tp1]=id_6_param3(Y,R,T);

cinf=(cp1*cp2-cm1^2)/(cp1+cp2-2*cm1);

v=(cinf-cm1)/(cp1-cinf);
xi=-log(v)/sqrt(pi^2+log(v)^2);
Km=cinf/(Kc*(A-cinf));
K=Km*Kc;
tau=(tm1-tp1)*sqrt(1-xi^2)/pi;
gamma1=-0.6143; gamma2=0.1247; delta=0.3866;
alpha=2*xi*tau*(1+K)/(delta+gamma1*K);
beta=-inv(delta+gamma1*K);

A1=beta^2*gamma2*K+beta*delta;
B1=2*gamma2*K*alpha*beta+alpha*delta;
C1=gamma2*K*alpha^2-tau^2*(1+K);

tm=(-B1+sqrt(B1^2-4*A1*C1))/(2*A1);
dm=alpha+beta*tm;
function [Km, tm, dm]=MF(Y,R,T,Kc,Ti,fc)

global d xi tau

[cp1,cp2,A,cm1,tm1,tp1]=id_6_param3(Y,R,T);

css=(cp1*cp2-cm1^2)/(cp1+cp2-2*cm1);

K=css/A;
p=-1/(2*pi)*log((cp2-css)/(cp1-css));
xi=sqrt(p^2/(1+p^2));
tau=(tp2-tp1)*sqrt(1-xi^2)/(2*pi);

Sc=css*ones(N,1)-Y; Sc=sum(Sc,1)*0.001;
d=Sc/css-2*xi*tau;
wu=0;
wc=fsolve(@fun,wu,optimset('Display','off'));
Kp=Ti/(Kc*Sc)*css;
M=K/sqrt((1-tau^2*wc^2)^2+(2*tau*xi*wc)^2);
tp=sqrt((Kc*Kp)^2*(1+Ti^2*wc^2)-M^2*Ti^2*wc^2)/(M*wc^2*Ti);
dp=1/wc*(atan(wc*Ti)+atan(1/(tau*wc)));
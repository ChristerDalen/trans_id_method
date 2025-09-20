function [kp, tau, to]=FSL(Y,R,T,kc,fc)


%if(fc~=-1)
%    [b,a]=butter(1, fc);
%    Y=filtfilt(b,a,Y); 
%end

[yp1,yp2,dR,ym1,tm1,tp1,tp2]=id_6_param3(Y,R,T);
yinf=(yp1*yp2-ym1^2)/(yp1+yp2-2*ym1);
ir = find(Y==yinf); tr=T(ir);
tr=1.694;
y0=0; OS=(yp1-yinf)/(yinf-y0);
wd=2*pi/(tp2-tp1);
K=(yinf-y0)/dR;
phi=(tp2-tr)*wd; xi=cos(phi);
wn=wd/sin(phi);
a=2*xi*wn; b=wn^2;
kp=K/(kc*(1-K));
tau=(-(-(1+kc*kp)*a)+sqrt((-(1+kc*kp)*a)^2-4*b*(1-kc^2*kp^2)))/(2*b);
A=-(1+kc*kp)/tau;
to=-2*A/b;
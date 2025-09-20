function [dyp,dyp2,dys,dyu,tu,tp,tp2]=id_6_param3(Y,R,T)

% time vector formulation
[dyp,ip]=max(Y);
tp = T(ip);
dys=max(R);
[dyu,iu]=min(Y(ip:end));
tu = T(ip+iu-1);
[dyp2,ip2]=max(Y(iu+iu-1:end));
tp2=T(ip2+2*iu-1);
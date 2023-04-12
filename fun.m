function f=fun(wu)
global d xi tau

f=-d*wu-atan((2*xi*tau*wu)/sqrt(1-tau^2*wu^2))+pi;

end


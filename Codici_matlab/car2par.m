function [a, e, i, OM, om, th] = car2par(rr, vv, mu)
% passaggio da cartesiane a parametriche: conoscendo il vettore posizione,
% il vettore velocità e la costante gravitazionale ricaviamo i parametri
% orbitali

r = norm(rr);
v = norm(vv);

a = ((2/r) - (v^2/mu))^(-1); % semiasse maggiore

hh = cross(rr,vv);
h = norm(hh);

ee = cross(vv,hh)/mu - rr./r;
e = norm(ee); % eccentricità

i = acos(hh(3)/h); % inclinazione

k = [0, 0, 1]';
N = cross(k,hh)/norm(cross(k,hh));
if (N(2) < 0)
    OM = 2*pi - acos(N(1)); % ascensione retta del nodo ascendente
else
    OM = acos(N(1));
end

if (ee(3) < 0)
    om = 2*pi - acos(dot(N,ee)/e); % anomalia del pericentro
else
    om = acos(dot(N,ee)/e);
end

vr = (vv'*rr)/r;
if (vr < 0)
    th = 2*pi - acos (dot(rr,ee)/(r*e)); % coordinata angolare
else
    th = acos (dot(rr,ee)/(r*e));
end

end
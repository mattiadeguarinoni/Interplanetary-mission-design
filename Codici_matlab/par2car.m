function [rr, vv] = par2car(a, e, i, OM, om, th, mu)
% passaggio da parametriche a cartesiane: conoscendo i parametri orbitali e
% la costante gravitazionale ricaviamo il vettore posizione e il vettore
% velocità

a = a(1); e = e(1); i = i(1); OM = OM(1); om = om(1); th = th(1);

R_OM = [cos(OM), sin(OM), 0; -sin(OM), cos(OM), 0; 0, 0, 1];
R_i = [1, 0, 0; 0, cos(i), sin(i); 0, -sin(i), cos(i)];
R_om = [cos(om), sin(om), 0; -sin(om), cos(om), 0; 0, 0, 1];

T = R_om*R_i*R_OM; % matrice di rotazione che ci consente di passare
% dal sistema geocentrico equatoriale al sistema perifocale

p = a*(1 - e^2);
r = p/(1 + e*cos(th));

% calcoliamo il vettore di stato nel sistema perifocale
r_pf = [r*cos(th), r*sin(th), 0]';
v_pf = [-sqrt(mu/p)*sin(th), sqrt(mu/p)*(e + cos(th)), 0]';

% passiamo al sistema geocentrico equatoriale
rr = T'*r_pf;
vv = T'*v_pf;

rr = rr(:)';
vv = vv(:)';

end
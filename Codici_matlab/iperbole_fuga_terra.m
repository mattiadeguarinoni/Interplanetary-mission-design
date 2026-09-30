%% iperbole di fuga dalla Terra

clear
clc

mu_t = 398600.4;
mu_s = 132712440018;

% parametri dell'orbita di trasferimento
a_trasf = 1.8101*10^(8);
e_trasf = 0.1648;
i = 0.0657;
OM_trasf = 2.4247;
om_trasf = 3.2638;
th_trasf = mod(-0.1227, 2*pi);

% parametri dell'orbita eliocentrica della Terra
a_terra = 1.4946 * 10^(8);
e_terra = 0.016;
i_terra = 9.1920*10^(-5);
OM_terra = 2.7847;
om_terra = 5.2643;
th_terra = 3.8; % trovato da scenario 2 come val_theta_i

% parametri dell'orbita di parcheggio intorno alla Terra
r_p = [-6017.4706  -5030.3532  4564.9360]';
v_p = [4.671  -6.799  -1.335]';
[a_p, e_p, i_p, OM_p, om_p, ~] = car2par(r_p, v_p, mu_t);

[r_terra, v_terra] = par2car(a_terra, e_terra, i_terra, OM_terra, om_terra, th_terra, mu_s);
[r_t, v_t] = par2car(a_trasf, e_trasf, i, OM_trasf, om_trasf, th_trasf, mu_s);

v_inf = norm(v_t-v_terra);
a_h = -mu_t / (v_inf^2);

% caratterizzazione dell'orbita di parcheggio intorno alla Terra
% posizione e velocità del pericentro
[r_p, v_p] = par2car(a_p, e_p, i, OM_p, om_p, 0, mu_t);

rp_h = norm(r_p);
r_oc = rp_h;
v_oc = sqrt(mu_t/r_oc);

% eccentricità
e_h = 1 - (rp_h/a_h);

% velocità al pericentro dell'iperbole
v_fuga = sqrt(2*mu_t/rp_h);
v_ph = sqrt(v_inf^2 + v_fuga^2);

% parametro di impatto
imp = -a_h*sqrt(e_h^2 - 1);

% angolo di deflessione
delta_h = 2 * asin(1/e_h);

% deltaV di arrivo
delta_V = abs(v_oc - v_ph);

% tempo di trasferimento
m_terra = 5.972*10^(24);
m_sole = 1.989*10^(30);
r_SOI = a_terra*((m_terra/m_sole)^(2/5));
p = abs(a_h)*(e_h^2 - 1);
theta1 = 0;
theta2 = acos((p/r_SOI -1)/e_h);

F1 = 2*atanh(sqrt((e_h-1)/(e_h+1))*tan(theta1/2));
F2 = 2*atanh(sqrt((e_h-1)/(e_h+1))*tan(theta2/2));

delta_t = sqrt(abs(a_h)^3/mu_t)*((e_h*sinh(F2)-F2)-(e_h*sinh(F1)-F1));
%% iperbole di arrivo all'asteroide

clear
clc

M_a = 1.5067*10^(11);
d_a = 0.52;
G = 6.674*10^(-20);
mu_a = G*M_a;
mu_s = 132712440018;
Delta_vec = [];
Delta_v_vec = [];

% parametri dell'orbita di trasferimento
a_t = 1.8101*10^(8);
e_t = 0.1648;
i = 0.0657;
Om_t = 2.4247;
om_t = 3.2638;
th_t = 2.9492;

% parametri dell'orbita eliocentrica dell'asteroide
a_a  = 165228155.43;
e_a  = 0.276277;
i_a = deg2rad(9.85);
Om_a = deg2rad(136.43);
om_a = deg2rad(186.57);
th_a = 3; % trovato da scenario 2 come val_theta_f

[r_a, v_a] = par2car(a_a, e_a, i_a, Om_a, om_a, th_a, mu_s);
[r_t, v_t] = par2car(a_t, e_t, i, Om_t, om_t, th_t, mu_s);

v_inf = norm(v_a - v_t);
a_h = -mu_a / (v_inf^2);

minimo = inf;

for rp_h = d_a/2 + [0, 0.25, 0.5, 0.750, 1]

    % caratterizzazione dell'orbita di arrivo intorno all'asteroide
    % posizione e velocità
    r_oc = rp_h;
    v_oc = sqrt(mu_a/r_oc);

    % eccentricità
    e_h = 1 - (rp_h/a_h);

    % velocità al pericentro dell'iperbole
    v_fuga = sqrt(2*mu_a/rp_h);
    v_ph = sqrt(v_inf^2 + v_fuga^2);

    % parametro di impatto
    imp = -a_h*sqrt(e_h^2 - 1);

    Delta_vec = [Delta_vec imp];
    % angolo di deflessione
    delta_h = 2 * asin(1/e_h);

    % deltaV di arrivo
    delta_V = abs(v_oc - v_ph);

    Delta_v_vec = [Delta_v_vec delta_V];

    if delta_V < minimo
        minimo = delta_V;
        e_h_best = e_h;
        imp_best = imp;
        rp_h_best = rp_h;
        delta_h_best = delta_h;
    end

    if rp_h == d_a/2 + 0.5
        delta_V_val = delta_V;
        e_h_val = e_h;
        imp_val = imp;
        rp_h_val = rp_h;
        delta_h_val = delta_h;
    end

end

rp_h_val
a_h
e_h_val
imp_val
delta_h_val
delta_V_val

% tempo di trasferimento
m_a = 1.5067*10^(11);
m_sole = 1.989*10^(30);
r_SOI = a_a*((m_a/m_sole)^(2/5));
p = abs(a_h) * (e_h_val^2 - 1);
theta1 = mod(acos((p/r_SOI - 1)/e_h_val), 2*pi);
theta2 = 2*pi;

F1 = mod(2*atanh(sqrt((e_h_val-1)/(e_h_val+1)) * tan(theta1/2)), 2*pi);
F2 = mod(2*atanh(sqrt((e_h_val-1)/(e_h_val+1)) * tan(theta2/2)), 2*pi);

delta_t = sqrt(abs(a_h)^3/mu_a)*((e_h_val*sinh(F2)-F2)-(e_h_val*sinh(F1)-F1));

delta_t

%%
r_p_h = d_a/2 + [0, 0.25, 0.5, 0.750, 1];
figure('Name','parametro di impatto')
grid on
hold on
plot(r_p_h-d_a/2,Delta_vec,'o','Markerfacecolor', 'b', 'MarkerSize',5);
xlabel('r_p-R_m [km]')
ylabel('\Delta [km]')
title('Parametro di impatto \Delta')

%%
figure('Name','DeltaV')
grid on 
hold on
plot(r_p_h-d_a/1,Delta_v_vec,'o','Markerfacecolor','b','MarkerSize',5);
xlabel('r_p-R_m [km]')
ylabel('\DeltaV [km/s]')
title('Costo di cattura \Delta V')
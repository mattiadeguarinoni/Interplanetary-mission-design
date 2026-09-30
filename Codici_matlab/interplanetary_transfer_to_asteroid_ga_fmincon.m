%% interplanetary transfer to asteroid ga fmincon

%% Trasferimento interplanetario con ga ed fmincon

clear
clc

mu_s = 132712440018;

% Orbita iniziale (Terra)
a_i  = 1.4946*10^(8);
e_i  = 0.016;
i_i  = 9.1920*10^(-5);
OM_i = 2.7847;
om_i = 5.2643;

% Orbita finale
a_f  = 165228155.43;
e_f  = 0.276277;
i_f  = 9.85; % gradi
Om_f = 136.43; % gradi
om_f = 186.57; % gradi

% Funzione obiettivo
fitnessFcn = @(x) objective_transfer(x, a_i, e_i, i_i, OM_i, om_i, ...
                                        a_f, e_f, i_f, Om_f, om_f, mu_s);

% Numero di variabili [theta_i, theta_f, om_t]
nvars = 3;

% limiti [theta_i, theta_f, om_t]
lb = [0,    0,    0   ];
ub = [2*pi, 2*pi, 2*pi];

options = optimoptions('ga', ...
    'PopulationSize',    300, ...   
    'MaxGenerations',    500, ...   
    'EliteCount',        15,  ...   
    'CrossoverFraction', 0.8, ...   
    'FunctionTolerance', 1e-6, ...  
    'PlotFcn',           @gaplotbestf, ...  
    'Display',           'iter');   

[x_opt, DeltaV_min, exitflag, output] = ga(fitnessFcn, nvars, ...
                                            [], [], [], [], ...
                                            lb, ub, [], options);

% Risultati
fprintf('\n===== RISULTATI =====\n');
fprintf('theta_i ottimale : %.4f rad (%.2f°)\n', x_opt(1), rad2deg(x_opt(1)));
fprintf('theta_f ottimale : %.4f rad (%.2f°)\n', x_opt(2), rad2deg(x_opt(2)));
fprintf('om_t ottimale    : %.4f rad (%.2f°)\n', x_opt(3), rad2deg(x_opt(3)));
fprintf('DeltaV minimo    : %.4f km/s\n', DeltaV_min);
fprintf('Generazioni      : %d\n', output.generations);

% TOF con la soluzione ottimale
[r_i_opt, v_i_opt] = par2car(a_i, e_i, i_i, OM_i, om_i, x_opt(1), mu_s);
[r_f_opt, v_f_opt] = par2car(a_f, e_f, deg2rad(i_f), deg2rad(Om_f), ...
                              deg2rad(om_f), x_opt(2), mu_s);

theta1_t_opt = atan2(r_i_opt(2)/norm(r_i_opt), r_i_opt(1)/norm(r_i_opt));
theta2_t_opt = atan2(r_f_opt(2)/norm(r_f_opt), r_f_opt(1)/norm(r_f_opt));
r_i_opt_n = norm(r_i_opt);
r_f_opt_n = norm(r_f_opt);
e_t_opt = (1 - r_i_opt_n/r_f_opt_n) / (r_i_opt_n/r_f_opt_n * cos(theta1_t_opt) - cos(theta2_t_opt));
a_t_opt = r_i_opt_n * (1 + e_t_opt*cos(theta1_t_opt)) / (1 - e_t_opt^2);
TOF = TOF(a_t_opt, e_t_opt, theta1_t_opt, theta2_t_opt, mu_s);

fprintf('Optimal Time of flight   : %.4f s\n',TOF);

options_local = optimoptions('fmincon', 'Display', 'iter', ...
                              'Algorithm', 'sqp');
x_refined = fmincon(fitnessFcn, x_opt, [], [], [], [], lb, ub, ...
                    [], options_local);

fprintf('DeltaV raffinato: %.6f km/s\n', fitnessFcn(x_refined));


function DeltaV = objective_transfer(x, a_i, e_i, i_i, OM_i, om_i, ...
                                      a_f, e_f, i_f, Om_f, om_f, mu_s)

theta_i = x(1);
theta_f = x(2);
om_t    = x(3);

% Costanti
R_sole  = 696700;
margine = 45000000;

% Orbita iniziale
[r_i, v_i] = par2car(a_i, e_i, i_i, OM_i, om_i, theta_i, mu_s);

% Orbita finale
[r_f, v_f] = par2car(a_f, e_f, deg2rad(i_f), deg2rad(Om_f), ...
                     deg2rad(om_f), theta_f, mu_s);

r1 = r_i';
r2 = r_f';


h_t  = cross(r1, r2) / norm(cross(r1, r2));
i_t  = acos(h_t(3));
k    = [0 0 1]';
N_t  = cross(k, h_t) / norm(cross(k, h_t));

if N_t(2) >= 0
    OM_t = acos(N_t(1));
else
    OM_t = 2*pi - acos(N_t(1));
end

R_OM = [cos(OM_t) sin(OM_t) 0; -sin(OM_t) cos(OM_t) 0; 0 0 1];
R_i = [1 0 0; 0 cos(i_t) sin(i_t); 0 -sin(i_t) cos(i_t)];
R_om = [cos(om_t)  sin(om_t) 0; -sin(om_t) cos(om_t) 0; 0 0 1];
T = R_om * R_i * R_OM;

r1_pf = T * r1;
r2_pf = T * r2;

theta1_t = atan2(r1_pf(2)/norm(r1), r1_pf(1)/norm(r1));
theta2_t = atan2(r2_pf(2)/norm(r2), r2_pf(1)/norm(r2));

r1_n = norm(r1);
r2_n = norm(r2);

e_t = (1 - r1_n/r2_n) / (r1_n/r2_n * cos(theta1_t) - cos(theta2_t));
a_t = r1_n * (1 + e_t*cos(theta1_t)) / (1 - e_t^2);
rp_t = a_t * (1 - e_t);

if e_t <= 0 || e_t >= 1 || rp_t <= (R_sole + margine)
    DeltaV = 1e6;
    return
end

% Calcolo DeltaV
[~, vv1_t] = par2car(a_t, e_t, i_t, OM_t, om_t, theta1_t, mu_s);
[~, vv2_t] = par2car(a_t, e_t, i_t, OM_t, om_t, theta2_t, mu_s);

DeltaV = norm(vv1_t - v_i) + norm(v_f - vv2_t);
end


%% Minimizzazione TOF con GA + fmincon

clear; 
clc;

mu_s = 132712440018;

% Orbita iniziale (Terra)
a_i  = 1.4946e8;
e_i  = 0.016;
i_i  = 9.1920e-5;
OM_i = 2.7847;
om_i = 5.2643;

% Orbita finale (Asteroide)
a_f  = 165228155.43;
e_f  = 0.276277;
i_f  = 9.85;
Om_f = 136.43;
om_f = 186.57;


fitnessFcn_TOF = @(x) objective_TOF(x, a_i, e_i, i_i, OM_i, om_i, ...
                                       a_f, e_f, i_f, Om_f, om_f, mu_s);

nvars = 3;
lb = [0,    0,    0   ];
ub = [2*pi, 2*pi, 2*pi];

options_ga = optimoptions('ga', ...
    'PopulationSize',    300, ...
    'MaxGenerations',    500, ...
    'EliteCount',        15,  ...
    'CrossoverFraction', 0.8, ...
    'FunctionTolerance', 1e-6, ...
    'PlotFcn',           @gaplotbestf, ...
    'Display',           'iter');

% GA minimizza TOF
[x_opt_TOF, TOF_min, exitflag, output] = ga(fitnessFcn_TOF, nvars, ...
                                             [], [], [], [], ...
                                             lb, ub, [], options_ga);

fprintf('\n===== RISULTATI GA =====\n');
fprintf('theta_i : %.4f rad (%.2f deg)\n', x_opt_TOF(1), rad2deg(x_opt_TOF(1)));
fprintf('theta_f : %.4f rad (%.2f deg)\n', x_opt_TOF(2), rad2deg(x_opt_TOF(2)));
fprintf('om_t    : %.4f rad (%.2f deg)\n', x_opt_TOF(3), rad2deg(x_opt_TOF(3)));
fprintf('TOF min : %.2f s (%.2f giorni)\n', TOF_min, TOF_min/86400);

options_local = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
x_refined_TOF = fmincon(fitnessFcn_TOF, x_opt_TOF, [], [], [], [], ...
                         lb, ub, [], options_local);

[TOF_finale, DeltaV_associato] = objective_TOF(x_refined_TOF, a_i, e_i, ...
                                  i_i, OM_i, om_i, a_f, e_f, i_f, Om_f, om_f, mu_s);

fprintf('\n===== RISULTATI FINALI =====\n');
fprintf('TOF ottimale     : %.2f s (%.2f giorni)\n', TOF_finale, TOF_finale/86400);
fprintf('DeltaV associato : %.4f km/s\n', DeltaV_associato);


% funzione obiettivo di TOF
function [TOF, DeltaV] = objective_TOF(x, a_i, e_i, i_i, OM_i, om_i, ...
                                        a_f, e_f, i_f, Om_f, om_f, mu_s)

theta_i = x(1);
theta_f = x(2);
om_t    = x(3);

R_sole  = 696700;
margine = 45000000;

[r_i, v_i] = par2car(a_i, e_i, i_i, OM_i, om_i, theta_i, mu_s);
[r_f, v_f] = par2car(a_f, e_f, deg2rad(i_f), deg2rad(Om_f), ...
                     deg2rad(om_f), theta_f, mu_s);

r1 = r_i';
r2 = r_f';

h_t = cross(r1, r2) / norm(cross(r1, r2));
i_t = acos(h_t(3));
k   = [0 0 1]';
N_t = cross(k, h_t) / norm(cross(k, h_t));

if N_t(2) >= 0
    OM_t = acos(N_t(1));
else
    OM_t = 2*pi - acos(N_t(1));
end

R_OM = [cos(OM_t)  sin(OM_t) 0; -sin(OM_t) cos(OM_t) 0; 0 0 1];
R_i  = [1 0 0; 0 cos(i_t) sin(i_t); 0 -sin(i_t) cos(i_t)];
R_om = [cos(om_t)  sin(om_t) 0; -sin(om_t) cos(om_t) 0; 0 0 1];
T    = R_om * R_i * R_OM;

r1_pf = T * r1;
r2_pf = T * r2;

theta1_t = atan2(r1_pf(2)/norm(r1), r1_pf(1)/norm(r1));
theta2_t = atan2(r2_pf(2)/norm(r2), r2_pf(1)/norm(r2));

r1_n = norm(r1);
r2_n = norm(r2);

e_t  = (1 - r1_n/r2_n) / (r1_n/r2_n * cos(theta1_t) - cos(theta2_t));
a_t  = r1_n * (1 + e_t*cos(theta1_t)) / (1 - e_t^2);
rp_t = a_t * (1 - e_t);

if e_t <= 0 || e_t >= 1 || rp_t <= (R_sole + margine)
    TOF    = 1e12;
    DeltaV = 1e6;
    return
end

TOF = TOF(a_t, e_t, theta1_t, theta2_t, mu_s);

[~, vv1_t] = par2car(a_t, e_t, i_t, OM_t, om_t, theta1_t, mu_s);
[~, vv2_t] = par2car(a_t, e_t, i_t, OM_t, om_t, theta2_t, mu_s);
DeltaV = norm(vv1_t - v_i) + norm(v_f - vv2_t);

end
%% reaching parking orbit - minimizziamo la velocità

% in questo scenario abbiamo progettato un trasferimento orbitale che ci
% permette di minimizzare la differenza di velocità da imporre --> se ci
% spostiamo su un'orbita di trasferimento molto grande e facciamo una o
% entrambe le manovre di cambio piano e cambio anomalia del pericentro
% mentre orbitiamo su di essa, allora la velocità sarà più bassa e le
% manovre saranno molto meno costose

%% manovre all'apocentro dell'orbita grande

clear
clc

% orbita iniziale
a_i = 24400;
e_i = 0.728300;
i_i = 0.104700;
OM_i = 1.040000;
om_i = 3.107000;
theta_i = 1.232000;
mu = 398600.4;

% orbita finale
r_f = [-6017.470600 -5030.353200 4564.936000];
v_f = [4.671000 -6.799000 -1.335000];
[a_f,e_f,i_f,OM_f,om_f,theta_f] = car2par(r_f,v_f,mu);

Delta_OM = OM_f-OM_i;
Delta_i = i_f-i_i;
% DeltaOM > 0; Deltai > 0

alpha = acos(cos(i_i)*cos(i_f) + sin(i_i)*sin(i_f)*cos(Delta_OM));
cos_ui = (-cos(i_f) + cos(alpha)*cos(i_i))/(sin(alpha)*sin(i_i));
cos_uf = (cos(i_i) - cos(alpha)*cos(i_f))/(sin(alpha)*sin(i_f));
sin_ui = sin(Delta_OM)/sin(alpha)*sin(i_f);
sin_uf = sin(Delta_OM)/sin(alpha)*sin(i_i);

u_i = atan2(sin_ui,cos_ui);
u_f = atan2(sin_uf,cos_uf);
theta_Man_cp_biellittico = pi;
om_i_Man_biellittico = u_i - theta_Man_cp_biellittico;
om_i_Man_biellittico = mod(om_i_Man_biellittico,2*pi);

[deltaV_cw_1, theta_i, theta_f] = changePericenterArg(a_i,e_i,om_i, ...
    om_i_Man_biellittico,mu);

%% manovre in un punto qualsiasi dell'orbita grande

clear
clc

% orbita iniziale
a_i = 24400;
e_i = 0.728300;
i_i = 0.104700;
OM_i = 1.040000;
om_i = 3.107000;
theta_i = 1.232000;
mu = 398600.4;

% orbita finale
r_f = [-6017.470600 -5030.353200 4564.936000]';
v_f = [4.671000 -6.799000 -1.335000]';
[a_f,e_f,i_f,OM_f,om_f,theta_f] = car2par(r_f,v_f,mu);
ra_t = 2000000;

Delta_t1 = TOF(a_i, e_i, theta_i, 0, mu);
r_p_t1 = a_i*(1-e_i);
r_a_t1 = ra_t;
r_p_t2 = a_f*(1-e_f); 
e_t1 = (ra_t - r_p_t1)/(ra_t + r_p_t1);
e_t2 = (ra_t - r_p_t2)/(ra_t + r_p_t2);

% primo impulso
a_t1 = (r_a_t1+r_p_t1)/2;
Delta_V1 = sqrt(mu)*sqrt(2/r_p_t1  - 1/a_t1) - sqrt(mu)*sqrt(2/r_p_t1 - 1/a_i);
Delta_t2 = TOF(a_t1, e_t1, 0, pi, mu);
Delta_V1;
Delta_t2;

% secondo impulso 
[DeltaV_cp,om_cp,theta_cp] = changeOrbitalPlane(a_t1,e_t1,i_i,OM_i,om_i,i_f,OM_f,mu);
theta_cp = mod(theta_cp, 2*pi);
DeltaV_cp;

[DeltaV_cw, theta_i_cw, theta_f_cw] = changePericenterArg(a_t1,e_t1,om_cp,om_f,mu);
theta_i_cw = mod(theta_i_cw,2*pi);
theta_f_cw = mod(theta_f_cw,2*pi);

DeltaV_cw;

Delta_t3 = TOF(a_t1, e_t1, theta_cp, theta_i_cw(1), mu);
Delta_t4 = TOF(a_t1, e_t1, theta_f_cw(1), pi, mu);

a_t2 = (r_a_t1 + r_p_t2)/2;
Delta_V2 = sqrt(mu)*sqrt(2/r_a_t1  - 1/a_t2) - sqrt(mu)*sqrt(2/r_a_t1 - 1/a_t1);

Delta_t5 = TOF(a_t2, e_t2, pi, 0, mu);
% terzo impulso
Delta_V3 = sqrt(mu)*sqrt(2/r_p_t2 - 1/a_f) - sqrt(mu)*sqrt(2/r_p_t2 - 1/a_t2);

% sommatorie
Delta_V_tot = abs(Delta_V1)+abs(DeltaV_cp)+abs(DeltaV_cw)+abs(Delta_V2) +abs(Delta_V3);
Delta_t_tot = Delta_t1+Delta_t2+Delta_t3+Delta_t4 + Delta_t5;

%% Cambio solo del piano su biellittico

clear
clc

% orbita iniziale
a_i = 24400;
e_i = 0.728300;
i_i = 0.104700;
OM_i = 1.040000;
om_i = 3.107000;
theta_i = 1.232000;
mu = 398600.4;

% orbita finale
r_f = [-6017.470600 -5030.353200 4564.936000];
v_f = [4.671000 -6.799000 -1.335000];
[a_f,e_f,i_f,OM_f,om_f,theta_f] = car2par(r_f,v_f,mu);

minimo = 10;
for ra_t = 35000 : 1 : 45000
Delta_t1 = TOF(a_i, e_i, theta_i, 0, mu);
r_p_t1 = a_i*(1-e_i);
r_a_t1 = ra_t;
r_p_t2 = a_f*(1-e_f); 
e_t1 = (ra_t - r_p_t1)/(ra_t + r_p_t1);
e_t2 = (ra_t - r_p_t2)/(ra_t + r_p_t2);

% primo impulso
a_t1 = (r_a_t1+r_p_t1)/2;
Delta_V1 = sqrt(mu)*sqrt(2/r_p_t1  - 1/a_t1) - sqrt(mu)*sqrt(2/r_p_t1 - 1/a_i);

% cambio piano e secondo impulso
[DeltaV_cp,om_cp,theta_cp] = changeOrbitalPlane(a_t1,e_t1,i_i,OM_i,om_i,i_f,OM_f,mu);
theta_cp = mod(theta_cp, 2*pi);
DeltaV_cp;
Delta_t2 = TOF(a_t1, e_t1, 0, theta_cp, mu);
Delta_t3 = TOF(a_t1, e_t1, theta_cp, pi, mu);

a_t2 = (r_a_t1 + r_p_t2)/2;
Delta_V2 = sqrt(mu)*sqrt(2/r_a_t1  - 1/a_t2) - sqrt(mu)*sqrt(2/r_a_t1 - 1/a_t1);

Delta_t4 = TOF(a_t2, e_t2, pi, 0, mu);

% terzo impulso
Delta_V3 = sqrt(mu)*sqrt(2/r_p_t2 - 1/a_f) - sqrt(mu)*sqrt(2/r_p_t2 - 1/a_t2);

% cambio argomento pericentro su orbita finale
[DeltaV_cw, theta_i_cw, theta_f_cw] = changePericenterArg(a_f,e_f,om_cp,om_f,mu);
theta_i_cw = mod(theta_i_cw,2*pi);
theta_f_cw = mod(theta_f_cw,2*pi);

Delta_t5 = TOF(a_f, e_f, pi, theta_i_cw(1), mu);
Delta_t6 = TOF(a_f, e_f, theta_f_cw(1), theta_f, mu);

DeltaV_cw;

% sommatorie
Delta_V_tot = abs(Delta_V1)+abs(DeltaV_cp)+abs(DeltaV_cw)+abs(Delta_V2) +abs(Delta_V3);
Delta_t_tot = Delta_t1+Delta_t2+Delta_t3+Delta_t4 + Delta_t5 + Delta_t6;

if Delta_V_tot < minimo
    minimo = Delta_V_tot;   
    optimal_ra_t = ra_t;    
end

end

minimo;
optimal_ra_t;

%% al posto del biellittico il bitangente

clear
clc

% orbita iniziale
a_i = 24400;
e_i = 0.728300;
i_i = 0.104700;
OM_i = 1.040000;
om_i = 3.107000;
theta_i = 1.232000;
mu = 398600.4;
minimo = 10;

% orbita finale
r_f = [-6017.470600 -5030.353200 4564.936000]';
v_f = [4.671000 -6.799000 -1.335000]';
[a_f,e_f,i_f,OM_f,om_f,theta_f] = car2par(r_f,v_f,mu);

for a_t = 25000 : 100 : 1500000
    for e_t = 0.1 : 0.001 : 0.99

% orbita grande
Delta_t0 = TOF(a_i,e_i,theta_i,pi,mu);
[Delta_V1, Delta_V2, Delta_t1] = bitangentTransfer(a_i,e_i,a_t,e_t,'ap',mu);

% cambio piano
[DeltaV_cp,om_cp,theta_cp] = changeOrbitalPlane(a_t,e_t,i_i,OM_i,om_i,i_f,OM_f,mu);
theta_cp = mod(theta_cp, 2*pi);
DeltaV_cp;
Delta_t2 = TOF(a_t,e_t,0,theta_cp,mu);

% cambio argomento pericentro
[DeltaV_cw, theta_i_cw, theta_f_cw] = changePericenterArg(a_t,e_t,om_cp,om_f,mu);
theta_i_cw = mod(theta_i_cw,2*pi);
theta_f_cw = mod(theta_f_cw,2*pi);

Delta_t3 = TOF(a_t,e_t,theta_cp,theta_i_cw(1),mu);
Delta_t4 = TOF(a_t,e_t,theta_f_cw(1),pi,mu);

% secondo bitangente per tornare sull'orbita finale
[Delta_V3, Delta_V4, Delta_t5] = bitangentTransfer(a_t,e_t,a_f,e_f,'ap',mu);

Delta_V_tot = abs(Delta_V1)+abs(Delta_V2)+abs(Delta_V3) + abs(DeltaV_cp)+...
                                       abs(DeltaV_cw) + abs(Delta_V4);
Delta_t_tot = Delta_t0 + Delta_t1 + Delta_t2 + Delta_t3 + Delta_t4 + Delta_t5;

if Delta_V_tot < minimo
    minimo = Delta_V_tot;
    optimal_a_t = a_t;
    optimal_e_t = e_t;
    tempo = Delta_t_tot;
end

    end
end

minimo
optimal_a_t
optimal_e_t
tempo
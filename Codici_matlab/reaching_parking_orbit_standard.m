%% reaching parking orbit - trasferimento standard

% in questo scenario abbiamo seguito il procedimento standard per il
% trasferimento orbitale

clear
clc

% costante gravitazionale
mu = 398600.4;

% anomalie vere
th0 = 0;
thf = 2*pi;
dth = 0.01;

% parametri orbitali dell'orbita iniziale
a_in = 24400.00;
e_in = 0.728300;
i_in = 0.104700;
OM_in = 1.040000;
om_in = 3.107000;
th_in = 1.232000;

% calcoliamo le coordinate cartesiane dell'orbita iniziale
[rr_in, vv_in] = par2car(a_in, e_in, i_in, OM_in, om_in, th_in, mu);

% plottiamo l'orbita iniziale
plotOrbit(a_in, e_in, i_in, OM_in, om_in, th0, thf, dth, mu)

% coordinate cartesiane dell'orbita finale
rr_f = [-6017.470600, -5030.353200, 4564.936000]';
vv_f = [4.671000, -6.799000, -1.335000]';

% calcoliamo i parametri orbitali dell'orbita finale
[a_f, e_f, i_f, OM_f, om_f, th_f] = car2par(rr_f, vv_f, mu);

% plottiamo l'orbita finale
hold on
plotOrbit(a_f, e_f, i_f, OM_f, om_f, th0, thf, dth, mu)

% cambio piano
[DeltaV_cp, om_2, th_cp] = changeOrbitalPlane(a_in, e_in, i_in, OM_in, om_in, i_f, OM_f, mu);
[Deltat_1] = TOF(a_in, e_in, th_in, th_cp, mu);

hold on
plotOrbit(a_in, e_in, i_f, OM_f, om_2, th0, thf, dth, mu)

% cambio anomalia del pericentro
% possiamo fare la manovra di cambio anomalia del pericentro in due punti
% diversi --> th_cp è circa 240 gradi, thi_cw(1) è circa 349 gradi,
% thi_cw(2) è circa 169 gradi --> scegliamo il punto di manovra più vicino
% per minimizzare il Deltat, ovvero thi_cw(1), mentre il DeltaV non è
% influenzato da questa scelta
[DeltaV_cw, thi_cw, thf_cw] = changePericenterArg(a_in, e_in, om_2, om_f, mu);
[Deltat_2] = TOF(a_in, e_in, th_cp, thi_cw(1), mu);

hold on
plotOrbit(a_in, e_in, i_f, OM_f, om_f, th0, thf, dth, mu)
legend('Orbita iniziale','Orbita finale','Orbita cambio piano','Orbita cambio anomalia percentro')

% trasferimento bitangente
% le due orbite hanno i vettori eccentricità paralleli e con lo stesso
% verso --> possiamo fare un trasferimento dal pericentro all'apocentro o
% dall'apocentro al pericentro --> thf_cw(1) è circa 11 gradi, arriviamo
% prima all'apocentro, facciamo un trasferimento ap --> questo tipo di
% trasferimento è più conveniente anche perchè il nostro th_f è 2pi, quindi
% arrivando al pericentro dell'orbita finale siamo già nel punto desiderato
% e il Deltat_5 sarà nullo --> sia in termini di tempo sia in termini di
% velocità è meno costoso fare un trasferimento bitangente dall'apocentro
% al pericentro
[Deltat_3] = TOF(a_in, e_in, thf_cw(1), pi, mu);
[DeltaV_pa1, DeltaV_pa2, Deltat_4] = bitangentTransfer(a_in, e_in, a_f, e_f, 'ap', mu);

% costo totale in termini di tempo
Deltat_tot = Deltat_1 + Deltat_2 + Deltat_3 + Deltat_4;

% costo totale in termini di velocità
DeltaV_tot = abs(DeltaV_cp) + abs(DeltaV_cw) + abs(DeltaV_pa1) + abs(DeltaV_pa2);
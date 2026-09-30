function [DeltaV, thi, thf] = changePericenterArg(a, e, omi, omf, mu)
% cambio dell'anomalia del pericentro: imponiamo la geometria delle orbite
% (semiasse maggiore ed eccentricità), l'anomalia del pericentro iniziale,
% l'anomalia del pericentro finale e la costante gravitazionale -->
% troviamo le coordinate angolari iniziali e finali dei due punti in cui
% fare manovra e il costo della manovra in termini di velocità

DeltaOm = omf - omi;

thi_1 = mod(DeltaOm/2, 2*pi);
thi_2 = mod(pi + DeltaOm/2, 2*pi);
thf_1 = mod(2*pi - DeltaOm/2, 2*pi);
thf_2 = mod(pi - DeltaOm/2, 2*pi);

thi = [thi_1, thi_2]';
thf = [thf_1, thf_2]';

p = a*(1 - e^2);
DeltaV = 2*sqrt(mu/p)*e*sin(DeltaOm/2);

end
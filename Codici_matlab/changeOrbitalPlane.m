function [DeltaV, om_f, th] = changeOrbitalPlane(a, e, i_i, OM_i, om_i, i_f, OM_f, mu)
% cambio di piano: imponiamo la geometria delle orbite, l'inclinazione
% iniziale e finale, l'ascensione retta iniziale e finale, l'anomalia del
% pericentro iniziale e la costante gravitazionale --> troviamo il costo
% della manovra in termini di velocità, l'anomalia del pericentro finale e
% la coordinata angolare

DeltaOM = OM_f - OM_i;
Deltai = i_f - i_i;

if (DeltaOM > 0 && Deltai > 0)

    alpha = acos(cos(i_i)*cos(i_f) + sin(i_i)*sin(i_f)*cos(DeltaOM));
    cos_ui = (-cos(i_f) + cos(alpha)*cos(i_i)) / (sin(alpha)*sin(i_i));
    cos_uf = (cos(i_i) - cos(alpha)*cos(i_f)) / (sin(alpha)*sin(i_f));
    sin_ui = (sin(DeltaOM)/sin(alpha))*sin(i_f);
    sin_uf = (sin(DeltaOM)/sin(alpha))*sin(i_i);
    u_i = atan2(sin_ui, cos_ui);
    u_f = atan2(sin_uf, cos_uf);

    th = mod(u_i - om_i, 2*pi);
    % th = u_i - om_i;
    om_f = mod(u_f - th, 2*pi);
    % om_f = u_f - th;
    p = a*(1 - e^2);
    DeltaV = 2*sqrt(mu/p)*(1 + e*cos(th))*sin(alpha/2);

elseif (DeltaOM > 0 && Deltai < 0)

    alpha = acos(cos(i_i)*cos(i_f) + sin(i_i)*sin(i_f)*cos(DeltaOM));
    cos_ui = (cos(i_f) - cos(alpha)*cos(i_i)) / (sin(alpha)*sin(i_i));
    cos_uf = (-cos(i_i) + cos(alpha)*cos(i_f)) / (sin(alpha)*sin(i_f));
    sin_ui = (sin(DeltaOM)/sin(alpha))*sin(i_f);
    sin_uf = (sin(DeltaOM)/sin(alpha))*sin(i_i);
    u_i = atan2(sin_ui, cos_ui);
    u_f = atan2(sin_uf, cos_uf);

    th = mod(2*pi - u_i - om_i, 2*pi);
    % th = 2*pi - u_i - om_i;
    om_f = mod(2*pi - u_f - th, 2*pi);
    % om_f = 2*pi - u_f - th;
    p = a*(1 - e^2);
    DeltaV = 2*sqrt(mu/p)*(1 + e*cos(th))*sin(alpha/2);

elseif (DeltaOM < 0 && Deltai > 0)

    alpha = acos(cos(i_i)*cos(i_f) + sin(i_i)*sin(i_f)*cos(abs(DeltaOM)));
    cos_ui = (-cos(i_f) + cos(alpha)*cos(i_i)) / (sin(alpha)*sin(i_i));
    cos_uf = (cos(i_i) - cos(alpha)*cos(i_f)) / (sin(alpha)*sin(i_f));
    sin_ui = (sin(abs(DeltaOM))/sin(alpha))*sin(i_f);
    sin_uf = (sin(abs(DeltaOM))/sin(alpha))*sin(i_i);
    u_i = atan2(sin_ui, cos_ui);
    u_f = atan2(sin_uf, cos_uf);

    th = mod(2*pi - u_i - om_i, 2*pi);
    % th = 2*pi - u_i - om_i;
    om_f = mod(2*pi - u_f - th, 2*pi);
    % om_f = 2*pi - u_f - th;
    p = a*(1 - e^2);
    DeltaV = 2*sqrt(mu/p)*(1 + e*cos(th))*sin(alpha/2);

elseif (DeltaOM < 0 && Deltai < 0)

    alpha = acos(cos(i_i)*cos(i_f) + sin(i_i)*sin(i_f)*cos(abs(DeltaOM)));
    cos_ui = (cos(i_f) - cos(alpha)*cos(i_i)) / (sin(alpha)*sin(i_i));
    cos_uf = (-cos(i_i) + cos(alpha)*cos(i_f)) / (sin(alpha)*sin(i_f));
    sin_ui = (sin(abs(DeltaOM))/sin(alpha))*sin(i_f);
    sin_uf = (sin(abs(DeltaOM))/sin(alpha))*sin(i_i);
    u_i = atan2(sin_ui, cos_ui);
    u_f = atan2(sin_uf, cos_uf);

    th = mod(u_i - om_i, 2*pi);
    % th = u_i - om_i;
    om_f = mod(u_f - th, 2*pi);
    % om_f = u_f - th;
    p = a*(1 - e^2);
    DeltaV = 2*sqrt(mu/p)*(1 + e*cos(th))*sin(alpha/2);

end

end
function plotOrbit(a, e, i, OM, om, th0, thf, dth, mu)
% plot 3D dell'orbita: conoscendo i parametri orbitali, l'anomalia vera
% iniziale, l'anomalia vera finale e la step size plottiamo l'orbita
% tridimensionale punto per punto

anomalie_vere = th0:dth:thf;
vet_stato = zeros(3, length(anomalie_vere));

for j = 1:length(anomalie_vere)
    [rr_i, ~] = par2car(a, e, i, OM, om, anomalie_vere(j), mu);
    vet_stato(:,j) = rr_i;
end

plot3(vet_stato(1,:), vet_stato(2,:), vet_stato(3,:))
xlabel("Rx")
ylabel("Ry")
zlabel("Rz")
title("3D Orbit Plot")
grid on

end
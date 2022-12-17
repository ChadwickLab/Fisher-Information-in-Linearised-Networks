
function [rE, rI] = SimulateNetwork_Euler(network, inputs, Nt, dt)

NE = network.cells.NE;
NI = network.cells.NI;

rE = zeros([NE, Nt]);
rI = zeros([NI, Nt]);

IE_FF = inputs.IE_FF;
II_FF = inputs.II_FF;

noise = inputs.noise;

JEE = network.connectivity.JEE;
JEI = network.connectivity.JEI;
JIE = network.connectivity.JIE;
JII = network.connectivity.JII;

tauE = network.cells.tauE;
tauI = network.cells.tauI;
gamma = network.cells.gamma;

    
    noiseE = noise * mean(IE_FF);
    noiseI = noiseE / 2;

    for t=1:Nt

        eta_E_FF = noiseE * randn(NE,1);
        eta_I_FF = noiseI * randn(NI,1);

        inputE = JEE * rE(:,t) - JEI * rI(:,t) + IE_FF  + eta_E_FF;
        inputI = JIE * rE(:,t) - JII * rI(:,t) + II_FF  + eta_I_FF;

        drE = (-rE(:,t) + (inputE .* (inputE > 0 )).^gamma) * dt / tauE;
        drI = (-rI(:,t) + (inputI .* (inputI > 0)).^gamma) * dt / tauI;

        rE(:,t+1) = rE(:,t) + drE;
        rI(:,t+1) = rI(:,t) + drI;

    end


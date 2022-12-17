function network = create_network(kEE,kEI,kIE,kII, JEE_mean, JEI_mean,JIE_mean, JII_mean)

% number of cells in each class

NE = 1000;
NI = NE / 5;

% Time constants

tauE = 10;
tauI = tauE / 2;

gamma = 2;

% Connectivity 

theta_pE = linspace(0, 2*pi, NE+1);
theta_pE = theta_pE(1:(end-1)); 
theta_pI = linspace(0, 2*pi, NI+1);
theta_pI = theta_pI(1:(end-1));

JEE = zeros(NE);
JIE = zeros([NI, NE]);
JEI = zeros([NE, NI]);
JII = zeros(NI);



JEE_max = JEE_mean / besseli(0,abs(kEE));

for i=1:NE

    JEE(i,:) = JEE_max * exp(kEE * cos(theta_pE(i) - theta_pE));
    
end


JEI_max = JEI_mean / besseli(0, abs(kEI));  
JIE_max = JIE_mean / besseli(0, abs(kIE));
JII_max = JII_mean / besseli(0, abs(kII));

for i=1:NE

    JIE(:,i) = JIE_max * exp(kIE * cos(theta_pE(i) - theta_pI));
    JEI(i,:) = JEI_max * exp(kEI * cos(theta_pE(i) - theta_pI));
    
end

for i=1:NI
    
    JII(i,:) = JII_max * exp(kII * cos(theta_pI(i) - theta_pI));

end


% create struct

network = struct;
network.connectivity = struct;
network.connectivity.JEE = JEE;
network.connectivity.JEI = JEI;
network.connectivity.JIE = JIE;
network.connectivity.JII = JII;

network.cells = struct;
network.cells.NE = NE;
network.cells.NI = NI;
network.cells.tauE = tauE;
network.cells.tauI = tauI;
network.cells.gamma = gamma;

end
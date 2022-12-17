%% This script finds a fixed point of the supralinear stabilised network and computes the eigenvectors/eigenvalues of the Jacobian

%% set parameters

% Euler (to initialise FixedPointFinder)
Nt = 10000;  % number of timesteps 
dt = 0.01;   % time step
noise = 0; % noise in Euler method simulation (find initialisation for fixed point using noise-free simulation)

% inputs
theta_s = pi;  % stimulus orientation
kE_FF = 0.5;   % feedforward E tuning
kI_FF = 0.0;   % feedforward I tuning
IE_FF_area = 0.5;  % feedforward E strength
II_FF_area = 0;  % feedforward I strength

% weights
JEE_mean = 0.019;  % EE strength
JEI_mean = 0.04;   
JIE_mean = 0.04;
JII_mean = JEE_mean * 1.1;      
kEE = 2;  % EE tuning
kIE = 0.1; 
kEI = 0.4; 
kII = 0.0; 

%% create network

network = create_network(kEE,kEI,kIE,kII, JEE_mean, JEI_mean, JIE_mean, JII_mean);
NE = network.cells.NE;
NI = network.cells.NI;
inputs  = create_inputs(theta_s, noise, kE_FF, kI_FF, IE_FF_area, II_FF_area, network);

%% Find fixed point and linearise
  
[rE, rI] = SimulateNetwork_Euler(network, inputs, Nt, dt);
R0 = [mean(rE(:,(Nt/2):end),2)', mean(rI(:,(Nt/2):end),2)']; % initial guess for fixed point             
FixedPointFinder; % find a fixed point
FixedPoint = max(0,rmin'); % correct for any negative values (due to numerical precision)
Phip = diag(2 * FixedPoint.^(1/2)); % Phi'  
Jtilde = (W * Phip - eye(NE+NI)) * inv(T); % Jacobian after change of variables

%% find eigenmodes of Jtilde

[Vleft,D] = eig(Jtilde');  
Vleft = Vleft';
Evals = diag(D);        

%% get input statistics   
     
inputs.noise = 2; % noise for computing covariance matrix
Inp = ([inputs.IE_FF .* (- kE_FF * sin(inputs.theta_pE - theta_s))'; inputs.II_FF .* (- kI_FF * sin(inputs.theta_pI - theta_s))']);  % derivative of input tuning curves    
CovInp =  diag([inputs.noise * mean(inputs.IE_FF) * ones(NE,1); inputs.noise/2 * mean(inputs.IE_FF) * ones(NI,1)]);  % covariance of input
lindisc_input =  inv(CovInp) * Inp ; 
lindisc_input = lindisc_input / norm(lindisc_input);   % linear discriminant of inputs      

%% mode SNRs and time constants          

SNRmode = (Vleft * Inp).^2 ./ diag(Vleft * CovInp * Vleft') / (Inp' * pinv(CovInp) * Inp);  % normalised mode input SNR
taumode = -1./real(diag(D));        % mode time constants
                      
%% compute stationary state response information
        
J = (Phip * W - eye(NE+NI)); % Jacobian of network output
Sigma = lyap((inv(T) * J), inv(T) * Phip * CovInp *inv(T) * Phip);   % covariance of network output
rp = -inv(J) * Phip * Inp;  % tuning curve slopes
InfOut =  rp' * pinv(Sigma) * rp  / (Inp' * pinv(CovInp) * Inp);  % normalised output information

%% Plot top N modes (left eigenvectors, real and imaginary parts)

N=5;
for i=1:N
subplot(N,1,i)
hold on
plot(real(Vleft(i,:)), 'linewidth', 3)
plot(imag(Vleft(i,:)), 'linewidth', 3)
xlabel('Neuron')
ylabel('Eigenvector')
title(strcat('eigenvalue=', num2str(D(i,i))))
legend('real', 'imag')
set(gca,'fontsize', 18)
end

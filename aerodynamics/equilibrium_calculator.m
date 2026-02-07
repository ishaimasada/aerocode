clc
clear

global A; % Integration Coefficients from Prabhu & Erickson

A = zeros(5, 7, 11);
 
A(:,:,1) = [
            3.77037, -0.289522e-2, 0.953322e-05, -0.924699e-08, 0.301919e-11, -18.8598, 3.69335;
            2.8969178, 0.23736545e-2, -0.1491710e-05, 0.4660339e-09, -0.5394517e-13, 82.240429, 7.5019384;
            2.8421112, 0.13220561e-2, -0.3391585e-06, 0.4465218e-10, -0.2291482e-14, 0.5835039e+03, 8.6255038;
            5.6821087, -0.75300588e-3, 0.2298007e-06, -0.2395592e-10, 0.8004847e-15, -0.2470037e+04, -9.8740054;
            -0.27258966, 0.20115141e-2, -0.2454771e-06, 0.1202520e-10, -0.2138986e-15, 0.7611774e+04, 31.631739;
           ];

A(:,:,2) = [
            3.46226, 0.000582024, -0.305255e-05, 0.622801e-08, -0.337560e-11, 0.879517, 3.21926;
            2.7022403, 0.0019443393, -8.930004e-06, 0.1973921e-09, -0.1696781e-13, 0.2018022e+03, 7.2040846;
            3.9143506, 0.00031537108, -0.5648104e-07, 0.3601230e-11, 0.5235943e-16, -0.5355196e+03, 0.021671998;
            1.2657466, 0.0018269791, -0.3758390e-06, 0.3303372e-10, -0.9365111e-15, 0.3142696e+04, 17.943292;
            27.715943, -0.0074173486, 0.8239597e-06, -0.3528530e-10, 0.4967120e-15, -0.5694269e+05, -0.1740287e+03;
           ];

A(:,:,3) = [
            3.21671, -0.378227e-2, 0.847468e-05, -0.886582e-08, 0.353656e-11, 0.296485e+05, 1.85264;
            2.6045368, -0.17235464e-3, 0.1157414e-06, -0.3641786e-10, 0.4601154e-14, 0.2972896e+05, 4.5865262;
            2.8101683, -0.2903991e-3, 0.9083335e-07, -0.9942782e-11, 0.3770436e-15, 0.2953661e+05, 3.2536467;
            1.9209266, 0.21776554e-3, -0.1828841e-07, 0.4905096e-12, 0.2850730e-17, 0.3078342e+05, 9.2748629;
            1.9209266, 0.21776554e-3, -0.1828841e-07, 0.4905096e-12, 0.2850730e-17, 0.3078342e+05, 9.2748629;
           ];

A(:,:,4) = [
            4.20644, -0.450984e-2, 0.1055739e-04, -0.8591938e-08, 0.2404710e-11, 0.108890e+05, 2.3137934;
            2.7543778, 0.23093284e-2, -0.1282336e-05, 0.3404352e-09, -0.3480754e-13, 0.1113432e+05, 9.0789671;
            3.8015418, 0.49857546e-3, -0.1253132e-06, 0.1409389e-10, -0.4482019e-15, 0.1066557e+05, 3.161939;
            4.9133164, 0.6175527e-06, -0.5522238e-07, 0.1168649e-10, -0.5464249e-15, 0.8845378e+04, -4.57869;
            20.456854, -0.61498079e-2, 0.8691427e-06, -0.5087751e-10, 0.1062416e-14, -0.2295524e+05, -0.1156196e+03;
           ];

A(:,:,5) = [
            2.5, 0, 0, 0, 0, 0.566267e+05, 4.18073;
            2.50751, -0.247979e-04, 0.296415e-07, -0.152881e-10, 0.289137e-14, 0.566249e+05, 4.14319;
            2.63761, -0.873733e-05, -0.647727e-07, 0.234734e-10, -0.173962e-14, 0.564523e+05, 3.22321;
            3.37206, -0.885546e-03, -0.252033e-06, -0.231879e-10, 0.704714e-15, 0.849338e+05, -1.05640;
            -10.2056, 0.428313e-2, -0.493107e-06, 0.248511e-10, -0.468699e-15, 0.849338e+05, 96.4037;
           ];

A(:,:,6) = [
            3.56593, -0.247878e-03, -0.522843e-06, 0.302282e-08, -0.191127e-11, 0.118414e+06, 3.59357;
            2.71526, 0.189981e-2, -0.853462e-06, 0.183903e-09, -0.153285e-13, 0.118621e+06, 7.93691;
            3.9311, 0.306608e-03, -0.567059e-07, 0.455735e-11, -0.777620e-16, 0.117863e+06, 0.701982;
            1.06562, 0.17918e-2, -0.330778e-06, 0.250527e-10, -0.545401e-15, 0.122142e+06, 20.3283;
            31.3770, -0.781707e-2, 0.754371e-06, -0.247341e-10, 0.159524e-15, 0.481245e+05, -0.202185e+03;
           ];

A(:,:,7) = [
            2.5, 0, 0, 0, 0, -0.042499, -11.7339;
            2.5, 0, 0, 0, 0, -0.042499, -11.7339;
            2.5, 0, 0, 0, 0, -0.042499, -11.7339;
            2.5, 0, 0, 0, 0, -0.042499, -11.7339;
            2.5, 0, 0, 0, 0, -0.042499, -11.7339;
            ];
    
A(:,:,8) = [
            3.08170, -0.377375e-2, 0.967715e-05, -0.110995e-07, 0.472486e-11, 0.225328e+06, 2.43159;
            2.5347, -0.550674e-04, 0.3644e-07, -0.123132e-10, 0.192247e-14, 0.225395e+06, 4.78956;
            2.80056, -0.295664e-03, 0.966487e-07, -0.110536e-10, 0.441452e-15, 0.225185e+06, 3.13053;
            2.00477, 0.191897e-03, -0.154991e-07, 0.411659e-12, 0.262876e-17, 0.226226e+06, 8.46350;
            2.05984, 0.188754e-03, -0.180755e-07, 0.802843e-12, -0.130896e-16, 0.226028e+06, 8.02544;
           ];

A(:,:,9) = [
            2.5, 0, 0, 0, 0, 0.187713e+06, 4.38020;
            2.50134, -0.404313e-05, 0.432715e-08, -0.195277e-11, 0.316328e-15, 0.187713e+06, 4.37342;
            2.22589, 0.296501e-03, -0.115611e-06, 0.184612e-10, -0.892492e-15, 0.187912e+06, 6.05761;
            4.1211, -0.867862e-03, 0.153375e-06, -0.922601e-11, 0.17842e-15, 0.185437e+06, -6.63889;
            1.73495, -0.807107e-04, 0.598063e-07, -0.459985e-11, 0.102952e-15, 0.191074e+06, 10.7919;
           ];

 A(:,:,10) = [
              2.5, 0, 0, 0, 0, -0.042499, 4.36650;
              2.5, 0, 0, 0, 0, -0.042499, 4.36650;
              2.5, 0, 0, 0, 0, -0.042499, 4.36650;
              2.63035, -0.773934e-04, 0.172197e-07, -0.170282e-11, 0.631988e-16, -0.17548e+03, 3.48901;
              -9.1424, 0.385082e-02, -0.468004e-06, 0.245069e-10, -0.456469e-15, 0.277519e+05, 89.4618;
             ];
 
 A(:,:,11) = [
              2.68188, -0.207203e-2, 0.747860e-05, -0.876720e-08, 0.344842e-11, 0.182871e+06, 5.07646;
              2.42206, 0.863470e-03, -0.784594e-06, 0.266740e-09, -0.318506e-13, 0.182853e+06, 5.92363;
              2.55884, 0.897362e-04, -0.500916e-07, 0.882982e-11, -0.523719e-15, 0.183014e+06, 5.55313;
              2.50142, 0.294466e-04, -0.804717e-08, 0.775358e-12, -0.257809e-16, 0.183252e+06, 6.0544;
              2.57061, -0.206211e-04, 0.297860e-08, -0.211019e-12, 0.574768e-17, 0.183223e+06, 5.6165;
             ];

global sigma_O2; global sigma_N2; global sigma_Ar;

sigma_O2 = 14.48;
sigma_N2 = 53.96;
sigma_Ar = 0.321;

% Function to get the (partial pressure) equilibrium constants "Kp"
function K = get_equilibrium_constants(T, rho)
    global A;
    R = 8.314;
    P0 = 101325;
    G_RT = zeros(1, 11);
    delta_N = [1, 0, 1, -1, -1, -1, -1];
    
    % Define the temperature regime
    switch true
        case (200 <= T && T < 800), k = 1;
        case (800 <= T && T < 3000), k = 2;
        case (3000 <= T && T < 6000), k = 3;
        case (6000 <= T && T < 10000), k = 4;
        case (10000 <= T && T <= 15000), k = 5;
        otherwise
            disp("Temperature out of range (200 - 15000)");
            return;
    end
 
    for i=1: length(G_RT)
        G_RT(i) = A(k, 1, i) * (1 - log(T)) - ...
                  A(k, 2, i) * (T/2) - ...
                  A(k, 3, i) * (T^2/6) - ...
                  A(k, 4, i) * (T^3/12) - ...
                  A(k, 5, i) * (T^4/20) + ...
                  A(k, 6, i) / T - ...
                  A(k, 7, i);
    end

    delta_G_RT = [
                  2 * G_RT(3) - G_RT(1), ...
                  2 * G_RT(4) - G_RT(1) - G_RT(2), ...
                  2 * G_RT(5) - G_RT(2), ...
                  G_RT(4) - G_RT(6) - G_RT(7), ...
                  G_RT(5) - G_RT(8) - G_RT(7), ...
                  G_RT(3) - G_RT(9) - G_RT(7), ...
                  G_RT(10) - G_RT(11) - G_RT(7)
                 ];
    
    Kp = exp(-delta_G_RT);
    K = Kp .* (P0 / (R * rho * T)).^delta_N;
end

% First Model for equilibrium composition calculator
function sigmas = model1(T, rho)
    global A;
    global sigma_O2; global sigma_N2; global sigma_Ar;

    if T > 1600, disp("Conditions beyond model applicability"); return; end

    K = get_equilibrium_constants(T, rho);

    sigma1_check = sigma_O2 / 2;
    sigma2_check = sigma_N2 / 2;
    sigma10 = sigma_Ar;

    sigma3 = sqrt(K(1) * sigma1_check);
    sigma5 = sqrt(K(3) * sigma2_check);
    sigma4 = sqrt(K(2) * sigma1_check * sigma2_check);

    sigma1 = (sigma_O2 - sigma3 - sigma4) / 2;
    sigma2 = (sigma_N2 - sigma4 - sigma5) / 2;

    error_threshold = 10^-3;
    error1 = abs(sigma1_check - sigma1);
    error2 = abs(sigma2_check - sigma2);

    while error1 > error_threshold || error2 > error_threshold
        sigma1_check = sigma1;
        sigma2_check = sigma2;
        
        sigma3 = sqrt(K(1) * sigma1_check);
        sigma5 = sqrt(K(3) * sigma2_check);
        sigma4 = sqrt(K(2) * sigma1_check * sigma2_check);

        sigma1 = (sigma_O2 - sigma3 - sigma4) / 2;
        sigma2 = (sigma_N2 - sigma4 - sigma5) / 2;
        error1 = abs(sigma1_check - sigma1);
        error2 = abs(sigma2_check - sigma2);
    end

    sigma7 = sqrt((sigma4/K(4)) + (sigma5/K(5)) + (sigma3/K(6)) + (sigma10/K(7)));
    sigma6 = sigma4 / (K(4) * sigma7);
    sigma8 = sigma5 / (K(5) * sigma7);
    sigma9 = sigma3 / (K(6) * sigma7);
    sigma11 = sigma10 / (K(7) * sigma7);
    
    sigmas = [sigma1, sigma2, sigma3, sigma4, sigma5, sigma6, sigma7, sigma8, sigma9, sigma10, sigma11];
    sigmas = log(sigmas ./ sum(sigmas));
end

function sigmas = model2(T, rho)
    global A;
    global sigma_O2; global sigma_N2; global sigma_Ar;

    function [sigmas, delta_sigma3] = get_sigmas()
        sigmaO2_check = sigma_O2 - sigma6 - sigma9;
        sigmaN2_check = sigma_N2 - sigma6 - sigma8;
        b0 = 2 * K(1) * sigmaO2_check^2;
        b1 = (K(3)*Ke - 4) * K(1) * sigmaO2_check;
        b2 = -8*sigmaO2_check + 2*K(1) + (sigmaO2_check - sigmaN2_check)*K(2) - K(1)*K(3)*Ke;
        b3 = 8 - K(2) - 2*K(3)*Ke;
        b4 = (8 - 2*K(2)) / K(1);
        F = dot([b0, b1, b2, b3, b4], [1, sigma3, sigma3^2, sigma3^3, sigma3^4]);
        Fprime = sum([1:4] .* [b1, b2, b3, b4] .* [1, sigma3, sigma3^2, sigma3^3]);
        delta_sigma3 = -F / Fprime;
        sigma3 = sigma3 + delta_sigma3;
        sigma5 = (sigma_O2 - (2*sigma3^2 / K(1)) - sigma3) / (Ke * sigma3);
        sigma4 = Ke * sigma3 * sigma5;
        sigma2 = (sigmaN2_check - sigma4 - sigma5) / 2;
        sigma1 = (sigmaO2_check - sigma4 - sigma3) / 2;
        sigma7 = sqrt((sigma4/K(4)) + (sigma5/K(5)) + (sigma3/K(6)) + (sigma10/K(7)));
        sigma6 = sigma4 / (K(4) * sigma7);
        sigma8 = sigma5 / (K(5) * sigma7);
        sigma9 = sigma3 / (K(6) * sigma7);
        sigma11 = sigma10 / (K(7) * sigma7);
        sigma10 = sigma_Ar - sigma11;
        sigmas = [sigma1, sigma2, sigma3, sigma4, sigma5, sigma6, sigma7, sigma8, sigma9, sigma10, sigma11];
    end

    error_threshold = 10^-3;
    K = get_equilibrium_constants(T, rho);
    Ke = sqrt(K(2) / (K(1) * K(3)));
    sigma6 = 0; sigma8 = 0; sigma9 = 0; sigma11 = 0;
    sigma10 = sigma_Ar;
    
    sigma1 = sigma_O2; % Assuming...
    sigma3 = (-(1+sqrt(K(1)/(sigma_N2*sigma1))) + sqrt((1+sqrt(K(1)/(sigma_N2*sigma1))^2 - 4*(2/K(1))*(-sigma_O2)))) / (4/K(1));
    if sigma3 < 0, sigma3 = (-(1+sqrt(K(1)/(sigma_N2*sigma1))) - sqrt((1+sqrt(K(1)/(sigma_N2*sigma1))^2 - 4*(2/K(1))*(-sigma_O2)))) / (4/K(1)); end % account for the plus/minus in the quadratic formula
    
    [sigmas, delta_sigma3] = get_sigmas();

    while delta_sigma3 > error_threshold, [sigmas, delta_sigma3] = get_sigmas(); end

    sigmas = log(sigmas ./ sum(sigmas));
end

function sigmas = model3(T, rho)
    global A;
    global sigma_O2; global sigma_N2; global sigma_Ar; %#ok<*GVMIS,*NUSED>

    function [sigmas, delta_sigma3] = get_sigmas(sigma3_check)
        sigmaN2_check = sigma_N2 - sigma4 - sigma6 - sigma8;
        b0 = -K(3) * sigmaN2_check / 2;
        b1 = K(3) / 2;
        b2 = 1;
        Ke = sqrt(K(2) / (K(1) * K(3)));
        sigma5 = (-b1 + sqrt(b1^2 - 4*b1*b0)) / (2*b2); if sigma5 < 0, sigma5 = (-b1 - sqrt(b1^2 - 4*b1*b0)) / (2*b2); end
        sigma2 = (sigmaN2_check - sigma5) / 2;
        sigma1 = sigma3_check^2 / K(1);
        sigma4 = Ke * sigma3_check * sigma5;
        sigma7 = sqrt((sigma3_check/K(6)) + (sigma4/K(4)) + (sigma5/K(5)) + (sigma10/K(7)));
        sigma6 = sigma4 / (K(4) * sigma7);
        sigma8 = sigma5 / (K(5) * sigma7);
        sigma9 = sigma3_check / (K(7) * sigma7);
        sigma11 = sigma10 / (K(7) * sigma7);
        sigma10 = sigma_Ar - sigma11;
        sigma3 = sigma_O2 - 2*sigma1 - sigma4 - sigma6 - sigma9;
        delta_sigma3 = abs(sigma3 - sigma3_check);
        sigmas = [sigma1, sigma2, sigma3, sigma4, sigma5, sigma6, sigma7, sigma8, sigma9, sigma10, sigma11];
    end
    
    K = get_equilibrium_constants(T, rho);
    sigma3_check = sigma_O2;
    sigma10 = sigma_Ar;
    sigma4 = 0; sigma6 = 0; sigma8 = 0;
    
    [sigmas, delta_sigma3] = get_sigmas(sigma3_check);
    error_threshold = 10^-3;
    
    while delta_sigma3 < error_threshold
        sigma3_check = sigmas(3);
        [sigmas, delta_sigma3] = get_sigmas(sigma3_check);
    end

    sigmas = log(sigmas ./ sum(sigmas));
end

function plot_sigmas(temperatures, sigmas, density, altitude)
    O2 = cellfun(@(sigmas) sigmas(1), sigmas);
    N2 = cellfun(@(sigmas) sigmas(2), sigmas);
    O = cellfun(@(sigmas) sigmas(3), sigmas);
    NO = cellfun(@(sigmas) sigmas(4), sigmas);
    N = cellfun(@(sigmas) sigmas(5), sigmas);
    NOplus = cellfun(@(sigmas) sigmas(6), sigmas);
    e = cellfun(@(sigmas) sigmas(7), sigmas);
    Nplus = cellfun(@(sigmas) sigmas(8), sigmas);
    Oplus = cellfun(@(sigmas) sigmas(9), sigmas);
    Ar = cellfun(@(sigmas) sigmas(10), sigmas);
    Arplus = cellfun(@(sigmas) sigmas(11), sigmas);

    figure

    plot(temperatures, O2)
    hold on
    plot(temperatures, N2)
    hold on
    plot(temperatures, O)
    hold on
    plot(temperatures, NO)
    hold on
    plot(temperatures, N)
    hold on
    plot(temperatures, NOplus)
    hold on
    plot(temperatures, e)
    hold on
    plot(temperatures, Nplus)
    hold on
    plot(temperatures, Oplus)
    hold on
    plot(temperatures, Ar)
    hold on
    plot(temperatures, Arplus)
    hold on

    legend("O2", "N2", "O", "NO", "N", "NOplus", "e", "Nplus", "Oplus", "Ar", "Arplus")
    title("Species mole fractions at an altitude of" + " " + altitude + " " + "[km]" + " " + "(density = " + density + " " + "[kg/m^3])")
    ylim([-5, 0])
end

% Driver code
densities = [0.41271, 0.08803, 0.01801];
altitudes = [10000, 20000, 30000];
model1_temperatures = linspace(200, 1600, 10^2);
model2_temperatures = linspace(1600, 3500, 10^2);
model3_temperatures = linspace(3500, 15000, 10^2);
temperatures = [model1_temperatures, model2_temperatures, model3_temperatures];

% Sorry for the copy-pasting. i didn't have the brain power to think of
% a more elegant solution.
for i=1: length(densities)
    model1_sigmas = arrayfun(@(temperature) model1(temperature, densities(i)), model1_temperatures, "UniformOutput", false);
    model2_sigmas = arrayfun(@(temperature) model2(temperature, densities(i)), model2_temperatures, "UniformOutput", false);
    model3_sigmas = arrayfun(@(temperature) model3(temperature, densities(i)), model3_temperatures, "UniformOutput", false);
    sigmas = [model1_sigmas, model2_sigmas, model3_sigmas];
    plot_sigmas(temperatures, sigmas, densities(i), altitudes(i))
end


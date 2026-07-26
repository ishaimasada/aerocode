import numpy
import os
import cantera
import matplotlib.pyplot as plt

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)


# Brute force iteration function
def iterate(function_name, LHS, guess=None, *args):
    error_threshold = 10**-4
    if guess == None: guess = error_threshold
    RHS = function_name(guess, *args)
    error = (RHS - LHS) / LHS
    while abs(error) > error_threshold:
        if error > 0: guess += error_threshold/2
        elif error < 0: guess -= error_threshold/2
        RHS = function_name(guess, *args)
        error = (RHS - LHS) / LHS
    return guess


# Bisection method (faster than brute force)
def bisection(function_name, LHS, guess_high, guess_low, trend, *args):
    error_threshold = 10**-4
    guess_mid = (guess_high + guess_low) / 2
    flow = function_name(guess_low, *args)
    fmid = function_name(guess_mid, *args)
    fhigh = function_name(guess_high, *args)
    match trend:
        case "increasing":
            if LHS < flow or LHS > fhigh:
                raise ValueError("Desired parameter is out of given bounds.")
            error = (LHS - fmid) / LHS
            while abs(error) > error_threshold:
                if fmid < LHS: 
                    guess_low = guess_mid
                elif fmid > LHS: 
                    guess_high = guess_mid
                guess_mid = (guess_high + guess_low) / 2
                fmid = function_name(guess_mid, *args)
                error = (LHS - fmid) / LHS
        case "decreasing":
            if LHS > flow or LHS < fhigh:
                raise ValueError("Desired parameter is out of given bounds.")
            error = (LHS - fmid) / LHS
            while abs(error) > error_threshold:
                if fmid < LHS: 
                    guess_high = guess_mid
                elif fmid > LHS: 
                    guess_low = guess_mid
                guess_mid = (guess_high + guess_low) / 2
                fmid = function_name(guess_mid, *args)
                error = (LHS - fmid) / LHS
    return guess_mid


# A.1
def isentropic(parameter, gamma, lookup_key="M"):
    def get_A_Astar(M):
        if M == 0: return numpy.inf #print("Mach of zero causes singularity. Returning infinity")
        Tt_T = 1 + ((gamma - 1)/2) * M**2
        A_Astar = ((gamma + 1) / 2)**((-gamma - 1) / (2 * (gamma - 1))) * (Tt_T**((gamma + 1) / (2 * (gamma - 1)))) / M
        return A_Astar
    def get_Tt_T(M):
        Tt_T = 1 + ((gamma - 1)/2) * M**2
        return Tt_T
    def get_Pt_P(M):
        Pt_P = (1 + ((gamma - 1)/2) * M**2)**(gamma / (gamma - 1))
        return Pt_P
    def get_rhot_rho(M):
        rhot_rho = (1 + ((gamma - 1)/2) * M**2)**(1 / (gamma - 1))
        return rhot_rho

    match lookup_key:
        case "M":
            M = parameter
            Tt_T = get_Tt_T(M)
            Pt_P = get_Pt_P(M)
            rhot_rho = get_rhot_rho(M)
            A_Astar = get_A_Astar(M)
        case "pressure ratio":
            Pt_P = parameter
            M = iterate(get_Pt_P, Pt_P)
            Tt_T = get_Tt_T(M)
            rhot_rho = get_rhot_rho(M)
            A_Astar = get_A_Astar(M)
        case "temperature ratio":
            Tt_T = parameter
            M = iterate(get_Tt_T, Tt_T)
            Pt_P = get_Pt_P(M)
            rhot_rho = get_rhot_rho(M)
            A_Astar = get_A_Astar(M)
        case "density ratio":
            rhot_rho = parameter
            M = iterate(get_rhot_rho, rhot_rho)
            Tt_T = get_Tt_T(M)
            Pt_P = get_Pt_P(M)
            A_Astar = get_A_Astar(M)
        case "area ratio":
            A_Astar = parameter
            if A_Astar < 1:
                raise ValueError("Area ratio cannot be less than 1.")
            M_subsonic = bisection(get_A_Astar, A_Astar, 1, 0.001, trend="decreasing")
            #M_supersonic = bisection(get_A_Astar, A_Astar, 10, 1, trend="increasing")
            Tt_T = get_Tt_T(M_subsonic)
            Pt_P = get_Pt_P(M_subsonic)
            rhot_rho = get_rhot_rho(M_subsonic)
            return [M_subsonic, Tt_T, Pt_P, rhot_rho, A_Astar]
    return [M, Tt_T, Pt_P, rhot_rho, A_Astar]

# A.2
def normal_shock(parameter, gamma, lookup_key="M"):
    def get_P2_P1(M0):
        P2_P1 = 1 + ((2*gamma) / (gamma + 1)) * (M0**2 - 1)
        return P2_P1
    def get_rho2_rho1(M0):
        rho2_rho1 = ((gamma + 1)*M0**2) / (2 + (gamma - 1)*M0**2)
        return rho2_rho1
    def get_T2_T1(M0):
        T2_T1 = (2*gamma*M0**2 - (gamma - 1)) * ((gamma - 1)*M0**2 + 2) / ((gamma + 1)**2 * M0**2)
        return T2_T1
    def get_Pt2_Pt1(M0):
        Pt2_Pt1 = ((gamma + 1) * M0**2 / ((gamma - 1)*M0**2 + 2))**(gamma / (gamma - 1)) * ((gamma + 1) / ((2*gamma*M0**2) - (gamma - 1)))**(1 / (gamma - 1))
        return Pt2_Pt1
    def get_Pt2_P1(M0):
        Pt2_P1 = (((((gamma + 1)**2)*M0**2) / (4*gamma*M0**2 - 2*(gamma - 1)))**(gamma / (gamma + 1))) * ((1 - gamma + 2*gamma*M0**2) / (gamma + 1))
        return Pt2_P1
    def get_M2(M0):
        M2 = numpy.sqrt(((gamma - 1) * M0**2 + 2) / (2 * gamma * M0**2 - (gamma - 1)))
        return M2

    match lookup_key:
        case "M":
            M0 = parameter
            P2_P1 = get_P2_P1(M0)
            rho2_rho1 = get_rho2_rho1(M0)
            T2_T1 = get_T2_T1(M0)
            Pt2_Pt1 = get_Pt2_Pt1(M0)
            Pt2_P1 = get_Pt2_P1(M0)
            M2 = get_M2(M0)
        case "static pressure ratio":
            M0 = iterate(get_P2_P1, parameter)
            P2_P1 = parameter
            rho2_rho1 = get_rho2_rho1(M0)
            T2_T1 = get_T2_T1(M0)
            Pt2_Pt1 = get_Pt2_Pt1(M0)
            Pt2_P1 = get_Pt2_P1(M0)
            M2 = get_M2(M0)
        case "density ratio":
            M0 = iterate(get_rho2_rho1, parameter)
            P2_P1 = get_P2_P1(M0)
            rho2_rho1 = parameter
            T2_T1 = get_T2_T1(M0)
            Pt2_Pt1 = get_Pt2_Pt1(M0)
            Pt2_P1 = get_Pt2_P1(M0)
            M2 = get_M2(M0)
        case "tempeature ratio":
            M0 = iterate(get_T2_T1, parameter)
            P2_P1 = get_P2_P1(M0)
            rho2_rho1 = get_rho2_rho1(M0)
            T2_T1 = parameter
            Pt2_Pt1 = get_Pt2_Pt1(M0)
            Pt2_P1 = get_Pt2_P1(M0)
            M2 = get_M2(M0)
        case "total to total ratio":
            M0 = iterate(get_Pt2_P1, parameter)
            P2_P1 = get_P2_P1(M0)
            rho2_rho1 = get_rho2_rho1(M0)
            T2_T1 = get_T2_T1(M0)
            Pt2_Pt1 = parameter
            Pt2_P1 = get_Pt2_P1(M0)
            M2 = get_M2(M0)
        case "total to static ratio":
            M0 = iterate(get_rho2_rho1, parameter)
            P2_P1 = get_P2_P1(M0)
            rho2_rho1 = get_rho2_rho1(M0)
            T2_T1 = get_T2_T1(M0)
            Pt2_Pt1 = get_Pt2_Pt1(M0)
            Pt2_P1 = parameter
            M2 = get_M2(M0)
        case "total to static ratio":
            M0 = iterate(get_M2, parameter)
            P2_P1 = get_P2_P1(M0)
            rho2_rho1 = get_rho2_rho1(M0)
            T2_T1 = get_T2_T1(M0)
            Pt2_Pt1 = get_Pt2_Pt1(M0)
            Pt2_P1 = get_Pt2_Pt1
            M2 = parameter

    return [M0, P2_P1, rho2_rho1, T2_T1, Pt2_Pt1, Pt2_P1, M2]

# Oblique shocks
def oblique_shock(M1, parameter, gamma, lookup_key="deflection angle"):
    ''' gives post wave properties given the deflection angle (theta) and incoming Mach number '''

    def theta_beta_mach(beta):
        RHS = 2 * (1 / numpy.tan(beta)) * ((M1**2 * numpy.sin(beta)**2 - 1) / ((M1**2 * (gamma + numpy.cos(2*beta)) + 2)))
        return RHS

    match lookup_key:
        case "deflection angle":
            # Angle assumed to be in radians
            theta = parameter
            LHS = numpy.tan(theta) # find beta using iteration function
            beta = numpy.atan(iterate(theta_beta_mach, LHS))
        case "wave angle":
            beta = parameter
            RHS = theta_beta_mach(beta)
            theta = numpy.atan(RHS)

    Mn1 = M1 * numpy.sin(beta)
    Mt2 = M1 * numpy.cos(beta)
    [Mn1, P2_P1, rho2_rho1, T2_T1, Pt2_Pt1, Pt2_P1, Mn2] = normal_shock(Mn1, gamma)
    M2 = Mn2 / numpy.sin(beta - theta)
    return [M2, theta, beta, Pt2_P1, P2_P1, T2_T1]


'''
RAYLEIGH FLOW FUNCTIONS
'''
def get_P_Pstar(M, gamma):
    P_Pstar = (1 + gamma) / (1 + gamma*M**2)
    return P_Pstar
def get_Pt_Ptstar(M, gamma):
    Pt_Ptstar = ((1 + gamma) / (1 + gamma*M**2)) * ((1 + ((gamma - 1)/2) * M**2) / ((gamma + 1) / 2))**(gamma / (gamma - 1))
    return Pt_Ptstar
def get_T_Tstar(M, gamma):
    T_Tstar = (M**2)*((1 + gamma) / (1 + gamma*M**2))**2
    return T_Tstar
def get_Tt_Ttstar(M, gamma):
    Tt_Ttstar = (M**2)*(((1 + gamma) / (1 + gamma*M**2))**2) * ((1 + ((gamma - 1)/2) * M**2) / ((gamma + 1) / gamma))
    return Tt_Ttstar
def get_rho_rhostar(M, gamma):
    rho_rhostar = ((1 + gamma*M**2) / (1 + gamma)) / M**2
    return rho_rhostar

def rayleigh(parameter, gamma, lookup_key="M"):
    match lookup_key:
        case "M":
            M = parameter
            P_Pstar = get_P_Pstar(parameter, gamma)
            Pt_Ptstar = get_Pt_Ptstar(parameter, gamma)
            T_Tstar = get_T_Tstar(parameter, gamma)
            Tt_Ttstar = get_Tt_Ttstar(parameter, gamma)
            rho_rhostar = get_rho_rhostar(parameter, gamma)
        case "static temperature ratio":
            M = iterate(get_T_Tstar, parameter)
            P_Pstar = get_P_Pstar(M, gamma)
            Pt_Ptstar = get_Pt_Ptstar(M, gamma)
            T_Tstar = parameter
            Tt_Ttstar = get_Tt_Ttstar(M, gamma)
            rho_rhostar = get_rho_rhostar(M, gamma)
        case "total temperature ratio":
            M = iterate(get_Tt_Ttstar, parameter)
            P_Pstar = get_P_Pstar(M, gamma)
            Pt_Ptstar = get_Pt_Ptstar(M, gamma)
            T_Tstar = get_T_Tstar(M, gamma)
            Tt_Ttstar = parameter
            rho_rhostar = get_rho_rhostar(M, gamma)
        case "static pressure ratio":
            M = iterate(get_P_Pstar, parameter)
            P_Pstar = parameter
            Pt_Ptstar = get_Pt_Ptstar(M, gamma)
            T_Tstar = get_P_Pstar(M, gamma)
            Tt_Ttstar = get_Tt_Ttstar(M, gamma)
            rho_rhostar = get_rho_rhostar(M, gamma)
        case "total pressure ratio":
            M = iterate(get_Pt_Ptstar, parameter)
            P_Pstar = get_P_Pstar(M, gamma)
            Pt_Ptstar = parameter
            T_Tstar = get_T_Tstar(M, gamma)
            Tt_Ttstar = get_Tt_Ttstar(M, gamma)
            rho_rhostar = get_rho_rhostar(M, gamma)
        case "density ratio":
            M = iterate(get_rho_rhostar, parameter)
            P_Pstar = get_P_Pstar(M, gamma)
            Pt_Ptstar = get_Pt_Ptstar(M, gamma)
            T_Tstar = get_T_Tstar(M, gamma)
            Tt_Ttstar = get_Tt_Ttstar(M, gamma)
            rho_rhostar = parameter

    return [M, P_Pstar, Pt_Ptstar, T_Tstar, Tt_Ttstar, rho_rhostar]

# A.4
def fanno(parameter, gamma, lookup_key="M"):
    def get_P_Pstar(M):
        P_Pstar = numpy.sqrt((1 + gamma) / (2 + (gamma - 1)*M**2)) / M
        return P_Pstar
    def get_Pt_Ptstar(M):
        Pt_Ptstar = (((2 + (gamma - 1)*M**2) / (1 + gamma))**((gamma + 1) / (2 * (gamma - 1)))) / M
        return Pt_Ptstar
    def get_T_Tstar(M):
        T_Tstar = (1 + gamma) / (2 + (gamma - 1)*M**2)
        return T_Tstar
    def get_rho_rhostar(M):
        rho_rhostar = numpy.sqrt((2 + (gamma - 1)*M**2) / (1 + gamma)) / M
        return rho_rhostar
    def get_length_term(M):
        length_term = (1 - M**2) / (gamma * M**2) + ((gamma + 1) / (2 * gamma)) * numpy.log(((gamma + 1) * M**2) / (2 + (gamma - 1)* M**2))
        return length_term

    match lookup_key:
        case "M":
            M = parameter
            P_Pstar = get_P_Pstar(parameter)
            Pt_Ptstar = get_Pt_Ptstar(parameter)
            T_Tstar = get_T_Tstar(parameter)
            rho_rhostar = get_rho_rhostar(parameter)
            length_term = get_length_term(parameter)
        case "temperature ratio":
            M = iterate(get_T_Tstar, parameter)
            P_Pstar = get_P_Pstar(M)
            Pt_Ptstar = get_Pt_Ptstar(M)
            T_Tstar = parameter
            rho_rhostar = get_rho_rhostar(M)
            length_term = get_length_term(M)
        case "static pressure ratio":
            M = iterate(get_P_Pstar, parameter)
            P_Pstar = parameter
            Pt_Ptstar = get_Pt_Ptstar(M)
            T_Tstar = get_T_Tstar(M)
            rho_rhostar = get_rho_rhostar(M)
            length_term = get_length_term(M)
        case "total pressure ratio":
            M = iterate(get_Pt_Ptstar, parameter)
            P_Pstar = get_P_Pstar(M)
            Pt_Ptstar = parameter
            T_Tstar = get_T_Tstar(M)
            rho_rhostar = get_rho_rhostar(M)
            length_term = get_length_term(M)
        case "density ratio":
            M = iterate(get_rho_rhostar, parameter)
            P_Pstar = get_P_Pstar(M)
            Pt_Ptstar = get_Pt_Ptstar(M)
            T_Tstar = get_T_Tstar(M)
            rho_rhostar = parameter
            length_term = get_length_term(M)
        case "length term":
            M = iterate(get_length_term, parameter)
            P_Pstar = get_P_Pstar(M)
            Pt_Ptstar = get_Pt_Ptstar(M)
            T_Tstar = get_T_Tstar(M)
            rho_rhostar = get_rho_rhostar(M)
            length_term = parameter

    return [M, P_Pstar, Pt_Ptstar, T_Tstar, rho_rhostar, length_term]

# A.5
def expansion_fan(M1, theta, gamma):
    def prandtl_meyer(M):
        nu = numpy.sqrt((gamma + 1)/(gamma - 1)) * numpy.atan(((gamma - 1)/(gamma + 1)) * (M**2 - 1)) - numpy.atan(numpy.sqrt(M**2 - 1))
        return nu


    theta = numpy.radians(theta)
    nu1 = prandtl_meyer(M1)
    nu2 = nu1 + theta

    # Iterate for nu2
    M2 = iterate(prandtl_meyer, nu2, guess=1)
    return M2, nu2

def shock_tube(T1, T4, P1, P4, gamma):
    def get_P4_P1(P2_P1):
        P4_P1 = P2_P1 * (1 - (((gamma - 1)*(a1/a4)*(P2_P1 - 1)) / numpy.sqrt(2*gamma*(2*gamma + (gamma + 1)*(P2_P1 - 1)))))**(-2*gamma / (gamma - 1))
        return P4_P1
    R = 287
    a1 = numpy.sqrt(gamma * R * T1)
    a4 = numpy.sqrt(gamma * R * T4)
    P4_P1 = P4 / P1
    P2_P1 = iterate(get_P4_P1, P4_P1)
    return P2_P1

def two_dimension_airofil(M_freestream, thetas):
    pass
    '''
    import numpy

    M_freestream = numpy.full(len(thetas), M_freestream)
    deflection_angles = numpy.array(thetas)

    Machs, nus = numpy.vectorize(expansion_fan)(M_freestream, deflection_angles)
    '''

def get_TM_velocities(theta, Vr, Vtheta, h):
    def get_Vr_double_prime(theta, Vr, Vtheta):
        return (Vr * (Vtheta**2) - ((gamma - 1) / 2) * (1 - Vr**2 - Vtheta**2) * (2*Vr + Vtheta * (1 / numpy.tan(theta)))) / (((gamma - 1) / 2) * (1 - Vr**2 - Vtheta**2) - Vtheta**2)
    
    # Runge-Kutta Coefficients
    K1 = Vtheta
    M1 = get_Vr_double_prime(theta, Vr, Vtheta)
    K2 = Vtheta + (h/2) * M1
    M2 = get_Vr_double_prime(theta + (h/2), Vr + K1*(h/2), Vtheta + M1*(h/2))
    K3 = Vtheta + M2 * (h/2)
    M3 = get_Vr_double_prime(theta + (h/2), Vr + K2*(h/2), Vtheta + M2*(h/2))
    K4 = Vtheta + M3 * h
    M4 = get_Vr_double_prime(theta + h, Vr + K3*h, Vtheta + M3*h)

    # Update variables
    Vtheta += (h/6) * (M1 + (2*M2) + (2*M3) + M4)
    Vr += (h/6) * (K1 + (2*K2) + (2*K3) + K4)

    return Vr, Vtheta

def TM_cone(theta_s, M0, gamma):
    # Runge-Kutta solution to Taylor-Maccoll equation for flowfield between shock and cone surface

    if theta_s == 0: return 0 # Check if the shock wave is zero
    [M1, deflection, _, _, _, T1_T0] = oblique_shock(M0, theta_s, gamma, lookup_key="wave angle")

    # Initial conditions
    V1 = numpy.sqrt(M1**2 / (2/(gamma - 1) + M1**2)) # Non-dimensionalize all velocities
    Vr = V1 * numpy.cos(theta_s - deflection)
    Vtheta = -V1 * numpy.sin(theta_s - deflection)
    theta = theta_s

    # Runge-Kutta Numerical Method
    h = -10**-5 # step size in theta (radians)
    Vr, Vtheta = get_TM_velocities(theta, Vr, Vtheta, h) # update velocities

    # Iterate until tangential velocity is zero
    while Vtheta < 0:
        Vr_updated, Vtheta_updated = get_TM_velocities(theta, Vr, Vtheta, h)

        # Linear interpolation for precise zeroing of Vtheta
        if Vtheta_updated >= 0:
            fraction = -Vtheta / (Vtheta_updated - Vtheta)
            return numpy.degrees(theta + fraction * h)

        theta += h
        Vr, Vtheta = Vr_updated, Vtheta_updated

    return theta

def busemann_inlet(M3, altitude, theta_s, throat_radius):
    # Busemann Inlet Design

    min_theta_s = numpy.deg2rad(12)
    gamma = 1.4
    R = 287
    tolerance = 10^-3
    h = 10^-4
    max_iterations = 10^5
    max_theta_s = numpy.asin(1 / M3)
    max_theta = numpy.pi
    Cp = gamma * R / (gamma - 1)

    '''
    %% CHECK CONE HALF ANGLE BOUNDS

    % Maximum
    if theta_s > max_theta_s
        fprintf("Shock angle exceeds the maximum --> decrease the shock half angle.\n")
        return
    end

    % Minimum
    if theta_s < min_theta_s
        fprintf("Shock angle is below the minimum --> increase the shock half angle.\n")
        return
    end


    %% DELTA ITERATION

    % Upstream Mach number and Deflection Angle
    Mn3 = M3 * sin(theta_s);
    Mn2 = sqrt((2 + Mn3^2 * (gamma - 1)) / (2 * Mn3^2 * gamma - (gamma - 1)));

    % Guess delta
    delta_guess = tolerance; % Guess deflection in radians
    beta = delta_guess + theta_s; % Wave angle in radians
    M2 = Mn2 / sin(beta);

    % Solve for new delta
    delta = atan((2 * cot(beta) * (M2^2 * sin(beta)^2 - 1) / (M2^2 * (gamma + cos(2*beta)) + 2))); % Deflection angle in radians
    error = abs(delta - delta_guess);
    iterations = 0;

    % Iterate delta (flow deflection) until it converges
    while error > tolerance && iterations < max_iterations
        iterations = iterations + 1;
        
        % Add or subtract the step size based on sign of error
        if delta_guess < delta
            delta_guess = delta_guess + h;
        elseif delta_guess > delta
            delta_guess = delta_guess - h;
        end

        % Update variables
        beta = delta_guess + theta_s;
        M2 = Mn2 / sin(beta);
        delta = atan((2 * cot(beta) * (M2^2 * sin(beta)^2 - 1) / (M2^2 * (gamma + cos(2*beta)) + 2))); % Deflection angle in radians
        error = abs(delta - delta_guess);
    end


    %% TAYLOR-MACCOLL INTEGRATION

    % Initial Conditions (pre-shock state --> "station" 2)
    V2 = sqrt(M2^2 / (2/(gamma-1) + M2^2)); % Non-dimensionalized velocity using Equation 10.16 from [1]
    Vr = V2 * cos(beta);
    Vtheta = - V2 * sin(beta);
    r = throat_radius / sin(theta_s);
    theta = theta_s;
    contour_x = [r * cos(theta)];
    contour_y = [r * sin(theta)];
    Vy = Vr*sin(theta) + Vtheta*cos(theta); % Cartesian y-velocity

    % March theta positively (CCW/upstream) until Vy is zero (force maximum
    % theta of 180deg)
    while abs(Vy) > tolerance && theta < max_theta

        % Solve Taylor-Maccoll for the velocity components
        [Vr_updated, Vtheta_updated] = get_TM_velocity(theta, Vr, Vtheta, gamma, h);

        % Handle the zero crossing with linear interpolation
        if Vtheta * Vtheta_updated < 0
            fraction = -Vtheta / (Vtheta_updated - Vtheta);
            theta_zero = theta + fraction * h;
            r_zero = r + fraction * (h * r * (Vr / Vtheta));
        
            % interpolate velocities
            Vr = Vr + fraction * (Vr_updated - Vr);
            Vtheta = Vtheta + fraction * (Vtheta_updated - Vtheta);
            contour_x(end+1) = r_zero * cos(theta_zero);
            contour_y(end+1) = r_zero * sin(theta_zero);
            break
        end

        % Update variables
        theta = theta + h;
        dr = h * r * (Vr / Vtheta);
        r = r + dr;
        Vr = Vr_updated;
        Vtheta = Vtheta_updated;

        % Convert to Cartesian and store contour coordinates
        contour_x(end+1) = r * cos(theta);
        contour_y(end+1) = r * sin(theta);
        Vy = Vr*sin(theta) + Vtheta*cos(theta);
    end


    %% REPORT RESULTS

    % Atmospheric Properties
    [~, Patm, Tatm] = standard_atm(altitude);

    % Area Contraction Ratio
    CR = (contour_y(end))^2 / throat_radius^2;

    % Freestream Mach number
    V1 = sqrt(Vr^2 + Vtheta^2);
    M1 = sqrt((2 * V1^2) /((gamma-1) * (1 - V1^2))); % Rearanging Eq 10.16 with Vmax=1 & V=V1

    % Inlet Loss
    Pt3_Pt1 = (((gamma + 1)*Mn2^2) / ((gamma - 1)*Mn2^2 + 2))^(gamma / (gamma - 1)) * ((gamma + 1) / ((2 * gamma * Mn2^2) - (gamma - 2)))^(1 / (gamma - 1));

    % Static Pressure Ratio
    P3_P2 = ((2 * gamma * Mn2^2) - (gamma - 1)) / (gamma + 1);
    Pt2_P2 = (1 + (gamma -1)/2 * M2^2)^(gamma / (gamma - 1));
    Pt1_P1 = (1 + (gamma -1)/2 * M1^2)^(gamma / (gamma - 1));
    P3_P1 = P3_P2 * (1/Pt2_P2) * Pt1_P1;

    % Static Temperature Ratio
    T3_T2 = (((2 * gamma * Mn2^2) - (gamma - 1)) * ((gamma - 1)*Mn2^2 + 2)) / ((gamma + 1)^2 * Mn2^2);
    Tt2_T2 = 1 + (gamma - 1)/2 * M2^2;
    Tt1_T1 = (1 + (gamma - 1)/2 * M1^2);
    T3_T1 = T3_T2 * (1/Tt2_T2) * Tt1_T1;


    %% CONTOUR

    figure;
    plot(contour_x, contour_y, 'b-', 'LineWidth', 2.5); 
    hold on;
    plot(contour_x, -contour_y, 'b-', 'LineWidth', 2.5);
    axis equal;
    xlabel('Axial coordinate (normalized)');
    ylabel('Radial coordinate (normalized)');
    title('Busemann Inlet Wall Contour');
    '''


def equilibrium_air(T):
    pass


'''
M3 = 2.27
altitude = 10000
theta_s = numpy.deg2rad(17.2)
throat_radius = 1
M = 10
gamma = 1.2
mu = numpy.rad2deg(numpy.asin(1 / M))
betas = numpy.linspace(mu + 0.1, 75, 100)
cone_angles = [TM_cone(beta, M, gamma) for beta in betas]

plt.plot(betas, cone_angles)
plt.show()
'''
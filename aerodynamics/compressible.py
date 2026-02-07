import numpy
import os
import cantera

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

from atmosphere import *

# General iteration function
def iterate(function_name, LHS, guess=None):
    error_threshold = 10**-4
    if guess == None:
        guess = error_threshold
    RHS = function_name(guess)
    error = abs(LHS - RHS) / LHS
    while error > error_threshold:
        guess += error_threshold/2
        RHS = function_name(guess)
        error = abs(LHS - RHS) / LHS
    return guess

def bisection(function_name, guess_max, guess_min):
    error_threshold = 10**-4
    max_iter=100
    fmax, fmin = function_name(guess_max), function_name(guess_min)
    if fmax * fmin > 0:
        raise ValueError("f(a) and f(b) must have opposite signs.")

    for _ in range(max_iter):
        guess_mid = 0.5 * (guess_max + guess_min)
        fmid = function_name(guess_mid)

        # Check convergence
        if abs(fmid) < error_threshold or 0.5 * (guess_min - guess_max) < error_threshold: return guess_mid

        # Narrow the interval
        if fmax * fmid < 0: b, fb = guess_mid, fmid
        else: a, fmax = guess_mid, fmid

    # If max iterations reached, return midpoint
    return 0.5 * (guess_max + guess_min)


# A.1
def isentropic(parameter, gamma, lookup_key="M"):
    def get_A_Astar(M):
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
            Tt_T = get_Tt_T(parameter)
            Pt_P = get_Pt_P(parameter)
            rhot_rho = get_rhot_rho(parameter)
            A_Astar = get_A_Astar(parameter)
        case "pressure ratio":
            M = iterate(Pt_P, parameter)
            Tt_T = get_Tt_T(M)
            Pt_P = parameter
            rhot_rho = get_rhot_rho(M)
            A_Astar = get_A_Astar(M)
        case "temperature ratio":
            M = iterate(Tt_T, parameter)
            Tt_T = parameter
            Pt_P = get_Pt_P(M)
            rhot_rho = get_rhot_rho(M)
            A_Astar = get_A_Astar(M)
        case "density ratio":
            M = iterate(rhot_rho, parameter)
            Tt_T = get_Tt_T(M)
            Pt_P = get_Pt_P(M)
            rhot_rho = parameter
            A_Astar = get_A_Astar(M)
        case "area ratio":
            M_subsonic, M_supersonic = iterate(A_Astar, parameter)
            Tt_T = get_Tt_T(M)
            Pt_P = get_Pt_P(M)
            rhot_rho = get_rhot_rho(M)
            A_Astar = parameter
            return [[M_subsonic, M_supersonic], Tt_T, Pt_P, rhot_rho, A_Astar]
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
        Pt2_Pt1 = ((gamma + 1) * M0**2 / ((gamma + 1)*M0**2 + 2))**(gamma / (gamma - 1)) * ((gamma + 1) / ((2*gamma*M0) - (gamma - 1)))**(1 / (gamma - 1))
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

# A.3
def rayleigh(parameter, gamma, lookup_key="M"):
    def get_P_Pstar(M):
        P_Pstar = (1 + gamma) / (1 + gamma*M**2)
        return P_Pstar
    def get_Pt_Ptstar(M):
        Pt_Ptstar = ((1 + gamma) / (1 + gamma*M**2)) * ((1 + ((gamma - 1)/2) * M**2) / ((gamma + 1) / 2))**(gamma / (gamma - 1))
        return Pt_Ptstar
    def get_T_Tstar(M):
        T_Tstar = (M**2)*((1 + gamma) / (1 + gamma*M**2))**2
        return T_Tstar
    def get_Tt_Ttstar(M):
        Tt_Ttstar = (M**2)*(((1 + gamma) / (1 + gamma*M**2))**2) * ((1 + ((gamma - 1)/2) * M**2) / ((gamma + 1) / gamma))
        return Tt_Ttstar
    def get_rho_rhostar(M):
        rho_rhostar = ((1 + gamma*M**2) / (1 + gamma)) / M**2
        return rho_rhostar

    match lookup_key:
        case "M":
            M = parameter
            P_Pstar = get_P_Pstar(parameter)
            Pt_Ptstar = get_Pt_Ptstar(parameter)
            T_Tstar = get_T_Tstar(parameter)
            Tt_Ttstar = get_Tt_Ttstar(parameter)
            rho_rhostar = get_rho_rhostar(parameter)
        case "static temperature ratio":
            M = iterate(get_T_Tstar, parameter)
            P_Pstar = get_P_Pstar(M)
            Pt_Ptstar = get_Pt_Ptstar(M)
            T_Tstar = parameter
            Tt_Ttstar = get_Tt_Ttstar(M)
            rho_rhostar = get_rho_rhostar(M)
        case "total temperature ratio":
            M = iterate(get_Tt_Ttstar, parameter)
            P_Pstar = get_P_Pstar(M)
            Pt_Ptstar = get_Pt_Ptstar(M)
            T_Tstar = get_T_Tstar(M)
            Tt_Ttstar = parameter
            rho_rhostar = get_rho_rhostar(M)
        case "static pressure ratio":
            M = iterate(get_P_Pstar, parameter)
            P_Pstar = parameter
            Pt_Ptstar = get_Pt_Ptstar(M)
            T_Tstar = get_P_Pstar(M)
            Tt_Ttstar = get_Tt_Ttstar(M)
            rho_rhostar = get_rho_rhostar(M)
        case "total pressure ratio":
            M = iterate(get_Pt_Ptstar, parameter)
            P_Pstar = get_P_Pstar(M)
            Pt_Ptstar = parameter
            T_Tstar = get_T_Tstar(M)
            Tt_Ttstar = get_Tt_Ttstar(M)
            rho_rhostar = get_rho_rhostar(M)
        case "density ratio":
            M = iterate(get_rho_rhostar, parameter)
            P_Pstar = get_P_Pstar(M)
            Pt_Ptstar = get_Pt_Ptstar(M)
            T_Tstar = get_T_Tstar(M)
            Tt_Ttstar = get_Tt_Ttstar(M)
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

def integrate_TM(theta_s, M0, gamma):
    # Runge-Kutta solution to Taylor-Maccoll equation for flowfield between shock and cone surface
    def get_Vr_double_prime(theta, Vr, Vtheta):
        return (Vr * (Vtheta**2) - ((gamma - 1) / 2) * (1 - Vr**2 - Vtheta**2) * (2*Vr + Vtheta * (1 / numpy.tan(theta)))) / (((gamma - 1) / 2) * (1 - Vr**2 - Vtheta**2) - Vtheta**2)
    
    def get_velocities(theta, Vr, Vtheta):
        # Runge-Kutta Coefficients
        K1 = Vtheta
        M1 = get_Vr_double_prime(theta, Vr, Vtheta)
        K2 = Vtheta + (1/2) * M1 * h
        M2 = get_Vr_double_prime(theta + (h/2), Vr + K1 * (h/2), Vtheta + M1 * (h/2))
        K3 = Vtheta + M2 * (h/2)
        M3 = get_Vr_double_prime(theta + (h/2), Vr + K2 * (h/2), Vtheta + M2 * (h/2))
        K4 = Vtheta + M3 * h
        M4 = get_Vr_double_prime(theta + h, Vr + K3 * h, Vtheta + M3 * h)

        # Update variables
        Vtheta += h * (M1 + (2*M2) + (2*M3) + M4) / 6
        Vr += h * (K1 + (2*K2) + (2*K3) + K4) / 6

        return Vr, Vtheta


    if theta_s == 0: return 0 # Check if the shock wave is zero
    [M1, deflection, _, _, _, T1_T0] = oblique_shock(M0, theta_s, gamma, lookup_key="wave angle")

    # Initial conditions
    V1 = numpy.sqrt(M1**2 / (2/(gamma - 1) + M1**2)) # Non-dimensionalize all velocities
    Vr = V1 * numpy.cos(theta_s - deflection)
    Vtheta = -V1 * numpy.sin(theta_s - deflection)
    theta = theta_s

    # Runge-Kutta Numerical Method
    h = -10**-5 # step size in theta (radians)
    Vr, Vtheta = get_velocities(theta, Vr, Vtheta) # update velocities

    # Iterate until tangential velocity is zero
    while Vtheta < 0:
        Vr_updated, Vtheta_updated = get_velocities(theta, Vr, Vtheta)

        # Linear interpolation for precise zeroing of Vtheta
        if Vtheta_updated >= 0:
            fraction = -Vtheta / (Vtheta_updated - Vtheta)
            return numpy.degrees(theta + fraction * h)

        theta += h
        Vr, Vtheta = Vr_updated, Vtheta_updated

def equilibrium_air(T):
    pass

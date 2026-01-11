import math
import os
import cantera

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

gamma = 1.4

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
def isentropic(parameter, lookup_key="M"):
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
def normal_shock(parameter, lookup_key="M"):
    def get_P2_P1(M1):
        P2_P1 = 1 + ((2*gamma) / (gamma + 1)) * (M1**2 - 1)
        return P2_P1
    def get_rho2_rho1(M1):
        rho2_rho1 = ((gamma + 1)*M1**2) / (2 + (gamma - 1)*M1**2)
        return rho2_rho1
    def get_T2_T1(M1):
        T2_T1 = (1 + ((2*gamma) / (gamma + 1)) * (M1**2 - 1)) * ((2 + (gamma - 1)*M1**2 / (gamma + 1)*M1**2))
        return T2_T1
    def get_Pt2_Pt1(M1):
        Pt2_Pt1 = ((gamma + 1) * M1**2 / ((gamma + 1)*M1**2 + 2))**(gamma / (gamma - 1)) * ((gamma + 1) / ((2*gamma*M1) - (gamma - 1)))**(1 / (gamma - 1))
        return Pt2_Pt1
    def get_Pt2_P1(M1):
        Pt2_P1 = (((((gamma + 1)**2)*M1**2) / (4*gamma*M1**2 - 2*(gamma - 1)))**(gamma / (gamma + 1))) * ((1 - gamma + 2*gamma*M1**2) / (gamma + 1))
        return Pt2_P1
    def get_M2(M1):
        M2 = (1 + ((gamma - 1)/2)*M1**2) / (gamma*M1**2 - (gamma - 1)/2)
        return M2

    match lookup_key:
        case "M":
            M1 = parameter
            P2_P1 = get_P2_P1(M1)
            rho2_rho1 = get_rho2_rho1(M1)
            T2_T1 = get_T2_T1(M1)
            Pt2_Pt1 = get_Pt2_Pt1(M1)
            Pt2_P1 = get_Pt2_P1(M1)
            M2 = get_M2(M1)
        case "static pressure ratio":
            M1 = iterate(get_P2_P1, parameter)
            P2_P1 = parameter
            rho2_rho1 = get_rho2_rho1(M1)
            T2_T1 = get_T2_T1(M1)
            Pt2_Pt1 = get_Pt2_Pt1(M1)
            Pt2_P1 = get_Pt2_P1(M1)
            M2 = get_M2(M1)
        case "density ratio":
            M1 = iterate(get_rho2_rho1, parameter)
            P2_P1 = get_P2_P1(M1)
            rho2_rho1 = parameter
            T2_T1 = get_T2_T1(M1)
            Pt2_Pt1 = get_Pt2_Pt1(M1)
            Pt2_P1 = get_Pt2_P1(M1)
            M2 = get_M2(M1)
        case "tempeature ratio":
            M1 = iterate(get_T2_T1, parameter)
            P2_P1 = get_P2_P1(M1)
            rho2_rho1 = get_rho2_rho1(M1)
            T2_T1 = parameter
            Pt2_Pt1 = get_Pt2_Pt1(M1)
            Pt2_P1 = get_Pt2_P1(M1)
            M2 = get_M2(M1)
        case "total to total ratio":
            M1 = iterate(get_Pt2_P1, parameter)
            P2_P1 = get_P2_P1(M1)
            rho2_rho1 = get_rho2_rho1(M1)
            T2_T1 = get_T2_T1(M1)
            Pt2_Pt1 = parameter
            Pt2_P1 = get_Pt2_P1(M1)
            M2 = get_M2(M1)
        case "total to static ratio":
            M1 = iterate(get_rho2_rho1, parameter)
            P2_P1 = get_P2_P1(M1)
            rho2_rho1 = get_rho2_rho1(M1)
            T2_T1 = get_T2_T1(M1)
            Pt2_Pt1 = get_Pt2_Pt1(M1)
            Pt2_P1 = parameter
            M2 = get_M2(M1)
        case "total to static ratio":
            M1 = iterate(get_M2, parameter)
            P2_P1 = get_P2_P1(M1)
            rho2_rho1 = get_rho2_rho1(M1)
            T2_T1 = get_T2_T1(M1)
            Pt2_Pt1 = get_Pt2_Pt1(M1)
            Pt2_P1 = get_Pt2_Pt1
            M2 = parameter

    return [M1, P2_P1, rho2_rho1, T2_T1, Pt2_Pt1, Pt2_P1, M2]

# Oblique shocks
def oblique_shock(M1, theta):
    '''
    gives post wave properties given the deflection angle (theta) and incoming Mach number
    '''
    def theta_beta_mach(beta):
        RHS = 2 * (1 / math.tan(beta)) * ((M1**2 * math.sin(beta)**2 - 1) / ((M1**2 * (gamma + math.cos(2*beta)) + 2)))
        return RHS

    theta = math.radians(theta)

    # find beta using iteration function
    LHS = math.tan(theta)
    beta = math.atan(iterate(theta_beta_mach, LHS))

    Mn1 = M1 * math.sin(beta)
    [Mn1, P2_P1, rho2_rho1, T2_T1, Pt2_Pt1, Pt2_P1, Mn2] = normal_shock(Mn1)
    M2 = Mn2 / math.sin(beta - theta)
    return [M2, P2_P1]

# A.3
def rayleigh(parameter, lookup_key="M"):
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
def fanno(parameter, lookup_key="M"):
    def get_P_Pstar(M):
        P_Pstar = math.sqrt((1 + gamma) / (2 + (gamma - 1)*M**2)) / M
        return P_Pstar
    def get_Pt_Ptstar(M):
        Pt_Ptstar = (((2 + (gamma - 1)*M**2) / (1 + gamma))**((gamma + 1) / (2 * (gamma - 1)))) / M
        return Pt_Ptstar
    def get_T_Tstar(M):
        T_Tstar = (1 + gamma) / (2 + (gamma - 1)*M**2)
        return T_Tstar
    def get_rho_rhostar(M):
        rho_rhostar = math.sqrt((2 + (gamma - 1)*M**2) / (1 + gamma)) / M
        return rho_rhostar
    def get_length_term(M):
        length_term = (1 - M**2) / (gamma * M**2) + ((gamma + 1) / (2 * gamma)) * math.log(((gamma + 1) * M**2) / (2 + (gamma - 1)* M**2))
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
def expansion_fan(M1, theta):
    print(M1, theta)
    def prandtl_meyer(M):
        nu = math.sqrt((gamma + 1)/(gamma - 1)) * math.atan(((gamma - 1)/(gamma + 1)) * (M**2 - 1)) - math.atan(math.sqrt(M**2 - 1))
        return nu
    gamma = 1.4
    theta = math.radians(theta)
    nu1 = prandtl_meyer(M1)
    nu2 = nu1 + theta

    # Iterate for nu2
    M2 = iterate(prandtl_meyer, nu2, guess=1)
    return M2, nu2

def shock_tube(T1, T4, P1, P4):
    def get_P4_P1(P2_P1):
        P4_P1 = P2_P1 * (1 - (((gamma - 1)*(a1/a4)*(P2_P1 - 1)) / math.sqrt(2*gamma*(2*gamma + (gamma + 1)*(P2_P1 - 1)))))**(-2*gamma / (gamma - 1))
        return P4_P1
    gamma = 1.4 # assuming gamma 1 == gamma 2
    R = 287
    a1 = math.sqrt(gamma * R * T1)
    a4 = math.sqrt(gamma * R * T4)
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

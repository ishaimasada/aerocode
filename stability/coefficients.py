'''
Calculator for stability derivatives
'''
import math

def get_lift_curve_slope(A, M, tan_half_lambda):
    '''
    Emperical formula for 3D lift curve slope
    A: weighting aspect ratio
    M: Mach number
    K: Constant (usually 1)
    tan_half_lambda: Tangent of the sweep at the half chord
    '''
    # Check for supersonic or subsonic
    if M < 1: beta = math.sqrt(1 - M**2)
    else: beta = math.sqrt(M**2 - 1)
    K = 1 # Commonly assumed value
    CL_alpha = 2 * math.pi * A / (2 + math.sqrt((A**2 * beta**2 / K**2) * (1 + (tan_half_lambda**2)) + 4))

    return CL_alpha

# GEOMETRY INPUTS
#{altitude ~= 0ft #}
weight = 0.8 * 9.8 # N
Ma = 0
alpha = math.radians(10) # degrees to radians
rho = 1.225 # kg/m^3
Ta = 300 # K
mu = 1.716 * 10**-5 * ((Ta / 273.15)**3/2) * ((273.15 + 110.4) / (Ta + 110.4))
V = 1 # m/sec
q = (1/2) * rho * V**2
Re = rho * V * 0.1 / mu
R_bar = 8.314 # J/mol*K
gamma = 1.4
bw = 0.631 # m
crw = 0.1 # m
ctw = 0 # m
taper_ratio = ctw/crw
lambda_w = math.radians(10) # Leading edge sweep in degrees
gamma_w = math.radians(0) # Dihedral angle in degrees
CLo = 0.2
CDo = 0.07
h = 0.1349
h1 = 0.1349
h2 = 0.1349
hh = 0 # no tail
lh = 0 # no tail
lf = 0.1 # m
SBS = 0.1349 * lf # Projected body side area
K = 1
e = 1
lCG = 0.5 * lf
xCG = lCG / lf
num_sections = 1
xLE = lf - (13.1 + 3.5)
y_mac = (2 * bw) / (3 * math.pi)
xac = y_mac * math.tan(lambda_w)
zw = -1.1 
h_max = 5.1
zv = 0 # no vertical tail
lv = 0 # no vertical tail
wing_planform = "elliptical"
fuselage = "blended"
horizontal_tail = "none"
vertical_tail = "conventional"

match wing_planform.lower():
    case "elliptical":
        a = 0.09823837 # m
        b = 0.63134878 # m
        Sw = math.pi * a * b # m^2
        Aw = bw**2 / Sw
        c_bar = ((8*b**3) / (2 * Sw) - (8 * b**5) / (24 * a * Sw)) # from equation of ellipse
        tan_half_lambda_w = math.tan(lambda_w) - (4 * (0.5) * (1 - taper_ratio) / (Aw * (1 + taper_ratio)))
        quarter_lambda = math.atan(math.tan(lambda_w) - ((1 - taper_ratio) / (Aw * (1 + taper_ratio))))
        Cl_alpha = 2 * math.pi # 2D lift (rad**-1)
        CL_alpha = Cl_alpha / (1 + (57.3*Cl_alpha) / (math.pi * e * Aw)) # 3D Lift coefficient
        CL = CL_alpha * math.degrees(alpha)
        L = q * Sw * CL
        CDCL = 2*CL / (math.pi * e * Aw)
        CDalpha = CL_alpha * math.pi / 180 * CDCL # (rad**-1)
        CD = CDalpha * math.radians(alpha) + CDo
        D = q * Sw * CD
        # wing dihedral contribution to directional stability
        cnbeta_gamma = ((2 * gamma_w) / (Sw * bw)) * (CL - CDalpha) * (b**3 / (6 * a**4))
        # wing sweep contribution to directional stability
        local_alpha = alpha * (1/math.cos(quarter_lambda))
        cnbeta_lambda = ((2 * math.tan(quarter_lambda)) / (Sw * bw)) * ((CL_alpha * local_alpha**2 * (1/math.cos(quarter_lambda))**2) + (2 * CDo * math.cos(quarter_lambda)) + (CDalpha * local_alpha * (1/math.cos(quarter_lambda)))) * (2/3) * ((1 - (b**2 / (4 * a**2))**(3/2)) - 1)
        # combined wing sweep and dihedral contributions to lateral stability
        clbeta_lambda_gamma = (((4/3) * b * a**2 * (gamma_w + alpha * math.tan(lambda_w))) / (Sw * bw)) * CL_alpha * ((1 - (b**2 / (4 * a**2))**(3/2)) - 1)

    case "tapered":
        Sw = (bw / 2) * crw * (1 + taper_ratio)
        Aw = bw**2 / Sw
        c_bar = (2 / 3) * crw * ((1 + taper_ratio + taper_ratio**2) / (1 + taper_ratio)) # m


match fuselage.lower():
    case "conventional":
        '''
        # Multhop's Method for Fuselage Contribution
        bf = [1, 2, 2, 2, 1]
        bf_2 = bf.^2
        delta_x = lf / num_sections
        x1 = []
        for i in range(0, num_sections):
            # Forward sections
            x1[i] = xLE - ((2*i - 1) * delta_x / 2)
            
            # Aft Sections
            if x1[i] <= 0:
                n = i - num_sections/2
                x1[i] = (2*n - 1) * (delta_x / 2)
        x1_cre = x1 ./ crw
        cmalpha_f = 0.399116 # from Excel file
        '''
    case "flying wing":
        cmalpha_f = 0
    case "blended":
        bf = 0.1
        cmalpha_f = (math.pi / (2 * Sw * c_bar)) * bf**2 * lf


match horizontal_tail.lower():
    case "none":
        crt = ctt = bt = lt = lambda_tLE = lambda_t = St = At = lt = 0
        downwash = 2 * CL_alpha / (math.pi * Aw)
        cmalpha_t = 0
    case "conventional":
        crt = 1 # m
        ctt = 1 # m
        bt = 1 # m
        lt = 1 # m
        lambda_tLE = 1 # m
        lambda_t = crt / ctt
        St = bt/2 * (ctt + crt)
        At = bt**2 / St
        lt = lf - lCG - 3/4 * (crt)
        # Downwash Model
        KA = 1/Aw - (1 / (1 + Aw**1.7))
        Klamba = (10 - 3*taper_ratio) / 7
        KH = (1 - (hh/bw)) / math.pow(2 * lh / bw, 1/3)
        downwash_M0 = 4.44 * (KA * Klamba * KH * math.sqrt(math.cos(quarter_lambda)))**1.19
        CL_alpha_M0 = get_lift_curve_slope(Aw, 0, tan_half_lambda_w)
        downwash = downwash_M0 * CL_alpha / CL_alpha_M0
        # Horizontal Tail Contribution
        tan_half_lambda_t = math.tan(lambda_tLE) - (4 * (0.5) * (1 - lambda_t) / (At * (1 + lambda_t)))
        CL_alpha_t = get_lift_curve_slope(At, Ma, tan_half_lambda_t)
        eta = 1
        V_bar = (St / Sw) * (lt / crw)
        cmalpha_t = -1 * CL_alpha_t * eta * V_bar


match vertical_tail.lower():
    case "none":
        crv = ctv = bv = lv = lambda_vLE = lambda_v = Sv = Av = lv = 0
        clbeta_v = 0
    case "conventional":
        crv = 0.07673321 # m
        ctv = 0.02181554 # m
        bv = 0.07285 # m
        lv = 0.5 * lf + (crv/4) # m
        Sv = bv/2 * (ctv + crv)
        Av = bv**2 / Sv
        lv = lf - lCG - 3/4 * (crv)
        r1 = 0.05
        lambda_v = crv / ctv
        lambda_vLE = math.radians(51.842773) # m
        # Lateral Stability
        # Vertical Tail Contribution
        Aeff = Av # Change this later because it requires graphs
        tan_half_lambda_v = math.tan(lambda_vLE) - (4 * (0.5) * (1 - lambda_v) / (Aeff * (1 + lambda_v)))
        CL_alpha_v = get_lift_curve_slope(Av, Ma, tan_half_lambda_v)
        quarter_lambda_v = math.atan(math.tan(lambda_vLE) - ((1 - lambda_v) / (Aeff * (1 + lambda_v))))
        sidewash_term = 0.724 + 3.06*(Sv / Sw) / (1 + math.cos(quarter_lambda_v)) + 0.4*(zw / h_max) + 0.009*Aeff

        emperical_parameter = bv / (2*r1)
        if 0 <= emperical_parameter <= 2:
            k = 0.75
        elif 2 < emperical_parameter < 3.5:
            graph_slope = (1 - 0.75) / (3.5 - 2)
            k = emperical_parameter * graph_slope + 0.75
        elif emperical_parameter > 3.5:
            k = 1
        clbeta_v = -k * CL_alpha_v * sidewash_term * (Sv / Sw) * ((zv * math.cos(alpha) - lv * math.sin(alpha)) / bw)

# wing contribution to longitudinal stability
xa = xCG - xac
cmalpha_w = CL_alpha * xa 
#No = (CL_alpha_t*(1 - downwash)*eta*V_bar - cmalpha_f + CLalpha*xac) / CLalpha


# Fuselage contribution to directional stability
Kn = 0.0015 # FROM GRAPH
KRl = 1.6 # FROM GRAPH
cnbeta_f = -1 * Kn * KRl * (SBS/Sw) * (lf/bw)

# Total longitudinal stability
cmalpha_total = cmalpha_f + cmalpha_w + cmalpha_t 
# Total directional stability
cnbeta_total = cnbeta_f + cnbeta_gamma + cnbeta_lambda
# Total lateral stability
clbeta_total = clbeta_lambda_gamma + clbeta_v
# Vstall

print(f"lift: {L}")
print(f"weight: {weight}")
print(f"L/W: {L/weight}")
print(f"drag: {D}")
print(f"L/D: {CL/CD}")
print(f"longitudinal stability: {cmalpha_total}")
print(f"directional stability: {cnbeta_total}")
print(f"lateral stability: {clbeta_total}")
import numpy as np
import math
np.seterr(divide='ignore', invalid='ignore') # to ignore the warning about /zero division
# CalcI with numpy matrices like Benoit
def CalcI(Cel, Eo, MAC):
    """
    @Cel shape: (1x92)
    @Eo shape: (1)
    @MAC shape: (92x92)
    """
    theta = 40 # degree according to Michel
    eps = 0.001
    lnJ = 0
    M = 0

    """ Eel for K levels """

    Eel = [13.60, 24.60, 54.80, 111.00, 188.00, 283.80, 401.60, 532.00, 685.40, 866.90, 1072.10, 1305.00, 1559.60, 1838.90, 2145.50, 2472.00, 2822.40,
    3202.90, 3607.40, 4038.10, 4492.80, 4966.40, 5465.10, 5989.20, 6539.00, 7112.00, 7708.90, 8332.80, 8978.90, 9658.60, 10367.10, 11103.10, 11866.70,
    12657.80, 13473.70, 14325.60, 15199.70, 16104.60, 17038.40, 17997.60, 18985.60, 19999.50, 21044.00, 22117.20, 23219.90, 24350.30, 25514.00, 26711.20,
    27939.90, 29200.10, 30491.20, 31813.80, 33169.40, 34561.40, 35984.60, 37440.60, 38924.60, 40443.00, 41990.60, 43568.90, 45184.00, 46834.20, 48519.00,
    50239.10, 51995.70, 53788.50, 55617.70, 57485.50, 59389.60, 61332.30, 63313.80, 65350.80, 67416.40, 69525.00, 71676.40, 73870.80, 76111.00, 78394.80,
    80724.90, 83102.30, 85530.40, 88004.50, 90525.90, 93105.00, 95729.90, 98404.00, 101137.00, 103921.90, 106755.30, 109650.90, 112601.40, 115606.10]
    Eel = np.array(Eel)/1000 # convert to keV
    """ Uo """
    Uo = Eo/Eel
    Uorep = np.reshape(np.tile(np.array(Uo), 3), (3, 92)) # checked, tile instead of repeat
    """ M - for sample ## SCALAR"""
    # print(Cel[np.where(Cel > eps)])
    # print(Z_np[np.where(Cel > eps)])
    # print(A_np[np.where(Cel > eps)])

    # same as benoit after multiplying by 100
    M = sum(Cel[0][MeasuredEl] * Z_np[0][MeasuredEl] / A_np[0][MeasuredEl])

    # print("M: ", M)

    """ J Calculation - same value like benoit """
    # Ji like benoit
    Ji = pow(10, -3) * Z_np * (10.04 + 8.25 * np.exp(-Z_np / 11.22))  # in keV

    lnJi = np.log(Ji)[0][MeasuredEl] # should be a scalar
    # J is Javg Benoit code
    # print("sum", sum(Cel[0][MeasuredEl] * (Z_np[0][MeasuredEl] / A_np[0][MeasuredEl]) * lnJi))
    # print('M', M)
    #lnJ same as benoit
    lnJ = sum(Cel[0][MeasuredEl] * (Z_np[0][MeasuredEl] / A_np[0][MeasuredEl]) * lnJi)/M
    J = math.exp(lnJ)  # in keV

    # print("Ji", Ji)
    # print("LnJi", lnJi)
    # print("J", J)

    """ V """
    V = Eo / J
    Vo = V  #### IS THIS TRUE???? ####

    """ m for K lines """
    m = np.exp(-pow(Z_np / 5, 2)) * 0.12 + 0.86
    # m_values = {'K': m_k,
    #             'L': 0.82,
    #             'M': 0.78}

    """ f(V) """
    f_V = 0
    Dk = np.array([6.6 * pow(10, -6), 1.12 * pow(10, -5) * (1.35 - 0.45 * pow(J, 2)), 2.2 * pow(10, -6) / J])
    Dkrep = np.reshape(np.repeat(np.array(Dk), 92), (3,92))

    Pk = [0.78, 0.1, -(0.5 - 0.25 * J)]
    Pkrep = np.reshape(np.repeat(np.array(Pk), 92), (3, 92))

    T = 1+Pkrep-m # similar to result of Benoit

    """ 1/S """
    # print(" (Uo / (Vo * M))",  (Uo / (Vo * M))) # checked
    # print(" ((Vo / Uorep)**Pkrep)",  ((Vo / Uorep)**Pkrep)) # checked
    # print(" (T * pow(Uorep, T) * np.log(Uorep)",  T * pow(Uorep, T) * np.log(Uorep)) # checked
    # print("pow(Uorep, T) + 1",  pow(Uorep, T) + 1) # checked
    # 1/S checked
    inverse_S = (Uo / (Vo * M)) * sum(Dkrep * ((Vo / Uorep)**Pkrep) * (T * pow(Uorep, T) * np.log(Uorep) - pow(Uorep, T) + 1) / pow(T, 2))
    # print('Inverse S', inverse_S)
    """ QlAEoPAP - Ionization Cross section """
    cross_section = np.log(Uo) / (pow(Uo, m) * Eel**2) # similar to result of benoit

    """ R """
    # Z_b_bar # scalar
    Z_b_bar = pow(sum(Cel[0][MeasuredEl]*Z_np[0][MeasuredEl]**0.5), 2) # agreement with benoit
    # eta bar
    eta_bar = 1.75 * pow(10, -3) * Z_b_bar + 0.37 * (1 - math.exp(-0.015 * pow(Z_b_bar, 1.3))) # in agreement with benoit
    # print('eta_bar', eta_bar)
    # W_bar
    W_bar = 0.595 + (eta_bar / 3.7) + (eta_bar**4.55) # in agreement with benoit
    # print('W_bar', W_bar)
    # J(Uo)
    J_Uo = 1 + Uo * (np.log(Uo) - 1)
    # print('J_Uo',J_Uo)
    # G(Uo)
    q = (2 * W_bar - 1) / (1 - W_bar)
    G_Uo = (Uo - 1 - (1 - 1 / pow(Uo, 1 + q)) / (1 + q)) / ((2 + q) * J_Uo)
    R = 1 - eta_bar * W_bar * (1 - G_Uo) # checked

    """ F area """
    # print("inverse S", inverse_S) # checked
    # print("cross-section", cross_section) # checked
    # print("R", R) # checked
    F = (R * inverse_S) / cross_section # checked
    # print('F', F)
    """ Phi(rhoz) """
    r = 2 - 2.3 * eta_bar

    # Phio
    phi0 = 1 + 3.3 * (1 - 1 / pow(Uo, r)) * pow(eta_bar, 1.2)
    lnZ_bar_n = sum(Cel[0][MeasuredEl]*np.log(Z_np[0][MeasuredEl])) # scalar
    Z_bar_n = np.exp(lnZ_bar_n) # scalar
    # Qo in agreement with benoit
    Qo = 1 - 0.535 * math.exp(-pow(21 / Z_bar_n, 1 / 2)) - 2.5 * pow(10, -4) * pow(Z_bar_n / 20,3.5)  # page 60 why 1/2 instead of 1.2!!!!!!!!!!!!!!!!!!!!!!!!!
    # Z bar
    Z_bar = sum(Cel[0][MeasuredEl]*Z_np[0][MeasuredEl])
    # b
    b = 40/Z_bar
    """ Q """
    Q = Qo + (1 - Qo) * np.exp(-(Uo - 1) / b) # agreement with benoit
    """ Ro and Rx """
    h = pow(Z_bar, 0.45)
    D = 1 + 1 / pow(Uo, h)
    Eelrep = np.reshape(np.tile(np.array(Eel), 3), (3, 92)) # tile instead of repeat to solve the problem
    Eorep = np.reshape(np.repeat(np.array(Eo), 92*3), (3, 92))
    Javgrep = np.reshape(np.repeat(np.array(J), 92*3), (3, 92))
    #print("Javgrep", Javgrep)
    """ Ro """
    #print("E0rep",Eorep)
    # print("1/M", 1/M) # checked (like benoit)
    # print("Javgrep**(1-Pkrep)", Javgrep**(1-Pkrep)) # checked
    # print("Dkrep", Dkrep) #checked
    # print("Eorep**(1+Pkrep)",Eorep**(1+Pkrep)) # checked
    # print("Eelrep**(1+Pkrep)",Eelrep**(1+Pkrep)) # checked
    # print("Eelrep", Eelrep) # checkd

    Ro = (1/M) * sum( (Javgrep**(1-Pkrep))*Dkrep*(Eorep**(1+Pkrep) - Eelrep**(1+Pkrep)) / (1+Pkrep) )

    """ Rx same as benoit """

    Rx = Q * D * Ro

    """ Rm """
    # Rm depth of the maximum of the distribution
    # print("Uo-1", Uo-1)
    # if (Uo - 1) > 1:
    ok = (Uo-1) > 1
    Rm = np.ones((92,))
    G1 = 0.11 + 0.41 * math.exp(-pow(Z_bar / 12.75, 0.75))
    G2 =  1 - np.exp(-(np.abs(Uo-1) ** 0.35) / 1.19) ## pay attenstion because Uo-1 can be negative but you have ignored this "abs"
    G3 = 1 - np.exp(-(Uo - 0.5) * pow(Z_bar, 0.4) / 4)
    # print('Z_bar', Z_bar)
    # print('G1', G1)
    # print('G2', G2)
    # print('G3', G3)

    Rm[ok] = G1 * G2[ok] * G3[ok] * Rx[ok] # checked

    """ Rc """

    # d
    # print('Rx - Rm', Rx - Rm)
    # print('(F - (phi0 * Rx) / 3)', (F - (phi0 * Rx) / 3))
    d = (Rx - Rm) * ((F - (phi0 * Rx) / 3) * ((Rx - Rm) * F - phi0 * Rx * (Rm + Rx / 3)))

    # checked
    Rc = (3 / 2)
    Rc *= (((F - phi0 * Rx / 3) / phi0) - (pow(d, 1 / 2)) / (phi0 * (Rx - Rm)))

    """ A1, B1 and A2 """
    A1 = phi0 / (Rm * (Rc - Rx * ((Rc / Rm) - 1))) # checked
    B1 = phi0 - A1 * Rm * Rm # checked
    A2 = A1 * (Rc - Rm) / (Rc - Rx) # checked
    # print("A1", A1)
    # print("B1", B1)
    # print("A2", A2)
    """ X """
    # print("MAC", MAC)
    # print("MAC * Cel", np.matmul(np.transpose(MAC), np.transpose(Cel[0])))
    MAC = np.transpose(MAC)
    MAC[MAC < 0] = 0 # cleaning
    Cel[0][Cel[0] < 0] = 0
    np.nan_to_num(Cel[0], copy=False, nan=0, posinf=0, neginf=0)
    Chi = (MAC @ Cel[0].conj().transpose() * 1/math.sin(math.radians(theta))).conj().T # checked
    #print('Chi', Chi) # checked

    """ First part of the integral to get intensity from Zero to Rc eq.21 """
    Int1atRc = np.exp(-Chi * Rc)* (-B1 * (Chi** 2) - A1 * (2 + 2 * Chi * (Rc - Rm) + (Chi**2)* (Rc - Rm)**2))/(Chi**3)
    Int1at0 = (-B1*(Chi**2) - A1*(2 - 2 *Chi* Rm + (Chi**2) * (Rm**2)))/ (Chi**3)
    Int1 = Int1atRc - Int1at0
    #print('Int1', Int1.shape) # checked

    """ Second part of the integral to get intensity from Rc to Rx """
    Int2atRx = -(2 * A2 * np.exp(-Chi* Rx))/ (Chi**3) #%(2 * A2. * exp(-Chi. * Rc)). / (Chi. ^ 3);
    Int2atRc = (A2* np.exp(-Chi* Rc)* (-2 - 2* Chi*(Rc - Rx) - (Chi** 2)* (Rc - Rx)**2))/ (Chi**3)
    Int2 = Int2atRx - Int2atRc # checked

    # remove the -inf & inf & nan
    np.nan_to_num(Int1, copy=False, nan=0, posinf=0, neginf=0)
    np.nan_to_num(Int2, copy=False, nan=0, posinf=0, neginf=0)


    condition = np.reshape(MeasuredEl, (1,92))
    Intensity = np.zeros((1,92))
    Intensity[condition] = np.reshape(Cel[0],(1,92))[condition]*(Int1[condition] + Int2[condition])
    # print("Intensity", Intensity) # checked
    return Intensity
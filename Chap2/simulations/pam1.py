# Define the PAM solver for the ODE systems

import numpy as np

def monod(S, KS):
    """
    Returns limiting term for substrate S

    Input: S: Concentration of substrate
    Input: KS: Half saturation const for substrate S (same units)
    Output: Limiting Monod term

    This can be used for all limitation terms in the PAM1 or IWA family of models.
    """
    return S / (KS + S)


def pam_batch(y, t, u):
    """ Definition of the system of differential equations to integrate."""
    
    # Unpack the states
    # Assign variables for convenience of notation
    SS = y[0]
    SAC  = y[1]
    SIC = y[2]
    SH2 = y[3]
    SIN = y[4]
    SIP = y[5]
    SI = y[6]
    XPB = y[7]
    XS = y[8]
    XI = y[9]
    #SCAT = 0.001+y[2]*0.3
    
    # Define parameters of the model (basis days)
    fSSXS = 1.6382408312061000E-01
    fSAXS = 1.166839250294608E-01
    fICXS = 1.3039707398869100E-06
    fH2XS = 8.4424680871970500E-02
    fINXS = 0.011622
    fIPXS = 0.002835
    fSIXS = 0.1518209
    fXIXS = 0.4330922

    fICPHAC = 6.44841269841271e-6    # molHCO3-C/mgCOD
    fICPHSS = -1.242761443579e-6     # molHCO3-C/mgCOD

    fACCH = 0.6691         # mgCOD/mgCOD
    fH2CH = 1.0 - fACCH    # mgCOD/mgCOD

    fICAU = 1.0          # molHCO3-C/molHCO3-C
    fH2AU = 40320.0      # mgCOD/molHCO3-C

    fNB = 0.086          # mgNH3-N/mgCOD
    fPB = 0.015          # mgPO4-P/mgCOD

    kHYD = 7.09e-2       # d-1
    kDEC = 9.00e-2       # d-1
    kMAC = 2.375         # d-1
    kMPH = 1.435         # d-1
    #kMCH = 7.386e-2      # d-1
    kMCH = 7.386e-1      # d-1
    kMIC = 6.0e-6        # mol HCO3-C/mgCOD/d

    KSS = 0.524235146324262     # mg COD/L
    KSAC = 20.222562           # mg COD/L
    KSIC = 4.2e-4              # mol HCO3-C/L
    KSIN = 0.02                # mgNH3-N/L
    KSIP = 0.081               # mgPO4-P/L
    KIFA = 7850.0              # mg NH4-N/L
    KSH2 = 1.0                 # mg COD/L

    YPBPH = 1.0                # mgCOD/mgCOD
    YPBCH = 0.680705638006362  # mgCOD/mgCOD
    YPBAU = 40320.00           # mgCOD/molHCO3-C

    KSE = 8.76                 # W/m2
    SE = u                     # W/m2
    fICDEC = -1.98412698412703E-07 # mol C-HCO3 / mgCOD
    fINDEC = 0.058                 # mgNH3-N / mgCOD
    fIPDEC = 0.01                  # mgPO4-P / mgCOD
    
    # Constitutive equations
    rHYD = kHYD * XS
    rDEC = kDEC * XPB
    
    IIN = SIN/(KSIN + SIN)
    IIP = SIP/(KSIP + SIP)
    IFA = KIFA/(KIFA + SIN)
    IE = monod(SE, 8.76)
    ICS = SAC / (SS + SAC)     # ACT inhibition due to SS
    ICAC = SS / (SS + SAC)     # PHT inhibition due to SAC

    
    rACT = kMAC * XPB * IFA * IIN * IIP * IE * (SAC/(KSAC + SAC)) * ICS
    rPHT = kMPH * XPB * IFA * IIN * IIP * IE * (SS/(KSAC + SAC)) * ICAC
    rCHE = kMCH * XPB * IFA * IIN * IIP * SS/(KSS + SS)
    rAUT = kMIC * XPB * IFA * IIN * IIP * IE * SIC/(SIC + KSIC)*SH2/(SH2 + KSH2);
    
    
    n = len(y)      # 10: number of states
    dydt = np.empty((n))
    
    dydt[0] = fSSXS * rHYD \
            - rPHT \
            - rCHE
    
    dydt[1] = fSAXS * rHYD \
            - rACT \
            + (1 - YPBCH)*fACCH*rCHE
        
    dydt[2] = fICXS * rHYD \
            + fICPHAC * rACT \
            + fICPHSS * rPHT \
            - fICAU * rAUT \
            + fICDEC * rDEC

    dydt[3] = fH2XS * rHYD \
            + (1 - YPBCH) * fH2CH * rCHE \
            - fH2AU * rAUT

    dydt[4] = fINXS*rHYD \
            - fNB*YPBPH*rACT \
            - fNB*YPBPH*rPHT \
            + fINDEC*rDEC \
            - fNB*YPBCH*rCHE \
            - fNB*YPBAU*rAUT 
            

    dydt[5] = fIPXS*rHYD \
            - fPB*YPBPH*rACT \
            - fPB*YPBPH*rPHT \
            + fIPDEC*rDEC \
            - fPB*YPBCH*rCHE \
            - fPB*YPBAU*rAUT 
            

    dydt[6] = fSIXS*rHYD
    
    dydt[7] = YPBPH*rACT \
            + YPBPH*rPHT \
            - rDEC \
            + YPBCH*rCHE \
            + YPBAU*rAUT 
            

    dydt[8] = rDEC - rHYD

    dydt[9] = fXIXS * rHYD

    return dydt

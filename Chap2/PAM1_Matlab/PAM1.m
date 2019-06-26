% Function file for the kinetic model.

% ------------------------------------------------------------------------
% 23-11-2016 Author Ed Barry UQ e.barry@uq.edu.au 23/11/2016 for Paper:
%-------------------------------------------------------------------------

%% Function of PAnMBR
function [dydt,tout,pHout] = PAM1(t ,y, tu, u ,V)
% Organise the input so that it can handle the time stepping.
    u = interp1(tu, u, t);

    persistent pH_list t_list
    if t==0
        clear pHlist tlist
        pH_list(1) = 7;
    end
    t_list = [t_list;t];

    %Assignment of components to output matrix
    % Soluble Components
    SS = y(1);             % mgCOD/L
    SAC = y(2);            % mgCOD/L
    SIC = y(3);            % mol HCO3-C/L
    SH2 = y(4);            % mgCOD/L
    SIN = y(5);            % mgNH3-N/L
    SIP = y(6);            % mgPO4-P/L
    SI = y(7);             % mgCOD/L
    SCAT = 0.001+y(3)*0.3;   % cations - M
    %SCAT = 0.003;          % Linked to Scat
    % Particulate Components
    XPB = y(8);            % mgCOD/L
    XS = y(9);             % mgCOD/L
    XI = y(10);            % mgCOD/L
%%  Parameters

    % Hydrolysis
    fSSXS = 1.6382408312061000E-01;        % mgCOD/mgCOD
    fSAXS = 1.166839250294608E-01;         % mgCOD/mgCOD
    fICXS = 1.3039707398869100E-06;        % molHCO3-C/mgCOD
    fH2XS = 8.4424680871970500E-02;        % mgCOD/mgCOD
    fINXS = 1.1622464400041300E-02;        % mgNH3-N/mgCOD
    fIPXS = 2.0754400714359400E-03;        % mgPO4-P/mgCOD
    fSIXS = 1.5182085409209700E-01;        % mgCOD/mgCOD
    fXIXS = 4.3309113162071400E-01;        % mgCOD/mgCOD

    % Photoheterotrophic uptake stoichiometry
    fICPHAC = 6.44841269841271e-6;    % molHCO3-C/mgCOD
    fICPHSS = -1.242761443579e-6;     % molHCO3-C/mgCOD

    % Chemoheterotrophic uptake stoichiometry
    fACCH = 0.6691;         % mgCOD/mgCOD
    fH2CH = 1.0 - fACCH;    % mgCOD/mgCOD

    % Autotrophic Uptake Stoichiometry
    fICAU = 1.0;          % molHCO3-C/molHCO3-C
    fH2AU = 40320.0;      % mgCOD/molHCO3-C

    % Nutrient uptake stoichiometry
    fNB = 0.086;          % mgNH3-N/mgCOD
    fPB = 0.015;          % mgPO4-P/mgCOD

    % Kinetic parameters - Uptake rates
    kHYD = 7.09e-2;       % d-1
    kDEC = 9.00e-2;       % d-1
    kMAC = 2.375;         % d-1
    kMPH = 1.435;         % d-1
    kMCH = 7.386e-2;      % d-1
    kMIC = 6.0e-6;        % mol HCO3-C/mgCOD/d

    % Half saturation and inhibition constants
    KSS = 0.524235146324262;     % mg COD/L
    KSAC = 20.222562;            % mg COD/L
    KSIC = 4.2e-4;               % mol HCO3-C/L
    KSIN = 0.02;                 % mgNH3-N/L
    KSIP = 0.081;                % mgPO4-P/L
    KIFA = 7850.0;               % mg NH4-N/L
    KSH2 = 1.0;                  % mg COD/L

    % Biomass Yields
    YPBPH = 1.0;                 % mgCOD/mgCOD
    YPBCH = 0.680705638006362;   % mgCOD/mgCOD
    YPBAU = 40320.00;            % mgCOD/molHCO3-C

    % Irradiance Parameters
    KSE = 88.0;                  % W/m2
    SE = 176.0;                  % W/m2

    % Decay Parameters
    fICDEC = -1.98412698412703E-07; % mol C-HCO3 / mgCOD / d
    fINDEC = 0.058;                 % mgNH3-N / mgCOD / d
    fIPDEC = 0.01;                  % mgPO4-P / mgCOD / d

    % Solve the pH
    % Electroactive compounds
    y_pH = [SS,SAC,SIC,SIN,SIP,SCAT];
    pH = pHsolve(y_pH,10^-(pH_list(end)));
    pH_list = [pH_list;pH];

%%
    % --------------------- RATE EQUATIONS ----------------------------
    % HYDROLYSIS/FERMENTATION
    rHYD = kHYD * XS;

    % Biomass Decay
    rDEC = kDEC * XPB;

    % INHIBITION EQUATIONS
    % Inorg N Limitation
    IIN = SIN/(KSIN + SIN);

    % Inorg P Limitation
    IIP = SIP/(SIP + KSIP);

    % Ammonium Inhibition
    IFA = KIFA/(KIFA + SIN);

    % Light Inhibition
    IE = SE/(SE + KSE);

    % Competitive inhibition org substrates
    ICs = SAC / (SS + SAC);     % ACT inhibition due to SS
    ICac = SS / (SS + SAC);     % PHT inhibition due to SAC



    % ACETATE UPTAKE
    rACT = kMAC * XPB * IFA * IIN * IIP * IE * (SAC/(KSAC + SAC)) * ICs;%

    % PHOTOHETEROTROPHIC UPTAKE
    rPHT = kMPH * XPB * IFA * IIN * IIP * IE * (SS/(KSAC + SAC)) * ICac;

    % CHEMOHETEROTROPHIC UPTAKE
    rCHE = kMCH * XPB * IFA * IIN * IIP * SS/(KSS + SS);

    % AUTOTROPHIC UPTAKE
    rAUT = kMIC * XPB * IFA * IIN * IIP * IE * ...
        (SIC/(SIC + KSIC))*(SH2/(SH2 + KSH2));

    %%
    % STATE VARIABLES
    dydt(1) = (u(1) - SS)*u(11)/V...
        + fSSXS * rHYD ...
        - rPHT ...
        - rCHE;                             % SS

    dydt(2) = (u(2) - SAC)*u(11)./V ...
        + fSAXS * rHYD ...
        - rACT ...
        + (1 - YPBCH) * fACCH * rCHE;       % SAC

    dydt(3) = (u(3) - SIC)*u(11)/V ...
        + fICXS * rHYD ...
        + fICPHAC * rACT ...
        + fICPHSS * rPHT ...
        - fICAU * rAUT ...
        + fICDEC*rDEC;                             % SIC

    dydt(4) = (u(4) - SH2)*u(11)/V ...
        + fH2XS * rHYD ...
        + (1 - YPBCH) * fH2CH * rCHE...
        - fH2AU * rAUT;                                         % SH2

    dydt(5) = (u(5) - SIN)*u(11)/V ...
        + fINXS*rHYD ...
        - fNB*YPBPH*rACT ...
        - fNB*YPBPH*rPHT ...
        - fNB*YPBCH*rCHE ...
        - fNB*YPBAU*rAUT ...
        + fINDEC*rDEC;         % SIN

    dydt(6) = (u(6) - SIP)*u(11)/V ...
        + fIPXS*rHYD ...
        - fPB*YPBPH*rACT ...
        - fPB*YPBPH*rPHT ...
        - fPB*YPBCH*rCHE ...
        - fPB*YPBAU*rAUT ...
        + fIPDEC*rDEC;          % SIP

    dydt(7) = (u(7) - SI)*u(11)/V ...
        + fSIXS*rHYD;    % SI

    dydt(8) = (u(8) - XPB)*u(11)/V + XPB*u(12)/V ...
        + YPBPH*rACT ...
        + YPBPH*rPHT ...
        + YPBCH*rCHE ...
        + YPBAU*rAUT ...
        - rDEC;                           % XPB

    dydt(9) = (u(9) - XS)*u(11)/V + XS.*u(12)/V ...
        - rHYD ...
        + rDEC;               % XS

    dydt(10) = (u(10)- XI)*u(11)/V + XI*u(12)/V ...
        + fXIXS * rHYD;            % XI
    dydt = dydt';
   if nargout > 1
     pHout = pH_list(2:(end-1),:);
     tout = t_list(1:(end-1),:);
   end
t
end

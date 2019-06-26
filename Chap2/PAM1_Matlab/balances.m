%% Mass Balances

feed2 = feed;
feed2(:,2) = feed(:,2) + SAC_ADDED(:);
feed2 = interp1(tu, feed2, t);
QinVec2 = interp1(tu, QinVec, t);
QoutVec2 = interp1(tu, QoutVec, t);
% CALL IN PARAMETERS
pars;

%% DYNAMIC BALANCE CHECK
% PROCESSES
HYD1 = kHYD .* PAMStateMatrix(:,9);
DEC1 = kDEC .* PAMStateMatrix(:,8);
IIN1 = PAMStateMatrix(:,5)./(KSIN + PAMStateMatrix(:,5));
IIP1 = PAMStateMatrix(:,6)./(PAMStateMatrix(:,6) + KSIP);
IFA1 = KIFA./(KIFA + PAMStateMatrix(:,5));
IE1 = SE/(KSE + SE);
ACT1 = kMAC .* PAMStateMatrix(:,8) .* IFA1 .* IIN1 .* IIP1 .* ...
    IE1 .*(PAMStateMatrix(:,2)./(KSAC + PAMStateMatrix(:,2)));
PHT1 = kMPH .* PAMStateMatrix(:,8) .* IFA1 .* IIN1 .* IIP1 .* ...
    IE1 .* (PAMStateMatrix(:,1)./(KSS + PAMStateMatrix(:,1)));
CHE1 = kMCH .* PAMStateMatrix(:,8) .* IFA1 .* IIN1 .* IIP1 .* ...
    PAMStateMatrix(:,1)./(KSS + PAMStateMatrix(:,1));
AUT1 = kMIC.*PAMStateMatrix(:,8).*IFA1.*IIN1.*IIP1.*IE1.* ...
    (PAMStateMatrix(:,3)./(PAMStateMatrix(:,3) + KSIC)).*(PAMStateMatrix(:,4)./(PAMStateMatrix(:,4) + KSH2));


% Rates of generation for each state variable
gen = zeros(size(PAMStateMatrix));
gen(:,1) = + fSSXS .* HYD1 - PHT1 - CHE1;

gen(:,2) = +  fSAXS .* HYD1 - ACT1 + (1 - YPBCH) .* fACCH .* CHE1;

gen(:,3) = fICXS .* HYD1 + fICPHAC.*ACT1 + fICPHSS.* PHT1 ...
    - fICAU.*AUT1 + fICDEC.*DEC1;

gen(:,4) = + fH2XS.*HYD1 + (1 - YPBCH).*fH2CH.*CHE1 - fH2AU.*AUT1;

gen(:,5) = + fINXS.*HYD1 - fNB.*YPBPH.*ACT1 - fNB.*YPBPH.*PHT1 ...
    - fNB.*YPBCH.*CHE1 - fNB.*YPBAU.*AUT1 + fINDEC.*DEC1;

gen(:,6) = + fIPXS.*HYD1 - fPB .* YPBPH .* ACT1 - fPB .* YPBPH ...
    .* PHT1 - fPB .* YPBCH .* CHE1 - fPB .* YPBAU .* AUT1 +  fIPDEC.*DEC1;

gen(:,7) = + fSIXS.* HYD1;
gen(:,8) = YPBPH .* ACT1 + YPBPH .* PHT1 + YPBCH .* CHE1 + YPBAU .* ...
    AUT1 - DEC1;
gen(:,9) = -HYD1 + DEC1;
gen(:,10) = fXIXS.* HYD1;


% FLOW I/O TERMS
flows = zeros(size(PAMStateMatrix));
for i = 1:7
    flows(:,i) = (feed2(:,i) - PAMStateMatrix(:,i)).*QinVec2./V;
end

for i = 8:10
    flows(:,i) = (feed2(:,i) - PAMStateMatrix(:,i)).*QinVec2./V + QoutVec2./V.*PAMStateMatrix(:,i);
end
    


% ACCUMULATION TERMS
accum = zeros(size(PAMStateMatrix));      % Preallocate the matrix

for i = 1:10
    accum(:,i) = flows(:,i) + gen(:,i);
end


%% BALANCES
epsilon = 1e-6;

% Fractions of TCOD, N, P, C in each state
TCOD_CONTENT = [1 1 0 1 0 0 1 1 1 1];
N_CONTENT = [0 0 0 0 1 0 0.028 0.086 0.028 0.028];
P_CONTENT = [0 0 0 0 0 1 0.005 0.015 0.005 0.005];
C_CONTENT = [2.35588258580082E-05	3.125E-05	1	0	0	0	2.5E-05	2.48015873015874E-05	2.5E-05	2.5E-05];

% Accumulation terms of TCOD, N, P and C
TCOD_accum = zeros(size(accum));
N_accum = zeros(size(accum));
P_accum = zeros(size(accum));
C_accum = zeros(size(accum));


TCOD_LHS = zeros(size(TCOD_accum(:,1)));
N_LHS = zeros(size(N_accum(:,1)));
P_LHS = zeros(size(P_accum(:,1)));
C_LHS = zeros(size(C_accum(:,1)));

TCOD_gen = zeros(size(gen));
N_gen = zeros(size(gen));
C_gen = zeros(size(gen));
P_gen = zeros(size(gen));

TCOD_flows = zeros(size(flows));
N_flows = zeros(size(flows));
P_flows = zeros(size(flows));
C_flows = zeros(size(flows));

TCOD_RHS = zeros(size(TCOD_accum(:,1)));
N_RHS = zeros(size(N_accum(:,1)));
P_RHS = zeros(size(P_accum(:,1)));
C_RHS = zeros(size(C_accum(:,1)));

for i = 1:10
    TCOD_accum(:,i) = TCOD_CONTENT(i) .* accum(:,i);
    TCOD_gen(:,i) = TCOD_CONTENT(i).* gen(:,i);
    TCOD_flows(:,i) = TCOD_CONTENT(i) .* flows(:,i);
    
    N_accum(:,i) = N_CONTENT(i) .* accum(:,i);
    N_gen(:,i) = N_CONTENT(i) .* gen(:,i);
    N_flows(:,i) = N_CONTENT(i) .* flows(:,i);
    
    C_accum(:,i) = C_CONTENT(i) .* accum(:,i);
    C_gen(:,i) = C_CONTENT(i) .* gen(:,i);
    C_flows(:,i) = C_CONTENT(i) .* flows(:,i);
    
    P_accum(:,i) = P_CONTENT(i) .* accum(:,i);
    P_gen(:,i) = P_CONTENT(i) .* gen(:,i);
    P_flows(:,i) = P_CONTENT(i) .* flows(:,i);
    
end

for i = 1:length(TCOD_accum)
    TCOD_LHS(i) = sum(TCOD_accum(i,:));
    TCOD_RHS(i) = sum(TCOD_flows(i,:)) + sum(TCOD_gen(i,:));

    N_LHS(i) = sum(N_accum(i,:));
    N_RHS(i) = sum(N_flows(i,:)) + sum(N_gen(i,:));
    
    P_LHS(i) = sum(P_accum(i,:));
    P_RHS(i) = sum(P_flows(i,:)) + sum(P_gen(i,:));
    
    C_LHS(i) = sum(C_accum(i,:));
    C_RHS(i) = sum(C_flows(i,:)) + sum(C_gen(i,:));
end

TCOD_Error = abs(TCOD_LHS - TCOD_RHS) ./ TCOD_LHS;
N_Error = abs(N_LHS - N_RHS) ./ N_LHS;
P_Error = abs(P_LHS - P_RHS) ./ P_LHS;
C_Error = abs(C_LHS - C_RHS) ./ C_LHS;


if max(TCOD_Error) > epsilon
    disp('TCOD doesn''t balance');
end

if max(N_Error) > epsilon
    disp('N doesn''t balance');
end

if max(P_Error) > epsilon
    disp('P doesn''t balance');
end

if max(C_Error) > epsilon
    disp('C doesn''t balance');
end


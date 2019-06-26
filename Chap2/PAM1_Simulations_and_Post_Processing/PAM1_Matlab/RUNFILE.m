%% RUNFILE for the PAM1 function
% This file can be changed to suit your requirements. Currently, the
% function suits a membrane photo bio CSTR, but will be expanded to general
% reactor types in the future
% Author Ed Barry UQ e.barry@uq.edu.au 23/11/2016 for Paper:
clear; clc;
global V
%% USER INPUT
% Read the input data and organise. Prompt whether the user wants a
% clarifier, and what HRT and SRT are required
data = csvread('dynamicInfluent.csv');   % Raw influent file included
tu = data(:,1); % Time (t) to which raw input (u) is mapped



%----------------------------------------------
%% Assign variable names to influent data
% Conditions preamble. Previously scripted as user prompts - can put these
% back in if required.
HRT = 12;
thetaS = 3;
V = mean(data(:,12))*HRT/24;

% If clarifier is used, the solids content will be reduced by 50% as per
% Tchobanoglous, Ch 5, Section 5-7. Assume same volume of water as
% 'non-clarifier' case.
str = 'N';              % Do we want a clarifier
solidsRemoved = zeros(length(tu),3);
if str == 'Y'
    removalFrac = normrnd(0.6,0.03,length(tu),1);
    for i = 9:11
        solidsRemoved(:,i-8) = data(:,i).*removalFrac;
        data(:,i) = data(:,i).*(1 - removalFrac);
    end
end

Ss0 = data(:,2);
Sac0 = data(:,3);
Sic0 = data(:,4);
Sh20 = data(:,5);
Sin0 = data(:,6);
Sip0 = data(:,7);
Si0 = data(:,8);
XPB0 = data(:,9);
Xs0 = data(:,10);
Xi0 = data(:,11);

% Flowrates
QinVec = data(:,12);
QslVec = V/thetaS;
QoutVec = QinVec - QslVec;       % Solids are removed from the effluent

feed = [Ss0, Sac0, Sic0, Sh20, Sin0, ...
    Sip0, Si0, XPB0, Xs0, Xi0, QinVec, QoutVec];


%% ADDITION OF SAC
% Determine how much SAC is needed to be added to the system
% We will show the system running steadily for 150 days then will
% increase the SAC such that the limiting reagent (N or P) be consumed.
% This will be done for a period of 150 days, then SAC will no longer
% be added to the system for the remainder of the period.

SAC_ADDED = zeros(length(tu), 1);
for j = 1:length(tu) - 1
    if (j+1 > floor(300/0.0104) && j+1 <= floor(450/0.0104))
        if Sin0(j) < 8.6/1.5*Sip0(j)
            SAC_ADDED(j) = 1.2*(Sin0(j)*100/8.6 - (Ss0(j) + Sac0(j)));
        else
            SAC_ADDED(j) = 1.2*(Sip0(j)*100/1.5 - (Ss0(j) + Sac0(j)));
        end
    else
        SAC_ADDED(j) = 0.0;
    end
end

% Pass this array of dynamic inputs into ODE call and the function will 
% take care of mapping it to the time array.
uOrig = feed;
uOrig(:,2) = feed(:,2) + SAC_ADDED;
tspace = [0 600];

%% 
tic
[t,y] = ode15s(@(t,y)PAM1(t,y,tu, uOrig, V), tspace, uOrig(1,1:10));
toc


%% Save files
PAMStateMatrix = y;
[yy,ts,pHs]=PAM1(t(1),y(1,:),tu,uOrig, V);
%can filter out failed iterations or just leave in. Leave in for now.
% Find start time
for i = 2:length(t)
    if t(i-1) < 150 && t(i+1) > 150
        break
    end
end
stateResults = PAMStateMatrix(i:end,:);


if str == 'Y'
    stateInput = [feed, SAC_ADDED, ...
    solidsRemoved, removalFrac, ...
    QinVec, QoutVec];
else
    stateInput = [feed, SAC_ADDED, QinVec, QoutVec];
end
stateInput = interp1(tu, stateInput, t);
stateInputs = stateInput(i:end,:);
%%
csvwrite(sprintf('stateOutputsRaw%iS%iH%s.csv', thetaS, HRT, str), stateResults);
csvwrite(sprintf('stateInputsRaw%iS%iH%s.csv', thetaS, HRT, str), stateInputs);
csvwrite('t.csv',t);
csvwrite('V.csv', V);
csvwrite('pH_vals.csv', [ts, pHs]);
%% Balances
balances;           % Check to see if model balanced over C, COD, N, P

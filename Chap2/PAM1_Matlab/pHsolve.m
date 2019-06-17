% This is a simple pH solver for domestic wastewater built for the PAM1.
% It considers two pairs per electroactive compound, or no pairs in the
% case of cations.
% Pairs are H2PO4/HPO4, HCO3/CO2, Ac-/Ac, NH3/NH4, and the VFA contribution
% from propionate.

% ------------------------------------------------------------------------
% 28-2-2017 Author Damien Batstone Contact Ed Barry UQ e.barry@uq.edu.au 
%  Input file     y = [SS,SAC,SIC,SIN,SIP,SCAT];
%-------------------------------------------------------------------------

function pH = pHsolve(y,z0)

%Define charge densities for COD compounds
CD_SAC = 1/64/1e3;
CD_SS = 1/74/2/1e3; %assume 50% propionate

%Set KA values Everything at 298 for now. Use van't Hoff equation to
%correct if required
Ka_SS = 10^-4.88;
Ka_AC = 10^-4.76;
Ka_CO2 = 10^-6.35;
Ka_NH4 = 10^-9.25;
Ka_H2PO4 = 10^-7.21;
Ka_H2O = 10^-14;

%define zerofun for H+
zerofun = @(SH) y(6) + SH*y(4)/(Ka_NH4+SH)/14000 + SH ...             %cations
    - Ka_SS*CD_SS*y(1)/(Ka_SS+SH)- Ka_AC*CD_SAC*y(2)/(Ka_AC+SH)... %SS, SAC       
    - Ka_CO2*y(3)/(Ka_CO2+SH)- Ka_H2PO4*y(4)/(Ka_H2PO4+SH)/82000-Ka_H2O/SH; %CO2, PO4, H2O (changed col 43 to 6 instead of 4)

%calculate initial guess and slope
SH_soln = fzero(zerofun,z0);
pH=-log10(SH_soln);
disp(pH)
end

% z1=z0*1.1;
% g0=zerofun(z0)
% g1=zerofun(z1)
% slope = (g1-g0)/(z1-z0);
% error = abs((z0-z1)/z1)
% input('')
% while (error>0.001)
%     slope = (g1-g0)/(z1-z0)
%     z0=z1
%     g0=g1
%     z1=z0-slope*g0
%     g1=zerofun(z1);
%     error = abs((z0-z1)/z1);
%     input('')
% end
% SH_out = z1;
% pH = -log10(SH_out)
% end
%    




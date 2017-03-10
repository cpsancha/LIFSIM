%% Load data to define landing gear properties

% Definition of geometrical parameters, as found in LIBIS_Report_Submission
% All dimensions in meters
LD.LandingGear.Xnose = 0.35; 
LD.LandingGear.Xmain = 1.45; 
LD.LandingGear.Track = 1.3;
LD.LandingGear.hcg   = 0.414;

% Using the previous parameters, vectors of position for all 3
% contact points of the landing gear are defined. Axis system is the same as 
% in motorPosition definition: Body axis direction but with origin at the nose.

LD.LandingGear.A.Position = [-LD.LandingGear.Xnose   0   LD.LandingGear.hcg];
LD.LandingGear.B.Position = [-LD.LandingGear.Xmain   LD.LandingGear.Track/2   LD.LandingGear.hcg];
LD.LandingGear.C.Position = [-LD.LandingGear.Xmain   -LD.LandingGear.Track/2  LD.LandingGear.hcg];

% Point A is contact point at the nose. B and C are the rear contact
% points, B has a positive y component in body axes (under the right
% semiwing); C has a negative y component (under the left semiwing).
% All dimensions in meters:
function T = DH_homog(theta, d, a, alpha)
%DH_HOMOG  Homogeneous transformation (i-1)T(i) from the four D-H parameters.
%
%   T = DH_homog(theta, d, a, alpha)
%
%   Inputs : theta, alpha in DEGREES (same units as the Puma 560 D-H table)
%            d, a in any length unit (T(1:3,4) comes out in that unit)
%   Output : 4x4 matrix relating coordinate system i to coordinate system i-1
%
%   D-H convention (Salcudean notes, Slide 2.8) -- product of 4 screw motions:
%       (i-1)T(i) = [e^(theta k x) 0;0' 1] [I d k;0' 1] [I a i;0' 1] [e^(alpha i x) 0;0' 1]
%                 = [ e^(theta k x) e^(alpha i x)    d k + a e^(theta k x) i ]
%                   [            0'                              1           ]
%   Multiplied out:
%       [ c_th  -s_th*c_al   s_th*s_al   a*c_th ]
%       [ s_th   c_th*c_al  -c_th*s_al   a*s_th ]
%       [  0       s_al        c_al        d    ]
%       [  0        0           0          1    ]
%
%   cosd/sind are used so that multiples of 90 deg give exact 0/1 entries.

ct = cosd(theta);  st = sind(theta);     % joint-angle terms
ca = cosd(alpha);  sa = sind(alpha);     % link-twist terms

T = [ ct, -st*ca,  st*sa, a*ct;
      st,  ct*ca, -ct*sa, a*st;
       0,     sa,     ca,    d;
       0,      0,      0,    1 ];
end

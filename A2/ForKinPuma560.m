function [T06, J, O] = ForKinPuma560(q, doPlot)
%FORKINPUMA560  Forward kinematics and manipulator Jacobian of the Puma 560.
%
%   [T06, J, O] = ForKinPuma560(q)          compute only
%   [T06, J, O] = ForKinPuma560(q, true)    also plot the link origins (Q2e)
%
%   Input : q      6x1 joint vector [theta1 ... theta6] in DEGREES
%           doPlot (optional) true -> plot3 of origins o0..o6 marked with '*'
%   Output: T06  4x4 homogeneous transformation 0T6 (base -> gripper), mm
%           J    6x6 manipulator Jacobian in BASE coordinates (C0):
%                   [ v6 ; w6 ] = J * qdot,   v6 in mm/s, w6 & qdot in rad/s
%           O    3x7 matrix, column i+1 = coordinates of origin o_i in C0 (mm)
%
%   D-H table: Assgt #1 Q3 (positive joint rotation right-handed about k_(i-1)).
%   Calls DH_homog.m (Q2b).

if nargin < 2, doPlot = false; end
q = q(:);

%        theta offset   d (mm)   a (mm)   alpha (deg)
DH = [        0,          0,       0,      -90;     % link 1
              0,          0,     431.8,    180;     % link 2
             90,       -149.1,    20.3,     90;     % link 3  (theta3 + 90)
              0,        433.1,     0,       90;     % link 4
              0,          0,       0,      -90;     % link 5
              0,         60.0,     0,        0 ];   % link 6

% ---- Forward kinematics: 0Ti = 0T(i-1) * (i-1)Ti  (Slide 2.4) -----------
T = eye(4);
O = zeros(3,7);          % origins  o0 ... o6
Z = zeros(3,6);          % joint axes k0 ... k5 (in C0 coordinates)
for i = 1:6
    Z(:,i) = T(1:3,3);   % k_(i-1): third column of 0T(i-1) = joint i axis
    O(:,i) = T(1:3,4);   % o_(i-1)
    T = T * DH_homog(q(i) + DH(i,1), DH(i,2), DH(i,3), DH(i,4));
end
O(:,7) = T(1:3,4);       % o6 = gripper centre
T06 = T;

% ---- Manipulator Jacobian (all joints revolute) --------------------------
%   column i = [ k_(i-1) x (o6 - o_(i-1)) ;  k_(i-1) ]
J = zeros(6,6);
for i = 1:6
    J(1:3,i) = cross(Z(:,i), O(:,7) - O(:,i));   % linear-velocity part
    J(4:6,i) = Z(:,i);                           % angular-velocity part
end

% ---- Q2(e): plot the link origins ----------------------------------------
if doPlot
    figure;
    plot3(O(1,:), O(2,:), O(3,:), '-*', 'LineWidth', 1.5, 'MarkerSize', 10);
    hold on; grid on; axis equal;
    for i = 1:7
        text(O(1,i), O(2,i), O(3,i), sprintf('  o_%d', i-1));
    end
    % gripper approach direction k6 for reference
    quiver3(O(1,7), O(2,7), O(3,7), 100*T06(1,3), 100*T06(2,3), 100*T06(3,3), 'r');
    xlabel('x_0 (mm)'); ylabel('y_0 (mm)'); zlabel('z_0 (mm)');
    title(sprintf('Puma 560 link origins, q = [%g %g %g %g %g %g] deg', q));
    view(30, 25); hold off;
end
end

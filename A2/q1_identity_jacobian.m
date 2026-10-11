%%Name: Utkarsh Gautam Kumar
% Student Number: 88300108
% ELEC442 Assignment #2 -- Question 1: manipulator with identity Jacobian
%  P-P-P-R-R-R robot: 3 prismatic joints along i0, j0, k0 (Cartesian gantry)

clear; clc; format short g;

skew = @(s) [0 -s(3) s(2); s(3) 0 -s(1); -s(2) s(1) 0];   % s x  (Slide 1.15)
ii = [1;0;0];  jj = [0;1;0];  kk = [0;0;1];
p_home = [0.4; 0.3; 0.5];          % constant gantry offset of o6 at home (m) -- any value works

% Jacobian for joint vector q = [d1 d2 d3 (m); th4 th5 th6 (rad)]
%   C6 = C0 e^(th4 i x) e^(th5 j x) e^(th6 k x)  and  o6 = o0 + C0 (p_home + [d1;d2;d3])
%   Wrist axes in base coords: k3 = i,  k4 = e^(th4 i x) j,  k5 = e^(th4 i x) e^(th5 j x) k
%   Every wrist axis passes through o6, so (o6 - o_(i-1)) is parallel to k_(i-1)
%   and the linear-velocity part of columns 4-6 is zero.
Jac = @(q) [ eye(3),       zeros(3,3);
             zeros(3,3),   [ ii, ...
                             expm(q(4)*skew(ii))*jj, ...
                             expm(q(4)*skew(ii))*expm(q(5)*skew(jj))*kk ] ];

% Gripper pose (for completeness): 0T6
T06 = @(q) [ expm(q(4)*skew(ii))*expm(q(5)*skew(jj))*expm(q(6)*skew(kk)), p_home + q(1:3);
             zeros(1,3), 1 ];

%% (a) Home position: all joint variables zero -> J must be I
q0 = zeros(6,1);
disp('J(home) ='); disp(Jac(q0));
fprintf('J(home) == I ? %d\n', norm(Jac(q0) - eye(6)) < 1e-12);
disp('0T6(home) = (C6 = C0, o6 at the wrist centre)'); disp(T06(q0));

%% (c) Singularities: hand result det(J) = cos(theta5)
fprintf('\n(c) det(J) vs cos(theta5) for random joint vectors:\n');
rng(1);
for n = 1:4
    q = [randn(3,1); 2*pi*rand(3,1)];
    fprintf('   det(J) = %+8.5f   cos(th5) = %+8.5f\n', det(Jac(q)), cos(q(5)));
end
q_sing = [0.1; -0.2; 0.3; 0.7; pi/2; 1.2];          % theta5 = 90 deg
fprintf('   theta5 = 90 deg -> det(J) = %g, rank(J) = %d  (k3 and k5 both along i0)\n', ...
        det(Jac(q_sing)), rank(Jac(q_sing)));

%% (d) J = I whenever theta4 = theta5 = 0 (mod 360 deg); d1, d2, d3, theta6 free
q_I = [0.25; -0.1; 0.4; 2*pi; -2*pi; 1.3];
fprintf('\n(d) d = [0.25 -0.1 0.4] m, th4 = 360, th5 = -360, th6 = 74.5 deg -> J == I ? %d\n', ...
        norm(Jac(q_I) - eye(6)) < 1e-12);
q_notI = [0; 0; 0; pi; pi; 0];                       % th4 = th5 = 180 deg (not multiples of 360)
fprintf('    th4 = th5 = 180 deg                             -> J == I ? %d\n', ...
        norm(Jac(q_notI) - eye(6)) < 1e-12);

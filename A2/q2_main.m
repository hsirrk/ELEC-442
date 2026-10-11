%%Name: Utkarsh Gautam Kumar
% Student Number: 88300108
% ELEC442 Assignment #2 -- Question 2(c)-(e): Puma 560
%  Computes 0T6 and the manipulator Jacobian J for qA, qB, qC, checks the
%  hand-derived determinant (Q2a), computes the 3 manipulability metrics
%  (Q2d) and plots the link origins (Q2e).
%
%  Files needed in this folder: DH_homog.m, ForKinPuma560.m
%  Run:  q2_main      (writes the session to q2_session.txt via diary)

clear; clc; close all; format short g;
if exist('q2_session.txt', 'file'), delete('q2_session.txt'); end
diary q2_session.txt

% Joint vector sets from the assignment (degrees)
Q = [ 0    0    45;      % theta1
      0    0   -45;      % theta2
      0  -90    45;      % theta3
      0    0     0;      % theta4
      0    0   -30;      % theta5
      0  180    90 ];    % theta6
names = {'qA', 'qB', 'qC'};

% Link lengths (mm) for the closed-form determinant from Q2(a)
a2 = 431.8;  a3 = 20.3;  d4 = 433.1;
detJ_hand = @(q) -a2 * sind(q(5)) ...                                  % wrist term
                 * (a3*cosd(q(3)) + d4*sind(q(3))) ...                 % elbow term
                 * (a2*cosd(q(2)) + a3*sind(q(2)-q(3)) + d4*cosd(q(2)-q(3)));  % shoulder term

% J mixes mm/s (rows 1-3) and rad/s (rows 4-6). To show how much the
% metrics depend on the length unit, also evaluate them with J in metres.
S_m = diag([1e-3 1e-3 1e-3 1 1 1]);

for n = 1:3
    q = Q(:,n);
    [T06, J] = ForKinPuma560(q, true);              % (c) + plot for (e)
    chop = @(M) M .* (abs(M) > 1e-9);               % hide round-off like 1e-14

    fprintf('\n================ %s = [%g %g %g %g %g %g] deg ================\n', names{n}, q);
    disp('0T6 (mm) ='); disp(chop(T06));
    disp('J (rows 1-3: mm/s per rad/s, rows 4-6: rad/s per rad/s) ='); disp(chop(J));
    fprintf('det(J) numeric = %.6g ,  hand formula (Q2a) = %.6g\n', det(J), detJ_hand(q));

    % (d) manipulability metrics (Slide 3.20)
    for u = 1:2
        if u == 1, Ju = J;       unit = 'mm';
        else,      Ju = S_m*J;   unit = 'm ';
        end
        sv    = svd(Ju);                     % singular values, largest first
        w     = sqrt(det(Ju*Ju'));           % Yoshikawa manipulability = prod(sv)
        sigm  = sv(end);                     % minimum singular value
        kinv  = sv(end)/sv(1);               % inverse condition number (0 = singular, 1 = isotropic)
        fprintf('[J in %s]  w = sqrt(det(JJ^T)) = %-11.4g sigma_min = %-11.4g 1/cond = %-11.4g (cond = %.4g)\n', ...
                unit, w, sigm, kinv, sv(1)/sv(end));
    end
end
diary off

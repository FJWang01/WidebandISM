function [ R ] = RZRYRZrad( alpha, beta, gamma )
% RZRYRZDEG Generates the 3x3 'zyz' matrix corresponding the Euler angles
% arguments in deg
% gamma is the first rotation about the z-axis
% beta is the second rotation about the y-axis
% alpha is the third rotation about the z-axis

% alpha, beta and gamma are in radian

% Adapted from http://users.cecs.anu.edu.au/~rod/resources/p-spharm.html
%%
Rz1=[cos(gamma) -sin(gamma) 0; sin(gamma) cos(gamma) 0; 0 0 1];
Ry=[cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)];
Rz2=[cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];
R=Rz2*Ry*Rz1;
end
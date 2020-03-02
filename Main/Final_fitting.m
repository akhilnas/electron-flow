function [ Difference ] = Final_fitting( fitting_factor,OutGenerate, ab_initio,Vbetterguess,A)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Density_Gradient = scsdg(fitting_factor,OutGenerate,Vbetterguess,A);
Difference = norm(ab_initio  - Density_Gradient )*1e06;



end


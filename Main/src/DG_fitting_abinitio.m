function [ Difference ] = DG_fitting_abinitio( fitting_factor,OutGenerate,n_charge_ref,Vbetterguess,A)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Difference = norm(n_charge_ref  - scsdg(fitting_factor,OutGenerate,Vbetterguess,A) )*1e06;


end

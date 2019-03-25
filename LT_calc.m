function LTsig=LT_calc(varargin)
% Joshua Yeh
% 19/03/22
% 
%% Description
% This script calculates the areal chain density using the Lake Thomas
% theory.
% 
%% INPUT VARAIBLES
% 
% NAME PAIR ARGUMENTS: Argolight_gridfit(...'<fieldname>',<value>)
% 'E': Elastic modulus (Pa)
% 
% 'rho': density of the polymer (kg/m^3)
% 
% 'C': characteristic ratio
% 
% 'T': temperature (K)
% 
% 'm': molecular weight of the monomer in the chain (kg/mol)

narginchk(0,inf);
params=inputParser;
params.CaseSensitive=false;
params.addParameter('E',0.4e6);
params.addParameter('rho',1224);
params.addParameter('C',8.1);
params.addParameter('T',298);
params.addParameter('m',86.09e-3);
params.parse(varargin{:});

% Extract out values from parsed input
E=params.Results.E;
rho=params.Results.rho;
C=params.Results.C;
T=params.Results.T;
m=params.Results.m;
kb=1.38e-23;%Boltzmann's constant J/K
Na=6.022e23;%Avogadro's number
l=1.54e-10;%length of C-C bond (m)

LTsig=l*sqrt((E*rho*Na*C)/(6*m*T*kb));%areal chain density (#/m^2)

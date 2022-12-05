%***********************************************************************
%****************** NANOMOS 2.5 ****************************************
%***********************************************************************
function nanomos

close all;
clear all;
fprintf(1,'%s\n','*******************nanoMOS2.5**********************');
fprintf(1,'%s\n','COPYRIGHT, 2000: By Purdue Research Foundation');
fprintf(1,'%s\n','West Lafayette,Indiana 47907.');
fprintf(1,'%s\n', 'All Rights Reserved.');  
fprintf(1,'%s\n','Unless permission is granted, this material');
fprintf(1,'%s\n','not be copied, reproduced or coded for');
fprintf(1,'%s\n','reproduction by any electrical, mechanical ');
fprintf(1,'%s\n','or chemical processes, or combinations thereof');
fprintf(1,'%s\n','now known or later developed.');
fprintf(1,'%s\n','***************************************************');

%***********************************************************************
%********DEFINE GLOBAL VARIABLES FOR INPUT PARAMETERS*******************
%***********************************************************************

global transport_model t_vall max_subband;
global DG_flag fermi_flag ox_pnt_flag dummy_flag;
global Lsd Lg_top Lg_bot t_top t_bot t_si dx dy refine;
global dopslope_s dopslope_d overlap_s overlap_d;
global Vg1 Vg2 Vs Vd Vd_initial Vg_step Ng_step Vd_step Nd_step;
global phi_top phi_bot eps_top eps_bot bar_top bar_bot;
global eps_si m_t m_l Te;
global mu_low beta Vel_sat;
global N_sd N_body criterion_outer criterion_inner;
global Ng_step Nd_step
global Nx Ny Ntotal;
global filename dirname;
global N_dos Trans E;

%Energy transport parameters

global ELE_TAUW ELE_CQ;

% options for ploting

global plot_IV plot_Ec3d plot_Ne3d plot_Ecsub plot_Te plot_Nesub;
global plot_Ec_IV plot_Ne_IV;
%************************************************************************
%**************** PARSER READ THE INPUT FILE ****************************
%************************************************************************

readinput;

%************************************************************************
%*************CALL THE MAIN ROUTINE TO DO THE COMPUTATION ***************
%************************************************************************

ttt=cputime;
[Ie,Ie_sub,Te_sub,Ne_sub,E_sub,Ne,Ec,conve]=main;
t_total=cputime-ttt

%************************************************************************
%************************SAVE AND PLOT OUTPUT****************************
%************************************************************************

saveoutput(Ec,Ne,Ie,Ne_sub,E_sub,Te_sub,conve)

%************************************************************************
%*********************** THE END OF NANOMOS *****************************
%************************************************************************

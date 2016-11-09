

disp('%%% ****************************************************************************************')

check_times=false
check_times=true

dh=5.0;

scrivi_altezze=true;
%  scrivi_altezze=false;
scrivi_vel=true;
%  scrivi_vel=false;

tmax_ch1=nan;
t_limits_ch23=nan;
echo_min_pcorr = nan;
echo_min_pdist = nan;
echo_min_pdop2 = nan;
v_max_acc1 = nan;
v_max_acc2 = [-1000 1000];
v_max_acc3 = [-1000 1000];
echo_min_pdop2 = 150

%  parametri_analisi_v_dop;

perc_size=0.8;
perc_size_no = 0.3;
echo_min_pcorr = 900;

skip_ch=[ 1];
skip_ch=[ ];



#
# tmax_ch1=11000;
# v_max_acc1 = [ -2000 1000];
# t_limits_ch23=[3080 7800];
# echo_min_pdist = 800;
%  v_max_acc3 = [-1000 1000];
%  v_max_acc2 = [-1000 1000];
%  %  forza_h_min_mis = 13.5
%  %  ti=8200;
%  %  tf=8300;
%  %  gi=26;
%  %  gf=28;
%  %  g_dist=20
%  %  echo_min_pdop2 = 250
%
%  %  ti=4540;
%  %  tf=4630;
%  %  yi=100;
%  %  yf=140;
%

%%% -------------------------------------------

%  loadsi=true;
loadsi=false;


%%% -------------------------------------------

path_base = ''
path_in = 'SampleTest/';

dist_dop_rampa(2) = 1.7;   %% strato plexiglass molto sottile
dist_dop_rampa(3) = 1.2;   %% strato plexiglass spesso qualche mm



path_out=path_in

analizza_vel_dop;
analizza_dop;





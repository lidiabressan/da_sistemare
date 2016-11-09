%  Lidia Bressan
%  gennaio 2015
%
%  Legge la matrice dei dati dei DOP
%  separa l'echo dai dati di velocità (output)
%  calcola la velocità media del profilo vel_media
% -----------------------------------


1;

%  %  Per lettura dati e input

function [time_ms, vel_profiles, echo_profiles, x_profiles] = extract_data_dop(data, xdata)

    if(length(find(xdata==xdata(1)))==2)
        ii=find(xdata==xdata(1));
        ngates = ii(2)-1;
        x_profiles = xdata(1:ngates);
        echo_profiles = data(:,ngates+1:2*ngates);
        vel_profiles = data(:,1:ngates);
    else
        ngates = length(xdata);
        x_profiles = xdata;
        vel_profiles = data(:,1:ngates);
        echo_profiles = [];
    end

    time_ms = data(:,end-3);
    flow_rate_mlmin = data(:,end-2);
    n_trigger_seq = data(:,end-1);
    ch = data(:,end);
end




function filess = riordina_files(filess)     ## ok
    num =[];
    for ii=1:size(filess, 1)
        num = [num; trova_num(filess(ii,:))];
    end
    [n_sort, ii] = sort(num);
    filess = filess(ii, :);
end



function ff_analisi = det_ff_analisi(H_dop, selez_H)
    if(selez_H~=0)
        ff_analisi = find(H_dop>=selez_H-0.5 & H_dop<=selez_H+0.5);
    else
        ff_analisi = 1:size(H_dop,1);
    end
end



function num_file = trova_num(filee)
    ii=strfind(filee, '.ADD')-1;
    kk=strfind(filee, '_')(end)+1;
    num_file = str2num(filee(kk:ii));
end



function ch_g = det_dop(xgates)
%      max(xgates)
    if(max(xgates)>100.0) % mm - 1MHz - grosso
        ch_g=3;
    elseif(max(xgates)<50.0) % mm - 4MHz - sottile
        ch_g=2;
    elseif(max(xgates)<100.0) % mm - 1MHz - grosso
        ch_g=1;
    end
end



function x_pos = get_dop_position(ch_g, id_set)
    if(id_set<3)
        x_punta_dop(1) = (413.0+26.0)/100.0+2.005;   %% grosso - 1MHz
        x_punta_dop(2) = (413.0+26.0-1.5)/100.0+2.005;  %% sottile - 4MHz
        x_punta_dop(3) = 0.0;
    elseif(id_set==3)
        x_punta_dop(1) = 247.5/100.0+2.005;   %% sottile - 4MHz
        x_punta_dop(2) = 364.3/100.0+2.005;   %% sottile - 4MHz
    end
    x_pos = x_punta_dop(ch_g);
end









%  %  Pulitura dati



function v_profile_clean = clean_outliers(v_profile, kdev)
    %% togli spikes
    v_profile_clean = v_profile;
    N = length(v_profile_clean);
    non_out = [1:N];
    non_zero = find(v_profile_clean~=0.0);
    if(length(non_zero)>=0.6*N)
        buoni = non_zero;
    else
        buoni = non_out;
    end

    med_p = mean(v_profile_clean);
    std_p = std(v_profile_clean);
    outliers = [];



%      outliers2 = [0];
%      iqr_d = 100;
%      while(~isempty(outliers2) & iqr_d>0 & length(buoni)>1)
%          mediana=median(v_profile_clean(buoni));
%          iqr_d = 5.0*iqr(v_profile_clean(buoni));
%          qestr = quantile(v_profile_clean(buoni), [0.25, 0.75]);
%          outliers2 = find(v_profile_clean(buoni)<qestr(1)-iqr_d | v_profile_clean(buoni)>qestr(2)+iqr_d);
%          buoni = setdiff(buoni, buoni(outliers2));
%          outliers = [outliers outliers2];
%      end


    while (any(abs(v_profile_clean(buoni)-med_p)>kdev.*std_p))
        outliers = find(abs(v_profile_clean(buoni).-med_p)>kdev.*std_p);
        buoni = setdiff(buoni, buoni(outliers));
        med_p = mean(v_profile_clean( buoni));
        std_p = std(v_profile_clean( buoni));
    end

%      v_profile_clean(outliers)=0;
    v_profile_clean(outliers)=nan;
end




function vel_profile = clean_v_from_echo(echo_profile, echo0, echo0_std, vel_profile, ndev)
%  %  Pulisci in base ai dati di echo
%      ndev = 4.0;
    non_zero=find(echo_profile>0.0);
    xx=find(abs(echo_profile(non_zero).-echo0(non_zero))> ndev.*echo0_std(non_zero));
    sum2diff = sum((echo_profile.-echo0).^2);
%  %      if( length(xx)<=2 | sum2diff<ndev*sum(echo0_std.*echo0));
    if( length(xx)<=2 | sum2diff<ndev*sum(echo0_std));
        vel_profile = 0.0;
    end
%
end




function vel_profile = clean_vneg_vprof(vel_profile)
%  %  Pulisci in base al verso della corrente (nel mio caso sempre negativa)
    ngates = length(vel_profile);
    nonzero = ~(vel_profile==0.0);
    n_mis_ok = sum(nonzero);
    vel_pos = find(vel_profile>0);
    vel_neg = find(vel_profile<0);
    nmin_ok = n_mis_ok - floor(n_mis_ok*0.2);

    if(length(vel_neg)>= nmin_ok && nmin_ok>1)
        dati_out = find(vel_profile >0);
    elseif(length(vel_pos)==nmin_ok && nmin_ok>ngates*0.7)
        dati_out = find(vel_profile<0);
    else
        dati_out = [];
    end
    vel_profile(dati_out) = 0.0;
end



function profilo_velocity = profilo_vel_no0(vel)
    zeri=find(vel==0);
    profilo_velocity = vel;
    if(length(zeri)<length(vel))
        profilo_velocity(zeri)=nan;
    end
end





function [vel_profiles, profiles_vel] = clean_data_dop(vel_profiles, echo_profiles, echo0, echo0_std, xdata, tempo, ndev, check)
    ngates = size(vel_profiles, 2);
    profiles_vel = nan(size(vel_profiles));
    %%      togli gli zeri a echo0_std
    if(length(find(echo0_std==0))>0 & length(find(echo0_std==0))<length(echo0_std))
        estd_min = min(echo0_std(find(echo0_std>0)));
        echo0_std(find(echo0_std==0)) = estd_min;
    elseif(length(find(echo0_std==0))==length(echo0_std))
        echo0_std = min(echo0);
    end
    echo0_min = max(echo0.-ndev.*echo0_std, 0.0);
    echo0_max = echo0.+ndev.*echo0_std;

    if(check)
        vel_media = nan(size(vel_profiles, 1),1);
        scrsz = get(0,'ScreenSize');
        fig_i = figure(2, 'Position',[50 100 0.95*scrsz(3) 0.35*scrsz(4)]);
        clf;
        fig_check=check;
    else
        fig_check=false;
    end

%  %  Togli gli spikes
    for ii = 1:size(vel_profiles, 1)
        vel = vel_profiles(ii,:);
        echo = echo_profiles(ii,:);

        if(fig_check)
            hold off;
    %       plot echo
            subplot(1,3,2); title('echo intensity (V?)'); hold off;
            plot(echo, xdata, 'linewidth', 4); hold on;
            plot(echo0, xdata, '-r', 'linewidth', 0.5)
            plot(echo0_min, xdata, '-r', 'linewidth', 1)
            plot(echo0_max, xdata, '-r', 'linewidth', 1)
            grid on
            xlabel('echo')
            ylabel('x (mm)')
%              axis([0 1200])
    %  %    velocity profiles
            tempo_str=['tempo (ms) ' num2str(tempo(ii))]
            subplot(1,3,1); title('velocità (m/s)'); hold off;
            plot(vel./1000, xdata, 'linewidth', 4); hold on
            text(0.0, 5, tempo_str);
            axis([-2.0 2.0])
            grid on
            xlabel('v (m/s)');
            ylabel('x (mm)');
        end

%          vel = clean_v_from_echo(echo, echo0, echo0_std, vel, ndev);
        vel = clean_outliers(vel, kdev=ndev);

%  %  %          vel_profile = clean_vneg_vprof(vel_profile);
        vel_profiles(ii,:) = vel;

        profiles_vel(ii,:) = profilo_vel_no0(vel);

        if(fig_check)
            vel_media(ii) = calc_vm_profile(vel);
            plot(vel./1000, xdata, '2', 'linewidth', 1); hold on
            v_media=vel_media(ii);
            %%%% v media
            plot(vel_media(ii)./1000.*[1 1], [min(xdata) max(xdata)], '1', 'linewidth', 1); hold on
            %%%% v media
            subplot(1,3,3); title('velocità (m/s)'); hold off;
            plot(tempo, vel_media./1000)
            xlim([0 max(ceil(tempo(ii)/1000)*1000, 1)])

            n_prof = size(profiles_vel,2);
            colori = rainbow(n_prof);
            for kk=1:n_prof
%                  plot(tempo([1:ii]), profiles_vel(1:ii,kk),'-+', 'color', colori(kk,:), 'linewidth',2); hold on;
                plot([1:ii], profiles_vel(1:ii,kk),'-+', 'color', colori(kk,:), 'linewidth',2); hold on;
            end
            xlabel('ii dati')
            pause()
        end
        if(tempo(ii)<2000 |tempo(ii)>9000)
            fig_check=false;
        elseif(check)
            fig_check=true;
        end
    end
end



%  %  Analisi dei dati dop



function ii_arrivo = trova_t_arrivo_echo(echo_profiles, echo0, echo0_std, ndev)

    for ii=1:size(echo_profiles,1)
        non_zero=find(echo_profiles(ii,:)>0.0);
        xx=find(abs(echo_profiles(ii,non_zero).-echo0(non_zero))> ndev.*echo0_std(non_zero));
        if(length(xx)>3 && abs(echo_profiles(ii,non_zero).-echo0(non_zero))./echo0_std(non_zero)>10.0)
            if(ii_arrivo==0)
                ii_arrivo = ii;
                ii_check = ii_arrivo;
            elseif(ii_check-ii_arrivo==10)
                ii_arr = [ii_arr; ii_arrivo];
                break;
%                  ii_arrivo = 0;
            elseif(ii_arrivo>0 && ii-ii_check==1 && ii_check-ii_arrivo<10)
                ii_check = ii;
            end
        else
            ii_check = 0;
            ii_arrivo = 0;
        end
    end
end



function ii_arrivo = trova_t_arrivo_vel(vel_data, v_lim_mms)
    check_ii = find(sum(abs(vel_data)>v_lim_mms,2)>0.3*size(vel_data,2));
    check_ii = find(sum(abs(vel_data)>v_lim_mms,2)>10);
    if(~isempty(check_ii))
        check_ii = check_ii(1:min(4,end));
        if(max(diff(check_ii))==1)
            ii_arrivo = check_ii(1);
        else
            ii_arrivo = check_ii(min(2,length(check_ii)));
        end
    else
        ii_arrivo = size(vel_data,1);
        ii_arrivo = 1;
    end
end






%  %  Processing dei dati -- ricavare la superficie libera


function filtro_pdist = filtro_E_pdist(ee_profilo_ch, echo_min_pdist)
    %% filtri
    filtro_pdist = find(ee_profilo_ch<echo_min_pdist | isnan(ee_profilo_ch));
end



function x_prof_mass = costr_matrix_Emass(x_prof, ee_profilo_ch, echo_min_pdist)
        x_prof_mass = repmat(x_prof, size(ee_profilo_ch,1), 1);
        filtro_pdist = filtro_E_pdist(ee_profilo_ch, echo_min_pdist);
        x_prof_mass(filtro_pdist) = nan;
end




function [dist_min_u, dist_min] = ricava_h_rifless_E(x_prof, x_prof_mass, perc_size, perc_size_no)
        %% ricava le altezze della base
        dist_min = [];
        dist_min_u = [];

        ii_min30 = find(x_prof<30);
        l_gates=size(x_prof_mass,2);
        check_nn = sum(isnan(x_prof_mass),2);
        ii_lim = find(check_nn<l_gates)([1, end]);
        x_prof_mass_s =x_prof_mass(ii_lim(1):ii_lim(2), ii_min30);
        l_gates=size(x_prof_mass,2);
        l_tempi = size(x_prof_mass_s,1);
        ii=1;
        for ii_prof=ii_min30
            nonan = find(~isnan(x_prof_mass_s(:,ii_prof)));
            nonanp = find(~isnan(x_prof_mass(:,ii_prof+1)));
            if(length(nonan)>perc_size*l_tempi & length(nonanp)<perc_size_no*l_tempi)
                dist_min(ii) = x_prof(ii_prof);
                dist_min_u(ii) = x_prof(ii_prof+1);
                ii=ii+1;
            end
        end

%  dist_min

        %%%%% check
%          ii_centrali = round(size(x_prof_mass,1)/2).+[-5:5];
        insieme = [];
        cerca_nans = [];
        for ii_t=1:size(x_prof_mass_s,1)
            nans = find(isnan(x_prof_mass_s(ii_t,:)));
            cerca_nans = [cerca_nans; min(nans)];
%              nonan = find(~isnan(x_prof_mass_s(ii_t,:)));
%              insieme = [insieme; max(nonan)];
        end
%          ii_sel = ceil(median(insieme));
        ii_sel = ceil(median(cerca_nans));
        dist_min(ii) = x_prof(ii_sel);
        dist_min_u(ii) = x_prof(ii_sel+1);

%  dist_min


end










%  %  Processing dei dati -- calcolo della media




function [echo_vuoto_mean, echo_vuoto_std] = calc_echo_vuoto(echo_profiles, time_ms, t_max)
    echo_vuoto = echo_profiles(find(time_ms<=t_max),:);
    echo_vuoto_mean = mean(echo_vuoto, 1);
    echo_vuoto_std = std(echo_vuoto, 1, 1);
end



function times = det_new_tempi(tempo, passo_t)
    times = sort(unique(floor(tempo/passo_t)*passo_t)).+0.5*passo_t;
end



function x_profilo = det_new_gates(x_gates, passo_cm)
    x_profilo = sort(unique([floor(min(x_gates)*100.0)/100.0:passo_cm:ceil(max(x_gates*100.0))/100.0]));
end



function vel_masked = schermazeri_v_prof(ii_iniz, velocity)
    zeri = find(velocity(ii_iniz:end)==0);
    vel_masked = velocity;
    vel_masked(zeri+ii_iniz-1) = nan;
end



function vel_npos = filtra_val_pos(vel_data)
    vel_npos = vel_data;
    vel_npos(find(vel_npos>0)) = 0.0;
end



function vel_masked = schermazeri_vel(ii_iniz, velocity)
    if(size(velocity,2)==1)
        vel_masked = schermazeri_v_prof(ii_iniz, velocity);
    else
        vel_masked = velocity;
        vparz = velocity(ii_iniz:end,:);
        zeri = find(vparz==0);
        vparz(zeri) = nan;
        vel_masked(ii_iniz:end,:) = vparz;
    end
end



function v_media = calc_vm_profile(vel_profile)
    ngates = size(vel_profile, 2);
    nonzero = find(vel_profile~=0.0);
    if(length(nonzero)>ngates*0.1)
        v_media = mean(vel_profile(nonzero));
    else
        v_media = 0.0;
    end
end



function v_media = calc_vm_profile_nan(vel_profile)
    ngates = size(vel_profile, 2);
    nonnan = find(~isnan(vel_profile));
    if(~isempty(nonnan))
        v_media = mean(vel_profile(nonnan));
    else
        v_media = 0.0;
    end
end



function vel_media_st = calc_vm_smt(tempo, velocity, tempi_nuovi, k_smooth)
    if(k_smooth==0)
        vel_media_st = interp1(tempo, velocity, tempi_nuovi, 'cubic');
    else
        nonnan = find(~isnan(velocity));
        vel_media_st = csaps(tempo(nonnan), velocity(nonnan), k_smooth, tempi_nuovi);
    end
end




function v_medie_sm = calc_vm_sms(x_profiles, vel_data, x_profilo_nuovo)

    dx_medio=0.5*mean(diff(x_profiles));
    xlim = [x_profiles(1)-dx_medio x_profiles(1:end-1)+0.5.*diff(x_profiles) x_profiles(end)+dx_medio];

    ii_prof=1;
    v_medie_sm=nan(size(vel_data,1), length(x_profilo_nuovo));
    for cm = x_profilo_nuovo
        ii_prof=ii_prof+1;
        xcm = find(x_profiles>=xlim(ii_prof) & x_profiles<xlim(ii_prof+1));
        if(~isempty(xcm))
            for ii=1:size(vel_data,1)
                nonan=find(~isnan(vel_data(ii, xcm)));
                if(~isempty(nonan))
                    v_medie_sm(ii,ii_prof) = mean(vel_data(ii, xcm(nonan)));
                else
                    v_medie_sm(ii,ii_prof) = nan;
                end
            end
        end
    end
end


function v_media = calc_vm_profile_prova(vel_profile)
    ngates = size(vel_profile, 2);
    nonzero = find(vel_profile~=0.0);
%      vel_neg = find(vel_profiles(ii,:)<0.0);

    if(length(nonzero)>ngates*0.4)
        v_media = mean(vel_profile(nonzero));
%          elseif(~isempty(vel_neg))
%              vel_media(ii) = mean(vel_profiles(ii,vel_neg));
    else
        v_media = 0.0;
    end

%          n_mis_ok = sum(nonzero);
%          vel_pos=find(vel_profiles(ii,:)>0);
%          vel_neg=find(vel_profiles(ii,:)<0);
%          nmin_ok = n_mis_ok - floor(n_mis_ok*0.2);
%
%          if(length(vel_neg)>=nmin_ok && nmin_ok>1)
%              dati_buoni = find(vel_profiles(ii,:)<0);
%          elseif(length(vel_pos)==nmin_ok && nmin_ok>ngates*0.7)
%              dati_buoni = find(vel_profiles(ii,:)>0);
%          else
%              dati_buoni = [];
%          end
%
%          if(~isempty(dati_buoni))
%              vel_media(ii) = mean(vel_profiles(ii,dati_buoni));
%  %          elseif(nmin_ok<=1)
%  %              vel_media(ii) = 0.0;
%          elseif(sum(nonzero)>0)
%              vel_media(ii) = mean(vel_profiles(ii,nonzero));
%          else
%              vel_media(ii) = 0.0;
%          end

end



function vel_media = calc_vm_profiles(vel_profiles)
    vel_media = zeros(size(vel_profiles, 1),1);
    for ii=1:size(vel_profiles, 1)
        vel_media(ii) = calc_vm_profile(vel_profiles(ii,:));
     end
end




function v_media = clean_v_media(v_media)
    v_media(v_media>0) = 0.0;
end





%  %  Per figure


function plot_v_medie(tempo_ms, vel_media_ms, tempo_mst, vel_media_mst)

    c_vari = ['r', 'b', 'g', 'y', 'm', 'c', 'k', 'r', 'b', 'g', 'y', 'm', 'c', 'k'];
    figure(11); clf;
    hold on;
    plot(tempo_ms, vel_media_ms, ';ms;', 'color', c_vari(3), 'linewidth', 3)
    plot(tempo_mst, vel_media_mst, '-;mst;', 'color', c_vari(1), 'linewidth', 1)
%      plot(tempo_smts, vel_media_smts, ';smts;', 'color', c_vari(2), 'linewidth', 1)
    grid on
    legend()
    xlabel('tempo (ms)')
    ylabel('V media profilo (mm/s)')
end




function  colori=col_prof(x_prof)
    dh = 10.0;
    di=10.0;
    h = [di:dh:di+dh*6];
    if(x_prof<h(1))
        colori = [1 0 0];
    elseif(x_prof<h(2))
        colori = [0 1 0];
    elseif(x_prof<h(3))
        colori = [0 0 1];
    elseif(x_prof<h(4))
        colori = [1 0 1];
    elseif(x_prof<h(5))
        colori = [1 1 0];
        colori = [0 0 0];
    elseif(x_prof<h(6))
        colori = [0 1 1];
    else
        colori = [0 0 0];
    end
end




function fig_i = fig_profili_tempo(tempo, xdata, vel_data, v_media, echo_data, empty_echo, empty_echo_std)
    scrsz = get(0,'ScreenSize');
    fig_i = figure(2, 'Position',[50 100 0.7*scrsz(3) 0.35*scrsz(4)]);
    clf;
    for ii=1:size(vel_data,1)
        fig_i = fig_profili_t1(tempo(ii), xdata, vel_data(ii,:), v_media(ii), echo_data(ii,:), empty_echo, empty_echo_std)
        non_zero=find(echo_data(ii,:)>0.0);
        for ndev = [4.0 ]
            xx=find(abs(echo_data(ii,non_zero).-empty_echo(non_zero))> ndev.*empty_echo_std(non_zero));
            if( length(xx)>=2 )
%                  [ sum((abs(echo_data(ii,non_zero).-empty_echo(non_zero))> ndev.*empty_echo_std(non_zero))) ...
%                  max(abs(echo_data(ii,non_zero).-empty_echo(non_zero))./empty_echo_std(non_zero)) ]
%                  [xdata(xx)' empty_echo(xx)' empty_echo_std(xx)' echo_data(ii,xx)'  ]
                pause()
            end
        end
    end
end



function fig_i = fig_profili_t1(tempo, xdata, vel_data, v_media, echo_data, empty_echo, empty_echo_std)
    scrsz = get(0,'ScreenSize');
    fig_i = figure(2, 'Position',[50 100 0.7*scrsz(3) 0.35*scrsz(4)]);
    clf;

    tempo_str=['tempo (ms) ' num2str(tempo)];
    subplot(1,2,1); title('velocità (m/s)'); hold off;
    plot(vel_data./1000, xdata, 'linewidth', 4); hold on
    plot(v_media./1000.*[1 1 ], [min(xdata) max(xdata)], '-r', 'linewidth', 2)

    axis([-2.0 2.0])
    text(0.0, 5, tempo_str);
    grid on
    xlabel('v (m/s)');
    ylabel('x (mm)');

    empty_echo_min = min(empty_echo.-3.0.*empty_echo_std, 0.0);
    empty_echo_max = empty_echo.+3.0.*empty_echo_std;

    subplot(1,2,2); title('echo intensity (V?)'); hold off;
    plot(echo_data, xdata, 'linewidth', 4); hold on;
    plot(empty_echo, xdata, '-r', 'linewidth', 2)
    plot(empty_echo_min, xdata, '-r', 'linewidth', 1)
    plot(empty_echo_max, xdata, '-r', 'linewidth', 1)

    grid on
    xlabel('echo')
    ylabel('x (mm)')
    axis([0 1200])
end


function ax = plot_echo(ax, echo_data, xdata, empty_echo, empty_echo_min, empty_echo_max)
    ax = subplot(1,2,2); title('echo intensity (V?)');
    plot(echo_data, xdata, 'linewidth', 4); hold on;
    plot(empty_echo, xdata, '-r', 'linewidth', 2)
    plot(empty_echo_min, xdata, '-r', 'linewidth', 1)
    plot(empty_echo_max, xdata, '-r', 'linewidth', 1)

    grid on
    xlabel('echo')
    ylabel('x (mm)')
    axis([0 1200])
end



function plot_vel(xdata, vel_profile, lw)
    subplot(1,2,1); title('velocità (m/s)');
    plot(vel_profile./1000, xdata, 'linewidth', lw);

    axis([-2.0 2.0])
    grid on
    xlabel('v (m/s)');
    ylabel('x (mm)');
end



%  %  Per output


function scrivi_vel_m_prof(tempo_pm, vm_profilo, x_profilo, pathfile)
    tv_profiles=[[0.0 x_profilo]; tempo_pm vm_profilo];
    fid = fopen (pathfile, 'w');
    formato = '%10.3f  ';
    for ii=1:length(x_profilo)
        formato=[formato '   %12e'  ];
    end
    formato=[formato '\n'  ];
    for ii=1:size(tv_profiles,1)
        fprintf(fid, formato, tv_profiles(ii,:));
    end
    fclose(fid);
end



function scrivi_vm(tempo, vel_media, ch, filename, x_gauge)
    if(~isempty(ch))
        ch = tempo.*0.+ch;
        formato = '%10.3f   %12e   %d\n';
    else
        formato = '%10.3f   %12e\n';
    end
    t_velm = [tempo vel_media ch];
    [fid, msg] = fopen(filename, 'w');
    if(fid==-1)
        disp(msg)
    end
    fprintf(fid, ['# x_gauge ' num2str(x_gauge) '\n']);
    for ii=1:size(t_velm,1)
        fprintf(fid, formato, t_velm(ii,:));
    end
    fclose(fid);
end



function [arr1, arr2, arr3] = sort_a1(arr1, arr2, arr3)
    [arr_ord, ii] = sort(arr1);
    arr1 = arr1(ii);
    arr2 = arr2(ii);
    arr3 = arr3(ii);
end






function componi_vett_v



end




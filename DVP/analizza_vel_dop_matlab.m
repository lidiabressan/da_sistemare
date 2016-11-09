%  Lidia Bressan
%  gennaio 2015
%
%  Legge la matrice dei dati dei DOP
%  separa l'echo dai dati di velocità (output)
%  calcola la velocità media del profilo vel_media
% -----------------------------------


function funz = analizza_vel_dop_matlab()
    funz.extract_data_dop = @extract_data_dop;
    funz.det_dop = @det_dop;
    funz.trova_t_arrivo_vel = @trova_t_arrivo_vel;
    funz.calc_echo_vuoto = @calc_echo_vuoto;
    funz.costr_matrix_Emass = @costr_matrix_Emass;
    funz.ricava_h_rifless_E = @ricava_h_rifless_E;
    funz.col_prof = @col_prof;
    disp('loading')
end



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
        outliers = find(abs(v_profile_clean(buoni)-med_p)>kdev.*std_p);
        buoni = setdiff(buoni, buoni(outliers));
        med_p = mean(v_profile_clean( buoni));
        std_p = std(v_profile_clean( buoni));
    end

%      v_profile_clean(outliers)=0;
    v_profile_clean(outliers)=nan;
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
    echo0_min = max(echo0-ndev.*echo0_std, 0.0);
    echo0_max = echo0+ndev.*echo0_std;

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

        vel = clean_outliers(vel, ndev);

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
end










%  %  Processing dei dati -- calcolo della media




function [echo_vuoto_mean, echo_vuoto_std] = calc_echo_vuoto(echo_profiles, time_ms, t_max)
    echo_vuoto = echo_profiles(find(time_ms<=t_max),:);
    echo_vuoto_mean = mean(echo_vuoto, 1);
    echo_vuoto_std = std(echo_vuoto, 1, 1);
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



%  %  Per figure




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



%  Lidia Bressan
%  gennaio 2015
%
%  Legge e analizza i dati DOP di un esperimento in sequenza
%
%  --------------------------------


disp('%%% -------------   INIZIO PROGRAMMA   ------------- ')

dop_vtest = [2]
passo_t=0.05
passo_cm=1.0
stringaDaCercare='*ADD';
%  tmax_taglia = 30000;


passo_t = passo_t*1000;

%  %  controlla dati in input o definisci  default

if(~exist('path_in'))
    error('Determinare path_in !!')
end
if(~exist('path_out'))
    if(~isempty(strfind(path_in, 'Test/TestFatti')))
        path_out = strrep(path_in, 'Test/TestFatti', 'AnalisiDati')
    else
        path_out='';
        warning('Default path_out !!')
    end
end
if(~exist('check_ch'))
    check_ch=[];
end
if(~exist('check_times'))
    check_times=false;
end
if(~exist('ndev'))
    ndev=4;
end

%  per il calcolo dell'echo
if(~exist('calc_echo0'))
    calc_echo0=true;
end
if(~exist('t_max'))
    t_max=2000;
end
if(~exist('scrivi_echo0'))
    scrivi_echo0=false;
end
if(~exist('t_profiles'))
    loadsi=true;
end
if(~exist('smooth_iter'))
    smooth_iter=false;
end
%  per le figure dei risultati
if(~exist('figura'))
    figura=false;
end


%  %  programma
%  stringaDaCercare='*ADD';


%  %  carica le funzioni
analizza_vel_dop;



%  Leggi e raccogli tutti i dati di una prova
if(loadsi)
    id_set = det_setting(path_in);
    [filess, nfiles] = find_files(path_in, stringaDaCercare, opzioni=' -L ');
    filess
    if(nfiles==0)
        error('Nessun file trovato')
    elseif(id_set>=4 && exist('file_h_xdop'))
        h_xdop = load(file_h_xdop);
    else
        h_xdop = [];
    end

    [filess_out, nfiles] = find_files(path_out, "v_m_std_dp_rampa_dh*.dat");
    if(nfiles>0)
        disp('')
        warning("File di analisi dati già presenti. Rianalizzare?")
        cont = input('');
        if(~cont)
            exit;
        end
    end


    path_in
    filess = riordina_files(filess);
%      [size(filess,1)]
    filess

    [v_profiles, e_profiles, t_profiles, x_gates, ndati, all_ch] = read_files_dop_exp(path_in, filess, n_ch_max=3);

    for ch=all_ch
        x_profiles = x_gates{ch};
%         x_vel{ch} = x_gates{ch};
%         ch_g = det_dop(x_profiles);
%         x_vel{ch} = get_dop_position(ch_g, id_set).-x_profiles*10**-3;

        %%% Limita i dati a tmax_taglia
        t = t_profiles{ch};
        v = v_profiles{ch};
        e = e_profiles{ch};
%  %          ii_end = max(find(t<=tmax_taglia));
%  %          ndati(ch) = ii_end;
%  %          t_profiles{ch} = t(1:ii_end);
%  %          v_profiles{ch} = v(1:ii_end,:);
%  %          e_profiles{ch} = e(1:ii_end,:);
    end
    ndati=ndati(1:max(all_ch));
    disp(['ndati '  num2str(ndati)])
    disp('... Loading dei dati fatto')
end
if(loadsi)
    puliturasi=true;
end
loadsi=false;





if(id_set<3)
    ch_boulders=all_ch;
elseif(id_set==3 | id_set==4)
    ch_boulders=[2];
elseif(id_set>=5)
    ch_boulders=[2 3];

%%% setting del canale
    ang_assex_rampa = asin(1/10); %%% radianti
    ang_rampa_dop20 = 2*pi*110/360;
    ang_rampa_dop00 = pi/2;
    ang_assex_dop20 = ang_assex_rampa + ang_rampa_dop20;
    ang_assex_dop00 = ang_assex_rampa + ang_rampa_dop00;
    ang_prampa_dop(2) = 2*pi*20/360;
    ang_prampa_dop(3) = 2*pi;
    ang_dop00_dop20 = 2*pi*20/360;
end



%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  Processing dei dati
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %

pkg load statistics


%% Analisi per il primo canale e poi per i canali 2 e 3
for ch=[1 2]
ch
    %%% determino i tempi di arrivo
    vel_data = v_profiles{ch};
    echo_data = e_profiles{ch};
    t_ch = t_profiles{ch};
%  %   calcolo del tempo di arrivo per canale
    if(ch==2)
        v_lim_mms=20;
        for ch=[2 3]
            vel_data = v_profiles{ch};
            echo_data = e_profiles{ch};
            t_ch = t_profiles{ch};
            ii_arrivo_v = trova_t_arrivo_vel(vel_data, v_lim_mms);
            it_arrivo(ch,:) = [ ii_arrivo_v, t_ch(ii_arrivo_v)];
    %  %      ii_arrivo_e = trova_t_arrivo_echo(echo_data, empty_echo, empty_echo_std, ndev=3.0);
    %  %      it_arrivo_e = [ ii_arrivo_e, t_ch(ii_arrivo_v)];
        end
        it_arrivo([2 3],:)

    elseif(ch==1)
        v_lim_mms=100;
        ii_arrivo_v = trova_t_arrivo_vel(vel_data, v_lim_mms);
        it_arrivo(ch,:) = [ ii_arrivo_v, t_ch(ii_arrivo_v)];
%  %      ii_arrivo_e = trova_t_arrivo_echo(echo_data, empty_echo, empty_echo_std, ndev=3.0);
%  %      it_arrivo_e = [ ii_arrivo_e, t_ch(ii_arrivo_v)];
        disp(['calcolo dei tempi di arrivo = ' num2str(it_arrivo(ch,2))])
%  %          it_arrivo(ch,:)
    end

    if(~isnan(t_limits_ch23))
        for chx=[2 3]
            t_ch = t_profiles{chx};
            it_arrivo(chx,1)=find(t_ch>=t_limits_ch23(1))(1);
            it_arrivo(chx,2)=t_ch(find(t_ch>=t_limits_ch23(1))(1));
        end
    end


    disp('Fai le figure per il check dei dati')
    if(ch==1)
        vm_profilo_ch = v_profiles{ch};
        ee_profilo_ch = e_profiles{ch};
        times = t_profiles{ch};

        figure(1); clf; hold on; grid on
        plot(it_arrivo(ch,2).*[1 1], [-1000 1000], 'r')
        plot(times, mean(vm_profilo_ch,2), 'linewidth',2)
        plot(times, vm_profilo_ch)
    else
        ch=2
            vm_profilo_ch = v_profiles{ch};
            ee_profilo_ch = e_profiles{ch};
            times = t_profiles{ch};
            x_prof=x_gates{ch};

            yM=max(max(vm_profilo_ch));
            ym=min(min(vm_profilo_ch));
            figure(20); clf; hold on; grid on
            subplot(2,1,1); hold on; grid on
            plot(it_arrivo(2,2).*[1 1], [ym yM], 'r')
            plot(it_arrivo(3,2).*[1 1], [ym yM], 'g')
            plot(times, vm_profilo_ch)
            subplot(2,1,2); hold on; grid on
            imagesc(times([1,end]), x_prof([1,end]), vm_profilo_ch')
            colorbar;
            axis([times([1,end])' x_prof([1,end])])

            if(~check_times)
                figure(21); clf; hold on;
    %              cmap=colormap(jet(256));
    %              imshow(flipud(ee_profilo_ch'), [0 max(max(ee_profilo_ch))], cmap)
    %              image(times([1,end]), x_prof([1,end]), flipud(ee_profilo_ch'), [0 max(max(ee_profilo_ch))])
                imagesc(times([1,end]), x_prof([1,end]), ee_profilo_ch')
                title(["max = " num2str(max(max(ee_profilo_ch)))])
                colorbar;
                axis([times([1,end])' x_prof([1,end])])
            end

        ch=3
            vm_profilo_ch = v_profiles{ch};
            ee_profilo_ch = e_profiles{ch};
            times = t_profiles{ch};
            x_prof=x_gates{ch};

            yM=max(max(vm_profilo_ch));
            ym=min(min(vm_profilo_ch));
            figure(30); clf; hold on; grid on
            subplot(2,1,1); hold on; grid on
            plot(it_arrivo(2,2).*[1 1], [ym yM], 'r')
            plot(it_arrivo(3,2).*[1 1], [ym yM], 'g')
            plot(times, vm_profilo_ch)
            subplot(2,1,2); hold on; grid on
            imagesc(times([1,end]), x_prof([1,end]), vm_profilo_ch')
            colorbar;
            axis([times([1,end])' x_prof([1,end])])

            if(~check_times)
                figure(31); clf; hold on;
                imagesc(times([1,end]), x_prof([1,end]), ee_profilo_ch')
                title(["max = " num2str(max(max(ee_profilo_ch)))])
                colorbar;
                axis([times([1,end])' x_prof([1,end])])
            end

%          if(~ismember(3, skip_ch))
%              pause()
%          end
    end



    puliturasi = true;
    if(puliturasi)
        if(calc_echo0)
            if(ch==1)
                all_ch=1;
            else
                all_ch=[2 3];
            end
            for ch=all_ch
                echo_data = e_profiles{ch};
                t_ch = t_profiles{ch};
                [empty_echo, empty_echo_std] = calc_echo_vuoto(echo_data, t_ch, min(t_max, it_arrivo(ch,2)*0.9));
                echo_vuoto{ch} = empty_echo;
                echo_vuoto_std{ch} = empty_echo_std;
                if(scrivi_echo0)
                    save([pathout 'empty_echo_ch' num2str(ch) '.txt'], 'empty_echo', '-ascii');
                    save(['empty_echo_std_ch' num2str(ch) '.txt'], 'empty_echo_std', '-ascii');
                end
            end
        else
            for ch=all_ch
                empty_echo = load(['empty_echo_ch' num2str(ch) '.txt']);
                echo_vuoto{ch} = empty_echo;
                empty_echo_std = load(['empty_echo_std_ch' num2str(ch) '.txt']);
                echo_vuoto_std{ch} = empty_echo_std;
            end
        end

        %  %   Pulitura dei dati
        for ch=all_ch
            vel_data = v_profiles{ch};
            echo_data = e_profiles{ch};
            empty_echo = echo_vuoto{ch};
            empty_echo_std = echo_vuoto_std{ch};
            xdata=x_gates{ch};
            t_ch = t_profiles{ch};
            %% Pulitura dei dati di velocità
            if (~isempty(find(check_ch==ch)))
                check_dop=true;
            else
                check_dop=false;
            end
            ndev=3
            %% vel_data e velocity_profiles sono uguali
            %% solo che velocity_profiles ha dei nan al posto degli zeri.

            disp('no clean')
%              [vel_data, velocity_profiles] = clean_data_dop(vel_data, echo_data, empty_echo, empty_echo_std, xdata, t_ch, ndev, check_dop);

            v_profiles{ch} = vel_data;
%              v_profiles_p{ch} = velocity_profiles;
        end
        disp('... Processing - pulitura dei dati - fatto')
    end



ch

    if(ch==1)
        %%%  ---------------------------------------------------
        %%%   ANALISI CANALE 1
        %%%  ---------------------------------------------------
        vm_profilo_ch = v_profiles{ch};
        x_prof=x_gates{ch};
        times = t_profiles{ch};

        if(isnan(tmax_ch1))
            tmax_ch1=input("tempo massimo per questa prova?\n")
        end

        for chx=[2 3]
            echo_data = e_profiles{chx};
            vm_profilo_ch = v_profiles{chx};
            times = t_profiles{chx};

            sel_times_v = find(times<tmax_ch1);
            echo_data = echo_data(sel_times_v,:);
            vm_profilo_ch = vm_profilo_ch(sel_times_v,:);
            times = times(sel_times_v);

            e_profiles{chx} = echo_data;
            v_profiles{chx} = vm_profilo_ch;
            t_profiles{chx} = times;
        end

        vm_profilo_ch = v_profiles{ch};
        x_prof=x_gates{ch};
        times = t_profiles{ch};

        sel_times_v = find(times<tmax_ch1);

        vm_profilo_ch = vm_profilo_ch(sel_times_v,:);
        times = times(sel_times_v);
        v_profiles{ch} = vm_profilo_ch;
        t_profiles{ch} = times;
%          %% echo
%          echo_data = e_profiles{ch};
%          echo_data = echo_data(sel_times_v,:);
%          e_profiles{ch} = echo_data;


        if(length(v_max_acc1)<2)
            v_max_acc1 = 2000.*[-1 1]
            migliora=true
            while(migliora)

                vel_data = vm_profilo_ch;
                vel_data(find(vel_data>v_max_acc1(2))) = nan;
                vel_data(find(vel_data<v_max_acc1(1))) = nan;

                figure(11); clf; hold on;
                plot(times, nanmean(vel_data,2), 'b', 'linewidth',2)
%                  plot(it_arrivo(ch,2).*[1 1], [-1000 1000], 'r', 'linewidth',3)
                plot(times, vel_data)
                grid on
                vmax = input('filtra per v_max = ');
                vmin = input('filtra per v_min = ');
                if(vmax == v_max_acc1(2) && vmin ==v_max_acc1(1))
%                      migliora = input("migliora ancora? ");
                    migliora = false;
                end
                v_max_acc1 = [vmin vmax];
            end
        else
                vel_data = vm_profilo_ch;
                vel_data(find(vel_data>v_max_acc1(2))) = nan;
                vel_data(find(vel_data<v_max_acc1(1))) = nan;

                figure(11); clf; hold on;
                plot(times, nanmean(vel_data,2), 'b', 'linewidth',2)
%                  plot(it_arrivo(ch,2).*[1 1], [-1000 1000], 'r', 'linewidth',3)
                plot(times, vel_data)
                grid on

        end

        t_dati = find(times>it_arrivo(ch,2));
        vel_data_parz=vel_data(t_dati,:);
        vel_data_parz(find(vel_data_parz==0))=nan;
        vel_data(t_dati,:)=vel_data_parz;

        plot(times, nanmean(vel_data,2), 'g', 'linewidth',2)

        if(~ismember(1,skip_ch))

            v_save = [[0 x_prof];  times/1000 vel_data/1000 ];
            save([path_out 'v_profiles_ch1.dat'], 'v_save', '-ascii');

            v_save = [times/1000 nanmean(vel_data')'/1000 ...
                                            nanstd(vel_data')'/1000 ];
            save([path_out 'v_m_std_ch1.dat'], 'v_save', '-ascii');
            path_out
            disp('scritto ch 1!')
            pause()
        end
        disp('Analisi canale 1 fatto')
        close(1)
        close(11)
    else
        %%%  ---------------------------------------------------
        %%%   ANALISI CANALI 2 e 3
        %%%  ---------------------------------------------------

        t_limits_ch23
        for chx=[2 3]
            echo_data = e_profiles{chx};
            vm_profilo_ch = v_profiles{chx};
            times = t_profiles{chx};

            sel_times_v = find(times<=t_limits_ch23(2));
            times = times(sel_times_v);
            echo_data = echo_data(sel_times_v,:);
            vm_profilo_ch = vm_profilo_ch(sel_times_v,:);

%              sel_times_v = find(times<t_limits_ch23(1));
%              echo_data(sel_times_v,:)=0;
%              vm_profilo_ch(sel_times_v,:)=0;
            sel_times_v = find(times>=t_limits_ch23(1));
            times = times(sel_times_v);
            echo_data = echo_data(sel_times_v,:);
            vm_profilo_ch = vm_profilo_ch(sel_times_v,:);

            e_profiles{chx} = echo_data;
            v_profiles{chx} = vm_profilo_ch;
            t_profiles{chx} = times;
        end
    end
end




%%%  ---------------------------------------------------
%% CALCOLO DELLA SUPERFICIE LIBERA
%%%  ---------------------------------------------------
if(~check_times)
    close(21)
    close(31)
end
ch=3
vm_profilo_ch = v_profiles{ch};
ee_profilo_ch = e_profiles{ch};
x_prof=x_gates{ch};
times = t_profiles{ch};

%%% cerca echo_min_pdist
if(isnan(echo_min_pdist))
    echo_min_pdist = 800
    migliora=true
    while(migliora)
        vel_data = vm_profilo_ch;
        vel_data(find(vel_data>v_max_acc1(2))) = nan;
        vel_data(find(vel_data<v_max_acc1(1))) = nan;

        x_prof_mass_u = costr_matrix_Emass(x_prof(find(x_prof> 20)), ee_profilo_ch(:, find(x_prof> 20)), echo_min_pdist);
        x_prof_mass_d = costr_matrix_Emass(x_prof(find(x_prof<=20)), ee_profilo_ch(:, find(x_prof<=20)), echo_min_pcorr);
        x_prof_mass = [x_prof_mass_d x_prof_mass_u];

%
        figure(34); clf; hold on;
        subplot(1,2,1)
        plot(times, x_prof_mass, '+')

%                  figure(33); clf; hold on;
        subplot(1,2,2); hold on
        imagesc(times([1,end]), x_prof([1,end]), ee_profilo_ch')
%                  h=contour(times, x_prof, ee_profilo_ch', echo_min_pdist.*[0.9 1.0 1.1], 'linewidth', 3);
        title(["max = " num2str(max(max(ee_profilo_ch)))])
        colorbar;
        axis([t_limits_ch23 x_prof([1,end])])
        echo_min = input('echo_min_pdist = ');
        if(echo_min_pdist == echo_min)
%                      migliora = input("migliora ancora? ");
            migliora = false;
        end
        echo_min_pdist = echo_min;
    end
end

x_prof_mass_u = costr_matrix_Emass(x_prof(find(x_prof> 20)), ee_profilo_ch(:, find(x_prof> 20)), echo_min_pdist);
x_prof_mass_d = costr_matrix_Emass(x_prof(find(x_prof<=20)), ee_profilo_ch(:, find(x_prof<=20)), echo_min_pcorr);
x_prof_mass = [x_prof_mass_d x_prof_mass_u];


[dist_min_u, dist_min] = ricava_h_rifless_E(x_prof, x_prof_mass, perc_size, perc_size_no);
h_min_mis(ch) = max(dist_min_u);
if(exist('forza_h_min_mis'))
    h_min_mis(ch) = forza_h_min_mis;
end

h_min_mis(ch)


%          figure(33); clf; hold on;
figure(34); clf; hold on;
subplot(1,2,1)
plot(times, x_prof_mass, '+')
for ii=1:length(dist_min)
    line(times([1,end]), dist_min(ii).*[1 1], 'color', 'r', 'linewidth', 1.5)
    line(times([1,end]), dist_min_u(ii).*[1 1], 'color', 'b', 'linewidth', 1)
end
%                  figure(33); clf; hold on;
subplot(1,2,2); hold on
imagesc(times([1,end]), x_prof([1,end]), ee_profilo_ch');
%          h=contour(times, x_prof, ee_profilo_ch', echo_min_pdist.*[0.9 1.0 1.1], 'linewidth', 3);
title(["max = " num2str(max(max(ee_profilo_ch)))])
colorbar;
axis([t_limits_ch23 x_prof([1,end])])

disp('check h_min_mis fig 34?')
%  pause()

%% Cerco la superficie libera
x_prof_mass_dist = x_prof_mass;
x_prof_mass_dist(:,x_prof<h_min_mis(ch)) = nan;


%%% pulizia estrema per dist
if(exist('gi'))
    "pulizia estrema per dist"
    sel_t = find(times>=ti & times<=tf );
    sel_g = find(x_prof>=gi & x_prof<=gf );
    x_prof_mass_dist(sel_t, sel_g) = nan;
end

%  "pulizia estrema"
%  sel_t = find(times>=ti & times<=tf );
%  y_lim_v = interp1([ti tf], [yi yf], times(sel_t));
%  parz = vm_profilo_ch(sel_t, :);
%  parz(find(parz>=y_lim_v)) = nan;
%  vm_profilo_ch(sel_t, :) = parz;





dist_dop_max = [];
%% tolgo il minimo e il massimo
if(~exist('g_dist'))
    g_dist=4;
end
for ii_prof=1:size(x_prof_mass_dist,1)
    nonan=find(~isnan(x_prof_mass_dist(ii_prof,:)));
    if(length(nonan)>=1)
%                      if(length(nonan)>3)
%                          nonan =nonan(2:end-1);
%                      end
        if(isempty(dist_dop_max))  %% primo dato
            if(min(diff(x_prof_mass_dist(ii_prof,nonan)))>10)
                continue
            else
                h_med = median(x_prof_mass_dist(ii_prof,nonan));
            end
        else
            check = find(nonan>=ii_min_prec-g_dist & nonan<=ii_max_prec+g_dist);
            if(~isempty(check))
                h_med = median(x_prof_mass_dist(ii_prof,nonan(check)));
            else
                ii_prova=nonan(find(abs(nonan-ii_min_prec)==min(abs(nonan-ii_min_prec))));
                h_med = x_prof_mass_dist(ii_prof,ii_prova);
%                          h_med = median(x_prof_mass_dist(ii_prof,nonan))
                if(abs(h_med-mean(dist_dop_max(end,2:end)))>10)
                    continue
                end
            end
        end

        ii_med=find(x_prof_mass_dist(ii_prof,nonan)>=h_med)(1);
        salti = find(diff(nonan)>3);
        if(~isempty(find(salti<ii_med)))
            ii_min = salti(find(salti<ii_med)(end))+1;
        elseif(~isempty(find(salti==ii_med)))
            ii_min = ii_med;
        else
            ii_min = 1;
        end
        if(~isempty(find(salti>=ii_med)))
            ii_max = salti(find(salti>=ii_med)(1));
        else
            ii_max = length(nonan);
        end
        dati_h = [x_prof_mass_dist(ii_prof,nonan(ii_min)) h_med x_prof_mass_dist(ii_prof,nonan(ii_max))];
        dist_dop_max = [dist_dop_max; [times(ii_prof) dati_h]];

        ii_min_prec=nonan(ii_min);
        ii_max_prec=nonan(ii_max);
    end
end
%%% togli zigzag in mezzo
%  "togli zigzag in mezzo"
%  dt=mean(diff(times));
%  dx=mean(diff(x_prof));
%  h_max=sort(dist_dop_max(:,4))(end-2);
%  ii_max=find(dist_dop_max(:,4)>=h_max-0.1*dx);
%  ii_max=[ii_max(1):ii_max(end)];
%  diff_ii=diff(dist_dop_max(ii_max,1));
%  diff_M =diff(dist_dop_max(ii_max,4));
%  ii_1=find(diff_ii>dt*0.9 & abs(diff_M)>5);
%  ii_1=find(diff_ii>dt*0.9 | abs(diff_M)>5);
%  ii_del=ii_max([ii_1(1)+1:ii_1(end)]);
%  dist_dop_max(ii_del,:)=[];




dist_dop_max_s = interp1(dist_dop_max(:,1), dist_dop_max(:,2:end), times);
x_prof_massimi{ch} = x_prof_mass;
dist_dop_max_dch{ch} = dist_dop_max_s;

ii_max = find(dist_dop_max_s(:,1)==max(dist_dop_max_s(:,1)))(1);
ii_ends = find(abs(dist_dop_max_s(:,1)-h_min_mis(ch))<=mean(diff(x_prof))*0.3);
ii_ends = ii_ends(find(ii_ends>ii_max));
if(~isempty(ii_ends))
    dist_dop_max_s(ii_ends(2):end,:)=nan;
end
if(any(isnan(dist_dop_max_s)))
    ii_centr = find(dist_dop_max_s(:,2)==max(dist_dop_max_s(:,2)))(1);
    for cc=1:3
        nans = find(isnan(dist_dop_max_s(1:ii_centr,cc)));
        dist_dop_max_s(nans,cc) = min(dist_dop_max_s(1:ii_centr,cc));
        nans = find(isnan(dist_dop_max_s(ii_centr+1:end,cc)));
        dist_dop_max_s(nans+ii_centr,cc) = min(dist_dop_max_s(ii_centr+1:end,cc));
    end
end

figure(34); clf; hold on;
subplot(1,2,1); hold on;
plot(times, x_prof_mass, '+')
for ii=1:length(dist_min)
    line(times([1,end]), dist_min(ii).*[1 1], 'color', 'r')
    line(times([1,end]), dist_min_u(ii).*[1 1], 'color', 'b')
end
plot(times, dist_dop_max_s, 'linewidth', 2)
axis(t_limits_ch23 )
subplot(1,2,2); hold on
imagesc(times([1,end]), x_prof([1,end]), ee_profilo_ch');
%          h=contour(times, x_prof, ee_profilo_ch', echo_min_pdist.*[0.9 1.0 1.1], 'linewidth', 3);
title(["max = " num2str(max(max(ee_profilo_ch)))])
colorbar;
axis([t_limits_ch23 x_prof([1,end])])

if(scrivi_altezze)
    x_prof_mass_save = [times./1000 dist_dop_max_s];
    save([path_out 'd_from_echo_ch' num2str(ch) '_mm.dat'], 'x_prof_mass_save', '-ascii');

    x_prof_mass_save = [times./1000 (dist_dop_max_s-dist_dop_rampa(3)).*cos(ang_assex_rampa)];
    save([path_out 'h_from_echo_ch' num2str(ch) '_mm.dat'], 'x_prof_mass_save', '-ascii');
end
disp('scritte altezze (da ch 3)? ')
pause()


%%%  ---------------------------------------------------
%% PER CALCOLO DELLA VELOCITA' - CANALE 2 e 3
%%%  ---------------------------------------------------

%% CANALE 3
%% Pulisci i dati del dop
%%%  ---------------------------------------------------

maschera = ones(size(vm_profilo_ch,1), length(x_prof));
maschera(repmat(x_prof, size(vm_profilo_ch,1), 1)>dist_dop_max_s(:,end)+5)=0;
maschera_d{ch}=maschera;

vm_profilo_ch(:,find(x_prof<dist_dop_rampa(ch))) = nan;
ee_profilo_ch(:,find(x_prof<dist_dop_rampa(ch))) = nan;
vm_profilo_ch(:,find(x_prof<h_min_mis(ch))) = nan;

vm_profilo_ch(~maschera) = nan;

if(isnan(v_max_acc3))
    v_max_acc3 = 1000*[-1 1]
    migliora=true
    while(migliora)
        vel_data = vm_profilo_ch;
        vel_data(find(vel_data>v_max_acc3(2))) = nan;
        vel_data(find(vel_data<v_max_acc3(1))) = nan;

        figure(30); clf;
        subplot(2,1,1); hold on;
        plot(times, nanmean(vel_data,2), 'b', 'linewidth',2)
%                  plot(it_arrivo(2,2).*[1 1], [-1000 1000], 'r', 'linewidth',3)
%                  plot(it_arrivo(3,2).*[1 1], [-1000 1000], 'g', 'linewidth',3)
        ii_2=find(x_prof>=dist_dop_rampa(ch));
        for ii_prof=ii_2:1:size(vm_profilo_ch,2)
            h_rif = (x_prof(ii_prof)-dist_dop_rampa(ch))*cos(ang_prampa_dop(ch));
            colori = col_prof(h_rif);
            plot(times, vel_data(:,ii_prof),'-+', 'color', colori, 'linewidth',1)
        end
        grid on
        subplot(2,1,2);
        imagesc(times([1,end]), x_prof([1,end]), vm_profilo_ch')
        colorbar;
        axis([times([1,end])' x_prof([1,end])])

        vmax = input('filtra per v_max = ');
        vmin = input('filtra per v_min = ');
        if(vmax == v_max_acc3(2) && vmin ==v_max_acc3(1))
%                      migliora = input("migliora ancora? ");
            migliora=false;
        end
        v_max_acc3 = [vmin vmax];
    end
    filtro_mM = find(vel_data<v_max_acc3(1) | vel_data>v_max_acc3(2));
    vel_data = vm_profilo_ch;
    vel_data(filtro_mM) = nan;
else
    filtro_mM = find(vm_profilo_ch<v_max_acc3(1) | vm_profilo_ch>v_max_acc3(2));
    vel_data = vm_profilo_ch;
    vel_data(filtro_mM) = nan;

    figure(30); clf;
    subplot(2,1,1); hold on;
    ii_2=find(x_prof>=dist_dop_rampa(ch));
    for ii_prof=ii_2:1:size(vm_profilo_ch,2)
        h_rif = (x_prof(ii_prof)-dist_dop_rampa(ch))*cos(ang_prampa_dop(ch));
        colori = col_prof(h_rif);
        plot(times, vel_data(:,ii_prof),'-+', 'color', colori, 'linewidth',1)
    end
    plot(times, nanmean(vel_data,2), 'b', 'linewidth',2)
%              plot(it_arrivo(2,2).*[1 1], [-1000 1000], 'r', 'linewidth',3)
%              plot(it_arrivo(3,2).*[1 1], [-1000 1000], 'g', 'linewidth',3)
    grid on
    subplot(2,1,2);
    imagesc(times([1,end]), x_prof([1,end]), vm_profilo_ch')
    colorbar;
    axis([times([1,end])' x_prof([1,end])])
end

filtro_zero = find(vm_profilo_ch==0);
%%%  filtro_e = find(ee_profilo_ch<echo_min);
filtro_e = [];

%          vm_profilo_ch(filtro_zero) = nan;
vm_profilo_ch(filtro_mM) = nan;

vm_profilo_ch(filtro_e) = nan;
ee_profilo_ch(filtro_e) = nan;
%
v_profiles{ch} = vm_profilo_ch;



%  %          figure(34); clf; hold on;
%  %          subplot(1,2,1); hold on;
%  %          plot(times, x_prof_mass, '+')
%  %          for ii=1:length(dist_corr)
%  %              line(times([1,end]), dist_min(ii).*[1 1], 'color', 'r')
%  %              line(times([1,end]), dist_min_u(ii).*[1 1], 'color', 'b')
%  %          end
%  %          plot(times, dist_dop_max_s, 'linewidth', 2)
%  %          plot(times, 10+dist_dop_max_s(:,end), '-m', 'linewidth', 3)
%  %          axis(t_limits_ch23 )
%  %          subplot(1,2,2); hold on
%  %          plot(times, vm_profilo_ch)
%  %          plot(it_arrivo(2,2).*[1 1], [-500 500], 'r')
%  %          plot(it_arrivo(3,2).*[1 1], [-500 500], 'g')

figure(30); clf; hold on;
subplot(2,1,1); hold on;
plot(times, nanmean(vm_profilo_ch,2), 'b', 'linewidth', 2)
plot(times, dist_dop_max_s*5, 'm', 'linewidth', 2)
ii_2=find(x_prof>=dist_dop_rampa(ch));
for ii_prof=ii_2:1:size(vm_profilo_ch,2)
    h_rif = (x_prof(ii_prof)-dist_dop_rampa(ch))*cos(ang_prampa_dop(ch));
    colori = col_prof(h_rif);
    plot(times, vm_profilo_ch(:,ii_prof),'-+', 'color', colori, 'linewidth',1)
end
%          plot(it_arrivo(2,2).*[1 1], [-1000 1000], 'r', 'linewidth',3)
%          plot(it_arrivo(3,2).*[1 1], [-1000 1000], 'g', 'linewidth',3)
grid on
axis(t_limits_ch23)
subplot(2,1,2); hold on; grid on
imagesc(times([1,end]), x_prof([1,end]), vm_profilo_ch')
colorbar;
axis([times([1,end])' x_prof([1,end])])

if(scrivi_vel)
    v_save = [[0 x_prof]; [times/1000 vm_profilo_ch/1000]];
    save([path_out 'v_profiles_ch' num2str(ch) '.dat'], 'v_save', '-ascii');
    e_save = [[0 x_prof]; [times/1000 ee_profilo_ch/1000]];
    save([path_out 'e_profiles_ch' num2str(ch) '.dat'], 'e_save', '-ascii');
end

%  disp('scritte velocità ch 3? ')
%  pause()


times_tot = sort([t_profiles{2}; t_profiles{3}]);
dist_dop_max_int = interp1(times, dist_dop_max_s(:,2), times_tot);

%%%  ---------------------------------------------------
%% PER CALCOLO DELLA VELOCITA' - CANALE 2 e 3
%%%  ---------------------------------------------------

%% CANALE 2 - pulisci i dati del dop
%%%  ---------------------------------------------------
ch=2

vm_profilo_ch = v_profiles{ch};
ee_profilo_ch = e_profiles{ch};
x_prof=x_gates{ch};
times = t_profiles{ch};

x_prof_mass = costr_matrix_Emass(x_prof, ee_profilo_ch, echo_min_pdop2);
dist_dop_max_t2 = interp1(times_tot, dist_dop_max_int, times);




%% figura dei profili di velocità
figure(20); clf();
subplot(2,1,1); hold on; grid on
ii_2=find(x_prof>=dist_dop_rampa(ch));
for ii_prof=ii_2:1:size(vm_profilo_ch,2)
    h_rif = (x_prof(ii_prof)-dist_dop_rampa(ch))*cos(ang_prampa_dop(ch));
    colori = col_prof(h_rif);
    plot(times, vm_profilo_ch(:,ii_prof),'-+', 'color', colori, 'linewidth',2)
end
%          plot(it_arrivo(2,2).*[1 1], [-500 500], 'r')
%          plot(it_arrivo(3,2).*[1 1], [-500 500], 'g')
grid on
axis(t_limits_ch23)
xlabel('tempo (ms)')
ylabel('V (mm/s)')
subplot(2,1,2); hold on; grid on
imagesc(times([1,end]), x_prof([1,end]), vm_profilo_ch')
colorbar;
axis([times([1,end])' x_prof([1,end])])

h_proiettata_d2 = (dist_dop_max_t2-dist_dop_rampa(3))/cos(ang_dop00_dop20)+dist_dop_rampa(2);

figure(24); clf; hold on;
subplot(1,2,1); hold on;
plot(times, x_prof_mass, '+')
plot(times, h_proiettata_d2, '-g', 'linewidth', 2)
plot(times, dist_dop_max_t2, '-b', 'linewidth', 2)
axis(t_limits_ch23 )

subplot(1,2,2); hold on
imagesc(times([1,end]), x_prof([1,end]), ee_profilo_ch');
%          h=contour(times, x_prof, ee_profilo_ch', echo_min_pdop2.*[0.9 1.0 1.1], 'linewidth', 3);
title(["max = " num2str(max(max(ee_profilo_ch)))])
colorbar;
axis([t_limits_ch23 x_prof([1,end])])

maschera = ones(size(vm_profilo_ch,1), length(x_prof));
maschera(repmat(x_prof, size(vm_profilo_ch,1), 1)>h_proiettata_d2+5)=0;
maschera_d{ch}=maschera;
vm_profilo_ch(~maschera) = nan;

[dist_min_u, dist_min] = ricava_h_rifless_E(x_prof, x_prof_mass, perc_size, perc_size_no);
h_min_mis(ch) = max(dist_min_u);


figure(24); clf; hold on;
subplot(1,2,1); hold on;
for ii=1:length(dist_min)
    line(times([1,end]), dist_min(ii).*[1 1], 'color', 'r')
    line(times([1,end]), dist_min_u(ii).*[1 1], 'color', 'b')
end
plot(times, x_prof_mass, '+')
plot(times, h_proiettata_d2, '-g', 'linewidth', 2)
plot(times, dist_dop_max_t2, '-b', 'linewidth', 2)
axis(t_limits_ch23 )
%%%
subplot(1,2,2); hold on
imagesc(times([1,end]), x_prof([1,end]), ee_profilo_ch');
%          h=contour(times, x_prof, ee_profilo_ch', echo_min_pdop2.*[0.9 1.0 1.1], 'linewidth', 3);
title(["max = " num2str(max(max(ee_profilo_ch)))])
colorbar;
axis([t_limits_ch23 x_prof([1,end])])


x_prof_mass(:,find(x_prof<=h_min_mis(ch))) = nan;
x_prof_massimi{ch} = x_prof_mass;
vm_profilo_ch(:,find(x_prof<dist_dop_rampa(ch))) = nan;
ee_profilo_ch(:,find(x_prof<dist_dop_rampa(ch))) = nan;
vm_profilo_ch(:,find(x_prof<h_min_mis(ch))) = nan;


%  %  "pulizia estrema"
%  %  sel_t = find(times>=ti & times<=tf );
%  %  y_lim_v = interp1([ti tf], [yi yf], times(sel_t));
%  %  parz = vm_profilo_ch(sel_t, :);
%  %  parz(find(parz>=y_lim_v)) = nan;
%  %  vm_profilo_ch(sel_t, :) = parz;



%% figura dei profili di velocità
figure(20); clf();
subplot(2,1,1); hold on; grid on
ii_2=find(x_prof>=dist_dop_rampa(ch));
for ii_prof=ii_2:1:size(vm_profilo_ch,2)
    h_rif = (x_prof(ii_prof)-dist_dop_rampa(ch))*cos(ang_prampa_dop(ch));
    colori = col_prof(h_rif);
    plot(times, vm_profilo_ch(:,ii_prof),'-+', 'color', colori, 'linewidth',2)
end
if(exist('y_lim_v'))
    plot(times(sel_t), y_lim_v, '-r')
end
%          plot(it_arrivo(2,2).*[1 1], [-500 500], 'r')
%          plot(it_arrivo(3,2).*[1 1], [-500 500], 'g')
grid on
axis(t_limits_ch23)
xlabel('tempo (ms)')
ylabel('V (mm/s)')
subplot(2,1,2); hold on; grid on
imagesc(times([1,end]), x_prof([1,end]), vm_profilo_ch')
colorbar;
axis([times([1,end])' x_prof([1,end])])

if(isnan(v_max_acc2))
    disp("Definisci v_max_acc2!!\n")
    pause()
end

filtro_mM = find(vm_profilo_ch<v_max_acc2(1) | vm_profilo_ch>v_max_acc2(2));
filtro_zero = find(vm_profilo_ch==0);
%%%  filtro_e = find(ee_profilo_ch<echo_min);
filtro_e = [];

vm_profilo_ch(filtro_zero) = nan;
vm_profilo_ch(filtro_mM) = nan;

vm_profilo_ch(filtro_e) = nan;
ee_profilo_ch(filtro_e) = nan;

v_profiles{ch} = vm_profilo_ch;

figure(20); clf();
subplot(2,1,1); hold on; grid on
plot(times, nanmean(vm_profilo_ch,2), 'b', 'linewidth', 2)
plot(times, dist_dop_max_t2*10, 'm', 'linewidth', 2)
ii_2=find(x_prof>=dist_dop_rampa(ch));
for ii_prof=ii_2:1:size(vm_profilo_ch,2)
    h_rif = (x_prof(ii_prof)-dist_dop_rampa(ch))*cos(ang_prampa_dop(ch));
    colori = col_prof(h_rif);
    plot(times, vm_profilo_ch(:,ii_prof),'-+', 'color', colori, 'linewidth',1)
end
%          plot(it_arrivo(2,2).*[1 1], [-500 500], 'r')
%          plot(it_arrivo(3,2).*[1 1], [-500 500], 'g')
grid on
axis(t_limits_ch23)
xlabel('tempo (ms)')
ylabel('V (mm/s)')
subplot(2,1,2); hold on; grid on
imagesc(times([1,end]), x_prof([1,end]), vm_profilo_ch')
colorbar;
axis([times([1,end])' x_prof([1,end])])


if(scrivi_vel)
    v_save = [[0 x_prof]; [times/1000 vm_profilo_ch/1000]];
    save([path_out 'v_profiles_ch' num2str(ch) '.dat'], 'v_save', '-ascii');
    e_save = [[0 x_prof]; [times/1000 ee_profilo_ch/1000]];
    save([path_out 'e_profiles_ch' num2str(ch) '.dat'], 'e_save', '-ascii');
end


%  disp('scritte velocità ch 2? ')
%  pause()




%%%  ---------------------------------------------------
%%% COMPONI LE VELOCITÀ
%%%  ---------------------------------------------------

%      %% cambio di base
%      v_x = v_d00*cos(angolo_assex_dop00) + v_d20*cos(angolo_assex_dop20);
%      v_y = v_d00*sin(angolo_assex_dop00) + v_d20*sin(angolo_assex_dop20);

%      %% cambio di base
%      v_x = v_d00*cos(angolo_rampa_dop00) + v_d20*cos(ang_rampa_dop20);
%      v_y = v_d00*sin(angolo_rampa_dop00) + v_d20*sin(ang_rampa_dop20);

tch_20_rampa = t_profiles{2};
tch_00_rampa = t_profiles{3};

x_20 = (x_gates{2}-dist_dop_rampa(2)).*cos(ang_prampa_dop(2));
x_00 = (x_gates{3}-dist_dop_rampa(3)).*cos(ang_prampa_dop(3));

v_20_rampa = v_profiles{2};
v_00_rampa = v_profiles{3};

v_20_rampa(find(x_20<h_min_mis(2)))=nan;
v_00_rampa(find(x_00<h_min_mis(3)))=nan;

%%% sincronizza i tempi e aggiusta le lunghezze
ii = find(min(tch_20_rampa(2).-tch_00_rampa(1:3)));
if(ii==2)
    ii00=1;
    ii20=1;
elseif(ii==1)
    ii20=2;
    ii00=1;
elseif(ii==3)
    ii20=1;
    ii00=2;
end
numdati = min(length(tch_00_rampa)-ii00+1,length(tch_20_rampa)-ii20+1);
tch_00_rampa = tch_00_rampa(ii00:numdati+ii00-1,:);
tch_20_rampa = tch_20_rampa(ii20:numdati+ii20-1,:);
v_00_rampa = v_00_rampa(ii00:numdati+ii00-1,:);
v_20_rampa = v_20_rampa(ii20:numdati+ii20-1,:);


tempi_vett = mean([tch_00_rampa, tch_20_rampa],2);
%  %  tempi_vett = sort([tch_00_rampa; tch_20_rampa]);
v_00_rampa = interp1(tch_00_rampa, v_00_rampa, tempi_vett);
v_20_rampa = interp1(tch_20_rampa, v_20_rampa, tempi_vett);

%      %      sel = find(isnan(v_00_rampa(find(maschera))));
%      %      v_00_rampa(find(maschera(sel))) = 0.0;


% h_prof_darampa = [ceil(dist_dop_rampa(ch)/dh)*dh-dh/2:dh:ceil(max(dist_dop_max)/dh)*dh-dh/2]
h_prof_darampa = [dh/2:dh:ceil(max(dist_dop_max(:,2:end))/dh)*dh-dh/2];

v_00_rett = [];
v_20_rett = [];

for h=h_prof_darampa
    sel00 = find(x_00>=h-dh/2 & x_00<h+dh/2);
    parz = v_00_rampa(:,sel00);
    vmean00 = nanmean(parz,2);
    v_00_rett = [v_00_rett vmean00];

    sel20 = find(x_20>=h-dh/2 & x_20<h+dh/2);
    parz = v_20_rampa(:,sel20);
    vmean20 = nanmean(parz,2);
    v_20_rett = [v_20_rett vmean20];
end

[v_dir_rampa, v_perp_rampa] = det_vett_vel(v_00_rett, v_20_rett);

vm_dir = nanmean(v_dir_rampa,2);
vm_perp = nanmean(v_perp_rampa,2);

v_std_dir = nanstd(v_dir_rampa')';
v_std_perp = nanstd(v_perp_rampa')';



%  %              figure(1); clf();
%  %              for ii=sel_times'
%  %                  for ii_prof=1:length(h_prof_darampa)
%  %                      plot([0 v_dir_rampa(ii,ii_prof)],h_prof_darampa(ii_prof).+[0  v_perp_rampa(ii,ii_prof)])
%  %                      hold on
%  %                  end
%  %                  axis([-1000 1000 -500 500])
%  %      %              pause()
%  %                  hold off
%  %                  grid on
%  %              end


n_prof = length(h_prof_darampa);
colori = rainbow(n_prof);

%      h_leg=cellstr(num2str([0.0; h_prof_darampa']));
%      h_leg{1} = 'm';
h_leg=cellstr(num2str([ h_prof_darampa'; 0.0;]));
h_leg{end} = 'm';

figure(5); clf()
subplot(211); hold on
%          plot(tempi_vett/1000, v_dir_rampa/1000, '-o', 'linewidth',1)
plot(tempi_vett/1000, vm_dir/1000, '-1','linewidth',4)
plot(tempi_vett/1000, (vm_dir.+v_std_dir)./1000, '-3','linewidth',1)
plot(tempi_vett/1000, (vm_dir.-v_std_dir)./1000, '-3','linewidth',1)
grid on
ylabel('v dir rampa (m/s)')


subplot(212); hold on
%          plot(tempi_vett/1000, v_perp_rampa/1000, ['-o' ';' num2str(h_prof_darampa) ';'],'linewidth',1)
plot(tempi_vett/1000, vm_perp/1000, '-+2','linewidth',4)
plot(tempi_vett/1000, (vm_perp.+v_std_perp)/1000, '-3','linewidth',1)
plot(tempi_vett/1000, (vm_perp.-v_std_perp)/1000, '-3','linewidth',1)
xlabel('tempo (s)')
ylabel('v perp rampa (m/s)')
legend(h_leg)
grid on
%          print([path_out 'vel.png'],"-dpng")
%          pause()



%          figure(6); clf(); hold on
%          plot(tempi_vett/1000, v_std_dir./1000, '-1','linewidth',2)
%          plot(tempi_vett/1000, abs(vm_perp)./1000, '-3','linewidth',2)
%          title('errore (std) v dir e |vm perp|')
%          grid on
%
%          figure(7); clf;
%          ii=find(tempi_vett<4000);
%          plot(abs(vm_perp(ii)), v_std_dir(ii), 'x')
%          title('correlazione vm_perp v_std_dir')
%          grid on

%          figure(8); clf; hold on;
%          ii=find(tempi_vett<4000);
%          plot(tempi_vett/1000, v_std_dir./1000,  '-1','linewidth',2)
%          plot(tempi_vett/1000, v_std_perp./1000, '-3','linewidth',2)
%          title('errore (std) vm perp v dir')
%          grid on




v_save = [[0 h_prof_darampa  h_prof_darampa]; ...
        [tempi_vett/1000 v_dir_rampa/1000 v_perp_rampa/1000]];
save([path_out 'v_profiles_dp_rampa_dh' num2str(dh) '.dat'], 'v_save', '-ascii');

v_save = [tempi_vett/1000 vm_dir/1000   vm_perp/1000 ...
                    v_std_dir/1000   v_std_perp/1000 ];
save([path_out 'v_m_std_dp_rampa_dh' num2str(dh) '.dat'], 'v_save', '-ascii');






analisi_spettrale=false;
if(analisi_spettrale)
    %  %  %  Analisi spettrale fft
    for jj_dop=dop_vtest
        pkg load signal
        tempi = t_profiles{jj_dop};
        dt = mean(diff(tempi));
        vel_data = v_profiles{jj_dop};
        i_end = find(tempi<maxtempo)(end);
        vel_data_z = vel_data(it_arrivo(jj_dop, 1):i_end, :);
        vel_data_z(find(isnan(vel_data(it_arrivo(jj_dop, 1):i_end,:)))) = 0.0;
        spectra_aib = [];
        spectra_welsh_aib = [];
        for ii=1:size(vel_data_z,2)
            spectra_aib = [spectra_aib fft(vel_data_z(:,ii))];
            [spectra, freq] = pwelch(vel_data_z(:,ii), 128, 0.3);
            spectra_welsh_aib = [spectra_welsh_aib spectra  ];
        end
        per = dt./freq;
        freq = freq/dt;
        N=size(spectra_aib,1);
        freqfft = [1:N]./N./dt./2;
        perfft = 2.*dt.*N./[1:N];
    end
    disp('... Analisi spettrale fatta')
end












function [v_profiles, e_profiles, t_profiles, x_gates, ndati, all_ch] = read_files_dop_exp(path_in, filess, n_ch_max)
#     funz = analizza_vel_dop_matlab;
    t_ms=[];
    ndati=zeros(n_ch_max);


    t_profiles = {};
    v_profiles = {};
    e_profiles = {};
    for ch=1:n_ch_max
        v_profiles{ch} = [];
        e_profiles{ch} = [];
        t_profiles{ch} = [];
    end

    t_in=0.0;
    for ff=1:size(filess,1)
        filename = deblank(filess(ff,:));
        if(~isempty(strfind(filename, '.ADD')))
            id_set = 5;
            [DOP_data, all_ch, x_gates, Nb_gates, Nb_profiles, Nb_sequence] = read_dop(filename);
            nchannels = length(all_ch);

            all_ch

            for ch=all_ch
                data=DOP_data{ch};
                xdata=x_gates{ch};

#                 [time_ms, vel_profiles, echo_profiles, x_profiles] = funz.extract_data_dop(data, xdata);
#                 ch_g = funz.det_dop(x_profiles);
                [time_ms, vel_profiles, echo_profiles, x_profiles] = extract_data_dop(data, xdata);
                ch_g = det_dop(x_profiles);
                if(id_set==2 && ch_g==3)
                    all_ch = setdiff(all_ch, [ch]);
                    continue
                end
                x_gates{ch} = x_profiles;

                v_profiles{ch} = [v_profiles{ch}; vel_profiles ];
                e_profiles{ch} = [e_profiles{ch}; echo_profiles ];
                t_profiles{ch} = [t_profiles{ch}; time_ms+t_in];

                %% da trovare un dt più affidabile..
                dt(ch)=mean(diff(time_ms));

                t_ms = [t_ms; time_ms+t_in];
                ndati(ch) = ndati(ch) + size(vel_profiles,1);
            end
            t_in = max(t_ms)+min(dt);

        else
            [vel,echo,z,t]=BinDop(filename);
            v_profiles{1} = vel;
            e_profiles{1} = echo;
            t_profiles{1} = t;
            x_gates{1} = z;
            nchannels =1;
            ndati(1)=size(vel,1);
        end
    end
    ndati=ndati(1:max(all_ch));
end




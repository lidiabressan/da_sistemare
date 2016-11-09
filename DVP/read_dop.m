


function [data, channels, d_gates, Nb_gates, Nb_profiles, Nb_sequence] = read_dop(filename)

    testo=fileread(filename);
    lines=strsplit(testo, '\n');
    nrighe=max(size(lines));
    data={};
    n_ch_max=3;
    for ii=1:n_ch_max
        data{ii}=[];
    end

    gates=false;
    leggi_ch=false;

    c_lines=zeros(1, 6+n_ch_max);
    conta_ch=zeros(1, n_ch_max);
    nchannels=0;
    channels=[];

    for ii=1:nrighe
         if(findstr(lines{ii}, 'Channel'))
            ch_str_a=strsplit(lines{ii}, 'Channel');
            ch_str=ch_str_a{2};


            ch=str2num(ch_str);
            channels = union(channels, [ch]);
            nchannels=max(nchannels, ch);
            c_lines(6+ch)=c_lines(6+ch)+1;
        elseif(findstr(lines{ii}, 'Depth of gates [mm]'))
            gates=true;
            leggi_ch=false;
            c_lines(1)=c_lines(1)+1;
        elseif(findstr(lines{ii}, 'Data values'))
            leggi_ch=true;
            gates=false;
            conta_ch(ch)=0;
            c_lines(2)=c_lines(2)+1;
        elseif(any(isletter((lines{ii}))))
            c_lines(3)=c_lines(3)+1;
            continue
        elseif(length(deblank(lines{ii}))==0)
            c_lines(4)=c_lines(4)+1;
            continue
        elseif(gates)
            c_lines(5)=c_lines(5)+1;
            d_gates{ch}=str2num(lines{ii});
        elseif(leggi_ch)
            data{ch}=[data{ch}; str2num(lines{ii})];
            conta_ch(ch)=conta_ch(ch)+1;
            c_lines(6)=c_lines(6)+1;
        end
    end
    if (sum(c_lines)~=nrighe)
        [sum(c_lines) nrighe]
        erroreee
    end

    for ch=channels
        Nb_gates(ch)=size(d_gates{ch},2);
    end
    Nb_profiles = conta_ch;
    Nb_sequence = c_lines(7);

    c_lines=c_lines(1:6+nchannels);
    conta_ch=conta_ch(1:nchannels);

end



function id_setting = det_setting(pathfile)
    if(strfind(pathfile, '20150122'))
        id_setting=1;
    elseif(length(strfind(pathfile, '20150127'))>0 || ...
           length(strfind(pathfile, '20150128'))>0 || ...
           length(strfind(pathfile, '20150209'))>0 || ...
           length(strfind(pathfile, '20150211'))>0 )
        id_setting=2;
    elseif(length(strfind(pathfile, '201502'))>0  || ...
           length(strfind(pathfile, '2015030'))>0 || ...
           length(strfind(pathfile, '2015031'))>0   )
        id_setting=3;
    elseif(length(strfind(pathfile, '201503'))>0 || ...
           length(strfind(pathfile, '20150401'))>0 )
        id_setting=4;
    else
        id_setting=5;
    end
end


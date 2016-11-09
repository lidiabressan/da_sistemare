%  %  febbraio 2015
%  %  Lidia Bressan
%  %
%  %  It looks for files in path_in

function [filess, nfiles] = find_files(path_in, stringaDaCercare, opzioni1='', opzioni2='')
    if(isunix)
        [stato_bash,filesls]=system(['find ' opzioni1 path_in opzioni2 '  -name  "' stringaDaCercare '"' ]);
    elseif(ispc)
        [stato_dos, pwd]=system(['echo %cd%' ]);
        stato_dos
        [stato_dos, cd]=system(['cd ' path_in ]);
        stato_dos
        [stato_dos, filesls]=system(['dir ' path_in  stringaDaCercare ' /s' ]);
        stato_dos
        [stato_dos, cd]=system(['cd ' pwd ]);
        stato_dos
    end
    if(stato_bash==0)
        clear aa; aa=strfind(filesls,"\n");
        if(~isempty(aa))
            filess=filesls(1:aa(1)-1);
            for ii=2:length(aa)
                filess=[filess; filesls(aa(ii-1)+1:aa(ii)-1)];
            end
        else
            filess = deblank(filesls);
        end

        nfiles=size(filess,1);
    elseif(isunix)
        stato_bash
        filesls
        nfiles=0;
        filess=[];
    else
        stato_dos
        filess=[];
        nfiles=0;
    end
end

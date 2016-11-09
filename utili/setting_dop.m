
## acqua (e gel)
c1=1500;
rho1=1000;
Z1 =c1*rho1;




function angolo_critico = det_ang_critico(c1, c2)
    angolo_critico = asin(c1/c2);
end




function ang_rifr = angolo_rifrazione_rad (ang_inc, c1, c2)
    ang_rifr = asin(c2*sin(ang_inc)/c1);
    if(iscomplex(ang_rifr))
        ang_rifr=0.0;
    end
end

function D = coeff_rifrazione(c1, c2, rho1, rho2, ang_inc)
    ang_rifr = angolo_rifrazione_rad (ang_inc, c1, c2);
    if(ang_rifr==0)
        D=0.0;
    else
        z1=c1*rho1;
        z2=c2*rho2;
#         D1=4*z1*z2*cos(ang_inc)*cos(ang_rifr)/((z2*cos(ang_inc)+z1*cos(ang_rifr))^2);
        D1=4*z1*z2*cos(ang_inc)*cos(ang_inc)/((z2*cos(ang_inc)+z1*cos(ang_rifr))^2);
# #         D2=((z2-z1)/(z2+z1))^2;
        D=D1;
    end
end



function tab, tab_show = fai_tab_rifraz(angolo_incidenza_daT, c1, c2, rho1, rho2)
    tab = [];
    for ang_grad=angolo_incidenza_daT
        ang_rad=ang_grad*2*pi/360.0;
        in_plex = [ang_rad angolo_rifrazione_rad(ang_rad, c1, c2)    coeff_rifrazione(c1, c2, rho1, rho2, ang_rad)];
        in_canale = [in_plex(2)  angolo_rifrazione_rad(in_plex(2), c2, c1)   coeff_rifrazione(c2, c1, rho2, rho1, in_plex(2))    coeff_rifrazione(c2, c1, rho2, rho1, in_plex(2))*in_plex(3)   ];
        tab = [tab; [in_plex in_canale]];
    end

    tab =[ tab  sin(tab(:,1)).*100 tab(:,end).**2 ];
    tab_show=tab;
    tab_show(:, [1 2 4 5])=tab(:, [1 2 4 5]).*360/2.0/pi;
#     disp('     AI-pl     AR-pl      CR-pl       AI-ch      AR-ch      CR-ch   CR-ch*CR-pl    %vel     +++  CR-lettura = (CR-ch*CR-pl)^2 ');
end




# # z_gel = 1.58, 1.68
## Plexiglass
disp('Plexiglass')
c2 = 2750;
rho2 = 1185;
# z2=3.26;
angolo_incidenza_daT = [5:2.5:30];



angolo_critico = det_ang_critico(c1, c2)*360.0/2/pi

[tab, tab_show] = fai_tab_rifraz(angolo_incidenza_daT, c1, c2, rho1, rho2);


disp('     AI-pl     AR-pl      CR-pl       AI-ch      AR-ch      CR-ch   CR-ch*CR-pl    %vel     +++  CR-lettura = (CR-ch*CR-pl)^2 ')
tab_show





# # z_gel = 1.58, 1.68
## Plexiglass
disp('Glass')
c2 = 5900;
rho2 = 2200;
# z2=3.26;
angolo_incidenza_daT = [1:1:15];

angolo_critico = det_ang_critico(c1, c2)*360.0/2/pi

[tab, tab_show] = fai_tab_rifraz(angolo_incidenza_daT, c1, c2, rho1, rho2);

disp('     AI-pl     AR-pl      CR-pl       AI-ch      AR-ch      CR-ch   CR-ch*CR-pl    %vel     +++  CR-lettura = (CR-ch*CR-pl)^2 ')
tab_show


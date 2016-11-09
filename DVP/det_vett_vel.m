
function [v_tot_x, v_tot_y] = det_vett_vel(v00, v20)

    ang = 90+20;
    v20_y = v20.*sin(2*pi*ang/360);
    v20_x = v20.*cos(2*pi*ang/360);
    ca = v20_y./v20_x;
    ca_perp = -1./ca;
    v_tot_x = (v00-v20_y+ca_perp.*v20_x)./ca_perp;
    v_tot_y = v00;

%      figure(1); clf()
%      hold on
%      plot([0 0], [0 v00],'-+k', 'linewidth',2)
%      plot([0 v20_x], [0 v20_y],'-+k', 'linewidth',2)
%      plot([0 v_tot_x], [0 v_tot_y],'-+2', 'linewidth',2)
%      axis('equal')
%      pause()

end

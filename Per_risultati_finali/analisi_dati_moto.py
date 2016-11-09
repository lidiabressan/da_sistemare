
import numpy as np
import matplotlib.pyplot as plt

dati_g = ['all', 'Vuoto', 'b8c', 'b4c', 'b7c', 'b4_', ]
#dati_g = ['all', 'Vuoto', 'b8c', 'b4c',  ]
#dati_g = ['all', 'Vuoto',  'b7c', 'b4_', ]
#dati_g = ['all', 'Vuoto', 'b4_', 'b7c']
#dati_g = [ 'b4c', 'b8c']
##dati_g = [ 'b4_']
#dati_g = [ 'all']



CL = 0.178
Cm = 10**-1

path_parz='../EspBoulder/Ris_now'

dv = 0.01




"""
  Dati esperimento
"""

coeff_drag=dict()



pesi=dict()
pesi['b4_'] =  53.0
pesi['b7c'] =  61.0
pesi['b4c'] =  97.5
pesi['b8c'] = 150.0
boulders = pesi.keys()


coeff_drag=dict.fromkeys(boulders)
coeff_drag_3d=dict.fromkeys(boulders)
for bbx in coeff_drag.keys():
    coeff_drag[bbx] = dict.fromkeys(('0','45','90'))
    coeff_drag_3d[bbx] = dict.fromkeys(('0','45','90'))

    if bbx=='b4_' or bbx=='b7c':
        ## cubi
        coeff_drag_3d[bbx]['0'] = 1.05
        coeff_drag_3d[bbx]['45'] = 0.8
        coeff_drag_3d[bbx]['90'] = coeff_drag_3d[bbx]['0']
        coeff_drag[bbx]['0'] = 2.05
        coeff_drag[bbx]['45'] = 1.55
        coeff_drag[bbx]['90'] = coeff_drag[bbx]['0']
    else:
        ### prisma rettangolare
        coeff_drag[bbx]['0'] =  1.7  ## flusso contro lato lungo
        coeff_drag[bbx]['90'] = 2.5  ## flusso contro lato corto
        coeff_drag[bbx]['45'] = 1.7
        coeff_drag_3d[bbx]['0'] = 1.95
        coeff_drag_3d[bbx]['45'] = 1.95
        coeff_drag_3d[bbx]['90'] = coeff_drag_3d[bbx]['0']







col_pesi=dict()
col_pesi['b4_'] = 0.8   ## #   53.0  ## bianco
col_pesi['b7c'] = 0.6   ## #   61.0
col_pesi['b4c'] = 0.6   ## #   97.5
col_pesi['b8c'] = 0.3   ## #  150.0  ## nero


k_attrito=dict()
k_attrito['b4_'] = 0.52 ##  0.52+-0.05
k_attrito['b7c'] = 0.51 ##  0.51+-0.06
k_attrito['b4c'] = 0.65 ##  0.65+-0.06
k_attrito['b8c'] = 0.56 ##  0.56+-0.08


dim_boulders=dict()

dim_boulders['b4_'] = np.array([ 3.0,  3.1,  2.95])
dim_boulders['b7c'] = np.array([ 3.0,  2.8,  3.0 ])
dim_boulders['b4c'] = np.array([ 3.1,  5.45, 3.0 ])
dim_boulders['b8c'] = np.array([ 3.14, 5.92, 3.08])

###
teta=np.arcsin(1./10)
d_cm=2.0
###
H_cm=np.arange(0.0, 5.0, 0.2)  #*10**-2
##
H = H_cm*10**-2             ## cm --> m
d = d_cm*10**-2
##
g=9.81




print H.shape, H_cm.shape

from calc_v_teor import calc_k
from calc_v_teor import v_soglia_nand_mod, v_soglia_nand_mod4


k = calc_k(H, teta, d=0.02)




def check_nans(h_moto_m, v_moto_m):
    #print(" ***************** check nans")
    for ii in np.arange(0, h_moto_m.shape[0]):
        ##if not np.any(np.isnan(h_moto_m[ii])) == np.all(np.isnan(h_moto_m[ii])) :
            ##print ii, np.any(np.isnan(h_moto_m[ii])), np.all(np.isnan(h_moto_m[ii]))
            ##raise
        if not np.any(np.isnan(v_moto_m[ii,:])) == np.all(np.isnan(v_moto_m[ii,:])) :
            print ii, np.any(np.isnan(v_moto_m[ii,:])), np.all(np.isnan(v_moto_m[ii,:]))
            raise
        if not np.any(np.isnan(v_moto_m[ii,:])) == np.any(np.isnan(h_moto_m[ii])):
            print ii, test_list[ii], v_moto_m[ii,:], h_moto_m[ii]
            raise
    return





fig12=None
fig12=plt.figure(12)
asse12 = plt.subplot(1,1,1)


for ii, bbx in enumerate(dati_g):
    if bbx in boulders:
        ii_primo_b = ii
        break

## ## Propagazione e figure
conta=0
for bbx in dati_g:
    print bbx
    conta=conta+1


    H = np.loadtxt('{1}/dati_H_{0}.txt'.format(bbx, path_parz))

    k = calc_k(H, teta, d=0.02)
    H_cm = H*10


    t_percorrenza = np.loadtxt('{1}/dati_t_percorrenza_{0}.txt'.format(bbx, path_parz))
    v_prop_media = np.loadtxt('{1}/dati_v_media_{0}.txt'.format(bbx, path_parz))
    t_h_bore = np.loadtxt('{1}/dati_t_h_bore_{0}.txt'.format(bbx, path_parz))
    ## v massima e v massima media
    vmax_bore = np.loadtxt('{1}/dati_vmax_bore_{0}.txt'.format(bbx, path_parz))
    ## dati moto
    moto = np.loadtxt('{1}/dati_moto_{0}.txt'.format(bbx, path_parz))
    v_moto = np.loadtxt('{1}/dati_v_moto_{0}.txt'.format(bbx, path_parz))
    err_v_moto = np.loadtxt('{1}/dati_err_v_moto_{0}.txt'.format(bbx, path_parz))
    ## gia' corretti
    h_moto = np.loadtxt('{1}/dati_h_moto_{0}.txt'.format(bbx, path_parz))

    fid = open('{1}/dati_bb_oo_{0}.txt'.format(bbx, path_parz), 'r')
    bb_oo = np.array(filter(None, np.array(fid.read().split('\n'))))
    fid.close()

    blocco = np.array([xx.split('_')[0] for xx in bb_oo])
    orient = np.array([xx.split('_')[1] if 'Vuoto' not in xx else '-' for xx in bb_oo ])
    set_oo = list(set(orient))
    set_oo = sorted(set_oo)
    print set_oo

    fid = open('{1}/dati_test_{0}.txt'.format(bbx, path_parz), 'r')
    test_list = np.array(fid.read().split('\n'))
    fid.close()

    motosi = moto==1
    motono = moto==0
    motobo = np.isnan(moto)


    for oo in set_oo:
        sel_oo = orient==oo
        if np.any(sel_oo):
            print bbx, oo, np.count_nonzero(sel_oo), np.count_nonzero(np.logical_not(np.isnan(v_moto[sel_oo,0]))), np.count_nonzero(motosi[sel_oo])
        else:
            print bbx, oo, np.count_nonzero(sel_oo), np.count_nonzero(np.logical_not(np.isnan(v_moto[sel_oo,0]))), np.count_nonzero(motosi[sel_oo])





    if False:
        continue




    if 'all' in bbx and False:
        print 'Figure sui dati di propagazione'
        """
        Figura sui dati di propagazione:
         * sensori:
             - WL1, WL2, WL3
             - DOP1
         * dati di propagazione
             - velocita' di propagazione WL1-WL2, WL2-WL3, WL1-WL3
             - Tempo tra arrivo onda - arrivo onda di ritorno a WL3
             - Altezza bore media tra arrivo onda - arrivo onda di ritorno a WL3
             - velocita' massima DOP1
        """
        nonan = np.logical_and(np.logical_not(np.isnan(v_prop_media[:,1])), \
                               np.logical_not(np.isnan(vmax_bore[:,1])))
        #p_v = np.polyfit (v_prop_media[nonan,1], vmax_bore[nonan,0], 1)
        ##print p_v


        plt.figure(1)
        plt.subplot(311)
        #plt.plot(v_prop_media[motobo,1], t_h_bore[motobo,0], '+b', linewidth=1)
        #plt.plot(v_prop_media[motono,1], t_h_bore[motono,0], 'xg', linewidth=2)
        #plt.plot(v_prop_media[motosi,1], t_h_bore[motosi,0], 'or', linewidth=1)
        plt.plot(v_prop_media[:,1], t_h_bore[:,0], 'or', linewidth=1)
        plt.ylabel('t bore (m/s)')
        plt.grid(True)

        #dv = 0.15
        #v_prop = np.arange(1.2, 2.2+dv, dv)
        #dh_bore = 0.15
        #h_bore = np.arange(1.2, 2.0+dh_bore, dh_bore)
        #plt.vlines(v_prop, np.min(h_bore), np.max(h_bore), alpha=0.5)
        #plt.hlines(h_bore, np.min(v_prop), np.max(v_prop), alpha=0.5)
        plt.subplot(312)
        #plt.plot(v_prop_media[motobo,1], t_h_bore[motobo,1], '+b', linewidth=1)
        #plt.plot(v_prop_media[motono,1], t_h_bore[motono,1], 'xg', linewidth=2)
        #plt.plot(v_prop_media[motosi,1], t_h_bore[motosi,1], 'or', linewidth=1)
        plt.plot(v_prop_media[:,1], t_h_bore[:,1], 'or', linewidth=1)
        plt.ylabel('h bore (m/s)')
        plt.grid(True)
        plt.grid(True)

        plt.subplot(313)
        #plt.plot(v_prop_media[motobo,1], vmax_bore[motobo,0], '+b', linewidth=1)
        #plt.plot(v_prop_media[motono,1], vmax_bore[motono,0], 'xg', linewidth=2)
        #plt.plot(v_prop_media[motosi,1], vmax_bore[motosi,0], 'or', linewidth=1)
        plt.plot(v_prop_media[:,1], vmax_bore[:,0], 'or', linewidth=1)
        plt.xlabel('v propagazione (m/s)')
        plt.ylabel('v max bore (m/s)')
        plt.grid(True)
        plt.savefig('propagazione.png', dpi=300, format='png')


        ### Figure articolo

        plt.figure(7)
        plt.clf()
        plt.plot(H, v_prop_media[:,1], 'o', label='c')
        plt.ylabel('Velocity (m/s)')
        plt.xlabel('Basin water level (cm)')
        plt.grid(True)
        plt.savefig('H_c.png', dpi=300, format='png')

        plt.figure(7)
        plt.clf()
        plt.plot(H, vmax_bore[:,0], 'o', label="v_M")
        plt.ylabel('Velocity (m/s)')
        plt.xlabel('Basin water level (cm)')
        plt.grid(True)
        plt.savefig('H_vmax.png', dpi=300, format='png')

        plt.plot(H, v_prop_media[:,1], 'o', label='c')
        plt.legend(loc="lower right")
        plt.savefig('H_vmax_c.png', dpi=300, format='png')

        plt.figure(7)
        plt.clf()
        plt.plot(v_prop_media[:,1], t_h_bore[:,0], 'o')
        plt.ylabel('time (s)')
        plt.xlabel('Velocity c (m/s)')
        plt.grid(True)
        plt.savefig('c_t.png', dpi=300, format='png')

        plt.figure(7)
        plt.clf()
        plt.plot(v_prop_media[:,1], t_h_bore[:,1], 'o')
        plt.ylabel('Bore height (V)')
        plt.xlabel('Velocity c (m/s)')
        plt.grid(True)
        plt.savefig('c_h.png', dpi=300, format='png')


        p_v = np.polyfit (v_prop_media[nonan,1], vmax_bore[nonan,0], 1)
        retta = np.poly1d(p_v)
        #print(p_v)

        plt.figure(7)
        plt.clf()
        plt.hold(True)
        plt.plot(v_prop_media[:,1], vmax_bore[:,0], 'o')
        plt.plot(np.arange(1, 2.3, 0.2), retta(np.arange(1, 2.3, 0.2)), '-b')
        plt.ylabel('Flow velocity (m/s)')
        plt.xlabel('Velocity c (m/s)')
        plt.grid(True)
        plt.savefig('c_v.png', dpi=300, format='png')

        plt.figure(7)
        plt.clf()
        plt.plot(vmax_bore[:,0], vmax_bore[:,1], 'o', label='c')
        plt.ylabel('Velocity DOP media (m/s)')
        plt.xlabel('Velocity DOP picco (m/s)')
        plt.axis('equal')
        plt.xlim(0.5, 1.5)
        plt.grid(True)
        plt.savefig('vmax_vmax.png', dpi=300, format='png')


    attenzione = np.logical_or(vmax_bore[:,0]<0.2, t_h_bore[:,0]<3.0)


    if 'Vuoto' in bbx or 'all' in bbx:
        continue
    plt.close('all')



    """
    CHECKS
    """
    #print h_moto.shape
    check_h = h_moto[:,0]>=30
    attenzione = np.logical_or(attenzione, check_h)
    ##print test_list[attenzione]

    check_v = v_moto<=0.1
    attenzione = np.logical_or(attenzione, np.any(check_v,1))
    ##print test_list[attenzione]

    check_nans(h_moto[:,0], v_moto)


    check_nans(h_moto[:,0], v_moto)

    h_moto_m, err_h_m = h_moto[:,0], h_moto[:,1]
    v_moto_m = v_moto[:,1]
    #v_moto_m = np.amax(v_moto,1)
    oo_h_m= orient



    """
        FIGURE

    """




    motosi = np.logical_not(np.isnan(v_moto_m))

    colore = {'b4_':'r', 'b7c':'b', 'b8c':'g', 'b4c':'y'}
    #marker_c = {'0':'o', '45':'^', '90':'s'}
    e_color = {'0':'k', '45':'r', '90':'b'}
    colore_o = {'0':'w', '45':'w', '90':'w'}
    colore_o = {'0':'none', '45':'none', '90':'none'}


    a, b, c = dim_boulders[bbx]
    a, b = sorted([a,b])
    ldiag = 0.5*np.sqrt(2.)
    marker_c = { '0':[(-0.5*a,-0.5*b),(-0.5*a,0.5*b),(0.5*a,0.5*b), (0.5*a,-0.5*b),(-0.5*a,-0.5*b),], \
                '45': [(-ldiag*a,0.),((b-a)*ldiag, b*ldiag),(b*ldiag,(b-a)*ldiag), (0.0,-ldiag*a),(-ldiag*a,0.),], \
                '90':[(-0.5*b,-0.5*a),(-0.5*b,0.5*a),(0.5*b,0.5*a), (0.5*b,-0.5*a),(-0.5*b,-0.5*a),]}

    size_m = {'0': 20*a, '45': 20*a*np.sqrt(2.), '90': 20*a}
    size_m = {'0': 20*a, '45': 20*a*2., '90': 20*a}

    if abs(a-b)<0.3:
        ceoff=0.5
    else:
        ceoff=1.
    size_m = {'0': ceoff*20*a, '45': ceoff*20*a*2., '90': ceoff*20*a}

    figsize=(12.,6.)

    if True:

        dim_b = dim_boulders[bbx]*10**-2

        CD = coeff_drag[bbx][oo]





        fig=plt.figure(11)
        asse1 = plt.subplot(1,1,1)
        plt.suptitle('Incipient motion')
        asse1.hold(True)
        asse1.grid(True)
        ##asse1.axis([0, 50, 0, 1.7])
        asse1.axis([0, 1.3, 0, 1.7])
        asse1.set_xlabel("Flow depth h$_n$")
        asse1.set_ylabel("Flow velocity V (m/s)")
        asse1.legend()




        print a,b, c
        from calc_v_teor import  h_somm


        for oo in set_oo:
            sel = oo_h_m == oo
            CD = coeff_drag[bbx][oo]


            from calc_v_teor import H_di_k, k_di_y

            ang_phi = np.radians(float(oo))
            phi_B = ang_phi
            phi_A = 0.5*np.pi-ang_phi
            yA = a * np.sin(phi_A)
            yB = b * np.sin(phi_B)
            yL = yA + yB
            hn = H_di_k(c,teta,d)
            hn = c*10
            k4 = k_di_y(yL, teta)  + c
            hn = 10 * H_di_k(k4,teta,d)

            ##print 'hn',hn


            hs = h_somm(dim_boulders[bbx], teta, int(oo)*2*np.pi/360., d=2.)




            asse1.scatter(calc_k(h_moto_m[sel], teta, d)/hn, v_moto_m[sel], \
                s=size_m[oo], alpha=0.5, \
                c=colore_o[oo], edgecolor=e_color[oo], marker=marker_c[oo], label="{0} {1}".format(bbx, oo))

            hh_err = np.empty((np.count_nonzero(sel),2))
            hh_err[:,0] = h_moto_m[sel]-err_h_m[sel]
            hh_err[:,1] = h_moto_m[sel]+err_h_m[sel]

            hv_err = np.empty((np.count_nonzero(sel),2))
            hv_err[:,0] = v_moto_m[sel]
            hv_err[:,1] = v_moto_m[sel]
            asse1.plot(calc_k(hh_err.T, teta, d)/hn, hv_err.T, alpha=0.5, \
                c=e_color[oo], label='')

            vh_err = np.empty((np.count_nonzero(sel),2))
            vh_err[:,0] = h_moto_m[sel]
            vh_err[:,1] = h_moto_m[sel]

            vv_err = np.empty((np.count_nonzero(sel),2))
            vv_err[:,0] = v_moto_m[sel]
            vv_err[:,1] = np.maximum(v_moto[sel,0]+err_v_moto[sel,0], v_moto[sel,1]+err_v_moto[sel,1])

            asse1.plot(calc_k(vh_err.T, teta, d)/hn, vv_err.T, alpha=0.3, lw=1.2, \
                c=e_color[oo], label='')




            try:
                fig12
                fig14
                fig15
            except:
                fig12=plt.figure(12)
                asse12 = plt.subplot(1,1,1)
                fig14=plt.figure(14, figsize=figsize)
                asse14a = plt.subplot(1,2,1)
                asse14b = plt.subplot(1,2,2)
                fig15=plt.figure(15, figsize=figsize)
                asse15a = plt.subplot(1,2,1)
                asse15b = plt.subplot(1,2,2)

            ##if 'b4c' == bbx or  oo=='0':
            if oo=='0':
                label = "{0} {1}".format(bbx, oo)
                label = "W={: >3d} g".format(int(pesi[bbx]))
            else:
                label = ""


            massa_b_kg = pesi[bbx]/1000

            col = colore[bbx]
            col = (0.6+0.4*abs(97-(pesi[bbx]-53))/97.) *np.ones(3,)
            col = col_pesi[bbx] * np.ones(3,)

            ls = {'0': '-', '45':'-.', '90':'--'}

            asse12.scatter(calc_k(h_moto_m[sel], teta, d)/hn, v_moto_m[sel], \
                s=size_m[oo], alpha=0.8, \
                c=col, marker=marker_c[oo], label=label)

            #asse12.plot(hh_err.T, hv_err.T, alpha=0.5, \
                #c=colore[bbx], label='')

            #asse12.plot(vh_err.T, vv_err.T, alpha=0.3, lw=0.5, \
                #c=colore[bbx], label='')


            asse14a.scatter(calc_k(h_moto_m[sel], teta, d)/hn, v_moto_m[sel], \
                s=size_m[oo], alpha=0.8, \
                c=col, marker=marker_c[oo], label=label)

            if bbx=='b4_' or bbx=='b7c':
                asse15a.scatter(calc_k(h_moto_m[sel], teta, d)/hn, v_moto_m[sel], \
                s=size_m[oo], alpha=0.8, \
                c=col, marker=marker_c[oo], label=label)
            else:
                asse15b.scatter(calc_k(h_moto_m[sel], teta, d)/hn, v_moto_m[sel], \
                s=size_m[oo], alpha=0.8, \
                c=col, marker=marker_c[oo], label=label)


            dim_b = dim_boulders[bbx]*10**-2

            h_nand = np.arange(0.0, 10.0, 0.01) *10**-2

            CD = coeff_drag[bbx][oo]

            print
            print oo,bbx,CD, coeff_drag_3d[bbx][oo],

            vmin_Nand = v_soglia_nand_mod(dim_b, massa_b_kg, h_nand, oo, teta, d, k_attrito[bbx], CD=CD, CL=CL, dens_fluido=1000.0, CI=0.0, v_punto=0.0)

            vmin_Nand3d = v_soglia_nand_mod(dim_b, massa_b_kg, h_nand, oo, teta, d, k_attrito[bbx], CD=coeff_drag_3d[bbx][oo], CL=CL, dens_fluido=1000.0, CI=0.0, v_punto=0.0)


            if oo=='45' and (bbx=='b7c' or bbx=='b4_') and False:
                plt.figure()
                plt.hold(True)
                plt.plot(calc_k(h_nand*1000, teta, d)/hn, vmin_Nand)
                plt.plot(calc_k(h_nand*1000, teta, d)/hn, vmin_Nand3d)


                print vmin_Nand
                print  vmin_Nand3d
                raw_input('aaaaaaaaaaaaaaa')
                plt.show()

            vmin_Nand3d = np.where(vmin_Nand3d==np.nanmin(vmin_Nand3d), vmin_Nand3d, np.nan)

            print np.nanmin(vmin_Nand3d)
            print


            if True:
                asse1.plot(calc_k(h_nand*1000, teta, d)/hn, vmin_Nand, color=e_color[oo], lw=2, label='', alpha=0.5)

            if True:
                asse12.plot(calc_k(h_nand*1000, teta, d)/hn, vmin_Nand, ls[oo], color=col, lw=3, alpha=0.8)
                asse14a.plot(calc_k(h_nand*1000, teta, d)/hn, vmin_Nand, ls[oo], color=col, lw=3, alpha=0.8)
                asse14b.plot(calc_k(h_nand*1000, teta, d)/hn, vmin_Nand, ls[oo], color=col, lw=3, alpha=0.8)

                asse14a.plot(calc_k(h_nand*1000, teta, d)/hn, vmin_Nand3d, ls[oo], color=col, lw=3, alpha=0.8)
                asse14b.plot(calc_k(h_nand*1000, teta, d)/hn, vmin_Nand3d, ls[oo], color=col, lw=3, alpha=0.8)


                if bbx=='b4_' or bbx=='b7c':
                    asse15a.plot(calc_k(h_nand*1000, teta, d)/hn, vmin_Nand, ls[oo], color=col, lw=3, alpha=0.8)
                    asse15a.plot(calc_k(h_nand*1000, teta, d)/hn, vmin_Nand3d, ls[oo], color=col, lw=3, alpha=0.8)
                else:
                    asse15b.plot(calc_k(h_nand*1000, teta, d)/hn, vmin_Nand, ls[oo], color=col, lw=3, alpha=0.8)
                    asse15b.plot(calc_k(h_nand*1000, teta, d)/hn, vmin_Nand3d, ls[oo], color=col, lw=3, alpha=0.8)


            asse14b.scatter(calc_k(h_moto_m[sel], teta, d)/hn, vv_err[:,1], \
                s=size_m[oo], alpha=0.8, \
                c=col, marker=marker_c[oo], label=label)




            """

              figura 13

            """


            fig13=plt.figure(13, figsize=figsize)
            asse13a = plt.subplot(1,2,1)
            asse13b = plt.subplot(1,2,2)

            plt.suptitle('Incipient motion')
            asse13a.hold(True)
            asse13a.grid(True)
            asse13a.set_title('Average conditions')
            asse13a.axis([0, 1.3, 0, 1.7])
            asse13a.set_xlabel("Flow depth h$_n$")
            asse13a.set_ylabel("Flow velocity V (m/s)")
            #asse13a.legend()

            asse13b.hold(True)
            asse13b.grid(True)
            asse13b.set_title('Max conditions')
            asse13b.axis([0, 1.3, 0, 1.7])
            asse13b.set_xlabel("Flow depth h$_n$")
            asse13b.set_ylabel("Flow velocity V (m/s)")
            #asse13b.legend()




            asse13a.plot(calc_k(h_nand*1000, teta, d)/hn, vmin_Nand, ls[oo],color=e_color[oo], lw=2, label='', alpha=0.5)
            asse13b.plot(calc_k(h_nand*1000, teta, d)/hn, vmin_Nand, ls[oo],color=e_color[oo], lw=2, label='', alpha=0.5)

            asse13a.plot(calc_k(h_nand*1000, teta, d)/hn, vmin_Nand3d, ls[oo],color=e_color[oo], lw=2, label='', alpha=0.5)

            asse13b.plot(calc_k(h_nand*1000, teta, d)/hn, vmin_Nand3d, ls[oo],color=e_color[oo], lw=2, label='', alpha=0.5)


            asse13a.scatter(calc_k(h_moto_m[sel], teta, d)/hn, v_moto_m[sel], \
                s=size_m[oo], alpha=0.5, \
                c=colore_o[oo], edgecolor=e_color[oo], marker=marker_c[oo], label="{0} {1}".format(bbx, oo))

            asse13b.scatter(calc_k(h_moto_m[sel], teta, d)/hn, vv_err[:,1], \
                s=size_m[oo], alpha=0.5, \
                c=colore_o[oo], edgecolor=e_color[oo], marker=marker_c[oo], label="{0} {1}".format(bbx, oo))





        asse1.legend()
        fig.savefig('cond_medie_{0}.png'.format(bbx), dpi=300, format='png')
        print 'cond_medie_{0}.png'.format(bbx)

        fig13.savefig('cond_medie_{0}_inc.png'.format(bbx), dpi=300, format='png')
        ##plt.show()
        print 'cond_medie_{0}_inc.png'.format(bbx)




    asse12.axis([0, 1.5, 0, 1.7])
    asse12.set_xlabel("Flow depth h$_n$")
    asse12.set_ylabel("Flow velocity V (m/s)")
    asse12.hold(True)
    asse12.grid(True)
    asse12.legend()
    fig12.savefig('cond_medie_all.png', dpi=300, format='png')
    print 'cond_medie_all.png'
    #plt.show()

    asse14a.axis([0, 1.3, 0, 1.7])
    asse14a.set_xlabel("Flow depth h$_n$")
    asse14a.set_ylabel("Flow velocity V (m/s)")
    asse14a.hold(True)
    asse14a.grid(True)
    asse14a.set_title('Average conditions')
    asse14a.legend()

    asse14b.set_title('Max conditions')
    asse14b.axis([0, 1.3, 0, 1.7])
    asse14b.set_xlabel("Flow depth h$_n$")
    asse14b.set_ylabel("Flow velocity V (m/s)")
    asse14b.hold(True)
    asse14b.grid(True)

    fig14.savefig('cond_max_medie_all_inc.png', dpi=300, format='png')

    asse15a.axis([0, 1.3, 0, 1.7])
    asse15a.set_xlabel("Flow depth h$_n$")
    asse15a.set_ylabel("Flow velocity V (m/s)")
    asse15a.hold(True)
    asse15a.grid(True)
    asse15a.set_title('Average conditions: C1 and C2')
    asse15a.legend()

    asse15b.set_title('Average conditions: R1 and R2')
    asse15b.axis([0, 1.3, 0, 1.7])
    asse15b.set_xlabel("Flow depth h$_n$")
    asse15b.set_ylabel("Flow velocity V (m/s)")
    asse15b.hold(True)
    asse15b.grid(True)
    asse15b.legend()

    fig15.savefig('cond_medie_all_inc.png', dpi=300, format='png')
    print 'cond_medie_all_inc.png'



#plt.show()

#!/usr/bin/python

"""
dim_b:

A lunghezza
C altezza
B lunghezza
"""
import numpy as np
g=9.81 ## kg/m2
teta=np.arcsin(1./10)
d=0.02

import matplotlib.pyplot as plt

def h_somm(dim_b, teta, ang_phi, d):
    A, B, C = dim_b

    phi_B = ang_phi
    phi_A = 0.5*np.pi-ang_phi
    yL = A * np.sin(phi_A) + B * np.sin(phi_B)
    k1 = k_di_y(yL, teta)
    k_somm = k1 + C

    h_somm = ( k_somm + d*np.tan(teta) ) * np.cos(teta)
    return h_somm


"""

per calcolare le distanze

"""

def limita(v, vmin=None, vmax=None):
    if not vmax is None:
        v = np.minimum(v, vmax)
    if not vmin is None:
        v = np.maximum(v, vmin)
    return v



def calc_k(H, teta, d):
    return np.maximum(H/np.cos(teta)-d*np.tan(teta),0)



def H_di_k(k, teta, d):
    return np.maximum( k*np.cos(teta) + d*np.sin(teta), 0.)


def y_di_k(k, teta, ymin=None, ymax=None):
    y = k/np.tan(teta)
    y = limita(y, ymin, ymax)
    return y


def k_di_y(y, teta, kmin=None, kmax=None):
    k = y*np.tan(teta)
    k = limita(k, kmin, kmax)
    return k


def k_trasl_di_k(k, A, teta, kmin=None, kmax=None):
    y_tr = y_di_k(k , teta)-A
    return k_di_y(y_tr, teta, kmin, kmax)

def k_trasl_di_y(y, A, teta, kmin=None, kmax=None):
    return k_di_y(y-A, teta, kmin, kmax)

def y_trasl_di_y(y, C, teta, ymin=None, ymax=None):
    k_tr = k_di_y(y , teta)-C
    return y_di_k(k_tr, teta, ymin, ymax)

def y_trasl_di_k(k, C, teta, ymin=None, ymax=None):
    return y_di_k(k-C, teta, ymin, ymax)


def x_di_y(y, ang_phi, xmin=None, xmax=None):
    x = y/np.tan(ang_phi)
    x = limita(x, xmin, xmax)
    return x


def d_di_y(y, ang_phi, dmin=None, dmax=None):
    d = y/np.sin(ang_phi)
    d = limita(d, dmin, dmax)
    return d

def d_di_xy(x, y, dmin=None, dmax=None):
    d = np.sqrt(x*x + y*y)
    d = limita(d, dmin, dmax)
    return d



"""

per calcolare le aree

"""


def calc_area_frontale_bagn(dim_b, oo, k, teta, ang_phi):
    A, B, C = dim_b
    #print oo=='0',oo=='90',oo=='45',oo
    if oo=='0':
        area_frontale_bagn = B * limita(k, 0.0, C)
    elif oo=='90':
        area_frontale_bagn = A * limita(k, 0.0, C)
    elif oo=='45':
        ang_phi = float(oo)*np.pi/180.

        phi_B = ang_phi
        phi_A = 0.5*np.pi-ang_phi

        area_B = calc_area_frontale_bagn_phi(k, (B, C), teta, phi_B)
        area_A = calc_area_frontale_bagn_phi(k, (A, C), teta, phi_A)

        area_frontale_bagn = area_A * np.cos(phi_A) + area_B * np.cos(phi_B)

    return area_frontale_bagn



def calc_area_frontale_bagn_phi(k, dim_fc, teta, ang_phi):
    """
    Calcolo area frontale bagnata per un lato A
    con angolo teta con la perpendicolare al flusso.
    """
    A, C = dim_fc
    xA = A * np.cos(ang_phi)
    yA = A * np.sin(ang_phi)
    #
    k1 = k_di_y(yA, teta)
    k2 = C
    k3 = C + k1
    ##
    if teta!=0:
        if ang_phi==0.5*np.pi or ang_phi==0.0:
            return A * limita(k, 0.0, C)
        else:
            y = y_di_k(k, teta, ymin=0.0)
            d = d_di_y(y, ang_phi, dmin=0.0, dmax=A)

            area1M = 0.5 * A * k1
            area1 = 0.5* d * limita(k, 0.0, k1)
            lv = C - k_trasl_di_k(k, yA, teta, kmin=C-k1, kmax=C)
            lh = A - d_di_y(y_trasl_di_k(k, C, teta, ymin=0.0, ymax=yA), ang_phi, dmin=0.0, dmax=A)
            area3 = area1M - 0.5 * lv * lh
        #print area3
    else:
        area1 = 0.0
        area3 = 0.0
    area2 =  A * limita(k-k1, vmin=0.0, vmax=C-k1)

    area = area1 + area2 + area3

    return area



def calc_area_base_bagn(k, dim_fc, teta, oo):
    A, B = dim_fc
    ang_phi = np.radians(float(oo))
    if oo=='0':
        area_base_bagn = B * y_di_k(k, teta, ymin=0.0, ymax=A)
    elif oo=='90':
        area_base_bagn = A * y_di_k(k, teta, ymin=0.0, ymax=B)
    elif oo=='45':
        ang_phi = np.pi*0.25
        area_base_bagn = calc_area_base_bagn_phi(k, dim_fc, teta, ang_phi)
    return area_base_bagn



def calc_area_base_bagn_phi(k, dim_fc, teta, ang_phi):
    """
        A o           oooo  B
           o      oooo
    90-phi  o oooo   phi
    ---------------------
    """
    A, B = dim_fc
    y = y_di_k(k, teta, ymin=0.0)
    #
    phi_B = ang_phi
    phi_A = 0.5*np.pi-ang_phi
    #
    xB = B*np.cos(phi_B)
    yB = B*np.sin(phi_B)
    xA = A*np.cos(phi_A)
    yA = A*np.sin(phi_A)
    yL = yA + yB
    #
    y_tr = np.minimum(yA, yB)
    y_pl = np.maximum(yA, yB)
    #
    xb_A=x_di_y(y_tr, phi_A)
    xb_B=x_di_y(y_tr, phi_B)
    #
    tr_1A = 0.5 * x_di_y(y, phi_A, xmin=0.0, xmax=xb_A) * limita(y, 0.0, y_tr)
    tr_1B = 0.5 * x_di_y(y, phi_B, xmin=0.0, xmax=xb_B) * limita(y, 0.0, y_tr)
    #
    parall_2 = ( xb_A + xb_B ) * limita(y-y_tr, 0.0, y_pl-y_tr)
    #
    h_trap_M = yL-y_pl
    bM_trap_B = h_trap_M * np.tan(phi_B)
    bM_trap_A = h_trap_M * np.tan(phi_A)
    #
    h_tr_sup = limita(yL-y, 0.0, yL-y_pl)
    bm_A =  h_tr_sup * np.tan(phi_A)
    bm_B =  h_tr_sup * np.tan(phi_B)
    #
    trap_3A = 0.5* (bM_trap_A + bm_A) * limita(y-y_pl, 0.0, yL-y_pl)
    trap_3B = 0.5* (bM_trap_B + bm_B) * limita(y-y_pl, 0.0, yL-y_pl)
    #
    area = tr_1A + tr_1B + parall_2 + trap_3A + trap_3B
    #
    return area




"""

per calcolare i baricentri

"""
def calc_bar_area_af(k, dim_b, teta, ang_phi):
    A, B, C = dim_b

    phi_B = ang_phi
    phi_A = 0.5*np.pi-ang_phi

    bar_a, area_a, area_proiett_a = calc_bar_area_af_phi(k, (A, C), teta, phi_B)
    bar_a[:,0] = -bar_a[:,0]
    bar_b, area_b, area_proiett_b = calc_bar_area_af_phi(k, (B, C), teta, phi_A)

    #plt.figure()
    #plt.plot(area_a)
    #plt.plot(area_b)
    #plt.show()

    for ii in range(3):
        bar[:,ii] = (bar_a[:,ii] * area_a + bar_b[:,ii] * area_b) / (area_a + area_b)

    area_proiett_a = area_a * np.cos(phi_A)
    area_proiett_b = area_b * np.cos(phi_B)
    area_proiett = area_proiett_a + area_proiett_b
    return bar, area_proiett






def calc_bar_area_af_phi(k, dim_fc, teta, ang_phi):
    #
    BAR = np.zeros((k.shape[0], 3))
    #
    B, C = dim_fc

    phi_B = ang_phi
    xB = B*np.cos(phi_B)
    yB = B*np.sin(phi_B)

    if teta==0:
        BAR[:,0] = 0.5 * xB
        BAR[:,1] = 0.5 * yB
        BAR[:,2] = 0.5 * limita(k, 0.0, C)
        area = B * limita(k, 0.0, C)
    elif ang_phi==0.0:
        BAR[:,0] = 0.5 * B
        BAR[:,1] = 0.0
        BAR[:,2] = 0.5 * limita(k, 0.0, C)
        area = B * limita(k, 0.0, C)
    elif ang_phi==0.5*np.pi:
        BAR[:,0] = 0.0
        BAR[:,1] = 0.5 * B
        BAR[:,2] = 0.5 * limita(k, 0.0, C)
        area = B * limita(k, 0.0, C)
    else:
        #
        y = y_di_k(k, teta, ymin=0.0, ymax=yB)
        x_b = x_di_y(y, phi_B, xmin=0.0, xmax=xB)
        d_b = d_di_y(y, phi_B, dmin=0.0, dmax=B)
        k1_b = k_di_y(yB, teta)
        ##

        ## tr1: punti (0,0,0), (0,0,k), (x,y,0)
        x_BAR_tr1 = ( 0. + 0. + x_b )/3.
        y_BAR_tr1 = ( 0. + 0. + limita(y, vmin=0.0, vmax=yB) )/3.
        z_BAR_tr1 = ( 0. + 0. + limita(k, vmin=0.0, vmax=k1_b) )/3.
        area_tr1 = 0.5* d_b * limita(k, 0.0, k1_b)

        ## parall_2: punti (0,0,k1_b), (0,0,k), (xB,yB,0), (xB,yB,k-k1_b) == simmetrico
        x_BAR_parall2 = 0.5 * xB
        y_BAR_parall2 = 0.5 * yB
        z_BAR_parall2 = 0.5 * limita(k, vmin=0.0, vmax=C)
        area_pl2 =  B * limita(k-k1_b, vmin=0.0, vmax=C-k1_b)

        ## tr3: punti (0,0,C), (x_TD,y_TD,C), (x_TD,y_TD,z_TD)
        y_TD = y_di_k(k-C, teta, ymin=0.0, ymax=yB)
        x_TD = x_di_y(y_TD, phi_B, xmin=0.0, xmax=xB)
        z_TD = limita( C-k_di_y(y_TD, teta), C-k1_b, C)
        x_BAR_tr3 = ( 0. + x_TD + x_TD )/3.
        y_BAR_tr3 = ( 0. + y_TD + y_TD )/3.
        z_BAR_tr3 = ( C + C + z_TD )/3.
        area_tr3 = 0.5 * d_di_xy(x_TD, y_TD, 0.0, B) * (C - z_TD)

        ## parall_4: punti(x_TD,y_TD,C), (x_TD,y_TD,z_TD), (xB,yB,C-k1_b), (xB, yB, k-k1_b)
        x_BAR_parall4 = 0.5 * ( xB + x_TD )
        y_BAR_parall4 = 0.5 * ( yB + y_TD )
        #z_BAR_parall4 = 0.25 * ( C + z_TD + C-k1_b + limita(k-k1_b, C-k1_b, C) )
        z_BAR_parall4 = z_TD*0.0 + C -0.5 * k1_b
        area_pl4 = limita( B-d_di_xy(x_TD, y_TD), 0.0, B) * ( C - z_TD)

        area = area_tr1 + area_pl2 + area_tr3 + area_pl4
        #area1M = 0.5 * B * k1_b

        lim_h = k<=k1_b
        BAR[lim_h,0] = x_BAR_tr1[lim_h]
        BAR[lim_h,1] = y_BAR_tr1[lim_h]
        BAR[lim_h,2] = z_BAR_tr1[lim_h]

        lim_h = (k>k1_b) & (k<=C)
        area_sum = (area_tr1 + area_pl2) [lim_h]
        BAR[lim_h,0] = ( x_BAR_tr1 *area_tr1 + x_BAR_parall2 *area_pl2 )[lim_h] / area_sum
        BAR[lim_h,1] = ( y_BAR_tr1 *area_tr1 + y_BAR_parall2 *area_pl2 )[lim_h] / area_sum
        BAR[lim_h,2] = ( z_BAR_tr1 *area_tr1 + z_BAR_parall2 *area_pl2 )[lim_h] / area_sum

        lim_h = k>C
        area_sum = area_tr1[lim_h] + area_pl2[lim_h] + area_tr3[lim_h] + area_pl4[lim_h]
        BAR[lim_h,0] = ( x_BAR_tr1 *area_tr1 + x_BAR_parall2 *area_pl2 + x_BAR_tr3 *area_tr3 + x_BAR_parall4 *area_pl4)[lim_h] / area_sum
        BAR[lim_h,1] = ( y_BAR_tr1 *area_tr1 + y_BAR_parall2 *area_pl2 + y_BAR_tr3 *area_tr3 + y_BAR_parall4 *area_pl4)[lim_h] / area_sum
        BAR[lim_h,2] = ( z_BAR_tr1 *area_tr1 + z_BAR_parall2 *area_pl2 + z_BAR_tr3 *area_tr3 + z_BAR_parall4 *area_pl4)[lim_h] / area_sum

        ## tr

        ##plt.figure()
        ##plt.subplot(121)
        ##plt.plot([0,xB],[0,yB],'-r', lw=4)
        ##plt.plot(x_TD, y_TD, '-+')
        ##plt.plot(x_BAR_tr3, y_BAR_tr3, '-o')
        ##plt.subplot(122)
        ##plt.plot([0,B],[C,C-0.5*k1_b])
        ##plt.plot([B,0.5*B],[C-k1_b,C])
        ##plt.plot([0,B,B,0],[C,C,C-k1_b,C]) ## triangolo
        ##d = np.sqrt(x_BAR_tr3*x_BAR_tr3+y_BAR_tr3*y_BAR_tr3)
        ##plt.plot(d, z_BAR_tr3, '-x')

        ## parall

        ##plt.figure()
        ##plt.subplot(121)
        ##plt.plot(x_BAR_parall4, y_BAR_parall4, '-+')
        ##plt.subplot(122)
        ##plt.plot([0.5*B,B],[C-0.5*k1_b,C])
        ##plt.plot([0,B],[C,C-0.5*k1_b])
        ##plt.plot([B,0.5*B],[C-k1_b,C])
        ##plt.plot([0,B,B,0],[C,C,C-k1_b,C])
        ##d = np.sqrt(x_BAR_parall4*x_BAR_parall4+y_BAR_parall4*y_BAR_parall4)
        ##plt.plot(d, z_BAR_parall4, '-x')

        ##plt.show()

        ##raw_input('eheheheh')

    area_proiett = area * np.cos(phi_B)

    return BAR, area, area_proiett






def calc_bar_area_ab(k, dim_fc, teta, ang_phi):
    #
    BAR = np.zeros((k.shape[0], 3))
    #
    A, B = dim_fc
    #
    y = y_di_k(k, teta, ymin=0.0)

    if teta==0.0:
        BAR[:,0] = 0.0
        BAR[:,1] = 0.5 * B
        BAR[:,2] = 0.0
        area = A*B
    elif ang_phi==0.0:
        area = B * limita(y, 0.0, A)
        BAR[:,0] = 0.5 * B
        BAR[:,1] = 0.5 * limita(y, 0.0, A)
        BAR[:,2] = 0.0
    elif ang_phi==0.5*np.pi:
        area = A * limita(y, 0.0, B)
        BAR[:,0] = -0.5 * A
        BAR[:,1] = 0.5 * limita(y, 0.0, B)
        BAR[:,2] = 0.0
    else:

        """
            A o           oooo  B
            o      oooo
        90-phi  o oooo   phi
        ---------------------
        """
        phi_B = ang_phi
        phi_A = 0.5*np.pi-ang_phi
        #
        xB = B*np.cos(phi_B)
        yB = B*np.sin(phi_B)
        xA = A*np.cos(phi_A)
        yA = A*np.sin(phi_A)
        yL = yA + yB
        xL = xB - xA
        #
        y_pl = np.maximum(yA, yB)
        y_tr = np.minimum(yA, yB)
        #
        xb_A=x_di_y(y_tr, phi_A)
        xb_B=x_di_y(y_tr, phi_B)
        #
        #
        ## tr_1x: punti (0,0), (x,y), (0,y)
        xa = x_di_y(y, phi_A, xmin=0.0, xmax=xb_A)
        xb = x_di_y(y, phi_B, xmin=0.0, xmax=xb_B)
        x_BAR_tr1A = ( 0. + 0. - limita(xa, vmin=0.0, vmax=xb_A) )/3.
        y_BAR_tr1A = ( 0. +  2 * limita(y, vmin=0.0, vmax=y_tr) )/3.
        x_BAR_tr1B = ( 0. + 0. + limita(xb, vmin=0.0, vmax=xb_B) )/3.
        y_BAR_tr1B = ( 0. +  2 * limita(y, vmin=0.0, vmax=y_tr) )/3.


        ##plt.figure()
        ##plt.title(str(ang_phi*360/2/3.14)+'-qui')
        ###plt.plot(xb,'-+')
        ###plt.plot(xa,'-x')
        ##plt.plot([0,0.5*xb_B],[y_tr,0.5*y_tr])
        ##plt.plot([0,0.5*xb_B],[0,y_tr])
        ##plt.plot([0,-0.5*xb_A],[0,y_tr])
        ##plt.plot([0,-0.5*xb_A],[y_tr,0.5*y_tr])
        ##plt.plot([-xb_A,xb_B],[y_tr,y_tr], '-k')
        ##plt.plot([0,0],[0,yL], '-k')
        ##plt.plot([0,-xA,xL,xB,0],[0,yA,yL,yB,0], '-r', lw=2)
        ##plt.plot(x_BAR_tr1A, y_BAR_tr1A, '-+')
        ##plt.plot(x_BAR_tr1B, y_BAR_tr1B, '-+')
        ##plt.show()
        #
        tr_1A = 0.5 * xa * limita(y, 0.0, y_tr)
        tr_1B = 0.5 * xb * limita(y, 0.0, y_tr)
        ##
        ## parall_2: (xb_A, y_tr), (xb_B, y_tr), (xas, y), (xbs, y)
        if B>=A:
            xbs = x_di_y(y, phi_B, xmin=xb_B, xmax=xB)
            xas = xbs - (xb_A+xb_B)
            x_pl_B = xB
            x_pl_A = x_pl_B - (xb_A+xb_B)
        else:
            xas = - x_di_y(y, phi_A, xmin=xb_A, xmax=xA)
            xbs = xas + (xb_A+xb_B)
            x_pl_A = xA
            x_pl_B = x_pl_A - (xb_A+xb_B)
        x_BAR_parall2 = 0.25 * ( -xb_A + xb_B + xas + xbs)
        y_BAR_parall2 = 0.5 * (y_tr + limita(y, vmin=y_tr, vmax=y_pl) )
        #
        parall_2 = ( xb_A + xb_B ) * limita(y-y_tr, 0.0, y_pl-y_tr)
        #
        ##plt.figure()
        ###plt.title(str(ang_phi*360/2/3.14)+'-qui')
        ###plt.plot([0,-0.5*xb_A],[y_tr,0.5*y_tr])
        ##plt.plot([-xb_A,xb_B],[y_tr,y_tr], '-k')
        ##plt.plot([x_pl_A,x_pl_B],[y_pl,y_pl], '-k')
        ##plt.plot([0,0],[0,yL], '-k')
        ##plt.plot([0,-xA,xL,xB,0],[0,yA,yL,yB,0], '-r', lw=2)
        ###plt.plot(x_BAR_parall2,y_BAR_parall2, '-+')
        ####plt.figure()
        ####plt.plot(parall_2)
        ####plt.hlines(( xb_A + xb_B ) *(y_pl-y_tr), 0.0, parall_2.shape[0])
        ##plt.show()

        # trap_3:
        # yG = h/3 + (b+2a)/(a+b)
        #
        b_magg_trpz = xb_A + xb_B
        h_trap_max = yL-y_pl
        h_trap = limita(y-y_pl, 0.0, yL-y_pl)
        #
        bm_A =  h_trap / np.tan(phi_A)
        bm_B =  h_trap / np.tan(phi_B)
        b_min_trpz = b_magg_trpz - (bm_A + bm_B)
        #
        xm_bM =  x_pl_A +0.5 * b_magg_trpz
        xm_bm =  x_pl_A + bm_B + 0.5*b_min_trpz
        #
        trap_3 = 0.5 * (b_magg_trpz + b_min_trpz) * limita(y-y_pl, 0.0, yL-y_pl)
        #
        hG = limita(y-y_pl, 0.0, h_trap_max) * (b_magg_trpz + 2*b_min_trpz)/(b_magg_trpz + b_min_trpz)/3.
        xG = xm_bM + (xm_bm-xm_bM) * hG/h_trap
        x_BAR_tr3 = xG
        y_BAR_tr3 = y_pl + hG

        ###plt.figure()
        ###plt.plot(trap_3)
        ###plt.hlines( b_magg_trpz*0.5* h_trap_max, 0., y.shape[0])
        ##plt.plot(limita(y-y_pl, 0.0, yL-y_pl))
        ##plt.hlines(yL-y_pl , 0., y.shape[0])
        ##plt.hlines(h_trap_max , 0., y.shape[0])

        ##plt.figure()
        ##plt.title(str(ang_phi*360/2/3.14)+'-qui')
        ##plt.plot(xm_bM,y_pl, 'og')
        ##plt.plot(xm_bm,y_pl+h_trap, 'xg')
        ##plt.plot([x_pl_A,0.5*(xL+x_pl_B)],[y_pl,y_pl+0.5*(yL-y_pl)], '-k')
        ##plt.plot([x_pl_A,x_pl_B],[y_pl,y_pl], '-k')
        ##plt.plot([0,-0.5*xb_A],[y_tr,0.5*y_tr])
        ##plt.plot([xL, x_pl_A+0.5*b_magg_trpz],[yL,y_pl], '-k')
        ##plt.plot([-xb_A,xb_B],[y_tr,y_tr], '-k')
        ##plt.plot([xb_A,xb_B],[y_pl,y_pl], '-k')
        ##plt.plot([0,0],[0,yL], '-k')
        ##plt.plot([0,-xA,xL,xB,0],[0,yA,yL,yB,0], '-r', lw=2)
        ##plt.plot(x_BAR_tr3,y_BAR_tr3, '-o')

        ##plt.figure()
        ###plt.plot(x_BAR_tr3)
        ##plt.plot(y_BAR_tr3)
        ##plt.show()


        area = tr_1A + tr_1B + parall_2 + trap_3


        lim_h = y<=y_tr
        area_tr = (tr_1A + tr_1B)[lim_h]
        BAR[lim_h,0] = ( x_BAR_tr1A[lim_h] * tr_1A[lim_h] + x_BAR_tr1B[lim_h] * tr_1B[lim_h] ) / area_tr
        BAR[lim_h,1] = ( y_BAR_tr1A[lim_h] * tr_1A[lim_h] + y_BAR_tr1B[lim_h] * tr_1B[lim_h] ) / area_tr

        lim_h = (y>y_tr) & (y<y_pl)
        area_sum = ( tr_1A + tr_1B + parall_2 )[lim_h]
        BAR[lim_h,0] = ( x_BAR_tr1A * tr_1A + x_BAR_tr1B * tr_1B + x_BAR_parall2 * parall_2)[lim_h] / area_sum
        BAR[lim_h,1] = ( y_BAR_tr1A * tr_1A + y_BAR_tr1B * tr_1B + y_BAR_parall2 * parall_2)[lim_h] / area_sum

        lim_h = y>y_pl
        area_sum = ( tr_1A + tr_1B + parall_2 + trap_3 )[lim_h]
        BAR[lim_h,0] = ( x_BAR_tr1A * tr_1A + x_BAR_tr1B * tr_1B + x_BAR_parall2 * parall_2 + x_BAR_tr3 * trap_3)[lim_h] / area_sum
        BAR[lim_h,1] = ( y_BAR_tr1A * tr_1A + y_BAR_tr1B * tr_1B + y_BAR_parall2 * parall_2 + y_BAR_tr3 * trap_3)[lim_h] / area_sum

        BAR[:,2] = 0.0

        ##plt.figure()
        ##plt.title('qui')
        ###plt.plot( (x_BAR_tr1A * tr_1A + x_BAR_tr1B * tr_1B )/ (tr_1A + tr_1B))
        ###plt.plot( (y_BAR_tr1A * tr_1A + y_BAR_tr1B * tr_1B )/ (tr_1A + tr_1B))
        ###plt.plot(BAR[:,0])
        ###plt.plot(BAR[:,1])
        ###plt.hlines(0.5*( xb_A + xb_B )*(yL-y_pl), 0., BAR.shape[0])
        ###plt.plot(trap_3)
        ###plt.plot(x_BAR_tr3)
        ##plt.show()


    return BAR, area





"""

per calcolare i volumi

"""

def calc_fraz_vol_sommerso(dim_b, oo, H, teta, d):
    A, B, C = dim_b
    Vol_tot = A*B*C
    if oo=='0' or oo=='90':
        bp=H/np.sin(teta)
        bcd=np.maximum(bp-d, 0)
        b1=np.minimum(bcd, B)
        c1=np.tan(teta)*b1
        vol1=A*c1*b1/2
        ##
        cp=np.tan(teta)*bcd
        c2=np.minimum(np.maximum(cp-c1, 0), C-c1)
        vol2=A*B*c2
        ##
        fp=bp-d-C/np.tan(teta)
        fpp=np.minimum(np.maximum(fp, 0), B)
        vol3=A*np.tan(teta)*(2*B-fpp)*fpp/2
        ##
        ##vol3 = A*B*B*np.tan(teta)/2 - ((C*np.cos(teta)+(B+d)*np.sin(teta)-H)/np.sin(teta))**2
        ##
        Vol_somm = vol1+vol2+vol3
    else:
        k = calc_k(H, teta, d)
        A, B, C = dim_b
        if B<A:
            dim_b[0], dim_b[1] = dim_b[1], dim_b[0]
            A, B, C = dim_b
            if B<A: raise

        ## BASI triangolare avanti
        def V_b_triang_avanti(dim_fc, teta, k):
            A, C = dim_fc
            ## k<=k1
            k1 = np.tan(teta)*A*np.sqrt(2.0)/2.
            #Abase1 = 0.5 * k**2/np.tan(teta)
            ll1 = np.minimum(k, k1)
            pir1 = ll1**3/6./np.tan(teta)**2
            ## k<=k2
            k2 = C
            #Abase2 = 0.5 * A*np.sqrt(2)/2. * k1/np.sin(teta)
            #Abase2 = A**2/4./np.cos(teta)
            #h2 = (k-k1)*np.cos(teta)
            #h2 = np.minimum( np.maximum(h2, 0.), (C-k1)*np.cos(teta))
            #fettaV2 = Abase2*h2
            fettaV2 = A**2/4. * np.minimum( np.maximum(k-k1, 0.), k2-k1)
            ##k<=k3
            k3 = k2+k1
            ll = np.minimum( np.maximum(k3-k, 0), k3-k2)
            pir3 =k1*A**2/6. - ll*A**2/6.
            Vol = pir1+fettaV2+pir3
            return Vol


        def V_fetta_centr(dim_b, teta, k):
            A, B, C = dim_b
            if B<A:
                A, B = B, A
            #print A, B, C
            ## k<=k2
            y = k/np.tan(teta)
            y1 = A*np.sqrt(2.0)/2.
            y2 = B*np.sqrt(2.0)/2.
            k1 = np.tan(teta)*y1
            k2 = np.tan(teta)*y2
            km = (y2-y1)*np.tan(teta)
            ll = np.minimum( np.maximum(y-y1, 0.), (B-A)/np.sqrt(2.) )
            Abase1 = 2*y1*ll
            hk1 = np.minimum(np.maximum(k-k1, 0.), C)
            hk2 = np.minimum(np.maximum(k-k2, 0.), C-km)
            h_media = 0.5*(hk1+hk2)
            pez1 = Abase1*h_media
            ## k<=k3
            k3 = C+k1
            k4 = k3+km
            vol3_tot = np.tan(teta)*y1*(y2-y1)**2
            hasc =-(k-k4)
            hasc = np.maximum(np.minimum( -(k-k4), km), 0. )
            Vasc2 = 0.5*hasc**2/np.tan(teta)*2*y1
            pez2 = vol3_tot - Vasc2
            ##
            Vol = pez1+pez2
            return Vol

        def V_dietro(dim_b, teta, k):
            A, B, C = dim_b
            if B<A:
                A, B = B, A
                if B<A: raise
            y = k/np.tan(teta)
            ##
            y1 = A*np.sqrt(2.0)/2.
            y2 = B*np.sqrt(2.0)/2.
            yF=y1+y2
            ##
            k1 = y1 * np.tan(teta)
            k2 = y2 * np.tan(teta)
            kF = yF * np.tan(teta)
            k2U = k2+C
            kFU = kF+C
            ks = (yF-y2) * np.tan(teta)
            ks = kF-k2
            ## k<=kF
            h = np.minimum( np.maximum(k-k2, 0.), ks )
            c = np.minimum( np.maximum(y-y2, 0.), yF-y2 )
            a = 2*y1
            b=2*np.maximum(yF-y, 0.0)
            pez1 = 0.5*h*c*(a-(a-b)/3.)
            ## k<=k2U
            ll = np.minimum(np.maximum(k-kF, 0.), k2U-kF)
            pez2 = y1**2 *ll
            ## k<=kFU
            vol3_tot = y1**2 * ks /3.
            hasc = np.minimum(np.maximum( -(k-kFU), 0.), ks )
            ll = hasc/np.tan(teta)
            Vasc3 =  ll**2 * hasc /3.
            pez3 = vol3_tot - Vasc3
            ##
            Vol = pez1+pez2+pez3
            return Vol

        Vol_somm = 2*V_b_triang_avanti([A, C], teta, k) + \
                   V_fetta_centr(dim_b, teta, k) + \
                   V_dietro(dim_b, teta, k)

    ##
    fraz_vol_sommerso = Vol_somm/Vol_tot
    ##
    return fraz_vol_sommerso, Vol_tot



#"""
#Funzioni utili
#"""


def BAR_peso(dim_b, ang_phi):
    A, B, C = dim_b
    ##print A, B, C
    phi_B = ang_phi
    phi_A = 0.5*np.pi-ang_phi

    xB = B*np.cos(phi_B)
    yB = B*np.sin(phi_B)
    xA = A*np.cos(phi_A)
    yA = A*np.sin(phi_A)
    yL = yA + yB
    xL = xB - xA

    vertici = np.array([[0., 0., 0.], [0., 0., C],
                        [-xA,     yA,     0.],
                        [-xA,     yA,     C ],
                        [ xB,     yB,     0.],
                        [ xB,     yB,     C ],
                        [ xL,     yL,     0.],
                        [ xL,     yL,     C ],
                        ])
    BAR_peso = np.mean(vertici, axis=0)
    return BAR_peso




def calc_braccio_xyz(dir_forza, bar_forza, fulcro):
    print 'calc_braccio_xyz:verificaa?'

    fulcro = np.where(np.isnan(fulcro), bar_forza, fulcro)

    #### fulcro F
    ##F = np.array([ xL, yL, 0.0])
    ####def M_rot_assex(ang):
        ####return np.array([[1., 0., 0.], [0., np.cos(ang), -np.sin(ang)], [0., np.sin(ang), np.cos(ang)]])
    ####BAR_p = np.dot(M_rot_assex(-teta), BAR_peso)
    ####F_r = np.dot(M_rot_assex(-teta), F)

    if dir_forza==0 or dir_forza == (1,0,0):
        braccio = np.sqrt((bar_forza[:,1]-fulcro[:,1])**2.+(bar_forza[:,2]-fulcro[:,2])**2.)
    elif dir_forza==1 or dir_forza == (0,1,0):
        braccio = np.sqrt((bar_forza[:,0]-fulcro[:,0])**2.+(bar_forza[:,2]-fulcro[:,2])**2.)
    elif dir_forza==2 or dir_forza == (0,0,1):
        braccio = np.sqrt((bar_forza[:,0]-fulcro[:,0])**2.+(bar_forza[:,1]-fulcro[:,1])**2.)

    return braccio




def calc_braccio_peso(bar_peso, fulcro, teta):

    if len(fulcro.shape)==1:
        fulcro = np.tile(fulcro, [bar_peso.shape[0],1])

    fulcro = np.where(np.isnan(fulcro), bar_peso, fulcro)

    def M_rot_assex(ang):
        return np.array([[1., 0., 0.], [0., np.cos(ang), -np.sin(ang)], [0., np.sin(ang), np.cos(ang)]])

    fulcro_r = np.empty((0,3))
    for ff in fulcro:
        fulcro_r = np.vstack((fulcro_r, np.dot(M_rot_assex(-teta), ff) ))

    bar_peso_r = np.empty((0,3))
    for bar_ii in bar_peso:
        bar_peso_r = np.vstack((bar_peso_r, np.dot(M_rot_assex(-teta), bar_ii)))

    braccio_peso = np.sqrt((bar_peso_r[:,0]-fulcro_r[:,0])**2.+(bar_peso_r[:,1]-fulcro_r[:,1])**2.)

    return braccio_peso







def calc_peso(dens_solid, vol_solid, fraz_vol_sommerso=0.0, dens_fluido=1000):
    forza_peso = g*vol_solid*(dens_solid - dens_fluido*fraz_vol_sommerso)
    return forza_peso




def fulcro(oo, dim_b, ang_phi):
    A, B, C = dim_b
    ##print A, B, C
    ang_phi = float(oo)*np.pi/180.
    phi_B = ang_phi
    phi_A = 0.5*np.pi-ang_phi

    xB = B*np.cos(phi_B)
    yB = B*np.sin(phi_B)
    xA = A*np.cos(phi_A)
    yA = A*np.sin(phi_A)
    yL = yA + yB
    xL = xB - xA
    if oo=='0' or oo=='90':
        fulcro = np.array([np.nan,max(yA,yB),0.0])
    else:
        fulcro = np.array([xL,yL,0.0])
    return fulcro





def check_vol_bar(k, teta, ang_phi, dim_b):

    import solidi

    A, B, C = dim_b
    ##print A, B, C
    phi_B = ang_phi
    phi_A = 0.5*np.pi-ang_phi

    xB = B*np.cos(phi_B)
    yB = B*np.sin(phi_B)
    xA = A*np.cos(phi_A)
    yA = A*np.sin(phi_A)
    yL = yA + yB
    xL = xB - xA

    vertici = np.array([[0., 0., 0.], [0., 0., C],
                        [-xA,     yA,     0.],
                        [-xA,     yA,     C ],
                        [ xB,     yB,     0.],
                        [ xB,     yB,     C ],
                        [ xL,     yL,     0.],
                        [ xL,     yL,     C ],
                        ])

    spigoli = [[0,1],[0,2],[0,4],[4,6],[2,6],[6,7],[1,3],[2,3],[3,7],[4,5],[1,5],[5,7]]
    facce = [[0,2,3,1],[0,1,5,4],[5,4,6,7],[6,2,3,7],[1,3,7,5],[0,2,6,4]]

    block = solidi.solidi(vertici, spigoli, facce)
    vol_tot, bar_tot = block.volume_e_baricentro()

    normale_al_piano = np.array([0, -np.sin(teta), np.cos(teta)])

    dist = k*np.cos(teta)

    vol_tot_wet = np.empty((0,))
    vol_tot_dry = np.empty((0,))
    bar_tot_wet = np.empty((0,3))
    bar_tot_dry = np.empty((0,3))

    for dd in dist:
        piano_acqua = solidi.piano(vn=normale_al_piano, do=dd)
        sol_sommerso, sol_asciutto = block.taglio_con_piano(piano_acqua)

        #print sol_sommerso.vertex
        #print sol_asciutto.vertex

        if len(sol_sommerso.vertex)>0:
            vol_w, bar_w = sol_sommerso.volume_e_baricentro()
            vol_tot_wet = np.concatenate((vol_tot_wet, np.array([vol_w])))
            bar_tot_wet = np.vstack((bar_tot_wet, bar_w))
        else:
            vol_tot_wet = np.concatenate((vol_tot_wet, np.array([np.nan])))
            bar_tot_wet = np.vstack((bar_tot_wet, np.array([np.nan,np.nan,np.nan])))

        if len(sol_asciutto.vertex)>0:
            vol_d, bar_d = sol_asciutto.volume_e_baricentro()
            vol_tot_dry = np.concatenate((vol_tot_dry, np.array([vol_d])))
            bar_tot_dry = np.vstack((bar_tot_dry, bar_d))
        else:
            vol_tot_dry = np.concatenate((vol_tot_dry, np.array([np.nan])))
            bar_tot_dry = np.vstack((bar_tot_dry, np.array([np.nan,np.nan,np.nan])))




    return vol_tot_dry, vol_tot_wet, bar_tot_dry, bar_tot_wet, vol_tot, bar_tot









def v_soglia_equilibrio(dim_b, massa_b_kg, H, oo, teta, d, k_attrito, CD, CL, dens_fluido=1000.0, CI=0.0, v_punto=0.0):

    v1= v_soglia_momenti(dim_b, massa_b_kg, H, oo, teta, d, k_attrito, CD, CL, dens_fluido=1000.0)
    v3= v_soglia_f_lift(dim_b, massa_b_kg, H, oo, teta, d, k_attrito, CD, CL, dens_fluido=1000.0)
    v4= v_soglia_f_drag(dim_b, massa_b_kg, H, oo, teta, d, k_attrito, CD, CL, dens_fluido=1000.0) ##, CI=CI, v_punto=v_punto)

    ##import matplotlib.pyplot as plt
    ##plt.figure()
    ##plt.plot(v1, label='v1')
    ###plt.plot(v2, label='v2')
    ##plt.plot(v3, label='v3')
    ##plt.plot(v4, label='v4')
    ##plt.legend()
    ##plt.show()

    #v = np.nanmin(np.array([v1, v2, v3, v4]))
    v = np.min(np.array([v1, v3, v4]), axis=0)

    return v





def v_soglia_momenti(dim_b, massa_b_kg, H, oo, teta, d, k_attrito, CD, CL, dens_fluido=1000.0):
    """

    equilibrio momenti

    """
    k = calc_k(H, teta, d)
    A, B, C = dim_b
    peso_b_N = massa_b_kg*g
    ang_phi = float(oo)*np.pi/180.
    phi_B = ang_phi
    phi_A = 0.5*np.pi-ang_phi

    f_somm, Vol_tot = calc_fraz_vol_sommerso(dim_b=dim_b, oo=oo, H=H, teta=teta, d=d)




    vol_tot_dry, vol_tot_wet, bar_tot_dry, bar_tot_wet, vol_tot, bar_tot = check_vol_bar(k, teta, ang_phi, dim_b)

    peso = calc_peso(dens_solid=massa_b_kg/np.prod(dim_b), vol_solid=np.prod(dim_b), \
                     fraz_vol_sommerso=f_somm, dens_fluido=dens_fluido)

    bar_peso_check = BAR_peso(dim_b, ang_phi)

    dens_dry = massa_b_kg/np.prod(dim_b)
    dens_wet = massa_b_kg/np.prod(dim_b)-1000

    vol_tot_wet = np.where(np.isnan(vol_tot_wet), 0.0, vol_tot_wet)
    vol_tot_dry = np.where(np.isnan(vol_tot_dry), 0.0, vol_tot_dry)


    bar_tot_dry = np.where(np.isnan(bar_tot_dry), 0.0, bar_tot_dry)

    bar_peso = ( np.expand_dims(vol_tot_dry*dens_dry,axis=1) * bar_tot_dry + \
                 np.expand_dims(vol_tot_wet*dens_wet,axis=1) * bar_tot_wet ) / \
                 np.expand_dims(vol_tot_wet*dens_wet + vol_tot_dry*dens_dry,axis=1)

    braccio_peso = calc_braccio_peso(bar_peso, fulcro(oo, dim_b, ang_phi), teta)

    #plt.figure()
    #plt.plot(braccio_peso)
    #plt.show()


    dim_fc = A,B
    bar_lift, area_base_bagn = calc_bar_area_ab(k, dim_fc, teta, ang_phi)

    bar_drag_a, area_front_bagn_a, area_fr_pr_a = calc_bar_area_af_phi(k, (A,C), teta, phi_A)
    bar_drag_b, area_front_bagn_b, area_fr_pr_b = calc_bar_area_af_phi(k, (B,C), teta, phi_B)


    braccio_lift =  calc_braccio_xyz(dir_forza=2, bar_forza=bar_lift, fulcro=fulcro(oo, dim_b, ang_phi))
    braccio_drag_a =  calc_braccio_xyz(dir_forza=1, bar_forza=bar_drag_a, fulcro=fulcro(oo, dim_b, ang_phi))
    braccio_drag_b =  calc_braccio_xyz(dir_forza=1, bar_forza=bar_drag_b, fulcro=fulcro(oo, dim_b, ang_phi))

    ##plt.figure()
    ##plt.plot(k, bar_lift[:,0], '-r')
    ##plt.plot(k, bar_lift[:,1], '-b')
    ##plt.plot(k, bar_lift[:,2], '-g')
    ###plt.plot(braccio_drag_a, '-g')
    ###plt.plot(braccio_drag_b, '-b')
    ##plt.show()

    den = dens_fluido*(CD * braccio_drag_a * area_fr_pr_a + CD * braccio_drag_b * area_fr_pr_b + \
                       CL * braccio_lift * area_base_bagn)

    v2 = 2. * braccio_peso * peso/den

    v = np.sqrt(v2)

    return v







def v_soglia_f_lift(dim_b, massa_b_kg, H, oo, teta, d, k_attrito, CD, CL, dens_fluido=1000.0):

    ang_phi = float(oo)*np.pi/180.
    k = calc_k(H, teta, d)
    A, B, C = dim_b

    f_somm, Vol_tot = calc_fraz_vol_sommerso(dim_b=dim_b, oo=oo, H=H, teta=teta, d=d)
    peso = calc_peso(dens_solid=massa_b_kg/np.prod(dim_b), vol_solid=np.prod(dim_b), \
                     fraz_vol_sommerso=f_somm, dens_fluido=dens_fluido)

    dim_fc = A,B
    bar_lift, area_base_bagn = calc_bar_area_ab(k, dim_fc, teta, ang_phi)

    v2 = 2. * peso * np.cos(teta) / (dens_fluido * CL * area_base_bagn)
    v = np.sqrt(v2)

    return v



def v_soglia_f_drag(dim_b, massa_b_kg, H, oo, teta, d, k_attrito, CD, CL, dens_fluido=1000.0, CI=0., v_punto=0.):
    ang_phi = float(oo)*np.pi/180.
    k = calc_k(H, teta, d)
    A, B, C = dim_b
    phi_B = ang_phi
    phi_A = 0.5*np.pi-ang_phi

    f_somm, Vol_tot = calc_fraz_vol_sommerso(dim_b=dim_b, oo=oo, H=H, teta=teta, d=d)
    V_somm = np.prod(dim_b) * f_somm
    peso = calc_peso(dens_solid=massa_b_kg/np.prod(dim_b), vol_solid=np.prod(dim_b), \
                     fraz_vol_sommerso=f_somm, dens_fluido=dens_fluido)


    bar_drag_a, area_front_bagn_a, area_fr_pr_a = calc_bar_area_af_phi(k, (A,C), teta, phi_A)
    bar_drag_b, area_front_bagn_b, area_fr_pr_b = calc_bar_area_af_phi(k, (B,C), teta, phi_B)

    bar_lift, area_base_bagn = calc_bar_area_ab(k, (A,B), teta, ang_phi)

    den = dens_fluido*(CD * (area_fr_pr_a + area_fr_pr_b) + k_attrito * CL * area_base_bagn)

    v2 = 2.* peso * (k_attrito * np.cos(teta) + np.sin(teta))/den
    v = np.sqrt(v2)

    return v




















if __name__=='__main__':
    import matplotlib.pyplot as plt
    ##

    H_cm=np.arange(0.0, 5.0, 0.005)  #*10**-2
    #H_cm=np.arange(0.0, 150.0, 0.1)  #*10**-2


    dim_b_cm = np.array([3.0, 3.0, 3.0])
    dim_b_cm = np.array([2.8, 3.1, 3.0])
    dim_b_cm = np.array([2.8, 5.1, 3.0])

    dim_b_cm = np.array([2.0, 3.0, 2.0])

    peso_b_gf=53.0  # g-forza

    #dim_b_cm = np.array([3.0, 6.0, 3.0])
    #massa_b_g=120.0  #kg

    k_attrito = 0.5
    #k_attrito = 0.0
    k_drag = 1.0

    ###
    teta=np.arcsin(1./10)
    d_cm=2.0


    ### esempio
    dim_b=3*np.ones(3,)
    #dim_b[1] = 6
    print dim_b
    print
    print
    peso_b_gf = 50
    k_attrito = 0.5
    CD = 1.05
    CL = 0.178
    print teta, d

    ## Conversioni

    dim_b = dim_b_cm*10**-2     ## cm --> m
    H = H_cm*10**-2             ## cm --> m
    peso_b_kgf = peso_b_gf*10**-3    ##  g --> kg
    peso_b_N = peso_b_kgf*g
    massa_b_kg = peso_b_N/g
    massa_b_g = massa_b_kg*10**-3
    dens_solid = massa_b_kg/np.prod(dim_b)
    d = d_cm*10**-2
    dens_fluido=1000.0





    ##print peso_b_gf, massa_b_kg, peso_b_N
    ##raw_input('ciao')

    A, B, C = dim_b
    print 'A, B, C',A, B, C

    ## fine conversioni
    k = calc_k(H, teta, d)
    y = y_di_k(k, teta, ymin=0, ymax=None)

    dim_b = A, B, C
    dim_fc = A, C


    ##dim_fc = A, C
    ##diff_a = calc_area_frontale_bagn_phi(k, dim_fc, teta, ang_phi=np.pi*0.25)
    #check, diff_a = calc_bar_area_af(k, (A, C), teta, ang_phi=np.pi*0.25)
    #check = calc_area_frontale_bagn_phi(k, (A, C), teta, ang_phi=np.pi*0.25)



    ##oo='0'
    ##ang_phi = 0.0
    ##oo='90'
    ##ang_phi = np.pi*0.5
    #oo='45'
    #ang_phi = np.pi*0.25
    #phi_B = ang_phi
    #phi_A = 0.5*np.pi-ang_phi

    #xB = B*np.cos(phi_B)
    #yB = B*np.sin(phi_B)
    #xA = A*np.cos(phi_A)
    #yA = A*np.sin(phi_A)
    #yL = yA + yB
    #xL = xB - xA

    #hfig = A*B
    #hfig =  (A/np.sqrt(2.0))**2.
    #hfig = 0
    #hfig = C * A
    #hfig = 0.5 * A/np.sqrt(2.0)
    #hfig = 0.5 * C
    #hfig = k1 * A *0.5

    #x_area_base = np.array([0.0, -xA, xL, xB, 0.0])
    #y_area_base = np.array([0.0, yA, yL, yB, 0.0])

    #punto = fulcro(oo, dim_b, ang_phi)
    #print 'fulcro',punto

    #plt.figure()
    #plt.grid(True)
    ##plt.hlines(hfig, 0, k.shape[0])
    ##plt.plot(check,'x-')
    ##plt.plot(diff_a,'+-')
    ##plt.plot(area_ab-diff_a,'+-')


    #plt.plot([0.0, 0.5*(xB-xA)], [yA, yB])
    #plt.plot([-xA+(xB-xA), xB], [yB, yB])

    #plt.plot([-xA, xB], [yA, yB])
    #xB=xA
    #yB=yA
    #plt.plot([-0.5*xA, xB], [0.5*yA, yB])
    #plt.plot([-xA, 0.5*xB], [yA, 0.5*yB])
    #plt.plot([-xA, xB], [yA, yB])
    #plt.plot([xL, 0.], [yL, 0.])
    #plt.plot(x_area_base, y_area_base)
    #plt.plot(bar[:,0], bar[:,1],'+-')
    #plt.axis("equal")

    #plt.plot(punto[0], punto[1], '*')

    #plt.show()

    #xA = A * np.cos(ang_phi)
    #yA = A * np.sin(ang_phi)
    #
    print 'A, B, C',A, B, C


    #k0 = 0
    #k1 = k_di_y(yA, teta)
    #k2 = C
    #k3 = C + k1
    ##print k2, k3, H_di_k(k2,teta,d)

    print '++++++++++++++++++++++++++++++++++++++++++'
    print '++++++++++++++++++++++++++++++++++++++++++'

    for oo in ['0', '90', '45']:
        ang_phi = np.radians(float(oo))
        ##print dim_b, massa_b_kg, oo, teta, d, k_attrito, CD, CL, dens_fluido

        diff_a = calc_area_base_bagn(k, (A,B), teta, oo)

        bar, area_ab = calc_bar_area_ab(k, (A,B), teta, ang_phi=ang_phi)

        phi_B = ang_phi
        phi_A = 0.5*np.pi-ang_phi
        yL = A * np.sin(phi_A) + B * np.sin(phi_B)
        print oo, 'phi_A, phi_B',np.degrees(phi_A), np.degrees(phi_B)

        xA = A * np.cos(phi_A)
        yA = A * np.sin(phi_A)
        xB = B * np.cos(phi_B)
        yB = B * np.sin(phi_B)
        yL = yA + yB
        xL = xB - xA
        print 'xA, yA, xB, yB', xA, yA, xB, yB

        y_tr = min(yA,yB)
        x_trA = -x_di_y(y_tr, phi_A)
        x_trB = x_di_y(y_tr, phi_B)
        y_pll =  max(yA,yB)
        if B>A:
            x_pllB = x_di_y(y_pll, phi_B)
            x_pllA = x_pllB - (x_trB-x_trA)
            #print B, A, x_pllB, x_pllA, '-----------------'
        else:
            x_pllA = -x_di_y(y_pll, phi_A)
            x_pllB = x_pllA + (x_trB-x_trA)

        x_area_base = np.array([0.0, -xA, xL, xB, 0.0])
        y_area_base = np.array([0.0, yA, yL, yB, 0.0])
        ####

        k0 = 0
        k1_a = k_di_y(yA, teta)
        k1_b = k_di_y(yB, teta)
        k2 = C
        k3_a = C + k1_a
        k3_b = C + k1_b
        k4 = k_di_y(yL, teta)  + C
        k3=max(k3_a, k3_b)

        if False:

            xb_A=x_di_y(y_tr, phi_A)
            xb_B=x_di_y(y_tr, phi_B)



            plt.figure()
            plt.grid(True)
            plt.title(oo)
            plt.axis("equal")

            plt.plot(x_area_base, y_area_base)
            plt.plot([-xb_A, xb_B], [y_tr, y_tr], '*')        # base triangolo
            plt.plot([-xA, xB], [yA, yB])                     ## diagonale AB
            plt.plot([xL, 0.], [yL, 0.])                      ## diagonale LO
            plt.plot([0.0, 0.5*(xB-xA)], [0.0, 0.5*(yB+yA)])  ## mediana AB

            plt.plot([x_trA, x_trB], [y_tr, y_tr])            ## base triangolo
            plt.plot([x_pllA, x_pllB], [y_pll, y_pll])        ## base parall

            plt.plot([0.5*x_trA, x_trB], [0.5*y_tr, y_tr])           ## mediana triangolo
            plt.plot([0.0, x_trB], [0.0, y_tr])           ## mediana triangolo lato AO
            plt.plot([0.5*(x_trB+x_trA), 0.5*(x_pllB+x_pllA)], [y_tr, y_pll])           ## mediana parall


            #plt.plot([-xA, 0.5*xB], [yA, 0.5*yB])
            #plt.plot([-xA, xB], [yA, yB])
            plt.text(-xA, yA, 'A')
            plt.text(xB, yB, 'B')
            plt.text(xL, yL, 'L')
            plt.text(0., 0., 'O')

            print xA, xB, yA, yB

            plt.plot(bar[:,0], bar[:,1],'+-')
            ##plt.show()

            ### check area base bagnata
            ###plt.figure()
            ####plt.plot(area_ab)
            ####plt.plot(diff_a)
            ###plt.plot(diff_a-area_ab)




            plt.show()
            #plt.plot(punto[0], punto[1], '*')

            xA = A * np.cos(ang_phi)
            yA = A * np.sin(ang_phi)

            plt.show()


        if False:

            area_p_check = calc_area_frontale_bagn((A, B, C), oo, k, teta, ang_phi)
            check_area = calc_area_frontale_bagn_phi(k, (B,C), teta, phi_B)
            bar, area = calc_bar_area_af_phi(k, (B, C), teta, phi_B)

            bar_c, area_p = calc_bar_area_af(k, dim_b, teta, ang_phi)

            #plt.figure()
            #plt.plot(area-check_area)
            #plt.figure()
            #plt.plot(area_p-area_p_check)

            ### baricentro in h
            #plt.figure()
            #plt.subplot(121)
            #plt.grid(True)
            #plt.plot([0,xB],[0,yB], '-r', lw=2)
            #plt.plot(bar[:,0], bar[:,1], '-+')
            #plt.xlim(-xA*1.1, xB*1.1)
            #plt.axis('equal')
            #plt.xlabel('x')
            #plt.ylabel('y')

            ### baricentro in v
            #plt.subplot(122)
            #plt.grid(True)
            #plt.plot([0,B,B,0],[0,0,C,C], '-r', lw=2)
            #d = np.sqrt(bar[:,0]**2 + bar[:,1]**2)
            #plt.plot(d, bar[:,-1], '-+')
            #plt.xlim(0.0, B*1.1)
            #plt.xlabel('d')
            #plt.ylabel('k')
            #plt.axis('equal')


            plt.figure()
            plt.subplot(121)
            plt.grid(True)
            plt.plot([0,xB,xL,-xA,0],[0,yB,yL,yA,0], '-r', lw=2)
            plt.plot(bar_c[:,0], bar_c[:,1], '-+')
            plt.xlim(-xA*1.1, xB*1.1)
            plt.axis('equal')
            plt.xlabel('x')
            plt.ylabel('y')

            ## baricentro in v
            plt.subplot(122)
            plt.grid(True)
            plt.plot([0,0],[0,C], '-r', lw=2)
            plt.plot([-xA,xB,xB,-xA,-xA],[0,0,C,C,0], '-r', lw=2)
            #d = np.sqrt(bar_c[:,0]**2 + bar_c[:,1]**2)
            plt.plot(bar_c[:,0], bar[:,-1], '-+')
            plt.xlim(-xA*1.1, xB*1.1)
            plt.xlim(0.0, B*1.1)
            plt.xlabel('d')
            plt.ylabel('k')
            plt.axis('equal')




            plt.show()


        v1 = v_soglia_momenti(dim_b, massa_b_kg, H, oo, teta, d, k_attrito, CD, CL, dens_fluido=1000.0)
        v3 = v_soglia_f_lift(dim_b, massa_b_kg, H, oo, teta, d, k_attrito, CD, CL, dens_fluido=1000.0)
        v4 = v_soglia_f_drag(dim_b, massa_b_kg, H, oo, teta, d, k_attrito, CD, CL, dens_fluido=1000.0)

        #plt.figure()
        #plt.title(oo)
        #plt.vlines(H_di_k(k0,teta,d),0.0, 100.)
        #plt.vlines(H_di_k(k1_a,teta,d),0.0, 100.)
        #plt.vlines(H_di_k(k1_b,teta,d),0.0, 100.)
        #plt.vlines(H_di_k(C,teta,d),0.0, 9.)
        #plt.vlines(H_di_k(k4,teta,d),0.0, 8.)
        #plt.plot(H,v1,'-r', label='1: M')
        #plt.plot(H,v3,'-g', label='2: L')
        #plt.plot(H,v4,'-b', label='3: D')
        #plt.xlim()
        #plt.legend()

        plt.figure()
        plt.title(oo+'$^\circ$')
        plt.grid(True)
        plt.vlines(k0/dim_b[-1],0.0, 100., linestyles='--')
        plt.vlines(k1_a/dim_b[-1],0.0, 100., linestyles='--')
        plt.vlines(k1_b/dim_b[-1],0.0, 100., linestyles='--')
        plt.vlines(C/dim_b[-1],0.0, 9., linestyles='--')
        plt.vlines(k4/dim_b[-1],0.0, 8., linestyles='--')
        plt.plot(k/dim_b[-1],v1,'-r', label='1: M', lw=2)
        plt.plot(k/dim_b[-1],v3,'-g', label='2: L', lw=2)
        plt.plot(k/dim_b[-1],v4,'-b', label='3: D', lw=2)
        plt.xlim(0.0, 1.2)
        plt.ylim(0.0, 10.)
        plt.xlabel('Normalized h$_n$')
        plt.ylabel('Velocity (m/s)')
        plt.legend()


        print np.nanmax(v1[H_di_k(k,teta,d)>k4])-np.nanmin(v1[H_di_k(k,teta,d)>k4])
        print np.nanmax(v3[H_di_k(k,teta,d)>k4])-np.nanmin(v3[H_di_k(k,teta,d)>k4])
        print np.nanmax(v4[H_di_k(k,teta,d)>k4])-np.nanmin(v4[H_di_k(k,teta,d)>k4])


    plt.show()


    if False:
        #if True:
        if False:
            ## Controlla integrale area
            Afull1 = dim_b[0]*dim_b[-1]
            Afull2 = dim_b[1]*dim_b[-1]
            Afull3 = dim_b[0]*dim_b[-1] + dim_b[1]*dim_b[-1]
            Aeff1 = Afull1
            Aeff2 = Afull2
            Aeff3 = Afull3*np.sqrt(2.0)/2.0

            A = dict()
            for oo in ['0', '45', '90']:
                A[oo] = calc_area_impatto(dim_b, oo, k, teta)

            import matplotlib.pyplot as plt

            plt.hlines(Aeff1, H[0], H[-1], 'r')
            plt.hlines(Aeff2, H[0], H[-1], 'g')
            plt.hlines(Aeff3, H[0], H[-1], 'b')
            for oo in ['0', '45', '90']:
                plt.plot(H, A[oo])
            plt.grid(True)


        ##for oo in ['0', '45', '90']:
        for oo in [ '45']:

            f_somm, Vol_tot = calc_fraz_vol_sommerso(dim_b=dim_b_cm, oo=oo, H=H_cm, teta=teta, d=d_cm)
            peso = calc_peso(dens_solid=massa_b_kg/np.prod(dim_b), vol_solid=np.prod(dim_b), \
                            fraz_vol_sommerso=f_somm, dens_fluido=dens_fluido)

            if False:
            #if True:
                plt.figure()
                plt.plot(H*100, peso, 'g')
                plt.hlines(calc_peso(dens_solid=massa_b_kg/np.prod(dim_b), vol_solid=np.prod(dim_b), \
                                fraz_vol_sommerso=0.0), 0, 5, 'g')
                plt.hlines(calc_peso(dens_solid=massa_b_kg/np.prod(dim_b), vol_solid=np.prod(dim_b), \
                                fraz_vol_sommerso=1.0), 0, 5, 'r')
                plt.xlabel('H (cm)')
                plt.ylabel('peso (cm)')
                plt.grid(True)



        import matplotlib as mpl
        from mpl_toolkits.mplot3d import Axes3D

        #fig = plt.figure()
        #ax = fig.gca(projection='3d')
        #ax.scatter(bar[:,0], bar[:,1], bar[:,2])
        #plt.xlabel('x')
        #plt.ylabel('y')
        #plt.grid(True)


        plt.figure(2)
        plt.plot(bar[:,0]*np.sqrt(2.), bar[:,2], '-+')
        plt.plot(bar[:,1]*np.sqrt(2.), bar[:,2], '-x')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.grid(True)



        plt.show()

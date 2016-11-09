#!/usr/bin/python
# author: Andrea Bressan
#

import numpy as np


def tetraedro(vertici):
    """
    calcolo il volume e il baricentro
    """
    matrice = np.ones((4,4))
    matrice[:vertici.shape[0], :vertici.shape[1]] = vertici
    vol = np.abs(np.linalg.det(matrice))/6.0

    bar = np.sum(vertici, axis=0)/4.0

    return vol,bar



class piano(object):
    def __init__(self, vn, do):
        self.vn = vn   ##vettore_norm
        self.do = do   ## dist_origine



class solidi(object):
    def __init__(self, vertici, spigoli, facce):
        self.faces = facce                     ## una faccia e' lista ordinata di indici di vertici
        self.edges = spigoli                   ## array di coordinate
        self.vertex = vertici                  ## array di coordinate


    def taglio_con_piano(self,eq_piano):  ## piano = vettore normale al piano, distanza dall'origine

        vertici0 = np.empty((0,3))
        vertici1 = np.empty((0,3))
        p_scal_vertici = np.dot(self.vertex, eq_piano.vn)

        nuovi_ii = dict()
        nuovi_ii[0] = np.zeros(p_scal_vertici.shape, dtype=int)
        nuovi_ii[1] = np.zeros(p_scal_vertici.shape, dtype=int)
        nuovi_ii[0][1:] = np.cumsum(p_scal_vertici[:-1]<=eq_piano.do)
        nuovi_ii[1][1:] = np.cumsum(p_scal_vertici[:-1]>=eq_piano.do)

        #print 'indici'
        #print nuovi_ii[0]
        #print nuovi_ii[1]

        id_solid_vert = np.zeros(p_scal_vertici.shape, dtype=int)
        id_solid_vert[p_scal_vertici >eq_piano.do] = 1
        id_solid_vert[p_scal_vertici <eq_piano.do] = 0
        id_solid_vert[p_scal_vertici==eq_piano.do] = -1

        vertici0 = np.vstack((vertici0, self.vertex[p_scal_vertici<=eq_piano.do]))
        vertici1 = np.vstack((vertici1, self.vertex[p_scal_vertici>=eq_piano.do]))




        """
        sistemo le facce e gli spigoli
        """
        nuovi_vertici = {}

        """
        sistemo gli spigoli e in caso aggiungo un vertice
        """
        spigolo = {0: [], 1: []}

        for edge in self.edges:
            if id_solid_vert[edge[0]] == -1:
                if id_solid_vert[edge[1]] == -1:
                    spigolo[0].append( [ nuovi_ii[0][edge[0]], nuovi_ii[0][edge[1]]] )
                    spigolo[1].append( [ nuovi_ii[1][edge[0]], nuovi_ii[1][edge[1]]] )
                else:
                    dest=id_solid_vert[edge[1]]
                    spigolo[dest].append( [ nuovi_ii[dest][edge[0]], nuovi_ii[dest][edge[1]]] )
            elif id_solid_vert[edge[1]] == -1:
                    dest=id_solid_vert[edge[0]]
                    spigolo[dest].append( [ nuovi_ii[dest][edge[0]], nuovi_ii[dest][edge[1]]] )
            else:
                if id_solid_vert[edge[0]] == id_solid_vert[edge[1]]:
                    dest = id_solid_vert[edge[0]]
                    indici_nuovi = [ nuovi_ii[dest][edge[0]], nuovi_ii[dest][edge[1]]]
                    spigolo[dest].append(indici_nuovi)
                else:
                    den = p_scal_vertici[edge[0]] - p_scal_vertici[edge[1]]
                    p_inter = self.vertex[edge[0]] * (eq_piano.do-p_scal_vertici[edge[1]])/den  + \
                            self.vertex[edge[1]] * (p_scal_vertici[edge[0]]-eq_piano.do)/den

                    vertici0 = np.vstack((vertici0, p_inter))
                    vertici1 = np.vstack((vertici1, p_inter))

                    if id_solid_vert[edge[0]]==1:
                        spigolo[0].append([nuovi_ii[0][edge[1]], vertici0.shape[0]-1])
                        spigolo[1].append([nuovi_ii[1][edge[0]], vertici1.shape[0]-1])

                        #print 'bb',spigolo[0][-1], spigolo[1][-1]

                    else:
                        spigolo[0].append([nuovi_ii[0][edge[0]], vertici0.shape[0]-1])
                        spigolo[1].append([nuovi_ii[1][edge[1]], vertici1.shape[0]-1])

                    ##print 'bb',spigolo[0][-1], spigolo[1][-1]

                    iijj = tuple(sorted([edge[0], edge[1]]))
                    nuovi_vertici[iijj] = [vertici0.shape[0]-1, vertici1.shape[0]-1]



        faccia0 = []
        faccia1 = []

        spigoli_nuovi = {0: [], 1: []}

        for faccia in self.faces:
            facce = {0: [], 1: []}

            ff_p = id_solid_vert[faccia[-1]]
            ii_p = faccia[-1]

            spigolo_nuovo = {0: [], 1: []}


            for ii in faccia:
                ff_n = id_solid_vert[ii]
                if ff_n == -1:
                    facce[0].append(nuovi_ii[0][ii])
                    facce[1].append(nuovi_ii[1][ii])

                else:
                    if ff_n == ff_p or ff_p==-1:
                        facce[ff_n].append(nuovi_ii[ff_n][ii])
                    else:
                        ii_p_inter = nuovi_vertici[tuple(sorted([ii_p, ii]))]
                        facce[ff_n].append(ii_p_inter[ff_n])
                        facce[ff_n].append(nuovi_ii[ff_n][ii])
                        facce[ff_p].append(ii_p_inter[ff_p])

                        spigolo_nuovo[0].append(nuovi_ii[ff_n][ii])

                ff_p = ff_n
                ii_p = ii

            faccia0.append(facce[0])
            faccia1.append(facce[1])

        solido0 = solidi(vertici0, spigolo[0], faccia0 )
        solido1 = solidi(vertici1, spigolo[1], faccia1 )

        return solido0, solido1



    def volume_e_baricentro(self):

        vertici = np.empty((4,3))
        vertici[-1,:] = self.vertex[-1,:]

        volume_tot = 0.0
        bar_tot = np.zeros((3,))

        for faccia in self.faces:
            for ii,jj in zip(faccia[1:-1],faccia[2:]):
                vertici[0,:] = self.vertex[faccia[0],:]
                vertici[1,:] = self.vertex[ii,:]
                vertici[2,:] = self.vertex[jj,:]
                vol,bar = tetraedro(vertici)
                volume_tot += vol
                bar_tot += bar * vol
        bar_tot = bar_tot/volume_tot

        return volume_tot, bar_tot












    def print_sol(self):
        print 'vertici'
        print self.vertex
        print
        print 'spigoli'
        print self.edges
        print
        print 'facce'
        print self.faces





if __name__=='__main__':

    vertici = np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
    print


    spigoli = [[0,1],[0,2],[0,3], [1,2],[1,3], [2,3]]
    facce = [[0,1,2],[0,1,3],[0,2,3],[1,2,3]]


    rettangolo = solidi(vertici, spigoli, facce)

    rettangolo.print_sol()
    print

    pianoo = piano(vn=np.array([0,0,1]), do=0.5)

    solido1, solido2 = rettangolo.taglio_con_piano(pianoo)

    print 'solido1'

    solido1.print_sol()
    print 'solido2'

    solido2.print_sol()


    vol,bar = rettangolo.volume_e_baricentro()

    print solido1.volume_e_baricentro()[0] + solido2.volume_e_baricentro()[0] - vol
    print vol, solido1.volume_e_baricentro()[0] , solido2.volume_e_baricentro()[0] , vol

    print bar
    print solido1.volume_e_baricentro()[1] , solido2.volume_e_baricentro()[1]


    print (solido1.volume_e_baricentro()[1] *solido1.volume_e_baricentro()[0] + solido2.volume_e_baricentro()[1] *solido2.volume_e_baricentro()[0] )/vol


#!/usr/bin/env python
# ! ## File: dynSIS.py
# ! Module: use networkx graphs!
# ! ## See README.md for more information and use
# !-----------------------------------------------------------------------------
# ! SIS epidemic model algorithm based on the article "Simulation of Markovian epidemic models on large networks"
# ! Copyright (C) 2016 Wesley Cota, Silvio C. Ferreira
# ! 
# ! Please cite the above cited paper as reference to our code.
# ! 
# !    This program is free software: you can redistribute it and/or modify
# !    it under the terms of the GNU General Public License as published by
# !    the Free Software Foundation, either version 3 of the License, or
# !    (at your option) any later version.
# !
# !    This program is distributed in the hope that it will be useful,
# !    but WITHOUT ANY WARRANTY; without even the implied warranty of
# !    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# !    GNU General Public License for more details.
# !
# !    You should have received a copy of the GNU General Public License
# !    along with this program.  If not, see <http://www.gnu.org/licenses/>.
# !-----------------------------------------------------------------------------
# ! Author    : Wesley Cota
# ! Email     : wesley.cota@ufv.br
# ! Date      : October 2016
# !-----------------------------------------------------------------------------
# ! See README.md for more details
# ! This code is available at <https://github.com/wcota/dynSIS-networkx>
# ! For performance, see <https://github.com/wcota/dynSIS> (Fortran implementation)
# ! For pure Python, see <https://github.com/wcota/dynSIS-py>

import numpy as np
import random

def dyn_run(nw, fnOutput, dynp_sam, dynp_lb, dynp_tmax, dynp_pINI):
    net_N = nw.number_of_nodes()
    net_kmax = max(nw.degree().values())
    
    avg_rho = np.zeros(dynp_tmax, np.float64)
    avg_t = np.zeros(dynp_tmax, np.float64)
    avg_sam = np.zeros(dynp_tmax, np.int)

    dyn_ocp = [None]*net_N # any node type is allowed
    dyn_sig = { i : 0 for i in nw.nodes()}
    
    print('\nOk! Doing the samples...')
    
    dyn_dt_pos_max = 0
    for sam in range(1,dynp_sam+1):
        print('\nSample #', sam)
        
        print('|| (random) Initial condition...')
        
        dyn_sig = dict.fromkeys(dyn_sig, 0)
        dyn_ocp = [None]*net_N
        dyn_voc = 0
        dyn_sk = 0
        for ver in np.random.permutation(nw.nodes()):
            dyn_ocp[dyn_voc] = ver
            dyn_voc += 1
            dyn_sig[ver]= 1
            dyn_sk += nw.degree(ver)
            if dyn_voc == int(net_N*dynp_pINI):
                break
    
        # RUN, Forest, RUN!
        dyn_t = 0
        dyn_dt = 0.0
        dyn_dt_pos = 1
    
        print('|| Running dynamics...')
        while dyn_t <= dynp_tmax and dyn_voc > 0:
            # SIS II ALGORITHM
            dyn_p = 1.0*dyn_voc / (dyn_voc + 1.0*dynp_lb * dyn_sk)
        
            if np.random.uniform() < dyn_p: # Cure
                # Select a random occupied vertex and heal.
                pos_ocp = np.random.randint(0,dyn_voc)
                ver = dyn_ocp[pos_ocp]
                
                # Healed
                dyn_sig[ver]= 0
                dyn_sk -= nw.degree(ver)
                dyn_voc -= 1
                dyn_ocp[pos_ocp] = dyn_ocp[dyn_voc]
            else: # Try to infect
                # Select an infected vertex and accept with the probability.
                while True:
                    pos_ocp = np.random.randint(0,dyn_voc)
                    ver = dyn_ocp[pos_ocp]
                    if np.random.uniform() < 1.0*nw.degree(ver) / (1.0*net_kmax):
                        break
            
                # Select one of its neighbors
                ver = np.random.choice(nw.neighbors(ver))
            
                if dyn_sig[ver]== 0: # Infect!
                    dyn_sig[ver]= 1
                    dyn_sk += nw.degree(ver)
                    dyn_ocp[dyn_voc] = ver
                    dyn_voc += 1
    
            if dyn_voc > 0:
                dyn_dt = 1.0/(dyn_voc + 1.0*dynp_lb * (1.0*dyn_sk))
                dyn_t += dyn_dt
            
                while (dyn_t >= dyn_dt_pos): # Save data
                    dyn_dt_pos_max = max(dyn_dt_pos,dyn_dt_pos_max)
                    avg_rho[dyn_dt_pos - 1] += 1.0*dyn_voc/net_N
                    avg_t[dyn_dt_pos - 1] += dyn_t
                    avg_sam[dyn_dt_pos - 1] += 1
                    dyn_dt_pos += 1
    
        # Write output file
        flOutput = open(fnOutput, 'wt')
        print(  '#@ Number of nodes: '+str(net_N),
                '#@ Number of edges: '+str(2*nw.number_of_edges()),
                '#@ Samples: '+str(dynp_sam),
                '#! Infection rate (lambda): '+str(dynp_lb),
                '#! Maximum time steps: '+str(dynp_tmax),
                '#! Fraction of infected vertices (initial condition): '+str(dynp_pINI),
                sep='\n',
                file=flOutput)
            
        for dt_pos in range(0,dyn_dt_pos_max):
            print(1.0*avg_t[dt_pos]/avg_sam[dt_pos], 1.0*avg_rho[dt_pos]/(1.0*sam),
                    file=flOutput)
                
        flOutput.close()
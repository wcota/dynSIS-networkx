#!/usr/bin/env python
# ! ## File: example_read.py
# ! Example of using dynSIS with edges list file and NetworkX
# ! ## See README.md for more information and use
# !-----------------------------------------------------------------------------
# ! SIS epidemic model algorithm based on the article 
# !           "Optimized Gillespie algorithms for the simulation of Markovian
# !            epidemic processes on large and heterogeneous networks"
# ! Copyright (C) 2017 Wesley Cota, Silvio C. Ferreira
# ! 
# ! Please cite the above cited paper (available at <http://wesleycota.com/ ) as reference
# ! to our code.
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
# ! Date      : 27 Mar 2017
# ! Version   : 1.0
# !-----------------------------------------------------------------------------
# ! See README.md for more details
# ! This code is available at <https://github.com/wcota/dynSIS-networkx>
# ! For performance, see <https://github.com/wcota/dynSIS> (Fortran implementation)
# ! For pure Python, see <https://github.com/wcota/dynSIS-py>

import networkx as nx
import dynSIS
import sys

if len(sys.argv) < 3:
    print('You must enter input and output names as arguments!')
    exit()

fnInput = sys.argv[1]
fnOutput = sys.argv[2]

dynp_sam = int(input('How much dynamics samples? '))
dynp_lb = float(input('Value of infection rate lambda (mu is defined as equal to 1) '))
dynp_tmax = int(input('Maximum time steps (it stops if the absorbing state is reached) '))
dynp_pINI = float(input('Fraction of infected vertices on the network as initial condition (is random \
for each sample) '))

G = nx.Graph()

# Read the network
with open(fnInput, 'rt') as f:
    for i in f:
        li = i.strip().split(',')
        G.add_edge(int(li[0]),int(li[1]))
        
# Run dynamics
dynSIS.dyn_run(G, fnOutput, dynp_sam, dynp_lb, dynp_tmax, dynp_pINI)
        
print('')
print('Everything ok!',True)
print('Input file (edges list): '+ fnInput)
print('Output file: '+ fnOutput)
print('')
print('*****Algorithm used: Optimized Gillespie Algorithm for SIS (SIS-OGA, NetworkX)*****')
print('Codes available at <https://github.com/wcota/dynSIS>.')
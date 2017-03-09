# Optimized Gillespie algorithms for the efficient simulation of Markovian epidemic processes on large and heterogeneous networks: SIS-OGA

### [NetworkX](https://networkx.github.io/) implementation

## Versions

[Fortran implementation - for performance](https://github.com/wcota/dynSIS)

[Python implementation - learn and use](https://github.com/wcota/dynSIS-py)

[(this) NetworkX Python implementation - range of options](https://github.com/wcota/dynSIS-networkx)

## Synopsis

This code is a implementation of the SIS-OGA algorithm, as detailed in our paper (to be cited). It receives as input a NetworkX graph object and the dynamical parameters.

For performance, see https://github.com/wcota/dynSIS (Fortran implementation)

## Installation

Python 3 is required, and also the [NumPy](http://www.numpy.org/) and [NetworkX](https://networkx.github.io/) libraries.

## Use

Import:

```
import networkx as nx
import dynSIS
```
After defining the NetworkX graph, just call

```dynSIS.dyn_run(<nx.graph_Object>, <output_file_path>, <number_of_samples>, <infection_rate>, <maximum_time_steps>, <fraction_of_initial_infected_vertices>)```

where ``<output_file_path>`` will be written with the average density of infected vertices versus time.

## Examples

### ```example_karate.py```

This uses the Karate Club network generated by NetworkX. To run, just use

```python example_karate.py <output_file>```

and type the asked parameters.

### ```example_read.py```

You need to provide a file containing the list of edges (__in__ and __out__, two collumns). ID of the vertices must be enumerated sequentially as `1, 2, 3,..., N`, where `N` is the total number of vertices of the network. Here, we assume  __undirected__ and __unweighted__ networks without multiple neither self connections.

Consider, for example, a network with `N=5` vertices represented by:

```
1,2
1,3
2,4
2,5
3,4
```

Examples of datasets and their specifications are available at http://goo.gl/Bm3VsR.

To run, just use

```python example_read.py <edges_file> <output_file>```

and type the asked parameters.

## License

This code is under [GNU General Public License v3.0](http://choosealicense.com/licenses/gpl-3.0/).

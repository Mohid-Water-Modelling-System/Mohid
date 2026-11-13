# Parallel processing

The historical need in numerical models to reduce computational time became a priority to the Mohid development team as an operational hydrodynamic and water quality model to the Tagus Estuary, in Lisbon, Portugal, was implemented using the [Mohid Water](../../water/mohid-water.md) model full capabilities. Thus, parallel processing has been implemented in Mohid Water in 2003, by using [MPICH](http://www.mcs.anl.gov/mpi/mpich/), a free portable implementation of [MPI](http://www.mcs.anl.gov/mpi/), the standard for message-passing libraries.

Currently, and due to the use of the new Intel Fortran compiler both Mohid Water and Mohid Land have parallelization features using [OpenMP](http://www.openmp.org/).

## Parallel processing via MPI

[Mohid Water](../../water/mohid-water.md) uses MIP to distribute the workload with 2 different strategies that can be apllyed simultaneously or in separation: Nested Models and [Domain Decomposition](../../reference/domain-decomposition.md). A tree sutructure (linked list) of Mohid Water instances is created and a map of model dependencies so that each model has a list of models that it needs to communicate with.

The [Mohid Water](../../water/mohid-water.md) ability to run Nested Models was accomplished by creating a linked list of all the models and by attributing to each one a father-son identification, through which the models communicate. The first stage for introducing parallel processing in Mohid was to add the possibility of launching a process by each model to run, and then, using MPICH, establish communication between models. This enables each sub-model to run in a different processor (even if the processor belongs to a different computer, as long as it is in the same network) and in parallel, instead of running all in the same processor and each model having to wait for the others to perform their calculations. Currently Nested Models have a one-way communication, border conditions being passed from father to son. This results in a assyncronous comunication where the fathe model can calculate the next time step without wayting its sons. In most applications fathers are much lighter (faster) than their sons...

Parallel processing as it is presently implemented in Mohid, could not be achieved without [object oriented programming](object-oriented-programming.md) philosophy, as each model is an instance of [class Model](../../water/modules/README.md#module-model#module-model) and no changes, exception made to the implementation of the MPI communications calls needed to be added. Using this feature, computational speed was improved (varying from application to application), as now the whole model will take the same time as the slowest model to run plus the time to communicate with the other processes. Here, the network communication speed plays an important role, as it can become limiting. However, the amount of information passing between models, depending of course on the memory allocated for each model, has not yet proven to be big enough to make a 100 Mbps network connection time limiting.

[Domain Decomposition](../../reference/domain-decomposition.md) is an ongoing project, The present version decouples a domain into several subdomains that communicate among them (2-way) via MPI. See the follow report 

$$
[MPI\_DD\_Report](http://www.mohid.com/PublicData/Products/Manuals/Mohid_MPI_DomainDecomposition.pdf)
$$

 where a description of the MPI/Domain Decomposition functionality is described.

Find here information on how to [setup a MOHID simulation using MPI](../setup-a-mohid-simulation-using-mpi.md) and on [compiling Mohid with MPI](../build/compiling-mohid-with-mpi.md).

The problem with the current MPI implementation is it's complexity. MPI directives are spread throughout the code. This complexity is very hard to understand from the programming point of view. The proposed solution is to consolidate MPI calls in a single point using the [Actor Model](../../reference/actor-model.md) ([

$$
1
$$

](http://en.wikipedia.org/wiki/Actor_model)).

## Parallel processing via OpenMP

Parallel processing using OpenMP is currently being implemented in Mohid by defining directives **to optimize loops**. These directives are defined as comments in the code and therefore need special compilation options. See more on [compiling Mohid with OpenMP](../build/compiling-mohid-with-openmp.md). Without special compilation these instructions are not considered in normal Fortran compilation.

OpenMP parallel processing can be used in **multi core processors** present in the same computer. It cannot be used to parallel processing using processors located in several computers in a cluster.

Loops optimization is introduced in a first phase in loops referring to grid variables (grid indexes, k, j, i) located in the Modifier section of MOHID modules. Modifier loops are possibly used several times or involve a large resource allocation in MOHID simulations, hence these are locations with larger potential resource gains involved in parallelization.

In case of loops with several looping variables, the parallelized variable is chosen according with cost involved in the loop through this variable. E.g. in a 3D loop (k, j, i loop variables) if j dimension is much larger than the k (number of layers) or i then parallel processing is introduced in j variable, since the resource costs and the savings achieved with the parallelization are larger in this loop than in the others.

There is a [Basic OpenMP syntax](basic-openmp-syntax.md) overview here.

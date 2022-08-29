## Cellular Automata Retro Voxel Fluids

![alt text](https://github.com/wvds98/CellularAutomataVoxelPhysics/images/tsunami.png)

![alt text](https://github.com/wvds98/CellularAutomataVoxelPhysics/images/dambreak.png)

![alt text](https://github.com/wvds98/CellularAutomataVoxelPhysics/images/shallow.png)

## Description
This project implements a c++ Cellular automata based fluid simulation added onto the Retro Voxel Template https://github.com/jbikker/WrldTmpl8.
The relevant code can be found in the WaterWorld project, specifically the CAPE class (cape.h, cape.cpp), which implements the simulation. The scenario is initialized in
WaterWorld.cpp. It is particle-free simulation, where the state of the cells is directly transformed to voxels in the world.

See the accompanying Master's Thesis PDF for more information.

##How Does it Work?
The simulation follows along the lines of a traditional Eulerian fluid simulation, trying to solve or approximate the Navier-Stokes equations, but
does so only locally and restricts the movement of material to adjacent cells only. In this cellular automata style approach, values in the grid are updated
as a function of values in its direct neighborhood. The simulation has been implemented on the GPU, using OpenCL.

In order to save space and determine what cells need to be visited, the world is split into bricks of customizable size (currently 4x4x4).
A brick only exists if it is considered relevant to the simulation, it must either contain water, or be adjacent to a cell that does so. The simulation only 
processes bricks that are considered alive, and will automatically allocate new bricks where required as the water flows into new areas.

The simulation pipeline is as follows:

1. Create and destroy new bricks to store simulation data as is required by the current state of the simulation
2. Perform material advection: update the material in each cell as a function of it's neighbours and flow
3. Perform velocity self advection: self advect the momentum and determine new flow
4. Determine divergence: Determine and store the divergence: the amount that the water violates the mass constraints
5. Solve pressure: Iteratively determine a pressure value that will allow us to correct the flow to be mass conserving
6. Apply this correction, and finalize flow by enforcing some constraints.

##Configuring the simulation
An initial simulation scenario is provided and set up in waterworld.cpp. Parameters of the simulation can be tweaked in cape.h and cape.cl.

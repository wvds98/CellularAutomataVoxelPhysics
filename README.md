## Cellular Automata Retro Voxel Fluids

![tsunami](https://github.com/wvds98/CellularAutomataVoxelPhysics/blob/main/images/tsunami.png?raw=true)

![dambreak](https://github.com/wvds98/CellularAutomataVoxelPhysics/blob/main/images/dambreak.png?raw=true)

![shallow](https://github.com/wvds98/CellularAutomataVoxelPhysics/blob/main/images/shallow.png?raw=true)

<br/>

## Description
This project implements a c++ Cellular automata based fluid simulation added onto the Retro Voxel Template https://github.com/jbikker/WrldTmpl8. The simulation is designed
so that new rules modelling more effects, like heat flow, can be easily added as additional independent kernels. The relevant code can be found in the WaterWorld project,
specifically the CAPE class (cape.h, cape.cpp), which implements the simulation. The scenario is initialized in WaterWorld.cpp. It is particle-free simulation,
where the state of the cells is directly transformed to voxels in the world.

See the accompanying Master's Thesis PDF for more information.

## How Does it Work?
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

## Configuring the simulation
An initial simulation scenario is provided and set up in waterworld.cpp. Parameters of the simulation can be tweaked in cape.h and cape.cl.
Some parameters are duplicated in the cape.cl kernel and cape.h, mind that the gpu kernel uses those defined there.

As an additional tip, by mindful and try to box the fluid in a bit, or atleast use decent evaporation levels. Letting a little bit of water
spread into nothingness across a large plane can be costly performance wise.

### Renderer - Invisible Red voxels
For thesis scenario demonstration purposes I simply made red voxels be considered transparant in the renderer
to edit this, goto trace.cl line ~51:
v = brick0[v]; if (v && v != RED) 

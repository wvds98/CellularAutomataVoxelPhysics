#pragma once

namespace Tmpl8
{

	// MANY PARAMS NOT USED BY GPU IMPLEMENTATION, LOCALLY REDIFINED IN CL FILE
	// CAPE - Cellular Automata Physics Engine
    #define USECONCURRENCY 0
    #define PRESSURE_ITERATIONS 8
	#define DEBUG_MODE 1
	#define CAPE_BRICKDIM 4
    #define CAPE_GRIDWIDTH (MAPWIDTH / CAPE_BRICKDIM + 2)
	#define CAPE_GRIDHEIGHT (MAPHEIGHT / CAPE_BRICKDIM + 2)
	#define CAPE_GRIDDEPTH (MAPDEPTH / CAPE_BRICKDIM + 2)
    #define CAPE_GRIDSIZE (CAPE_GRIDWIDTH * CAPE_GRIDHEIGHT * CAPE_GRIDDEPTH)
	#define CAPE_BRICKSIZE (CAPE_BRICKDIM * CAPE_BRICKDIM * CAPE_BRICKDIM)
    #define BIX(x, y, z) ((x)+(y) * (CAPE_GRIDWIDTH) + (z) * (CAPE_GRIDHEIGHT) * (CAPE_GRIDDEPTH))
    #define MIN_BRICKMASS 0.001f
	#define VELOCITY_DAMPENING 0.1f 
    #define AL 0.5f //Advection limiter, prevents oscillations
    #define INVAL (1.0f/AL) 
    #define MINBRICKS 1000 //Minimum amount of bricks to reserve in memory
	#define SIMS_PER_RENDER 3

    #define EVAPORATION 0.001f //amount of mass to remove per cell per second (helps clean up low mass cells that are not visible.)
    #define MINRENDERMASS 0.001f //Dont render voxels below this much mass
    #define GRAVITYENABLED 1
    #define CELLSIZE 0.25f //In meters, essentially multiplier for gravity

	//0.333 is maximum 100% save speed, use up to 1.0f for less clamping and faster flowing water (less viscous), but requires that an appropriate
	//timestep is selected by the user, or simulation may blow up if local velocity becomes too high and negative mass is created
    #define MAXV 0.33f

	class CAPE
	{
	public:
		CAPE() {};
		~CAPE();

		World* world;

		void ExecuteGPUMaterialAdvection();

		void ExecuteGPUVelocityAdvection();

		void ExecuteGPUBrickUpdate();

		void ExecuteGPUWorldSetter();

		void CopyBuffersToDevice();

		void CopyBuffersFromDevice();

		void ExecuteGPUDivergence();

		void ExecuteGPUPressure();

		void ExecuteGPUPressureGradient();

		//Run CA physics for this frame
		void Tick(float deltaTime);
		void Initialise(World* w, uint updateRate);
		void SetMaterialBlock(uint x, uint y, uint z, uint width, uint height, uint depth, float amount, bool clear = true);
		void ConvertToVoxels();
		void ClearVoxels();
		void ClearColorForCell(uint x, uint y, uint z, float timeStep);
		void SetColorForCell(uint x, uint y, uint z, float timeStep);
		void AddMaterial(uint x, uint y, uint z, float amount);
		void ClearMaterial(uint x, uint y, uint z);

	private:
		//Configuration
		uint solveIterations = 1;
		float timeStep = 0.01;
		float brickTime = 0;
		float pipelineTime = 0;
		float worldsetTime = 0;
		float otherTime = 0;
		float advectTime = 0;
		float velupdateTime = 0;
		float divergenceupdatetime = 0;
		float pressuresolvetime = 0;
		float pressuregradienttime = 0;
		float buffertime = 0;
		float prevtime = 0;
		int updates = 0;

		bool reallocated = false;

		Timer timer;
		Timer timer2;
		float simulationTime = 0;
		int cellUpdates = 0;

		vector<float3> ga;

		uint GetBrickIDX(const uint x, const uint y, const uint z);
		float GetData(const uint x, const uint y, const uint z, float* data);
		void SetData(const uint x, const uint y, const uint z, float v, float* data);
		void AddData(const uint x, const uint y, const uint z, float v, float* data);
		void FreeBrick(uint i);
		void ScaleBricksBuffer(float scale);
		void ReallocBricks();
		void DoubleBricks();
		uint NewBrick(uint bx, uint by, uint bz);
		bool CheckCompressMemory();
		void UpdateBricks();

		//Brick addresses, or uint max if brick empty
		uint* grid;
		uint bricks_alive = 0;
		uint bricks_killed = 0;
		uint bricks_allocated = 0;
		uint bricks_reserved = 0;
		uint brick_update_jobs = 0;
		
		//Brick data
		float* brick_m;
		char* brick_static;
		vector<uint> trash;
		vector<uint> trash_addr;
		uint* brick_x;
		uint* brick_y;
		uint* brick_z;

		//Brick update
		uint* brick_jobs;
		uint* brick_a;
		uint* brick_oa;

		float* m_bricks;
		float* m0_bricks;
		float* p_bricks;
		float* p0_bricks;
		float* vx_bricks;
		float* vy_bricks;
		float* vz_bricks;
		float* vx0_bricks;
		float* vy0_bricks;
		float* vz0_bricks;

		//GPU
		Kernel* capeKernel;
		Kernel* materialAdvectionKernel;				//Sets particle values
		Kernel* velocityAdvectionKernel;				//Draws particles
		Kernel* divergenceKernel;
		Kernel* pressureSolverKernel;
		Kernel* pressureGradientKernel;
		Kernel* brickUpdateKernel;
		Kernel* worldSetKernel;

		//Kernel events
		cl_event brickUpdateEvent;
		cl_event materialAdvectionEvent;
		cl_event velocityAdvectionEvent;
		cl_event divergenceEvent;
		cl_event pressureSolverEvent;
		cl_event pressureGradientEvent;
		cl_event worldSetEvent;

		//Brick buffers
		Buffer* grid_buffer;
		Buffer* brick_static_buffer;
		Buffer* brick_m_buffer;
		Buffer* brick_x_buffer;
		Buffer* brick_y_buffer;
		Buffer* brick_z_buffer;

		Buffer* brick_a_buffer;
		Buffer* brick_oa_buffer;
		Buffer* brick_jobs_buffer;

		//Brick data buffers
		Buffer* m_bricks_buffer;
		Buffer* m0_bricks_buffer;
		Buffer* vx0_bricks_buffer;
		Buffer* vy0_bricks_buffer;
		Buffer* vz0_bricks_buffer;
		Buffer* vx_bricks_buffer;
		Buffer* vy_bricks_buffer;
		Buffer* vz_bricks_buffer;
		Buffer* p_bricks_buffer;
		Buffer* p0_bricks_buffer;

		float TotalDivergence();
		void PrintState();
		float MaterialChange(uint x, uint y, uint z, float vxl, float vxr, float vyl, float vyr, float vzl, float vzr);
		float IncomingMomentumX(uint x, uint y, uint z);
		float IncomingMomentumY(uint x, uint y, uint z);
		float IncomingMomentumZ(uint x, uint y, uint z);
		bool IsCellStatic(uint x, uint y, uint z);
		void CellVelocityUpdate(uint x, uint y, uint z, float timeStep);
		void SolvePressure(uint x, uint y, uint z, float timeStep);
		void CellDivergenceUpdate(uint x, uint y, uint z, float timeStep);
		void SetKernelArguments();
		void CreateBuffers(bool reallocate);
		float TotalMass();
		void MaterialAdvection(uint x, uint y, uint z, float timeStep);
		void PressureGradient(uint x, uint y, uint z, float timeStep);
		void RunOverAllBricks(CAPE* cape, void(CAPE::* func)(uint, uint, uint, float), float timeStep);
	};

}

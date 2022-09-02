#include "precomp.h"
#include "CAPE.h"
#include <ppl.h>

__forceinline uint CAPE::GetBrickIDX(const uint x, const uint y, const uint z)
{
	// calculate brick location in top-level grid
	const uint bx = x / CAPE_BRICKDIM;
	const uint by = y / CAPE_BRICKDIM;
	const uint bz = z / CAPE_BRICKDIM;
	if (bx >= CAPE_GRIDWIDTH || by >= CAPE_GRIDHEIGHT || bz >= CAPE_GRIDDEPTH) return UINT32_MAX;
	const uint brickIdx = BIX(bx,by,bz);
	return brickIdx;
}

__forceinline float CAPE::GetData(const uint x, const uint y, const uint z, float* data)
{
	// obtain brick reference from top-level grid if brick does not exist, return "default" value
	const uint bID = grid[GetBrickIDX(x, y, z)]; 
	if (bID == UINT32_MAX) return 0; 
	float* d = data + bID * CAPE_BRICKSIZE;// [bID] ;
	const uint lx = x & (CAPE_BRICKDIM - 1), ly = y & (CAPE_BRICKDIM - 1), lz = z & (CAPE_BRICKDIM - 1);
	return d[lx + ly * CAPE_BRICKDIM + lz * CAPE_BRICKDIM * CAPE_BRICKDIM];
}

__forceinline void CAPE::SetData(const uint x, const uint y, const uint z, float v, float* data)
{
	const uint bID = grid[GetBrickIDX(x, y, z)];
	if (bID == UINT32_MAX) return;
	float* d = data + bID * CAPE_BRICKSIZE;
	const uint lx = x & (CAPE_BRICKDIM - 1), ly = y & (CAPE_BRICKDIM - 1), lz = z & (CAPE_BRICKDIM - 1);
	uint cellIdx = lx + ly * CAPE_BRICKDIM + lz * CAPE_BRICKDIM * CAPE_BRICKDIM;
	d[cellIdx] = v;
}

__forceinline void CAPE::AddData(const uint x, const uint y, const uint z, float v, float* data)
{
	const uint bID = grid[GetBrickIDX(x, y, z)];
	if (bID == UINT32_MAX) return;
	float* d = data + bID * CAPE_BRICKSIZE;
	const uint lx = x & (CAPE_BRICKDIM - 1), ly = y & (CAPE_BRICKDIM - 1), lz = z & (CAPE_BRICKDIM - 1);
	uint cellIdx = lx + ly * CAPE_BRICKDIM + lz * CAPE_BRICKDIM * CAPE_BRICKDIM;
	d[cellIdx] += v;
}

//Reallocates bricks, doubles the array size, inits new memory to 0
void CAPE::ReallocBricks()
{
	uint obr = bricks_reserved;
	bricks_reserved = bricks_reserved * 2 + 1;
	uint ns = bricks_reserved * sizeof(float) * CAPE_BRICKSIZE;
	m_bricks = (float*)realloc(m_bricks, ns);
	m0_bricks = (float*)realloc(m0_bricks, ns);
	vx_bricks = (float*)realloc(vx_bricks, ns);
	vy_bricks = (float*)realloc(vy_bricks, ns);
	vz_bricks = (float*)realloc(vz_bricks, ns);
	vx0_bricks = (float*)realloc(vx0_bricks, ns);
	vy0_bricks = (float*)realloc(vy0_bricks, ns);
	vz0_bricks = (float*)realloc(vz0_bricks, ns);
	p0_bricks = (float*)realloc(p0_bricks, ns);
	p_bricks = (float*)realloc(p_bricks, ns);

	//realloc brick metadata
	ns = bricks_reserved * sizeof(float);
	brick_m = (float*)realloc(brick_m, ns);
	brick_x = (uint*)realloc(brick_x, ns);
	brick_y = (uint*)realloc(brick_y, ns);
	brick_z = (uint*)realloc(brick_z, ns);
	brick_a = (uint*)realloc(brick_a, ns);
	brick_oa = (uint*)realloc(brick_oa, ns);
	brick_jobs = (uint*)realloc(brick_jobs, ns);
	brick_static = (char*)realloc(brick_static, bricks_reserved);//1 byte per element

	//Set the new part of memory to default 0
	ns = (obr + 1) * sizeof(float) * CAPE_BRICKSIZE;
	uint is = obr * CAPE_BRICKSIZE;
	memset(m_bricks + is, 0.f, ns);
	memset(m0_bricks + is, 0.f, ns);
	memset(vx_bricks + is, 0.f, ns);
	memset(vy_bricks + is, 0.f, ns);
	memset(vz_bricks + is, 0.f, ns);
	memset(vx0_bricks + is, 0.f, ns);
	memset(vy0_bricks + is, 0.f, ns);
	memset(vz0_bricks + is, 0.f, ns);
	memset(p0_bricks + is, 0.f, ns);
	memset(p_bricks + is, 0.f, ns);

	//Set the metadata fields
	ns = (obr + 1) * sizeof(float);
	is = obr;
	memset(brick_m + is, 0.f, ns);
	memset(brick_x + is, UINT32_MAX, ns);
	memset(brick_y + is, UINT32_MAX, ns);
	memset(brick_z + is, UINT32_MAX, ns);
	memset(brick_a + is, UINT32_MAX, ns);
	memset(brick_oa + is, UINT32_MAX, ns);
	memset(brick_jobs + is, UINT32_MAX, ns);
	memset(brick_static + is, false, (obr + 1));
}

//Remove brick from active bricks buffer
void CAPE::FreeBrick(uint i)
{
	return; //Temporarily disable
	uint brick_addr = BIX(brick_x[i], brick_y[i], brick_z[i]);
	bricks_killed++;
	brick_x[i] = UINT32_MAX; //Mark as dead
	brick_m[i] = 0; //empty brick
	trash.push_back(i);
	grid[brick_addr] = UINT32_MAX;//Unset brick adress in lookup table
	brick_static[i] = true; //Marked as dead
	brick_oa[i] = brick_addr; //remember the address this brick had.
	brick_jobs[i] = 0; //mark brick to be freed
	bricks_alive--;
	memset(m0_bricks + i * CAPE_BRICKSIZE, 0, CAPE_BRICKSIZE * sizeof(float));
	memset(vx0_bricks + i * CAPE_BRICKSIZE, 0, CAPE_BRICKSIZE * sizeof(float));
	memset(vy0_bricks + i * CAPE_BRICKSIZE, 0, CAPE_BRICKSIZE * sizeof(float));
	memset(vz0_bricks + i * CAPE_BRICKSIZE, 0, CAPE_BRICKSIZE * sizeof(float));
}

//Add brick to active brick buffer and create memory or reuse a dead brick
uint CAPE::NewBrick(uint bx, uint by, uint bz)
{
	uint brick_addr = BIX(bx, by, bz);
	uint bidx = UINT32_MAX;
	if (trash.size() > 0 && false) //reuse brick
	{
		bidx = trash[trash.size() - 1];
		trash.pop_back();
		brick_a[bidx] = brick_addr; //record new adress for gpu
		brick_jobs[bidx] = 1;
		grid[brick_addr] = bidx; //store new adress
		brick_x[bidx] = bx; 
		brick_y[bidx] = by;
		brick_z[bidx] = bz;

		//make sure crucial data is initialised to default values
		brick_m[bidx] = 0;
		brick_static[bidx] = bx == 0 || bx >= CAPE_GRIDWIDTH || by == 0 || by >= CAPE_GRIDHEIGHT || bz == 0 || bz >= CAPE_GRIDDEPTH;
		memset(m0_bricks + bidx * CAPE_BRICKSIZE, 0, CAPE_BRICKSIZE * sizeof(float));
		memset(vx0_bricks + bidx * CAPE_BRICKSIZE, 0, CAPE_BRICKSIZE * sizeof(float));
		memset(vy0_bricks + bidx * CAPE_BRICKSIZE, 0, CAPE_BRICKSIZE * sizeof(float));
		memset(vz0_bricks + bidx * CAPE_BRICKSIZE, 0, CAPE_BRICKSIZE * sizeof(float));
	}
	else //insert at the end
	{
		//If we have ran out of reserved memory, reallocate
		if (bricks_allocated >= bricks_reserved)
		{
			//Make sure we have all relevant info
			if (updates > 1) CopyBuffersFromDevice();
			ReallocBricks();
			reallocated = true;
		}

		//Set default metadata values
		brick_m[bricks_allocated] = 0;
		brick_static[bricks_allocated] = bx == 0 || bx >= CAPE_GRIDWIDTH || by == 0 || by >= CAPE_GRIDHEIGHT || bz == 0 || bz >= CAPE_GRIDDEPTH;
		brick_jobs[bricks_allocated] = 2;
		brick_x[bricks_allocated] = bx;
		brick_y[bricks_allocated] = by;
		brick_z[bricks_allocated] = bz;
		brick_a[bricks_allocated] = brick_addr;
		grid[brick_addr] = bricks_allocated;

		//Set brick lookup in lookup table
		bricks_allocated++;
		bidx = bricks_allocated - 1;
	}
	bricks_alive++;
	return bidx;
}

//Initialises / deletes bricks no longer relevant in the simulation
//Is essentially a simple top level cellular automata, that checks every
//iteration if a brick is alive or dead, and will make neighbours alive if
//required, so that simulation data is always available, but empty regions do not
//take up any memory or require visits in update step.
void CAPE::UpdateBricks()
{
	reallocated = false;

	bool compressed = CheckCompressMemory();
	reallocated = reallocated || compressed;

	for (int i = 0; i < bricks_allocated; i++)
	{
		if (brick_x[i] != UINT32_MAX && !brick_static[i])
		{
			float bm = brick_m[i];

			//Get brick location
			uint bx = brick_x[i];
			uint by = brick_y[i];
			uint bz = brick_z[i];

			if (grid[BIX(bx, by, bz)] == UINT32_MAX)
				int x = 0;

			//Check if neighbours exist
			uint nbc = 0;
			if (bm < MIN_BRICKMASS)
			{
				//if all existing neighbours have near 0 mass remove this brick
				uint x0 = grid[BIX(bx - 1, by, bz)];
				uint x1 = grid[BIX(bx + 1, by, bz)];
				uint y0 = grid[BIX(bx, by - 1, bz)];
				uint y1 = grid[BIX(bx, by + 1, bz)];
				uint z0 = grid[BIX(bx, by, bz - 1)];
				uint z1 = grid[BIX(bx, by, bz + 1)];

				float mx0 = (x0 == UINT32_MAX ? 0 : brick_m[x0]);
				float mx1 = (x1 == UINT32_MAX ? 0 : brick_m[x1]);
				float my0 = (y0 == UINT32_MAX ? 0 : brick_m[y0]);
				float my1 = (y1 == UINT32_MAX ? 0 : brick_m[y1]);
				float mz0 = (z0 == UINT32_MAX ? 0 : brick_m[z0]);
				float mz1 = (z1 == UINT32_MAX ? 0 : brick_m[z1]);

				float nbm = mx0 + mx1 + my0 + my1 + mz0 + mz1;
				if (nbm < MIN_BRICKMASS)
					FreeBrick(i);
			}
			else
			{
				//If any neighbour is not loaded, create it
				if (grid[BIX(bx - 1, by, bz)] == UINT32_MAX) NewBrick(bx - 1, by, bz);
				if (grid[BIX(bx + 1, by, bz)] == UINT32_MAX) NewBrick(bx + 1, by, bz);
				if (grid[BIX(bx, by - 1, bz)] == UINT32_MAX) NewBrick(bx, by - 1, bz);
				if (grid[BIX(bx, by + 1, bz)] == UINT32_MAX) NewBrick(bx, by + 1, bz);
				if (grid[BIX(bx, by, bz - 1)] == UINT32_MAX) NewBrick(bx, by, bz - 1);
				if (grid[BIX(bx, by, bz + 1)] == UINT32_MAX) NewBrick(bx, by, bz + 1);
			}
		}
	}
}

//Check if we should compress memory (when we have many dead bricks not being recycled), e.g. most extreme case
//First and last block in buffer are alive, all others are dead. Helps prevent uncessary brick visits, also 
//at the same time also rescales reserved memory
bool CAPE::CheckCompressMemory()
{
	if (bricks_alive * 4 < bricks_allocated && bricks_allocated > MINBRICKS * 2)
	{
		if (updates > 1) CopyBuffersFromDevice();

		//Contract memory
		bricks_reserved = max((uint)MINBRICKS * 2, bricks_alive * 2);
		uint ns = bricks_reserved * CAPE_BRICKSIZE;

		//New buffers
		float* m_bricks_new = (float*)calloc(ns, sizeof(float));
		float* m0_bricks_new = (float*)calloc(ns, sizeof(float));
		float* vx_bricks_new = (float*)calloc(ns, sizeof(float));
		float* vy_bricks_new = (float*)calloc(ns, sizeof(float));
		float* vz_bricks_new = (float*)calloc(ns, sizeof(float));
		float* vx0_bricks_new = (float*)calloc(ns, sizeof(float));
		float* vy0_bricks_new = (float*)calloc(ns, sizeof(float));
		float* vz0_bricks_new = (float*)calloc(ns, sizeof(float));
		float* p0_bricks_new = (float*)calloc(ns, sizeof(float));
		float* p_bricks_new = (float*)calloc(ns, sizeof(float));

		//Brick metadata
		float* brick_m_new = (float*)calloc(bricks_reserved, sizeof(float));
		uint* brick_x_new = (uint*)calloc(bricks_reserved, sizeof(uint));
		uint* brick_y_new = (uint*)calloc(bricks_reserved, sizeof(uint));
		uint* brick_z_new = (uint*)calloc(bricks_reserved, sizeof(uint));
		uint* brick_a_new = (uint*)calloc(bricks_reserved, sizeof(uint));
		uint* brick_oa_new = (uint*)calloc(bricks_reserved, sizeof(uint));
		uint* brick_jobs_new = (uint*)calloc(bricks_reserved, sizeof(uint));
		char* brick_static_new = (char*)calloc(bricks_reserved, sizeof(char));

		memset(brick_m_new, 0.f, bricks_reserved * sizeof(float));
		memset(brick_x_new, UINT32_MAX, bricks_reserved * sizeof(float));
		memset(brick_y_new, UINT32_MAX, bricks_reserved * sizeof(float));
		memset(brick_z_new, UINT32_MAX, bricks_reserved * sizeof(float));
		memset(brick_a_new, UINT32_MAX, bricks_reserved * sizeof(float));
		memset(brick_oa_new, UINT32_MAX, bricks_reserved * sizeof(float));
		memset(brick_jobs_new, UINT32_MAX, bricks_reserved * sizeof(float));
		memset(brick_static_new, false, bricks_reserved);

		//Copy all alive bricks to new buffers, while maintaining new brick references
		uint brick_idx = 0;
		for (int i = 0; i < bricks_allocated; i++)
		{
			uint brick_addr = BIX(brick_x[i], brick_y[i], brick_z[i]);
			uint t = brick_x[i];
			if (grid[brick_addr] != UINT32_MAX)
			{
				grid[brick_addr] = brick_idx; //Update brick reference in grid 
				//Copy brick to new buffer
				memcpy(m_bricks_new + brick_idx * CAPE_BRICKSIZE, m_bricks + i * CAPE_BRICKSIZE, CAPE_BRICKSIZE * sizeof(float));
				memcpy(m0_bricks_new + brick_idx * CAPE_BRICKSIZE, m0_bricks + i * CAPE_BRICKSIZE, CAPE_BRICKSIZE * sizeof(float));
				memcpy(vx_bricks_new + brick_idx * CAPE_BRICKSIZE, vx_bricks + i * CAPE_BRICKSIZE, CAPE_BRICKSIZE * sizeof(float));
				memcpy(vy_bricks_new + brick_idx * CAPE_BRICKSIZE, vy_bricks + i * CAPE_BRICKSIZE, CAPE_BRICKSIZE * sizeof(float));
				memcpy(vz_bricks_new + brick_idx * CAPE_BRICKSIZE, vz_bricks + i * CAPE_BRICKSIZE, CAPE_BRICKSIZE * sizeof(float));
				memcpy(vx0_bricks_new + brick_idx * CAPE_BRICKSIZE, vx0_bricks + i * CAPE_BRICKSIZE, CAPE_BRICKSIZE * sizeof(float));
				memcpy(vy0_bricks_new + brick_idx * CAPE_BRICKSIZE, vy0_bricks + i * CAPE_BRICKSIZE, CAPE_BRICKSIZE * sizeof(float));
				memcpy(vz0_bricks_new + brick_idx * CAPE_BRICKSIZE, vz0_bricks + i * CAPE_BRICKSIZE, CAPE_BRICKSIZE * sizeof(float));
				memcpy(p0_bricks_new + brick_idx * CAPE_BRICKSIZE, p0_bricks + i * CAPE_BRICKSIZE, CAPE_BRICKSIZE * sizeof(float));
				memcpy(p_bricks_new + brick_idx * CAPE_BRICKSIZE, p_bricks + i * CAPE_BRICKSIZE, CAPE_BRICKSIZE * sizeof(float));

				//brick metadata
				brick_x_new[brick_idx] = brick_x[i];
				brick_y_new[brick_idx] = brick_y[i];
				brick_z_new[brick_idx] = brick_z[i];
				brick_a_new[brick_idx] = brick_a[i];
				brick_oa_new[brick_idx] = brick_oa[i];
				brick_jobs_new[brick_idx] = brick_jobs[i];
				brick_m_new[brick_idx] = brick_m[i];
				brick_static_new[brick_idx] = brick_static[i];
				brick_idx++;
			}
		}
		bricks_allocated = bricks_alive;

		//swap pointers and delete copy buffers
		Swap(m_bricks, m_bricks_new); free(m_bricks_new);
		Swap(m0_bricks, m0_bricks_new); free(m0_bricks_new);
		Swap(vx_bricks, vx_bricks_new); free(vx_bricks_new);
		Swap(vy_bricks, vy_bricks_new); free(vy_bricks_new);
		Swap(vz_bricks, vz_bricks_new); free(vz_bricks_new);
		Swap(vx0_bricks, vx0_bricks_new); free(vx0_bricks_new);
		Swap(vy0_bricks, vy0_bricks_new); free(vy0_bricks_new);
		Swap(vz0_bricks, vz0_bricks_new); free(vz0_bricks_new);
		Swap(p_bricks, p_bricks_new); free(p_bricks_new);
		Swap(p0_bricks, p0_bricks_new); free(p0_bricks_new);

		Swap(brick_x, brick_x_new); free(brick_x_new);
		Swap(brick_y, brick_y_new); free(brick_y_new);
		Swap(brick_z, brick_z_new); free(brick_z_new);
		Swap(brick_a, brick_a_new); free(brick_a_new);
		Swap(brick_oa, brick_oa_new); free(brick_oa_new);
		Swap(brick_jobs, brick_jobs_new); free(brick_jobs_new);
		Swap(brick_m, brick_m_new); free(brick_m_new);
		Swap(brick_static, brick_static_new); free(brick_static_new);
		return true;
	}
	return false;
}

CAPE::~CAPE()
{
	//Delete bricks
	free(m_bricks);
	free(m0_bricks);
	free(vx_bricks);
	free(vy_bricks);
	free(vz_bricks);
	free(vx0_bricks);
	free(vy0_bricks);
	free(vz0_bricks);
	free(p_bricks);
	free(p0_bricks);
	free(grid);
}

//Init memory and parameters
void CAPE::Initialise(World* w, uint updateRate)
{
	world = w;
	timer = Timer();
	timer2 = Timer();
	timeStep = 1.0f / updateRate;
	if (GRAVITYENABLED == 1)
		ga.push_back(float3(0, -9.81 * CELLSIZE, 0));

	//Initialise brick memory, reserve 20 bricks to start
	bricks_reserved = MINBRICKS * 100;
	uint is = bricks_reserved * CAPE_BRICKSIZE * sizeof(float);
	m_bricks = (float*)calloc(bricks_reserved * CAPE_BRICKSIZE, sizeof(float));
	m0_bricks = (float*)calloc(bricks_reserved * CAPE_BRICKSIZE, sizeof(float));
	vx_bricks = (float*)calloc(bricks_reserved * CAPE_BRICKSIZE, sizeof(float));
	vy_bricks = (float*)calloc(bricks_reserved * CAPE_BRICKSIZE, sizeof(float));
	vz_bricks = (float*)calloc(bricks_reserved * CAPE_BRICKSIZE, sizeof(float));
	vx0_bricks = (float*)calloc(bricks_reserved * CAPE_BRICKSIZE, sizeof(float));
	vy0_bricks = (float*)calloc(bricks_reserved * CAPE_BRICKSIZE, sizeof(float));
	vz0_bricks = (float*)calloc(bricks_reserved * CAPE_BRICKSIZE, sizeof(float));
	p_bricks = (float*)calloc(bricks_reserved * CAPE_BRICKSIZE, sizeof(float));
	p0_bricks = (float*)calloc(bricks_reserved * CAPE_BRICKSIZE, sizeof(float));

	brick_m = (float*)calloc(bricks_reserved, sizeof(float));
	brick_x = (uint*)calloc(bricks_reserved, sizeof(uint));
	brick_y = (uint*)calloc(bricks_reserved, sizeof(uint));
	brick_z = (uint*)calloc(bricks_reserved, sizeof(uint));
	brick_a = (uint*)calloc(bricks_reserved, sizeof(uint));
	brick_oa = (uint*)calloc(bricks_reserved, sizeof(uint));
	brick_jobs = (uint*)calloc(bricks_reserved, sizeof(uint));
	memset(brick_a, UINT32_MAX, bricks_reserved * sizeof(uint));
	memset(brick_x, UINT32_MAX, bricks_reserved * sizeof(uint));
	memset(brick_y, UINT32_MAX, bricks_reserved * sizeof(uint));
	memset(brick_z, UINT32_MAX, bricks_reserved * sizeof(uint));
	brick_static = (char*)calloc(bricks_reserved, sizeof(bool));

	//Top level grid pointing to active bricks
	grid = (uint*)(malloc(CAPE_GRIDSIZE * sizeof(uint)));
	if(grid != nullptr) memset(grid, UINT32_MAX, CAPE_GRIDSIZE * sizeof(uint));

	capeKernel = new Kernel("cl/cape.cl", "run");
	materialAdvectionKernel = new Kernel(capeKernel->GetProgram(), "materialAdvection");
	divergenceKernel = new Kernel(capeKernel->GetProgram(), "calculateDivergence");
	pressureSolverKernel = new Kernel(capeKernel->GetProgram(), "solvePressure");
	pressureGradientKernel = new Kernel(capeKernel->GetProgram(), "pressureGradient");
	velocityAdvectionKernel = new Kernel(capeKernel->GetProgram(), "velocityAdvection");
    brickUpdateKernel = new Kernel(capeKernel->GetProgram(), "brickUpdate");
	worldSetKernel = new Kernel(capeKernel->GetProgram(), "worldSetter");
}

void CAPE::ClearMaterial(uint x, uint y, uint z)
{
	if (grid[GetBrickIDX(x, y, z)] != UINT32_MAX)
		SetData(x, y, z, 0, m0_bricks);
}

void CAPE::AddMaterial(uint x, uint y, uint z, float amount)
{
	//Add a single brick offset, to account for the layer of ghost bricks around the grid
	//Similarly substract when indexing back to the worlds voxel space
	x += CAPE_BRICKDIM;
	y += CAPE_BRICKDIM;
	z += CAPE_BRICKDIM;
	uint brickIdx = GetBrickIDX(x, y, z);
	uint brick_addr = grid[brickIdx];
	if (brick_addr == UINT32_MAX)
	{
		const uint bx = x / CAPE_BRICKDIM;
		const uint by = y / CAPE_BRICKDIM;
		const uint bz = z / CAPE_BRICKDIM;
		brick_addr = NewBrick(bx, by, bz);
	}
	AddData(x, y, z, amount, m0_bricks);
	brick_m[brick_addr] += amount;
}

//Writes a block of material of width w, height h, depth z starting at location x, y, z
//Won't overwrite world voxel if clear not true
void CAPE::SetMaterialBlock(uint x0, uint y0, uint z0, uint w, uint h, uint d, float amount, bool clear)
{
	if (x0 + w > MAPWIDTH || y0 + h > MAPWIDTH || z0 + d > MAPWIDTH)
	{
		cout << "Block outside of grid" << endl;
		return;
	}
	for(uint x = x0; x < x0 + w; x++)
		for (uint y = y0; y < y0 + h; y++)
			for (uint z = z0; z < z0 + d; z++)
			{
				uint v = world->Get(x, y, z);
				if (clear)
				{
					world->Set(x, y, z, 0);
					ClearMaterial(x, y, z);
					AddMaterial(x, y, z, amount);
				}
				else if(v == 0)
					AddMaterial(x, y, z, amount);
			}
}

void CAPE::SetKernelArguments()
{
	materialAdvectionKernel->SetArgument(0, grid_buffer);
	materialAdvectionKernel->SetArgument(1, brick_x_buffer);
	materialAdvectionKernel->SetArgument(2, brick_y_buffer);
	materialAdvectionKernel->SetArgument(3, brick_z_buffer);
	materialAdvectionKernel->SetArgument(4, brick_static_buffer);
	materialAdvectionKernel->SetArgument(5, brick_m_buffer);

	materialAdvectionKernel->SetArgument(6, m_bricks_buffer);
	materialAdvectionKernel->SetArgument(7, m0_bricks_buffer);
	materialAdvectionKernel->SetArgument(8, vx0_bricks_buffer);
	materialAdvectionKernel->SetArgument(9, vy0_bricks_buffer);
	materialAdvectionKernel->SetArgument(10, vz0_bricks_buffer);
	materialAdvectionKernel->SetArgument(11, vx_bricks_buffer);
	materialAdvectionKernel->SetArgument(12, vy_bricks_buffer);
	materialAdvectionKernel->SetArgument(13, vz_bricks_buffer);
	materialAdvectionKernel->SetArgument(14, p0_bricks_buffer);
	cl_mem gm = world->GetGridMap();
	materialAdvectionKernel->SetArgument(15, &gm);
	materialAdvectionKernel->SetArgument(16, world->GetBrickBuffer());

	velocityAdvectionKernel->SetArgument(0, grid_buffer);
	velocityAdvectionKernel->SetArgument(1, brick_x_buffer);
	velocityAdvectionKernel->SetArgument(2, brick_y_buffer);
	velocityAdvectionKernel->SetArgument(3, brick_z_buffer);
	velocityAdvectionKernel->SetArgument(4, brick_static_buffer);
	velocityAdvectionKernel->SetArgument(5, brick_m_buffer);
	velocityAdvectionKernel->SetArgument(14, p0_bricks_buffer);

	divergenceKernel->SetArgument(0, grid_buffer);
	divergenceKernel->SetArgument(1, brick_x_buffer);
	divergenceKernel->SetArgument(2, brick_y_buffer);
	divergenceKernel->SetArgument(3, brick_z_buffer);
	divergenceKernel->SetArgument(4, brick_static_buffer);
}

//Creates OpenCL Buffers
//Delete and replace old buffers if reallocate = true
void CAPE::CreateBuffers(bool reallocate = false)
{
	if (reallocate)
	{
		clReleaseMemObject(grid_buffer->deviceBuffer);
		clReleaseMemObject(brick_static_buffer->deviceBuffer);
		clReleaseMemObject(brick_m_buffer->deviceBuffer);
		clReleaseMemObject(brick_x_buffer->deviceBuffer);
		clReleaseMemObject(brick_y_buffer->deviceBuffer);
		clReleaseMemObject(brick_z_buffer->deviceBuffer);
		clReleaseMemObject(brick_a_buffer->deviceBuffer);
		clReleaseMemObject(brick_oa_buffer->deviceBuffer);
		clReleaseMemObject(brick_jobs_buffer->deviceBuffer);

		clReleaseMemObject(m_bricks_buffer->deviceBuffer);
		clReleaseMemObject(m0_bricks_buffer->deviceBuffer);
		clReleaseMemObject(vx_bricks_buffer->deviceBuffer);
		clReleaseMemObject(vy_bricks_buffer->deviceBuffer);
		clReleaseMemObject(vz_bricks_buffer->deviceBuffer);
		clReleaseMemObject(vx0_bricks_buffer->deviceBuffer);
		clReleaseMemObject(vy0_bricks_buffer->deviceBuffer);
		clReleaseMemObject(vz0_bricks_buffer->deviceBuffer);
		clReleaseMemObject(p_bricks_buffer->deviceBuffer);
		clReleaseMemObject(p0_bricks_buffer->deviceBuffer);

		delete(grid_buffer);
		delete(brick_static_buffer);
		delete(brick_m_buffer);
		delete(brick_x_buffer);
		delete(brick_y_buffer);
		delete(brick_z_buffer);
		delete(brick_a_buffer);
		delete(brick_oa_buffer);
		delete(brick_jobs_buffer);

		delete(m_bricks_buffer);
		delete(m0_bricks_buffer);
		delete(vx_bricks_buffer);
		delete(vy_bricks_buffer);
		delete(vz_bricks_buffer);
		delete(vx0_bricks_buffer);
		delete(vy0_bricks_buffer);
		delete(vz0_bricks_buffer);
		delete(p_bricks_buffer);
		delete(p0_bricks_buffer);
		reallocated = false;
	}

	grid_buffer = new Buffer(CAPE_GRIDSIZE, Buffer::DEFAULT, grid);
	brick_x_buffer = new Buffer(bricks_reserved, Buffer::DEFAULT, brick_x);
	brick_y_buffer = new Buffer(bricks_reserved, Buffer::DEFAULT, brick_y);
	brick_z_buffer = new Buffer(bricks_reserved, Buffer::DEFAULT, brick_z);
	brick_a_buffer = new Buffer(bricks_reserved, Buffer::DEFAULT, brick_a, true);
	brick_oa_buffer = new Buffer(bricks_reserved, Buffer::DEFAULT, brick_oa, true);
	brick_jobs_buffer = new Buffer(bricks_reserved, Buffer::DEFAULT, brick_jobs, true);
	brick_m_buffer = new Buffer(bricks_reserved, Buffer::DEFAULT, brick_m, true);
	brick_static_buffer = new Buffer(bricks_reserved / 4, Buffer::DEFAULT, brick_static);

	//Brick data buffers
	uint s = bricks_reserved * CAPE_BRICKSIZE;
	m_bricks_buffer = new Buffer(s, Buffer::DEFAULT, m_bricks);
	m0_bricks_buffer = new Buffer(s, Buffer::DEFAULT, m0_bricks);
	vx0_bricks_buffer = new Buffer(s, Buffer::DEFAULT, vx0_bricks);
	vy0_bricks_buffer = new Buffer(s, Buffer::DEFAULT, vy0_bricks);
	vz0_bricks_buffer = new Buffer(s, Buffer::DEFAULT, vz0_bricks);
	vx_bricks_buffer = new Buffer(s, Buffer::DEFAULT, vx_bricks);
	vy_bricks_buffer = new Buffer(s, Buffer::DEFAULT, vy_bricks);
	vz_bricks_buffer = new Buffer(s, Buffer::DEFAULT, vz_bricks);
	p_bricks_buffer = new Buffer(s, Buffer::DEFAULT, p_bricks);
	p0_bricks_buffer = new Buffer(s, Buffer::DEFAULT, p0_bricks);

	CopyBuffersToDevice();
}

//Fully copy cpu buffers and copy them to gpu => later only if buffer size changes and we have implemented all simulation on gpu
void CAPE::CopyBuffersToDevice()
{
	grid_buffer->CopyToDevice(false);
	brick_x_buffer->CopyToDevice(false);
	brick_y_buffer->CopyToDevice(false);
	brick_z_buffer->CopyToDevice(false);
	brick_a_buffer->CopyToDevice(false);
	brick_oa_buffer->CopyToDevice(false);
	brick_jobs_buffer->CopyToDevice(false);
	brick_m_buffer->CopyToDevice(false);
	brick_static_buffer->CopyToDevice(false);

	m_bricks_buffer->CopyToDevice(false);
	m0_bricks_buffer->CopyToDevice(false);
	vx0_bricks_buffer->CopyToDevice(false);
	vy0_bricks_buffer->CopyToDevice(false);
	vz0_bricks_buffer->CopyToDevice(false);
	vx_bricks_buffer->CopyToDevice(false);
	vy_bricks_buffer->CopyToDevice(false);
	vz_bricks_buffer->CopyToDevice(false);
	p_bricks_buffer->CopyToDevice(false);
	p0_bricks_buffer->CopyToDevice(true);
}

void CAPE::CopyBuffersFromDevice()
{
	grid_buffer->CopyFromDevice(false);
	brick_x_buffer->CopyFromDevice(false);
	brick_y_buffer->CopyFromDevice(false);
	brick_z_buffer->CopyFromDevice(false);
	brick_a_buffer->CopyFromDevice(false);
	brick_oa_buffer->CopyFromDevice(false);
	brick_jobs_buffer->CopyFromDevice(false);
	brick_m_buffer->CopyFromDevice(false);
	brick_static_buffer->CopyFromDevice(false);

	m_bricks_buffer->CopyFromDevice(false);
	m0_bricks_buffer->CopyFromDevice(false);
	vx0_bricks_buffer->CopyFromDevice(false);
	vy0_bricks_buffer->CopyFromDevice(false);
	vz0_bricks_buffer->CopyFromDevice(false);
	vx_bricks_buffer->CopyFromDevice(false);
	vy_bricks_buffer->CopyFromDevice(false);
	vz_bricks_buffer->CopyFromDevice(false);
	p_bricks_buffer->CopyFromDevice(false);
	p0_bricks_buffer->CopyFromDevice(true);
}

//__global char* brick_static, __global float* m_bricks, __global float* m0_bricks, __global float* vx0_bricks, __global float* vy0_bricks,
//__global float* vz0_bricks, __global float* vx_bricks, __global float* vy_bricks, __global float* vz_bricks
void CAPE::ExecuteGPUDivergence()
{
	divergenceKernel->SetArgument(5, m_bricks_buffer);
	divergenceKernel->SetArgument(6, m0_bricks_buffer);
	divergenceKernel->SetArgument(7, vx0_bricks_buffer);
	divergenceKernel->SetArgument(8, vy0_bricks_buffer);
	divergenceKernel->SetArgument(9, vz0_bricks_buffer);
	divergenceKernel->SetArgument(10, vx_bricks_buffer);
	divergenceKernel->SetArgument(11, vy_bricks_buffer);
	divergenceKernel->SetArgument(12, vz_bricks_buffer);
	divergenceKernel->SetArgument(13, p0_bricks_buffer);
	cl_event divergenceDone;
	prevtime = timer.elapsed();
	divergenceKernel->Run(bricks_allocated * CAPE_BRICKSIZE, 0, 0, &divergenceEvent);
}

void CAPE::ExecuteGPUPressure()
{
	pressureSolverKernel->SetArgument(0, grid_buffer);
	pressureSolverKernel->SetArgument(1, brick_x_buffer);
	pressureSolverKernel->SetArgument(2, brick_y_buffer);
	pressureSolverKernel->SetArgument(3, brick_z_buffer);
	pressureSolverKernel->SetArgument(4, brick_static_buffer);
	pressureSolverKernel->SetArgument(5, m_bricks_buffer);
	pressureSolverKernel->SetArgument(6, m0_bricks_buffer);
	pressureSolverKernel->SetArgument(7, p_bricks_buffer);
	pressureSolverKernel->SetArgument(8, p0_bricks_buffer);
	cl_mem gm = world->GetGridMap();
	pressureSolverKernel->SetArgument(9, &gm);
	pressureSolverKernel->SetArgument(10, world->GetBrickBuffer());

	for (int i = 0; i < PRESSURE_ITERATIONS; i++)
	{
		pressureSolverKernel->Run(bricks_allocated * CAPE_BRICKSIZE, 0, 0, &pressureSolverEvent);

		//Swap buffers
		cl_mem tmp = p_bricks_buffer->deviceBuffer;
		p_bricks_buffer->deviceBuffer = p0_bricks_buffer->deviceBuffer;
		p0_bricks_buffer->deviceBuffer = tmp;

		pressureSolverKernel->SetArgument(7, p_bricks_buffer);
		pressureSolverKernel->SetArgument(8, p0_bricks_buffer);
	}
}

void CAPE::ExecuteGPUPressureGradient()
{
	pressureGradientKernel->SetArgument(0, grid_buffer);
	pressureGradientKernel->SetArgument(1, brick_x_buffer);
	pressureGradientKernel->SetArgument(2, brick_y_buffer);
	pressureGradientKernel->SetArgument(3, brick_z_buffer);
	pressureGradientKernel->SetArgument(4, brick_static_buffer);
	pressureGradientKernel->SetArgument(5, m_bricks_buffer);
	pressureGradientKernel->SetArgument(6, m0_bricks_buffer);
	pressureGradientKernel->SetArgument(7, vx0_bricks_buffer);
	pressureGradientKernel->SetArgument(8, vy0_bricks_buffer);
	pressureGradientKernel->SetArgument(9, vz0_bricks_buffer);
	pressureGradientKernel->SetArgument(10, vx_bricks_buffer);
	pressureGradientKernel->SetArgument(11, vy_bricks_buffer);
	pressureGradientKernel->SetArgument(12, vz_bricks_buffer);
	pressureGradientKernel->SetArgument(13, p0_bricks_buffer);
	cl_mem gm = world->GetGridMap();
	pressureGradientKernel->SetArgument(14, &gm);
	pressureGradientKernel->SetArgument(15, world->GetBrickBuffer());

	prevtime = timer.elapsed();
	pressureGradientKernel->Run(bricks_allocated * CAPE_BRICKSIZE, 0, 0, &pressureGradientEvent);

	//Swap velocity
	cl_mem tmp = vx_bricks_buffer->deviceBuffer;
	vx_bricks_buffer->deviceBuffer = vx0_bricks_buffer->deviceBuffer;
	vx0_bricks_buffer->deviceBuffer = tmp;

	tmp = vy_bricks_buffer->deviceBuffer;
	vy_bricks_buffer->deviceBuffer = vy0_bricks_buffer->deviceBuffer;
	vy0_bricks_buffer->deviceBuffer = tmp;

	tmp = vz_bricks_buffer->deviceBuffer;
	vz_bricks_buffer->deviceBuffer = vz0_bricks_buffer->deviceBuffer;
	vz0_bricks_buffer->deviceBuffer = tmp;
}

void CAPE::ExecuteGPUMaterialAdvection()
{
	materialAdvectionKernel->Run(bricks_allocated * CAPE_BRICKSIZE, 0, &brickUpdateEvent, &materialAdvectionEvent);

	//Swap buffers
	cl_mem tmp = m_bricks_buffer->deviceBuffer;
	m_bricks_buffer->deviceBuffer = m0_bricks_buffer->deviceBuffer;
	m0_bricks_buffer->deviceBuffer = tmp;

	tmp = vx_bricks_buffer->deviceBuffer;
	vx_bricks_buffer->deviceBuffer = vx0_bricks_buffer->deviceBuffer;
	vx0_bricks_buffer->deviceBuffer = tmp;

	tmp = vy_bricks_buffer->deviceBuffer;
	vy_bricks_buffer->deviceBuffer = vy0_bricks_buffer->deviceBuffer;
	vy0_bricks_buffer->deviceBuffer = tmp;

	tmp = vz_bricks_buffer->deviceBuffer;
	vz_bricks_buffer->deviceBuffer = vz0_bricks_buffer->deviceBuffer;
	vz0_bricks_buffer->deviceBuffer = tmp;
}

void CAPE::ExecuteGPUVelocityAdvection()
{

	velocityAdvectionKernel->SetArgument(6, m_bricks_buffer);
	velocityAdvectionKernel->SetArgument(7, m0_bricks_buffer);
	velocityAdvectionKernel->SetArgument(8, vx0_bricks_buffer);
	velocityAdvectionKernel->SetArgument(9, vy0_bricks_buffer);
	velocityAdvectionKernel->SetArgument(10, vz0_bricks_buffer);
	velocityAdvectionKernel->SetArgument(11, vx_bricks_buffer);
	velocityAdvectionKernel->SetArgument(12, vy_bricks_buffer);
	velocityAdvectionKernel->SetArgument(13, vz_bricks_buffer);

	cl_mem gm = world->GetGridMap();
	velocityAdvectionKernel->SetArgument(15, &gm);
	velocityAdvectionKernel->SetArgument(16, world->GetBrickBuffer());

	prevtime = timer.elapsed();
	velocityAdvectionKernel->Run(bricks_allocated * CAPE_BRICKSIZE, 0, 0, &velocityAdvectionEvent);
}

void CAPE::ExecuteGPUBrickUpdate()
{
	brickUpdateKernel->SetArgument(0, grid_buffer);
	brickUpdateKernel->SetArgument(1, brick_x_buffer);
	brickUpdateKernel->SetArgument(2, brick_y_buffer);
	brickUpdateKernel->SetArgument(3, brick_z_buffer);
	brickUpdateKernel->SetArgument(4, brick_a_buffer);
	brickUpdateKernel->SetArgument(5, brick_oa_buffer);
	brickUpdateKernel->SetArgument(6, brick_jobs_buffer);
	brickUpdateKernel->SetArgument(7, brick_m_buffer);
	brickUpdateKernel->SetArgument(8, brick_static_buffer);
	brickUpdateKernel->SetArgument(9, m0_bricks_buffer);
	brickUpdateKernel->SetArgument(10, vx0_bricks_buffer);
	brickUpdateKernel->SetArgument(11, vy0_bricks_buffer);
	brickUpdateKernel->SetArgument(12, vz0_bricks_buffer);
	brickUpdateKernel->Run(bricks_allocated * CAPE_BRICKSIZE, 0, 0, &brickUpdateEvent);
}

void CAPE::ExecuteGPUWorldSetter()
{
	worldSetKernel->SetArgument(0, grid_buffer);
	worldSetKernel->SetArgument(1, brick_x_buffer);
	worldSetKernel->SetArgument(2, brick_y_buffer);
	worldSetKernel->SetArgument(3, brick_z_buffer);
	worldSetKernel->SetArgument(4, brick_static_buffer);
	worldSetKernel->SetArgument(5, m0_bricks_buffer);
	cl_mem gm = world->GetGridMap();
	worldSetKernel->SetArgument(6, &gm);
	worldSetKernel->SetArgument(7, world->GetBrickBuffer());
	worldSetKernel->SetArgument(8, world->GetZeroesBuffer());
	worldSetKernel->Run(bricks_allocated * CAPE_BRICKSIZE, 0, 0, &worldSetEvent);
}

//Returns kernel execution time in milliseconds
double KernelExecutionTime(cl_event& e)
{
	cl_ulong time_start;
	cl_ulong time_end;
	clGetEventProfilingInfo(e, CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
	clGetEventProfilingInfo(e, CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);
	double nanoSeconds = time_end - time_start;
	return nanoSeconds / 1000000.0;
}

//Performs a single update - also maintains timers and calls print function
void CAPE::Tick(float deltaTime)
{
	clFinish(capeKernel->GetQueue());
	if (updates == 0) //init, move to initialise
		CreateBuffers(false);
	else
		otherTime += timer2.elapsed();

	for (int i = 0; i < SIMS_PER_RENDER; i++)
	{
		updates++;
		if (updates > SIMS_PER_RENDER)//I have no clue, but there is some sort of race condition if we try to run the simulation on the first visit
		{
			timer.reset();

			prevtime = timer.elapsed();
			UpdateBricks();
			brickTime += timer.elapsed() - prevtime;

			//Update buffers
			if (reallocated)
				CreateBuffers(true);
			brick_a_buffer->CopyToDevice(false);
			brick_oa_buffer->CopyToDevice(false);
			brick_jobs_buffer->CopyToDevice(false);
			SetKernelArguments();

			//Step 0: Update bricks
			ExecuteGPUBrickUpdate();

			//Step 1: Let all cells advect
			ExecuteGPUMaterialAdvection();

			//data in buffer is finalised, retrieve for future brickupdate (we can do this out of order.. and asynchronously.
			brick_m_buffer->CopyFromDevice(false);

			//Step 2: Let all cells transition their velocity and apply global acceleration
			ExecuteGPUVelocityAdvection();

			//Step 3: Determine expected divergence with new velocities
			ExecuteGPUDivergence();

			//Step 4: Solve for a pressure that will prevent divergence
			ExecuteGPUPressure();

			//Step 5: Apply pressure gradient to correct velocity
			ExecuteGPUPressureGradient();

			//Step 6: Convert sim data to voxels in render data
			ExecuteGPUWorldSetter();

			clFinish(capeKernel->GetQueue());

			for (int i = 0; i < bricks_allocated; i++) //Jobs done
				brick_jobs[i] = UINT32_MAX;
		}
	}

	if (updates > SIMS_PER_RENDER)
		PrintState();
}

//Retrieves amount of incoming momentum on the X-axis from tangential moving material across the y and z axis
float CAPE::IncomingMomentumX(uint x, uint y, uint z)
{
	float yl = GetData(x, y - 1, z, vx_bricks);
	float yr = GetData(x, y + 1, z, vx_bricks);
	//incoming for x across y axis
	float xym = -min(0.0f, GetData(x - 1, y + 1, z, vy_bricks)) * AL * GetData(x - 1, y + 1, z, m_bricks) * max(0.0f, yr) * AL //left above
		+ max(0.0f, GetData(x - 1, y, z, vy_bricks)) * AL * GetData(x - 1, y - 1, z, m_bricks) * max(0.0f, yl) * AL//left below
		+ -min(0.0f, GetData(x, y + 1, z, vy_bricks)) * AL * GetData(x, y + 1, z, m_bricks) * min(0.0f, yr) * AL //right above
		+ max(0.0f, GetData(x, y, z, vy_bricks)) * AL * GetData(x, y - 1, z, m_bricks) * min(0.0f, yl) * AL; //right below
	float zl = GetData(x, y, z - 1, vx_bricks);
	float zr = GetData(x, y, z + 1, vx_bricks);
	//incoming for x across z axis
	float xzm = -min(0.0f, GetData(x - 1, y, z + 1, vz_bricks)) * AL * GetData(x - 1, y, z + 1, m_bricks) * max(0.0f, zr) * AL //left above
		+ max(0.0f, GetData(x - 1, y, z, vz_bricks)) * AL * GetData(x - 1, y, z - 1, m_bricks) * max(0.0f, zl) * AL//left below
		+ -min(0.0f, GetData(x, y, z + 1, vz_bricks)) * AL * GetData(x, y, z + 1, m_bricks) * min(0.0f, zr) * AL //right above
		+ max(0.0f, GetData(x, y, z, vz_bricks)) * AL * GetData(x, y, z - 1, m_bricks) * min(0.0f, zl) * AL; //right below
	return xym + xzm;
}

//Retrieves amount of incoming momentum on the Y-axis from tangential moving material across the x and z axis
float CAPE::IncomingMomentumY(uint x, uint y, uint z)
{
	float xl = GetData(x - 1, y, z, vy_bricks);
	float xr = GetData(x + 1, y, z, vy_bricks);
	//incoming for x across y axis
	float yxm = -min(0.0f, GetData(x + 1, y - 1, z, vx_bricks)) * AL * GetData(x + 1, y - 1, z, m_bricks) * max(0.0f, xr) * AL //left above
		+ max(0.0f, GetData(x, y - 1, z, vx_bricks)) * AL * GetData(x - 1, y - 1, z, m_bricks) * max(0.0f, xl) * AL//left below
		+ -min(0.0f, GetData(x + 1, y, z, vx_bricks)) * AL * GetData(x + 1, y, z, m_bricks) * min(0.0f, xr) * AL //right above
		+ max(0.0f, GetData(x, y, z, vx_bricks)) * AL * GetData(x - 1, y, z, m_bricks) * min(0.0f, xl) * AL; //right below

	float zl = GetData(x, y, z - 1, vy_bricks);
	float zr = GetData(x, y, z + 1, vy_bricks);
	//incoming for x across z axis
	float yzm = -min(0.0f, GetData(x, y - 1, z + 1, vz_bricks)) * AL * GetData(x, y - 1, z + 1, m_bricks) * max(0.0f, zr) * AL //left above
		+ max(0.0f, GetData(x, y - 1, z, vz_bricks)) * AL * GetData(x, y - 1, z - 1, m_bricks) * max(0.0f, zl) * AL//left below
		+ -min(0.0f, GetData(x, y, z + 1, vz_bricks)) * AL * GetData(x, y, z + 1, m_bricks) * min(0.0f, zr) * AL //right above
		+ max(0.0f, GetData(x, y, z, vz_bricks)) * AL * GetData(x, y, z - 1, m_bricks) * min(0.0f, zl) * AL; //right below
	return yxm + yzm;
}

//Retrieves amount of incoming momentum on the Z-axis from tangential moving material across the x and y axis
float CAPE::IncomingMomentumZ(uint x, uint y, uint z)
{
	float xl = GetData(x - 1, y, z, vz_bricks);
	float xr = GetData(x + 1, y, z, vz_bricks);
	//incoming for x across y axis
	float zxm = -min(0.0f, GetData(x + 1, y, z - 1, vx_bricks)) * AL * GetData(x + 1, y, z - 1, m_bricks) * max(0.0f, xr) * AL //left above
		+ max(0.0f, GetData(x, y, z - 1, vx_bricks)) * AL * GetData(x - 1, y, z - 1, m_bricks) * max(0.0f, xl) * AL//left below
		+ -min(0.0f, GetData(x + 1, y, z, vx_bricks)) * AL * GetData(x + 1, y, z, m_bricks) * min(0.0f, xr) * AL //right above
		+ max(0.0f, GetData(x, y, z, vx_bricks)) * AL * GetData(x - 1, y, z, m_bricks) * min(0.0f, xl) * AL; //right below

	float yl = GetData(x, y - 1, z, vz_bricks);
	float yr = GetData(x, y + 1, z, vz_bricks);
	//incoming for x across z axis
	float zym = -min(0.0f, GetData(x, y + 1, z - 1, vy_bricks)) * AL * GetData(x, y + 1, z - 1, m_bricks) * max(0.0f, yr) * AL //left above
		+ max(0.0f, GetData(x, y, z - 1, vy_bricks)) * AL * GetData(x, y - 1, z - 1, m_bricks) * max(0.0f, yl) * AL//left below
		+ -min(0.0f, GetData(x, y + 1, z, vy_bricks)) * AL * GetData(x, y + 1, z, m_bricks) * min(0.0f, yr) * AL //right above
		+ max(0.0f, GetData(x, y, z, vy_bricks)) * AL * GetData(x, y - 1, z, m_bricks) * min(0.0f, yl) * AL; //right below
	return zxm + zym;
}

//Converts simulation data to actual voxels placed in the world
void CAPE::ConvertToVoxels()
{
	RunOverAllBricks(this, &CAPE::SetColorForCell, timeStep);
}

//Remove simulation data from world
void CAPE::ClearVoxels()
{
	RunOverAllBricks(this, &CAPE::ClearColorForCell, timeStep);
}


//Converts rgb value to the 4-4-4 12 bit rgb format used 
uint LerpToRGB(float r, float g, float b)
{
	r = clamp(r, 0.0f, 1.0f);
	g = clamp(g, 0.0f, 1.0f);
	b = clamp(b, 0.0f, 1.0f);
	return (uint)(15 * b) + ((uint)(15 * g) << 4) + ((uint)(15 * r) << 8);
}

//Converts cell state to a color and sets it into the world
void CAPE::ClearColorForCell(uint x, uint y, uint z, float timeStep)
{
	uint color = 0;
	float mass = GetData(x, y, z, m0_bricks);
	if (mass > 0.00001)
		world->Set(x - CAPE_BRICKDIM, y - CAPE_BRICKDIM, z - CAPE_BRICKDIM, 0);
}

//Converts cell state to a color and sets it into the world
void CAPE::SetColorForCell(uint x, uint y, uint z, float timeStep)
{
	uint color = 0;
	float mass = GetData(x,y,z, m0_bricks);
	if (mass > 0.0001)
	{
		//interpolate mass between 0 and - max mass for the 15 different colors.
		float mass = GetData(x, y, z, m0_bricks);
		float diffFromEmpty = 1 - clamp(1 - mass, 0.0f, 1.0f);
		if (diffFromEmpty < MINRENDERMASS)
			world->Set(x - CAPE_BRICKDIM, y - CAPE_BRICKDIM, z - CAPE_BRICKDIM, 0);
		else
		{
			float b = 1 - diffFromEmpty / 2.0f + 0.5f;
			float g = diffFromEmpty > 0.5f ? 1 - (diffFromEmpty - 0.5f) : 1;
			world->Set(x - CAPE_BRICKDIM, y - CAPE_BRICKDIM, z - CAPE_BRICKDIM, LerpToRGB(0, g, b));
		}	
	}
}

//Print simulation debug information
void CAPE::PrintState()
{
	system("cls");
	cout << "Timestep: " << timeStep << endl;
	cout << "brick update time " << brickTime / updates * 1000 << endl;
	cout << "device buffer update time " << buffertime / updates * 1000 << endl;
	printf("BrickUpdate OpenCl Execution time is: %0.3f milliseconds \n", KernelExecutionTime(brickUpdateEvent));
	printf("MaterialAdvection OpenCl Execution time is: %0.3f milliseconds \n", KernelExecutionTime(materialAdvectionEvent));
	printf("VelocityAdvection OpenCl Execution time is: %0.3f milliseconds \n", KernelExecutionTime(velocityAdvectionEvent));
	printf("Divergence OpenCl Execution time is: %0.3f milliseconds \n", KernelExecutionTime(divergenceEvent));
	printf("PressureSolver OpenCl Execution time is: %0.3f milliseconds \n", KernelExecutionTime(pressureSolverEvent) * 8); //assume latest pressure solve kernel is average, fine over time
	printf("PressureGradient OpenCl Execution time is: %0.3f milliseconds \n", KernelExecutionTime(pressureGradientEvent));
	printf("WorldSet OpenCl Execution time is: %0.3f milliseconds \n", KernelExecutionTime(worldSetEvent));
	float pipelinekerneltime = KernelExecutionTime(materialAdvectionEvent) + KernelExecutionTime(velocityAdvectionEvent) + KernelExecutionTime(divergenceEvent)
		+ KernelExecutionTime(pressureSolverEvent) * 8 + KernelExecutionTime(pressureGradientEvent) + KernelExecutionTime(worldSetEvent);
	cout << "pipelinekerneltime: " << pipelinekerneltime << " milliseconds" << endl;

	float potpipelinekernelUpdates = 1000 / (pipelinekerneltime);

	cout << "pipeline time " << pipelineTime / updates << endl;
	cout << "tick time: " << (otherTime + simulationTime) / updates * 1000 << endl;
	cout << "fps: " << 1 / ((otherTime + simulationTime) / updates) << endl;
	cout << "Update: " << updates << endl;
	cout << "Alive Bricks: " << bricks_alive << endl;
	cout << "Allocated Bricks: " << bricks_allocated << endl;
	cout << "Active Cells: " << bricks_alive * CAPE_BRICKSIZE << endl;
	cout << "Cells updated per second: " << bricks_alive * CAPE_BRICKSIZE * (1 / (simulationTime / updates)) << endl;
	cout << "potential pipeline kernel cells updated per second: " << bricks_alive * CAPE_BRICKSIZE * potpipelinekernelUpdates << endl;
	cout << "world set time: " << worldsetTime / (updates / SIMS_PER_RENDER) * 1000 << endl;
}

//Old CPU Sim Rules

//Update all bricks, if they are allocated, alive and non-static
//When using concurrency, assign a brick to each thread
void CAPE::RunOverAllBricks(CAPE* cape, void (CAPE::* func)(uint, uint, uint, float), float timeStep)
{
	const size_t ac = bricks_allocated;
#if USECONCURRENCY == 1
	Concurrency::parallel_for(size_t(0), ac, [&](size_t i)
		{
			if (brick_x[i] != UINT32_MAX && !brick_static[i])
			{
				uint xo = brick_x[i] * CAPE_BRICKDIM;
				uint yo = brick_y[i] * CAPE_BRICKDIM;
				uint zo = brick_z[i] * CAPE_BRICKDIM;
				for (int x = 0; x < CAPE_BRICKDIM; x++)
					for (int y = 0; y < CAPE_BRICKDIM; y++)
						for (int z = 0; z < CAPE_BRICKDIM; z++)
							(cape->*func)(xo + x, yo + y, zo + z, timeStep);
			}
		});
#else
	for (int i = 0; i < ac; i++)
	{
		if (brick_x[i] != UINT32_MAX && !brick_static[i])
		{
			uint xo = brick_x[i] * CAPE_BRICKDIM;
			uint yo = brick_y[i] * CAPE_BRICKDIM;
			uint zo = brick_z[i] * CAPE_BRICKDIM;
			for (int x = 0; x < CAPE_BRICKDIM; x++)
				for (int y = 0; y < CAPE_BRICKDIM; y++)
					for (int z = 0; z < CAPE_BRICKDIM; z++)
						(cape->*func)(xo + x, yo + y, zo + z, timeStep);
		}
	}
#endif
}

//Determine a suitable pressure value that prevents divergence, given 
//the cells divergence, and a best estimation of the pressure of neighbours
//Repeatedly calling this improves this estimation
void CAPE::SolvePressure(uint x, uint y, uint z, float timeStep)
{
	float div = GetData(x, y, z, m_bricks);
	float nbp = 0;
	float k = 0;
	nbp += max(0.0f, GetData(x + 1, y, z, p0_bricks));
	nbp += max(0.0f, GetData(x - 1, y, z, p0_bricks));
	nbp += max(0.0f, GetData(x, y + 1, z, p0_bricks));
	nbp += max(0.0f, GetData(x, y - 1, z, p0_bricks));
	nbp += max(0.0f, GetData(x, y, z + 1, p0_bricks));
	nbp += max(0.0f, GetData(x, y, z - 1, p0_bricks));

	k += (IsCellStatic(x + 1, y, z) == 0);
	k += (IsCellStatic(x - 1, y, z) == 0);
	k += (IsCellStatic(x, y + 1, z) == 0);
	k += (IsCellStatic(x, y - 1, z) == 0);
	k += (IsCellStatic(x, y, z + 1) == 0);
	k += (IsCellStatic(x, y, z - 1) == 0);

	float np = (div + nbp) / k;
	SetData(x, y, z, np, p_bricks);
}

//Determine pressure gradient and apply to correct velocities
void CAPE::PressureGradient(uint x, uint y, uint z, float timeStep)
{
	//solid cells do not experience pressure 
	if (IsCellStatic(x, y, z))
		return;

	float m = GetData(x, y, z, m0_bricks);
	float cellPressure = max(0.0f, GetData(x, y, z, p0_bricks));

	//No pressure gradient to solid cells, therefore simply set to equal pressure
	float ovxl = GetData(x, y, z, vx0_bricks);
	float ovyl = GetData(x, y, z, vy0_bricks);
	float ovzl = GetData(x, y, z, vz0_bricks);
	float leftPressure = max(0.0f, IsCellStatic(x - 1, y, z) ? cellPressure : GetData(x - 1, y, z, p0_bricks));
	float botPressure = max(0.0f, IsCellStatic(x, y - 1, z) ? cellPressure : GetData(x, y - 1, z, p0_bricks));
	float backPressure = max(0.0f, IsCellStatic(x, y, z - 1) ? cellPressure : GetData(x, y, z - 1, p0_bricks));

	//Determine acceleration from pressure gradient
	//Since we are also finalising the velocity for the current update here, we want to dampen
	//the original remaining velocity and clamp the velocity so we do not exceed the maximum here
	float dampen = 1 - VELOCITY_DAMPENING * timeStep;
	float vxa = leftPressure - cellPressure;
	float vya = botPressure - cellPressure;
	float vza = backPressure - cellPressure;
	float vxl = ovxl * dampen + vxa;
	float vyl = ovyl * dampen + vya;
	float vzl = ovzl * dampen + vza;
	SetData(x, y, z, clamp(vxl, -MAXV, MAXV), vx_bricks);
	SetData(x, y, z, clamp(vyl, -MAXV, MAXV), vy_bricks);
	SetData(x, y, z, clamp(vzl, -MAXV, MAXV), vz_bricks);
}

//Determine change in mass given current velocities
//Assume cells cannot go past maximum density (e.g 1), therefore clamp to 1.0
float CAPE::MaterialChange(uint x, uint y, uint z, float vxl, float vxr, float vyl, float vyr, float vzl, float vzr)
{
	float mat = GetData(x, y, z, m0_bricks);
	float matxl = GetData(x - 1, y, z, m0_bricks);
	float matyl = GetData(x, y - 1, z, m0_bricks);
	float matzl = GetData(x, y, z - 1, m0_bricks);
	float matxr = GetData(x + 1, y, z, m0_bricks);
	float matyr = GetData(x, y + 1, z, m0_bricks);
	float matzr = GetData(x, y, z + 1, m0_bricks);

	float dxl = vxl < 0 ? vxl * min(1.0f, mat) : vxl * min(1.0f, matxl);
	float dyl = vyl < 0 ? vyl * min(1.0f, mat) : vyl * min(1.0f, matyl);
	float dzl = vzl < 0 ? vzl * min(1.0f, mat) : vzl * min(1.0f, matzl);
	float dxr = vxr < 0 ? -vxr * min(1.0f, matxr) : -vxr * min(1.0f, mat);
	float dyr = vyr < 0 ? -vyr * min(1.0f, matyr) : -vyr * min(1.0f, mat);
	float dzr = vzr < 0 ? -vzr * min(1.0f, matzr) : -vzr * min(1.0f, mat);
	return dxl + dyl + dzl + dxr + dyr + dzr;
}

//Update material content of cell given velocities
void CAPE::MaterialAdvection(uint x, uint y, uint z, float timeStep)
{
	float vxl = GetData(x, y, z, vx0_bricks) * AL;
	float vyl = GetData(x, y, z, vy0_bricks) * AL;
	float vzl = GetData(x, y, z, vz0_bricks) * AL;
	float vxr = GetData(x + 1, y, z, vx0_bricks) * AL;
	float vyr = GetData(x, y + 1, z, vy0_bricks) * AL;
	float vzr = GetData(x, y, z + 1, vz0_bricks) * AL;

	//node focussed advection
	float mat = GetData(x, y, z, m0_bricks);
	float matxl = GetData(x - 1, y, z, m0_bricks);
	float matyl = GetData(x, y - 1, z, m0_bricks);
	float matzl = GetData(x, y, z - 1, m0_bricks);
	float matxr = GetData(x + 1, y, z, m0_bricks);
	float matyr = GetData(x, y + 1, z, m0_bricks);
	float matzr = GetData(x, y, z + 1, m0_bricks);

	//Material movement in every direction == momentum
	float dxl = vxl < 0 ? vxl * min(1.0f, mat) : vxl * min(1.0f, matxl);
	float dyl = vyl < 0 ? vyl * min(1.0f, mat) : vyl * min(1.0f, matyl);
	float dzl = vzl < 0 ? vzl * min(1.0f, mat) : vzl * min(1.0f, matzl);
	float dxr = vxr < 0 ? -vxr * min(1.0f, matxr) : -vxr * min(1.0f, mat);
	float dyr = vyr < 0 ? -vyr * min(1.0f, matyr) : -vyr * min(1.0f, mat);
	float dzr = vzr < 0 ? -vzr * min(1.0f, matzr) : -vzr * min(1.0f, mat);

	//update incoming matter into current cell and outgoing matter neighbour cells
	float nmat = mat + MaterialChange(x, y, z, vxl, vxr, vyl, vyr, vzl, vzr);
#if DEBUG_MODE == 1
	if (nmat < 0) cout << "negative mass! (timestep too high and or not using save MAXV?)" << endl;
#endif
	//Apply some evaporation and set new material
	nmat = max(0.0f, nmat - EVAPORATION * timeStep);
	SetData(x, y, z, nmat, m_bricks);

	//Remember outgoing momentum (useful for velocity update)
	//Store temporarily in unused velocity buffers
	float om = min(0.0f, dxl) + min(0.0f, dyl) + min(0.0f, dzl)
		+ min(0.0f, dxr) + min(0.0f, dyr) + min(0.0f, dzr);
	if (vxl < 0) SetData(x, y, z, vxl * om, vx_bricks);
	if (vyl < 0) SetData(x, y, z, vyl * om, vy_bricks);
	if (vzl < 0) SetData(x, y, z, vzl * om, vz_bricks);
	if (vxr > 0) SetData(x + 1, y, z, vxr * om, vx_bricks);
	if (vyr > 0) SetData(x, y + 1, z, vyr * om, vy_bricks);
	if (vzr > 0) SetData(x, y, z + 1, vzr * om, vz_bricks);

	//Add to total brick mass (for brick_updating step, assumes only 1 thread per brick active)
	uint brickIdx = GetBrickIDX(x, y, z);
	uint brick_addr = grid[brickIdx];
	brick_m[brick_addr] += nmat;
}

//Is the cell active in the simulation:
//False if it is a voxel existing in the world
//but no material is defined for it, or the brick
//that it is in is marked static.
bool CAPE::IsCellStatic(uint x, uint y, uint z)
{
	// calculate brick location in top-level grid
	uint brickIDX = GetBrickIDX(x, y, z);
	uint brick_addr = grid[brickIDX];
	if (brick_addr == UINT32_MAX)
		return true;

	bool static_brick = brick_static[brick_addr];
	float mass = GetData(x, y, z, m0_bricks);
	return (mass == 0 && world->Get(x - CAPE_BRICKDIM, y - CAPE_BRICKDIM, z - CAPE_BRICKDIM) != 0) || static_brick;
}

//Updates velocity by applying global acceleration and self-advection
//As with material advection, both sides of material are retrieved, even
//If not required (e.g. they get multiplied by 0), however since we are likely
//to visit a cell that does require this material, within the same brick, this
//should hopefully remain in cache, this allows us to keep a cell centric functional approach
void CAPE::CellVelocityUpdate(uint x, uint y, uint z, float timeStep)
{
	//Load velocity and mass data so we can process momentum
	float vxl = GetData(x - 1, y, z, vx_bricks) * AL;
	float vyl = GetData(x, y - 1, z, vy_bricks) * AL;
	float vzl = GetData(x, y, z - 1, vz_bricks) * AL;
	float vxc = GetData(x, y, z, vx_bricks) * AL;
	float vyc = GetData(x, y, z, vy_bricks) * AL;
	float vzc = GetData(x, y, z, vz_bricks) * AL;
	float vxr = GetData(x + 1, y, z, vx_bricks) * AL;
	float vyr = GetData(x, y + 1, z, vy_bricks) * AL;
	float vzr = GetData(x, y, z + 1, vz_bricks) * AL;

	float omat = GetData(x, y, z, m_bricks);
	float omatxl2 = GetData(x - 2, y, z, m_bricks);
	float omatyl2 = GetData(x, y - 2, z, m_bricks);
	float omatzl2 = GetData(x, y, z - 2, m_bricks);
	float omatxl = GetData(x - 1, y, z, m_bricks);
	float omatyl = GetData(x, y - 1, z, m_bricks);
	float omatzl = GetData(x, y, z - 1, m_bricks);
	float omatxr = GetData(x + 1, y, z, m_bricks);
	float omatyr = GetData(x, y + 1, z, m_bricks);
	float omatzr = GetData(x, y, z + 1, m_bricks);

	float mat = GetData(x, y, z, m0_bricks);
	float matxl = GetData(x - 1, y, z, m0_bricks);
	float matyl = GetData(x, y - 1, z, m0_bricks);
	float matzl = GetData(x, y, z - 1, m0_bricks);
	float matxr = GetData(x + 1, y, z, m0_bricks);
	float matyr = GetData(x, y + 1, z, m0_bricks);
	float matzr = GetData(x, y, z + 1, m0_bricks);

	float avxl = max(0.0f, vxl);
	float avxr = min(0.0f, vxr);
	float avxcl = max(0.0f, vxc);
	float avxcr = min(0.0f, vxc);
	float avyl = max(0.0f, vyl);
	float avyr = min(0.0f, vyr);
	float avycl = max(0.0f, vyc);
	float avycr = min(0.0f, vyc);
	float avzl = max(0.0f, vzl);
	float avzr = min(0.0f, vzr);
	float avzcl = max(0.0f, vzc);
	float avzcr = min(0.0f, vzc);

	//incoming colinear momentum across given axis                           
	float imx = (avxl * avxl * omatxl2) - (avxr * avxr * omatxr);
	float imy = (avyl * avyl * omatyl2) - (avyr * avyr * omatyr);
	float imz = (avzl * avzl * omatzl2) - (avzr * avzr * omatzr);

	//old momentum across axis
	float mx0 = avxcl * omatxl + avxcr * omat;
	float my0 = avycl * omatyl + avycr * omat;
	float mz0 = avzcl * omatzl + avzcr * omat;

	//outgoing momentum, stored temporarily during material advection step
	float omx = GetData(x, y, z, vx0_bricks);
	float omy = GetData(x, y, z, vy0_bricks);
	float omz = GetData(x, y, z, vz0_bricks);

	//Remaining momentum
	float rmx = mx0 + omx;
	float rmy = my0 + omy;
	float rmz = mz0 + omz;

	//New momentum = remainin momentum + incoming colinear momentum + incoming tangential momentum
	float vcx = rmx + imx + 0;// IncomingMomentumX(x, y, z);
	float vcy = rmy + imy + 0;// IncomingMomentumY(x, y, z);
	float vcz = rmz + imz + 0;// IncomingMomentumZ(x, y, z);

	//mass of source cell
	float massx = vcx > 0 ? matxl : mat;
	float massy = vcy > 0 ? matyl : mat;
	float massz = vcz > 0 ? matzl : mat;

	//New velocity *
	float nvcx = massx == 0 ? 0 : vcx / massx * INVAL;
	float nvcy = massy == 0 ? 0 : vcy / massy * INVAL;
	float nvcz = massz == 0 ? 0 : vcz / massz * INVAL;

	//Add global acceleration
	for (int f = 0; f < ga.size(); f++)
	{
		nvcx += timeStep * ga[f].x;
		nvcy += timeStep * ga[f].y;
		nvcz += timeStep * ga[f].z;
	}

	//Check if not involving static cell, then finalise a valid velocity by clamping and dampening
	bool staticc = IsCellStatic(x, y, z);
	bool staticx = staticc || IsCellStatic(x - 1, y, z);
	bool staticy = staticc || IsCellStatic(x, y - 1, z);
	bool staticz = staticc || IsCellStatic(x, y, z - 1);

	float nvx = !staticx * nvcx;
	float nvy = !staticy * nvcy;
	float nvz = !staticz * nvcz;
	SetData(x, y, z, nvx, vx0_bricks);
	SetData(x, y, z, nvy, vy0_bricks);
	SetData(x, y, z, nvz, vz0_bricks);
}

//Calculate and set divergence of a cell we need to solve pressure for
void CAPE::CellDivergenceUpdate(uint x, uint y, uint z, float timeStep)
{
	SetData(x, y, z, 0, m_bricks);
	SetData(x, y, z, 0, p0_bricks);
	float vxl = GetData(x, y, z, vx0_bricks) * AL;
	float vyl = GetData(x, y, z, vy0_bricks) * AL;
	float vzl = GetData(x, y, z, vz0_bricks) * AL;
	float vxr = GetData(x + 1, y, z, vx0_bricks) * AL;
	float vyr = GetData(x, y + 1, z, vy0_bricks) * AL;
	float vzr = GetData(x, y, z + 1, vz0_bricks) * AL;

	float mc = MaterialChange(x, y, z, vxl, vxr, vyl, vyr, vzl, vzr); //divergence = amount of mass over 1 and under 0
	float div = (GetData(x, y, z, m0_bricks) + mc - 1.0f) * 2.0f; //new mass, premultiply by 2 to prevent gradient division 
	//Store in temporarily unused material buffer
	SetData(x, y, z, div, m_bricks);
}


//Some functions that aggregate simulation data for debugging
float CAPE::TotalDivergence()
{
	float divergence = 0;
	for (int i = 0; i < bricks_allocated; i++)
	{
		if (brick_x[i] != UINT32_MAX)
		{
			uint xo = brick_x[i] * CAPE_BRICKDIM;
			uint yo = brick_y[i] * CAPE_BRICKDIM;
			uint zo = brick_z[i] * CAPE_BRICKDIM;
			for (int x = 0; x < CAPE_BRICKDIM; x++)
				for (int y = 0; y < CAPE_BRICKDIM; y++)
					for (int z = 0; z < CAPE_BRICKDIM; z++)
						divergence += max(0.0f, GetData(xo + x, yo + y, zo + z, m0_bricks) - 1);
		}
	}
	return divergence;
}

float CAPE::TotalMass()
{
	float totalMass = 0;
	for (int i = 0; i < bricks_allocated; i++)
	{
		if (brick_x[i] != UINT32_MAX)
		{
			uint xo = brick_x[i] * CAPE_BRICKDIM;
			uint yo = brick_y[i] * CAPE_BRICKDIM;
			uint zo = brick_z[i] * CAPE_BRICKDIM;
			for (int x = 0; x < CAPE_BRICKDIM; x++)
				for (int y = 0; y < CAPE_BRICKDIM; y++)
					for (int z = 0; z < CAPE_BRICKDIM; z++)
						totalMass += GetData(xo + x, yo + y, zo + z, m0_bricks);
		}
	}
	return totalMass;
}
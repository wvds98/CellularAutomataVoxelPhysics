#include "template/common.h"
#include "cl/trace.cl"
#include "cl/tools.cl"

//Required on 6900XT
#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable

//Sloppy copy for now, make this a proper include
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
#define VELOCITY_DAMPENING 0.01f 
#define AL 0.5f //Advection limiter, prevents oscillations
#define INVAL (1.0f/AL) 
#define MINBRICKS 50 //Minimum amount of bricks to reserve in memory
#define UINT32_MAX 0xffffffff

#define EVAPORATION 0.0f //amount of mass to remove per cell per second (helps clean up low mass cells that are not visible.)
#define MINRENDERMASS 0.001f //Dont render voxels below this much mass
#define GRAVITYENABLED 1
#define CELLSIZE 0.10f //In meters, essentially multiplier for gravity
#define timestep 0.01f

//0.333 is maximum 100% save speed, use up to 1.0f for less clamping and faster flowing water (less viscous), but requires that an appropriate
//timestep is selected by the user, or simulation may blow up if local velocity becomes too high and negative mass is created
#define MAXV 0.33f

__kernel void run(){}

inline uint GetBrickIDX(const uint x, const uint y, const uint z)
{
	// calculate brick location in top-level grid
	const uint bx = x / CAPE_BRICKDIM;
	const uint by = y / CAPE_BRICKDIM;
	const uint bz = z / CAPE_BRICKDIM;
	if (bx >= CAPE_GRIDWIDTH || by >= CAPE_GRIDHEIGHT || bz >= CAPE_GRIDDEPTH) return UINT32_MAX;
	const uint brickIdx = (bx)+(by) * (CAPE_GRIDWIDTH) + (bz) * (CAPE_GRIDHEIGHT) * (CAPE_GRIDDEPTH);
	return brickIdx;
}

inline float GetData(const uint x, const uint y, const uint z, __global uint* grid, __global float* data)
{
	// obtain brick reference from top-level grid if brick does not exist, return "default" value
	const uint bID = grid[GetBrickIDX(x, y, z)]; 
	if (bID == UINT32_MAX) return 0; 
	__global float* d = data + bID * CAPE_BRICKSIZE;// [bID] ;
	const uint lx = x & (CAPE_BRICKDIM - 1), ly = y & (CAPE_BRICKDIM - 1), lz = z & (CAPE_BRICKDIM - 1);
	return d[lx + ly * CAPE_BRICKDIM + lz * CAPE_BRICKDIM * CAPE_BRICKDIM];
}

inline void SetData(const uint x, const uint y, const uint z, float v, __global uint* grid, __global float* data)
{
    const uint bidx = GetBrickIDX(x, y, z);
	const uint bID = grid[bidx];
	if (bID == UINT32_MAX) return;
	global float* d = data + bID * CAPE_BRICKSIZE;
	const uint lx = x & (CAPE_BRICKDIM - 1), ly = y & (CAPE_BRICKDIM - 1), lz = z & (CAPE_BRICKDIM - 1);
	uint cellIdx = lx + ly * CAPE_BRICKDIM + lz * CAPE_BRICKDIM * CAPE_BRICKDIM;
	d[cellIdx] = v;
}

//Determine change in mass given current velocities
//Assume cells cannot go past maximum density (e.g 1), therefore clamp to 1.0
float MaterialChange(const uint x, const uint y, const uint z, const float vxl, const float vxr, const float vyl, const float vyr, const float vzl, const float vzr,
                     __global uint* grid, __global float* m0_bricks)
{
	float mat   = GetData(x, y, z, grid, m0_bricks);
	float matxl = GetData(x - 1, y, z, grid, m0_bricks); 
	float matyl = GetData(x, y - 1, z, grid, m0_bricks); 
	float matzl = GetData(x, y, z - 1, grid, m0_bricks);
	float matxr = GetData(x + 1, y, z, grid, m0_bricks); 
	float matyr = GetData(x, y + 1, z, grid, m0_bricks); 
	float matzr = GetData(x, y, z + 1, grid, m0_bricks);

	float dxl = vxl < 0 ? vxl * min(1.0f, mat) : vxl * min(1.0f, matxl);
	float dyl = vyl < 0 ? vyl * min(1.0f, mat) : vyl * min(1.0f, matyl);
	float dzl = vzl < 0 ? vzl * min(1.0f, mat) : vzl * min(1.0f, matzl);
	float dxr = vxr < 0 ? -vxr * min(1.0f, matxr) : -vxr * min(1.0f, mat);
	float dyr = vyr < 0 ? -vyr * min(1.0f, matyr) : -vyr * min(1.0f, mat);
	float dzr = vzr < 0 ? -vzr * min(1.0f, matzr) : -vzr * min(1.0f, mat);
	return dxl + dyl + dzl + dxr + dyr + dzr;
}

__kernel void calculateDivergence(__global uint* grid, __global uint* bx, __global uint* by, __global uint* bz, __global bool* brick_static, __global float* m_bricks, __global float* m0_bricks, __global float* vx0_bricks, __global float* vy0_bricks,
                  __global float* vz0_bricks, __global float* vx_bricks, __global float* vy_bricks, __global float* vz_bricks, __global float* p0_bricks)
{
	int id = get_global_id(0);

	//Get brick location
	const uint bidx = id / CAPE_BRICKSIZE;
	const uint bxv = bx[bidx];
	const uint byv = by[bidx];
	const uint bzv = bz[bidx];
	if(brick_static[bidx]) return;

	//inverse mapping to x y and z inside brick
	const uint ib = id - bidx * CAPE_BRICKSIZE;
	const uint lx = ib & (CAPE_BRICKDIM - 1); 
	const uint ly = (ib / CAPE_BRICKDIM) & (CAPE_BRICKDIM - 1);
	const uint lz = ib / (CAPE_BRICKDIM * CAPE_BRICKDIM);

	//global x y z
	const uint x = bxv * CAPE_BRICKDIM + lx;
	const uint y = byv * CAPE_BRICKDIM + ly;
	const uint z = bzv * CAPE_BRICKDIM + lz;

	SetData(x, y, z, 0, grid, m_bricks);
	SetData(x, y, z, 0, grid, p0_bricks);
	float vxl = GetData(x, y, z, grid, vx0_bricks) * AL;
	float vyl = GetData(x, y, z, grid, vy0_bricks) * AL;
	float vzl = GetData(x, y, z, grid, vz0_bricks) * AL;
	float vxr = GetData(x + 1, y, z, grid, vx0_bricks) * AL;
	float vyr = GetData(x, y + 1, z, grid, vy0_bricks) * AL;
	float vzr = GetData(x, y, z + 1, grid, vz0_bricks) * AL;

	float mc = MaterialChange(x, y, z, vxl, vxr, vyl, vyr, vzl, vzr, grid, m0_bricks); //divergence = amount of mass over 1 and under 0
	float div = (GetData(x, y, z, grid, m0_bricks) + mc - 1.0f) * 2.0f; //new mass, premultiply by 2 to prevent gradient division 
	//Store in temporarily unused material buffer
	SetData(x, y, z, div, grid, m_bricks);
}

#define BMSK (BRICKDIM - 1)
#define BRICKDIM 8
#define BRICKSIZE (BRICKDIM * BRICKDIM * BRICKDIM)
#define BDIM2 (BRICKDIM * BRICKDIM)
#define GRIDWIDTH 128
#define GRIDHEIGHT 128
#define GRIDDEPTH 128
uint worldGet(const uint x, const uint y, const uint z, __read_only image3d_t grid, __global unsigned short* brick)
{
	// calculate brick location in top-level grid
	const uint bx = (x / BRICKDIM) & (GRIDWIDTH - 1);
	const uint by = (y / BRICKDIM) & (GRIDHEIGHT - 1);
	const uint bz = (z / BRICKDIM) & (GRIDDEPTH - 1);
	const uint cellIdx = bx + bz * GRIDWIDTH + by * GRIDWIDTH * GRIDDEPTH;
	const uint g = read_imageui(grid, (int4)(bx, bz, by, 0)).x;
	if ((g & 1) == 0 /* this is currently an empty cell */) return 0;
	// calculate the position of the voxel inside the brick
	const uint lx = x & (BRICKDIM - 1), ly = y & (BRICKDIM - 1), lz = z & (BRICKDIM - 1);
	return brick[(g >> 1) * BRICKSIZE + lx + ly * BRICKDIM + lz * BRICKDIM * BRICKDIM];
}

void worldSet(const uint x, const uint y, const uint z, __write_only image3d_t grid, __global unsigned short* brick, __global uint* zeroes,  uint v)
{
	// calculate brick location in top-level grid
	const uint bx = (x / BRICKDIM) & (GRIDWIDTH - 1);
	const uint by = (y / BRICKDIM) & (GRIDHEIGHT - 1);
	const uint bz = (z / BRICKDIM) & (GRIDDEPTH - 1);
	const uint cellIdx = bx + bz * GRIDWIDTH + by * GRIDWIDTH * GRIDDEPTH; //index of brick

	// calculate the position of the voxel inside the brick
	const uint lx = x & (BRICKDIM - 1), ly = y & (BRICKDIM - 1), lz = z & (BRICKDIM - 1);
	uint voxelIdx = cellIdx * BRICKSIZE + lx + ly * BRICKDIM + lz * BRICKDIM * BRICKDIM;
	const uint cv = brick[voxelIdx]; //original value

	atomic_add(zeroes + cellIdx, (cv != 0 && v == 0) - (cv == 0 && v != 0)); //Apply the change in zeroes
	write_imageui(grid, (int4)(bx, bz, by, 0), (uint4)((cellIdx << 1) | (zeroes[cellIdx] < BRICKSIZE), 0,0,0)); //Update zero status
	brick[voxelIdx] = v;
}

bool IsCellStatic(const uint x, const uint y, const uint z, __global uint* grid, __global char* brick_static, __global float* m0_bricks, __read_only image3d_t worldGrid, __global const unsigned short* worldBricks)
{
	// calculate brick location in top-level grid
	uint brickIDX = GetBrickIDX(x, y, z);
	uint brick_addr = grid[brickIDX];
	if (brick_addr == UINT32_MAX)
		return true;

	bool static_brick = brick_static[brick_addr];
	float m = GetData(x, y, z, grid, m0_bricks);

	//Get world voxel
	bool ws = worldGet(x - CAPE_BRICKDIM,y - CAPE_BRICKDIM,z - CAPE_BRICKDIM, worldGrid, worldBricks) != 0;

	return (m == 0 && ws) || static_brick; 
}

__kernel void solvePressure(__global uint* grid, __global uint* bx, __global uint* by, __global uint* bz, __global bool* brick_static, __global float* m_bricks, __global float* m0_bricks, __global float* p_bricks, __global float* p0_bricks, __read_only image3d_t worldGrid, __global const unsigned short* worldBricks)
{
	int id = get_global_id(0);

	//Get brick location
	const uint bidx = id / CAPE_BRICKSIZE;
	const uint bxv = bx[bidx];
	const uint byv = by[bidx];
	const uint bzv = bz[bidx];
	if(brick_static[bidx]) return;

	//inverse mapping to x y and z inside brick
	const uint ib = id - bidx * CAPE_BRICKSIZE;
	const uint lx = ib & (CAPE_BRICKDIM - 1); 
	const uint ly = (ib / CAPE_BRICKDIM) & (CAPE_BRICKDIM - 1);
	const uint lz = ib / (CAPE_BRICKDIM * CAPE_BRICKDIM);

	//global x y z
	const uint x = bxv * CAPE_BRICKDIM + lx;
	const uint y = byv * CAPE_BRICKDIM + ly;
	const uint z = bzv * CAPE_BRICKDIM + lz;

	float div = GetData(x, y, z, grid, m_bricks);
	float nbp = 0;
	float k = 0;
	nbp += max(0.0f, GetData(x + 1, y, z, grid, p0_bricks));
	nbp += max(0.0f, GetData(x - 1, y, z, grid, p0_bricks));
	nbp += max(0.0f, GetData(x, y + 1, z, grid, p0_bricks));
	nbp += max(0.0f, GetData(x, y - 1, z, grid, p0_bricks));
	nbp += max(0.0f, GetData(x, y, z + 1, grid, p0_bricks));
	nbp += max(0.0f, GetData(x, y, z - 1, grid, p0_bricks));

	k += (IsCellStatic(x + 1, y, z, grid, brick_static, m0_bricks, worldGrid, worldBricks) == 0);
	k += (IsCellStatic(x - 1, y, z, grid, brick_static, m0_bricks, worldGrid, worldBricks) == 0);
	k += (IsCellStatic(x, y + 1, z, grid, brick_static, m0_bricks, worldGrid, worldBricks) == 0);
	k += (IsCellStatic(x, y - 1, z, grid, brick_static, m0_bricks, worldGrid, worldBricks) == 0);
	k += (IsCellStatic(x, y, z + 1, grid, brick_static, m0_bricks, worldGrid, worldBricks) == 0);
	k += (IsCellStatic(x, y, z - 1, grid, brick_static, m0_bricks, worldGrid, worldBricks) == 0);

	float np = (div + nbp) / k;
	SetData(x, y, z, np, grid, p_bricks);
}

__kernel void pressureGradient(__global uint* grid, __global uint* bx, __global uint* by, __global uint* bz, __global bool* brick_static, __global float* m_bricks, __global float* m0_bricks, __global float* vx0_bricks, __global float* vy0_bricks,
                  __global float* vz0_bricks, __global float* vx_bricks, __global float* vy_bricks, __global float* vz_bricks, __global float* p0_bricks, __read_only image3d_t worldGrid, __global const unsigned short* worldBricks)
{
	int id = get_global_id(0);

    //Get brick location
	const uint bidx = id / CAPE_BRICKSIZE;
	const uint bxv = bx[bidx];
	const uint byv = by[bidx];
	const uint bzv = bz[bidx];
	if(brick_static[bidx]) return;

	//inverse mapping to x y and z inside brick
	const uint ib = id - bidx * CAPE_BRICKSIZE;
	const uint lx = ib & (CAPE_BRICKDIM - 1); 
	const uint ly = (ib / CAPE_BRICKDIM) & (CAPE_BRICKDIM - 1);
	const uint lz = ib / (CAPE_BRICKDIM * CAPE_BRICKDIM);

	//global x y z
	const uint x = bxv * CAPE_BRICKDIM + lx;
	const uint y = byv * CAPE_BRICKDIM + ly;
	const uint z = bzv * CAPE_BRICKDIM + lz;

	//solid cells do not experience pressure 
	if (IsCellStatic(x, y, z, grid, brick_static, m0_bricks, worldGrid, worldBricks))
		return;

	float cellPressure = max(0.0f, GetData(x,y,z, grid, p0_bricks));

	//No pressure gradient to solid cells, therefore simply set to equal pressure
	float ovxl = GetData(x, y, z, grid, vx0_bricks);
	float ovyl = GetData(x, y, z, grid, vy0_bricks);
	float ovzl = GetData(x, y, z, grid, vz0_bricks);
	float leftPressure = max(0.0f, IsCellStatic(x - 1, y, z, grid, brick_static, m0_bricks, worldGrid, worldBricks) ? cellPressure : GetData(x - 1, y, z, grid, p0_bricks));
	float botPressure = max(0.0f, IsCellStatic(x, y - 1, z, grid, brick_static, m0_bricks, worldGrid, worldBricks) ? cellPressure : GetData(x, y - 1,z, grid, p0_bricks));
	float backPressure = max(0.0f, IsCellStatic(x, y, z - 1, grid, brick_static, m0_bricks, worldGrid, worldBricks) ? cellPressure : GetData(x, y, z - 1, grid, p0_bricks));

	//Determine acceleration from pressure gradient
	//Since we are also finalising the velocity for the current update here, we want to dampen
	//the original remaining velocity and clamp the velocity so we do not exceed the maximum here
	float dampen = 1 - VELOCITY_DAMPENING * timestep;
	float vxa = leftPressure - cellPressure;
	float vya = botPressure - cellPressure;
	float vza = backPressure - cellPressure;
	float vxl = ovxl * dampen + vxa;
	float vyl = ovyl * dampen + vya;
	float vzl = ovzl * dampen + vza;
	SetData(x, y, z, clamp(vxl, -MAXV, MAXV), grid, vx_bricks);
	SetData(x, y, z, clamp(vyl, -MAXV, MAXV), grid, vy_bricks);
	SetData(x, y, z, clamp(vzl, -MAXV, MAXV), grid, vz_bricks);
}

__kernel void materialAdvection(__global uint* grid, __global uint* bx, __global uint* by, __global uint* bz, __global bool* brick_static, __global float* brick_m, __global float* m_bricks, __global float* m0_bricks, __global float* vx0_bricks, __global float* vy0_bricks,
                                __global float* vz0_bricks, __global float* vx_bricks, __global float* vy_bricks, __global float* vz_bricks, __global float* p0_bricks, __read_only image3d_t worldGrid, __global const unsigned short* worldBricks)
{
	int id = get_global_id(0);

    //Get brick location
	const uint bidx = id / CAPE_BRICKSIZE;
	const uint bxv = bx[bidx];
	const uint byv = by[bidx];
	const uint bzv = bz[bidx];
	if(brick_static[bidx]) return;

	//inverse mapping to x y and z inside brick
	const uint ib = id - bidx * CAPE_BRICKSIZE;
	const uint lx = ib & (CAPE_BRICKDIM - 1); 
	const uint ly = (ib / CAPE_BRICKDIM) & (CAPE_BRICKDIM - 1);
	const uint lz = ib / (CAPE_BRICKDIM * CAPE_BRICKDIM);

	//global x y z
	const uint x = bxv * CAPE_BRICKDIM + lx;
	const uint y = byv * CAPE_BRICKDIM + ly;
	const uint z = bzv * CAPE_BRICKDIM + lz;

	float vxl = GetData(x, y, z, grid, vx0_bricks) * AL;
	float vyl = GetData(x, y, z, grid, vy0_bricks) * AL;
	float vzl = GetData(x, y, z, grid, vz0_bricks) * AL;
	float vxr = GetData(x + 1, y, z, grid, vx0_bricks) * AL;
	float vyr = GetData(x, y + 1, z, grid, vy0_bricks) * AL;
	float vzr = GetData(x, y, z + 1, grid, vz0_bricks) * AL;

	//node focussed advection
	float mat = GetData(x, y, z, grid, m0_bricks);
	float matxl = GetData(x - 1, y, z, grid, m0_bricks); 
	float matyl = GetData(x, y - 1, z, grid, m0_bricks); 
	float matzl = GetData(x, y, z - 1, grid, m0_bricks); 
	float matxr = GetData(x + 1, y, z, grid, m0_bricks); 
	float matyr = GetData(x, y + 1, z, grid, m0_bricks); 
	float matzr = GetData(x, y, z + 1, grid, m0_bricks); 

	//Material movement in every direction == momentum
	float dxl = vxl < 0 ? vxl * min(1.0f, mat) : vxl * min(1.0f, matxl);
	float dyl = vyl < 0 ? vyl * min(1.0f, mat) : vyl * min(1.0f, matyl);
	float dzl = vzl < 0 ? vzl * min(1.0f, mat) : vzl * min(1.0f, matzl);
	float dxr = vxr < 0 ? -vxr * min(1.0f, matxr) : -vxr * min(1.0f, mat);
	float dyr = vyr < 0 ? -vyr * min(1.0f, matyr) : -vyr * min(1.0f, mat);
	float dzr = vzr < 0 ? -vzr * min(1.0f, matzr) : -vzr * min(1.0f, mat);

	//update incoming matter into current cell and outgoing matter neighbour cells
	float nmat = mat + MaterialChange(x,y,z, vxl, vxr, vyl, vyr, vzl, vzr, grid, m0_bricks);
	//Apply some evaporation and set new material

	nmat = max(0.0f, nmat - EVAPORATION * timestep);
	SetData(x, y, z, nmat, grid, m_bricks);

	//Remember outgoing momentum (useful for velocity update)
	//Store temporarily in unused velocity buffers
	float om = min(0.0f, dxl) + min(0.0f, dyl) + min(0.0f, dzl)
		     + min(0.0f, dxr) + min(0.0f, dyr) + min(0.0f, dzr);
	if (vxl < 0) SetData(x, y, z, vxl * om, grid, vx_bricks);
	if (vyl < 0) SetData(x, y, z, vyl * om, grid, vy_bricks);
	if (vzl < 0) SetData(x, y, z, vzl * om, grid, vz_bricks);
	if (vxr > 0) SetData(x + 1, y, z, vxr * om, grid, vx_bricks);
	if (vyr > 0) SetData(x, y + 1, z, vyr * om, grid, vy_bricks);
	if (vzr > 0) SetData(x, y, z + 1, vzr * om, grid, vz_bricks);

	//mark brick as nonempty
	uint brickIdx = GetBrickIDX(x, y, z);
	uint brick_addr = grid[brickIdx];
	if(nmat > 0.001) brick_m[brick_addr] = 1;
}

//Retrieves amount of incoming momentum on the X-axis from tangential moving material across the y and z axis
inline float IncomingMomentumX(uint x, uint y, uint z, __global uint* grid, __global float* vx_bricks, __global float* vy_bricks, __global float* vz_bricks, __global float* m_bricks)
{
	float yl = GetData(x, y - 1, z, grid, vx_bricks);
	float yr = GetData(x, y + 1, z, grid, vx_bricks);
	//incoming for x across y axis
	float xym = -min(0.0f, GetData(x - 1, y + 1, z, grid, vy_bricks)) * AL * GetData(x - 1, y + 1, z, grid, m_bricks) * max(0.0f, yr) * AL //left above
		+ max(0.0f, GetData(x - 1, y, z, grid, vy_bricks)) * AL * GetData(x - 1, y - 1, z, grid, m_bricks) * max(0.0f, yl) * AL//left below
		+ -min(0.0f, GetData(x, y + 1, z, grid, vy_bricks)) * AL * GetData(x, y + 1, z, grid, m_bricks) * min(0.0f, yr) * AL //right above
		+ max(0.0f, GetData(x, y, z, grid, vy_bricks)) * AL * GetData(x, y - 1, z, grid, m_bricks) * min(0.0f, yl) * AL; //right below
	float zl = GetData(x, y, z - 1, grid, vx_bricks);
	float zr = GetData(x, y, z + 1, grid, vx_bricks);
	//incoming for x across z axis
	float xzm = -min(0.0f, GetData(x - 1, y, z + 1, grid, vz_bricks)) * AL * GetData(x - 1, y, z + 1, grid, m_bricks) * max(0.0f, zr) * AL //left above
		+ max(0.0f, GetData(x - 1, y, z, grid, vz_bricks)) * AL * GetData(x - 1, y, z - 1, grid, m_bricks) * max(0.0f, zl) * AL//left below
		+ -min(0.0f, GetData(x, y, z + 1, grid, vz_bricks)) * AL * GetData(x, y, z + 1, grid, m_bricks) * min(0.0f, zr) * AL //right above
		+ max(0.0f, GetData(x, y, z, grid, vz_bricks)) * AL * GetData(x, y, z - 1, grid, m_bricks) * min(0.0f, zl) * AL; //right below
	return xym + xzm;
}

//Retrieves amount of incoming momentum on the Y-axis from tangential moving material across the x and z axis
float IncomingMomentumY(uint x, uint y, uint z, __global uint* grid, __global float* vx_bricks, __global float* vy_bricks, __global float* vz_bricks, __global float* m_bricks)
{
	float xl = GetData(x - 1, y, z, grid, vy_bricks);
	float xr = GetData(x + 1, y, z, grid, vy_bricks);
	//incoming for x across y axis
	float yxm = -min(0.0f, GetData(x + 1, y - 1, z, grid, vx_bricks)) * AL * GetData(x + 1, y - 1, z, grid, m_bricks) * max(0.0f, xr) * AL //left above
		+ max(0.0f, GetData(x, y - 1, z, grid, vx_bricks)) * AL * GetData(x - 1, y - 1, z, grid, m_bricks) * max(0.0f, xl) * AL//left below
		+ -min(0.0f, GetData(x + 1, y, z, grid, vx_bricks)) * AL * GetData(x + 1, y, z, grid, m_bricks) * min(0.0f, xr) * AL //right above
		+ max(0.0f, GetData(x, y, z, grid, vx_bricks)) * AL * GetData(x - 1, y, z, grid, m_bricks) * min(0.0f, xl) * AL; //right below
	float zl = GetData(x, y, z - 1, grid, vy_bricks);
	float zr = GetData(x, y, z + 1, grid, vy_bricks);
	//incoming for x across z axis
	float yzm = -min(0.0f, GetData(x, y - 1, z + 1, grid, vz_bricks)) * AL * GetData(x, y - 1, z + 1, grid, m_bricks) * max(0.0f, zr) * AL //left above
		+ max(0.0f, GetData(x, y - 1, z, grid, vz_bricks)) * AL * GetData(x, y - 1, z - 1, grid, m_bricks) * max(0.0f, zl) * AL//left below
		+ -min(0.0f, GetData(x, y, z + 1, grid, vz_bricks)) * AL * GetData(x, y, z + 1, grid, m_bricks) * min(0.0f, zr) * AL //right above
		+ max(0.0f, GetData(x, y, z, grid, vz_bricks)) * AL * GetData(x, y, z - 1, grid, m_bricks) * min(0.0f, zl) * AL; //right below
	return yxm + yzm;
}

//Retrieves amount of incoming momentum on the Z-axis from tangential moving material across the x and y axis
float IncomingMomentumZ(uint x, uint y, uint z, __global uint* grid, __global float* vx_bricks, __global float* vy_bricks, __global float* vz_bricks, __global float* m_bricks)
{
	float xl = GetData(x - 1, y, z, grid, vz_bricks);
	float xr = GetData(x + 1, y, z, grid, vz_bricks);
	//incoming for x across y axis
	float zxm = -min(0.0f, GetData(x + 1, y, z - 1, grid, vx_bricks)) * AL * GetData(x + 1, y, z - 1, grid, m_bricks) * max(0.0f, xr) * AL //left above
		+ max(0.0f, GetData(x, y, z - 1, grid, vx_bricks)) * AL * GetData(x - 1, y, z - 1, grid, m_bricks) * max(0.0f, xl) * AL//left below
		+ -min(0.0f, GetData(x + 1, y, z, grid, vx_bricks)) * AL * GetData(x + 1, y, z, grid, m_bricks) * min(0.0f, xr) * AL //right above
		+ max(0.0f, GetData(x, y, z, grid, vx_bricks)) * AL * GetData(x - 1, y, z, grid, m_bricks) * min(0.0f, xl) * AL; //right below
	float yl = GetData(x, y - 1, z, grid, vz_bricks);
	float yr = GetData(x, y + 1, z, grid, vz_bricks);
	//incoming for x across z axis
	float zym = -min(0.0f, GetData(x, y + 1, z - 1, grid, vy_bricks)) * AL * GetData(x, y + 1, z - 1, grid, m_bricks) * max(0.0f, yr) * AL //left above
		+ max(0.0f, GetData(x, y, z - 1, grid, vy_bricks)) * AL * GetData(x, y - 1, z - 1, grid, m_bricks) * max(0.0f, yl) * AL//left below
		+ -min(0.0f, GetData(x, y + 1, z, grid, vy_bricks)) * AL * GetData(x, y + 1, z, grid, m_bricks) * min(0.0f, yr) * AL //right above
		+ max(0.0f, GetData(x, y, z, grid, vy_bricks)) * AL * GetData(x, y - 1, z, grid, m_bricks) * min(0.0f, yl) * AL; //right below
	return zxm + zym;
}

__kernel void velocityAdvection(__global uint* grid, __global uint* bx, __global uint* by, __global uint* bz, __global bool* brick_static, __global float* brick_m, __global float* m_bricks, __global float* m0_bricks, __global float* vx0_bricks, __global float* vy0_bricks,
                                __global float* vz0_bricks, __global float* vx_bricks, __global float* vy_bricks, __global float* vz_bricks, __global float* p0_bricks, __read_only image3d_t worldGrid, __global const unsigned short* worldBricks)
{
	int id = get_global_id(0);
    //Get brick location
	const uint bidx = id / CAPE_BRICKSIZE;
	const uint bxv = bx[bidx];
	const uint byv = by[bidx];
	const uint bzv = bz[bidx];
	if(brick_static[bidx]) return;

	//inverse mapping to x y and z inside brick
	const uint ib = id - bidx * CAPE_BRICKSIZE;
	const uint lx = ib & (CAPE_BRICKDIM - 1); 
	const uint ly = (ib / CAPE_BRICKDIM) & (CAPE_BRICKDIM - 1);
	const uint lz = ib / (CAPE_BRICKDIM * CAPE_BRICKDIM);

	//global x y z
	const uint x = bxv * CAPE_BRICKDIM + lx;
	const uint y = byv * CAPE_BRICKDIM + ly;
	const uint z = bzv * CAPE_BRICKDIM + lz;

	//Load velocity and mass data so we can process momentum
	float vxl = GetData(x - 1, y, z, grid, vx_bricks) * AL;
	float vyl = GetData(x, y - 1, z, grid, vy_bricks) * AL;
	float vzl = GetData(x, y, z - 1, grid, vz_bricks) * AL;
	float vxc = GetData(x, y, z, grid, vx_bricks) * AL;
	float vyc = GetData(x, y, z, grid, vy_bricks) * AL;
	float vzc = GetData(x, y, z, grid, vz_bricks) * AL;
	float vxr = GetData(x + 1, y, z, grid, vx_bricks) * AL;
	float vyr = GetData(x, y + 1, z, grid, vy_bricks) * AL;
	float vzr = GetData(x, y, z + 1, grid, vz_bricks) * AL;

	float omat = GetData(x, y, z, grid, m_bricks);
	float omatxl2 = GetData(x - 2, y, z, grid, m_bricks);
	float omatyl2 = GetData(x, y - 2, z, grid, m_bricks);
	float omatzl2 = GetData(x, y, z - 2, grid, m_bricks);
	float omatxl = GetData(x - 1, y, z, grid, m_bricks);
	float omatyl = GetData(x, y - 1, z, grid, m_bricks);
	float omatzl = GetData(x, y, z - 1, grid, m_bricks);
	float omatxr = GetData(x + 1, y, z, grid, m_bricks);
	float omatyr = GetData(x, y + 1, z, grid, m_bricks);
	float omatzr = GetData(x, y, z + 1, grid, m_bricks);

	float mat = GetData(x, y, z, grid, m0_bricks);
	float matxl = GetData(x - 1, y, z, grid, m0_bricks);
	float matyl = GetData(x, y - 1, z, grid, m0_bricks);
	float matzl = GetData(x, y, z - 1, grid, m0_bricks);
	float matxr = GetData(x + 1, y, z, grid, m0_bricks);
	float matyr = GetData(x, y + 1, z, grid, m0_bricks);
	float matzr = GetData(x, y, z + 1, grid, m0_bricks);

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
	float omx = GetData(x,y,z, grid, vx0_bricks);
	float omy = GetData(x, y, z, grid, vy0_bricks);
	float omz = GetData(x, y, z, grid, vz0_bricks);

	//Remaining momentum
	float rmx = mx0 + omx;
	float rmy = my0 + omy;
	float rmz = mz0 + omz;

	//New momentum = remainin momentum + incoming colinear momentum + incoming tangential momentum
	float vcx = rmx + imx + IncomingMomentumX(x, y, z, grid, vx_bricks, vy_bricks, vz_bricks, m_bricks);
	float vcy = rmy + imy + IncomingMomentumY(x, y, z, grid, vx_bricks, vy_bricks, vz_bricks, m_bricks);
	float vcz = rmz + imz + IncomingMomentumZ(x, y, z, grid, vx_bricks, vy_bricks, vz_bricks, m_bricks);

	//mass of source cell
	float massx = vcx > 0 ? matxl : mat;
	float massy = vcy > 0 ? matyl : mat;
	float massz = vcz > 0 ? matzl : mat;

	//New velocity *
	float nvcx = massx == 0 ? 0 : vcx / massx * INVAL;
	float nvcy = massy == 0 ? 0 : vcy / massy * INVAL;
	float nvcz = massz == 0 ? 0 : vcz / massz * INVAL;

	//Add global acceleration
	//for (int f = 0; f < ga.size(); f++)
	//{
		nvcx += timestep * 0.0f;//ga[f].x;
		nvcy += timestep * -9.81f * CELLSIZE;// ga[f].y;
		nvcz += timestep * 0.0f; //ga[f].z;
	//}

	//Check if not involving static cell, then finalise a valid velocity by clamping and dampening
	bool staticc = IsCellStatic(x,y,z, grid, brick_static, m0_bricks, worldGrid, worldBricks);
	bool staticx = staticc || IsCellStatic(x - 1, y, z, grid, brick_static, m0_bricks, worldGrid, worldBricks);
	bool staticy = staticc || IsCellStatic(x, y - 1, z, grid, brick_static, m0_bricks, worldGrid, worldBricks);
	bool staticz = staticc || IsCellStatic(x, y, z - 1, grid, brick_static, m0_bricks, worldGrid, worldBricks);

	float nvx = !staticx * nvcx;
	float nvy = !staticy * nvcy;
	float nvz = !staticz * nvcz;
	SetData(x, y, z, nvx, grid, vx0_bricks);
	SetData(x, y, z, nvy, grid, vy0_bricks);
	SetData(x, y, z, nvz, grid, vz0_bricks);
}

		//float4 f = (float4)(1.0f, 2.0f, 3.0f, 4.0f);
		//uchar4 uc = (uchar4)(0xFA, 0xFB, 0xFC, 0xFD);

		//printf("f4 = %2.2v4hlf\n", f);
		//printf("uc = %#v4hhx\n", uc);

//For now just run for every allocated brick, using brick jobs can be reduced to just running for bricks that actually need updating
//However cost is very low, so not a priority, also runs jobs for every cell if individual brick data needs to be set.
__kernel void brickUpdate(__global uint* grid, __global uint* brick_x, __global uint* brick_y, __global uint* brick_z, __global uint* brick_a, __global uint* brick_oa, __global uint* brick_jobs, __global float* brick_m, __global char* brick_static, __global float* m0_bricks, __global float* vx0_bricks, __global float* vy0_bricks, __global float* vz0_bricks)
{
	int id = get_global_id(0);
	uint brick_id = id / CAPE_BRICKSIZE;
	char s = brick_static[brick_id];
	uint j = brick_jobs[brick_id];
	uint brick_addr = brick_a[brick_id];

	//Reset default values for the upcoming update (for any brick)
	brick_m[brick_id] = 0;

	//inverse mapping to bx by and bz
	const uint bx = brick_addr % (CAPE_GRIDWIDTH); 
	const uint by = (brick_addr / CAPE_GRIDWIDTH) % (CAPE_GRIDHEIGHT);
	const uint bz = brick_addr / (CAPE_GRIDWIDTH * CAPE_GRIDHEIGHT);

	if(j == UINT32_MAX) return; //nothing to do here.

	brick_x[brick_id] = bx;
	brick_y[brick_id] = by;
	brick_z[brick_id] = bz;

	if(j == 0) //brick has been freed, make it static, clear values
	{ 
		grid[brick_addr] = UINT32_MAX;
		m0_bricks[id] = 0; //Clear brick cells of remaining small mass
		brick_x[brick_id] = UINT32_MAX;
		brick_static[brick_id] = true;
	} 
	else if(j == 1) //freed brick is recycled, assume it may have not been previously freed (freed and recycled in same update)
	{	
		uint oa = brick_oa[brick_id];
		grid[oa] = UINT32_MAX;
		m0_bricks[id] = 0; //Clear brick cells of remaining small mass
		brick_static[brick_id] = (bx == 0 || bx >= CAPE_GRIDWIDTH || by == 0 || by >= CAPE_GRIDHEIGHT || bz == 0 || bz >= CAPE_GRIDDEPTH);
		grid[brick_addr] = brick_id;
		vx0_bricks[id] = 0;
		vy0_bricks[id] = 0;
		vz0_bricks[id] = 0;

		//Reset any values in possibly existing neighbouring bricks
		//Work this into adjusted enqueue  --inverse neighbour mapping a good performance booster?
		uint bxo = bx * CAPE_BRICKDIM;
		uint byo = by * CAPE_BRICKDIM;
		uint bzo = bz * CAPE_BRICKDIM;
		for(int y = 0; y < CAPE_BRICKDIM; y++)
			for(int z = 0; z < CAPE_BRICKDIM; z++)
			{
				SetData((bx+1) * CAPE_BRICKDIM, byo + y, bzo + z, 0, grid, vx0_bricks);
				SetData((bx+1) * CAPE_BRICKDIM, byo + y, bzo + z, 0, grid, vy0_bricks);
				SetData((bx+1) * CAPE_BRICKDIM, byo + y, bzo + z, 0, grid, vz0_bricks);
			}

		for(int x = 0; x < CAPE_BRICKDIM; x++)
			for(int z = 0; z < CAPE_BRICKDIM; z++)
			{
				SetData(bxo + x,(by+1) * CAPE_BRICKDIM, bzo + z, 0, grid, vx0_bricks);
				SetData(bxo + x,(by+1) * CAPE_BRICKDIM, bzo + z, 0, grid, vy0_bricks);
				SetData(bxo + x,(by+1) * CAPE_BRICKDIM, bzo + z, 0, grid, vz0_bricks);
			}

		for(int x = 0; x < CAPE_BRICKDIM; x++)
			for(int y = 0; y < CAPE_BRICKDIM; y++)
			{
				SetData(bxo + x, byo + y, (bz+1) * CAPE_BRICKDIM, 0, grid, vx0_bricks);
				SetData(bxo + x, byo + y, (bz+1) * CAPE_BRICKDIM, 0, grid, vy0_bricks);
				SetData(bxo + x, byo + y, (bz+1) * CAPE_BRICKDIM, 0, grid, vz0_bricks);
			}
	} 
	if (j == 2)//a new brick allocated and initialised at end of list
	{	
		brick_static[brick_id] = (bx == 0 || bx >= CAPE_GRIDWIDTH || by == 0 || by >= CAPE_GRIDHEIGHT || bz == 0 || bz >= CAPE_GRIDDEPTH);
		grid[brick_addr] = brick_id;
	}
}

//Converts rgb value to the 4-4-4 12 bit rgb format used 
inline uint LerpToRGB(float r, float g, float b)
{
	r = clamp(r, 0.0f, 1.0f);
	g = clamp(g, 0.0f, 1.0f);
	b = clamp(b, 0.0f, 1.0f);
	return (uint)(15 * b) + ((uint)(15 * g) << 4) + ((uint)(15 * r) << 8);
}

__kernel void worldSetter(__global uint* grid, __global uint* bx, __global uint* by, __global uint* bz, __global bool* brick_static, __global float* m0_bricks, __write_only image3d_t worldGrid, __global const unsigned short* worldBricks, __global uint* zeroes)
{
	int id = get_global_id(0);

	//Get brick location
	const uint bidx = id / CAPE_BRICKSIZE;
	const uint bxv = bx[bidx];
	const uint byv = by[bidx];
	const uint bzv = bz[bidx];
	if(brick_static[bidx]) return;

	//inverse mapping to x y and z inside brick
	const uint ib = id - bidx * CAPE_BRICKSIZE;
	const uint lx = ib & (CAPE_BRICKDIM - 1); 
	const uint ly = (ib / CAPE_BRICKDIM) & (CAPE_BRICKDIM - 1);
	const uint lz = ib / (CAPE_BRICKDIM * CAPE_BRICKDIM);

	//global x y z
	const uint x = bxv * CAPE_BRICKDIM + lx;
	const uint y = byv * CAPE_BRICKDIM + ly;
	const uint z = bzv * CAPE_BRICKDIM + lz;

	uint color = 0;
	float mass = GetData(x,y,z, grid, m0_bricks);
	if (mass > 0.0001f)
	{
		//interpolate mass between 0 and - max mass for the 15 different colors.
		float diffFromEmpty = 1 - clamp(1 - mass, 0.0f, 1.0f);
		if (diffFromEmpty < MINRENDERMASS)
			worldSet(x - CAPE_BRICKDIM,y - CAPE_BRICKDIM,z - CAPE_BRICKDIM, worldGrid, worldBricks, zeroes, 0);
		else
		{
			float b = 1 - diffFromEmpty / 2.0f + 0.5f;
			float g = diffFromEmpty > 0.5f ? 1 - (diffFromEmpty - 0.5f) : 1;
			worldSet(x - CAPE_BRICKDIM,y - CAPE_BRICKDIM,z - CAPE_BRICKDIM, worldGrid, worldBricks, zeroes, LerpToRGB(0, g, b));
		}	
	}
}
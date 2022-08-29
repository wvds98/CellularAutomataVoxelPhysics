#include "precomp.h"
#include "waterworld.h"

Game* CreateGame() { return new WaterWorld(); }

void WaterWorld::SetStaticBlock(uint x0, uint y0, uint z0, uint w, uint h, uint d, uint v)
{
	if (x0 + w > MAPWIDTH || y0 + h > MAPWIDTH || z0 + d > MAPWIDTH)
	{
		cout << "Block outside of grid" << endl;
		return;
	}
	for (uint x = x0; x < x0 + w; x++)
		for (uint y = y0; y < y0 + h; y++)
			for (uint z = z0; z < z0 + d; z++)
				Plot(x, y, z, v);
}

//Little dam holding water with a hole in it
void WaterWorld::InitialiseDamHoleScenario()
{
	//boundary
	SetStaticBlock(500, 500, 500, 200, 100, 1, RED);
	SetStaticBlock(700, 500, 500, 1, 100, 100, RED);
	SetStaticBlock(500, 500, 600, 200, 100, 1, RED);
	SetStaticBlock(500, 500, 500, 1, 100, 100, RED);
	SetStaticBlock(500, 500, 500, 200, 1, 100, RED);

	//Obstacle
	SetStaticBlock(540, 500, 500, 3, 100, 100, WHITE);

	SetStaticBlock(540, 500, 500, 3, 10, 10, 0); //hole 1

	SetStaticBlock(540, 500, 500, 160, 10, 1, WHITE); //tube
	SetStaticBlock(540, 500, 510, 160, 10, 1, WHITE); //tube
	SetStaticBlock(540, 500, 500, 160, 1, 10, WHITE); //tube
	SetStaticBlock(540, 510, 500, 150, 1, 10, WHITE); //tube
	SetStaticBlock(700, 500, 500, 1, 5, 10, WHITE); //tube

	SetStaticBlock(540, 550, 575, 3, 8, 8, 0); //hole 2
	SetStaticBlock(540, 501, 535, 3, 12, 30, 0); //hole 3
	SetStaticBlock(543, 513, 535, 6, 1, 30, WHITE); //hole 3

	//water
	SetMaterialBlock(500, 500, 500, 40, 80, 100, 1, false);
}

//Little dam holding water with a hole in it
void WaterWorld::InitialiseWaterLevelScenario()
{
	//boundary
	SetStaticBlock(500, 500, 500, 200, 100, 1, RED);
	SetStaticBlock(700, 500, 500, 1, 100, 100, RED);
	SetStaticBlock(500, 500, 600, 200, 100, 1, RED);
	SetStaticBlock(500, 500, 500, 1, 100, 100, RED);
	SetStaticBlock(500, 500, 500, 200, 1, 100, RED);

	//Obstacle
	SetStaticBlock(540, 510, 500, 3, 90, 100, WHITE);
	SetStaticBlock(540, 510, 500, 150, 1, 100, WHITE);

	//water
	SetMaterialBlock(500, 500, 500, 40, 80, 100, 1, false);
}

//Spawns a wall of water that will collapse
void WaterWorld::InitialiseDamBreakScenario()
{
	//boundary
	SetStaticBlock(500, 500, 500, 200, 100, 1, RED);
	SetStaticBlock(700, 500, 500, 1, 100, 100, RED);
	SetStaticBlock(500, 500, 600, 200, 100, 1, RED);
	SetStaticBlock(500, 500, 500, 1, 100, 100, RED);
	SetStaticBlock(500, 500, 500, 200, 1, 100, RED);

	//Obstacle
	SetStaticBlock(660, 500, 540, 20, 20, 20, WHITE);

	//water
	SetMaterialBlock(500, 500, 500, 40, 30, 100, 1, false);
}

//Scenario used for evaluation: Drop a block 40x40x40 water into 100x100x100 cube
void WaterWorld::InitialiseWaterBlockDropScenario()
{
	//water
	SetMaterialBlock(515, 530, 515, 30, 30, 30, 1, false);

	//boundary
	SetStaticBlock(500, 500, 500, 60, 60, 1, RED);
	SetStaticBlock(560, 500, 500, 1, 60, 60, RED);
	SetStaticBlock(500, 500, 560, 60, 60, 1, RED);
	SetStaticBlock(500, 500, 500, 1, 60, 60, RED);
	SetStaticBlock(500, 500, 500, 60, 1, 60, RED);
}

//Scenario used for evaluation: Drop a block 40x40x40 water into 100x100x100 cube
void WaterWorld::InitialiseBuildingDropScenario()
{
	//water
	SetMaterialBlock(320, 500, 510, 40, 10, 40, 1, false);

	ship = LoadSprite("assets/flyingapts.vx"), corvette = LoadSprite("assets/corvette.vx");
	StampSpriteTo(ship, 200, 250, 400);
	int c[15] = { 100, 280, 400, 500, 320, 350, 600, 480, 250, 390, 350, 530, 290, 310, 30 };
	//for (int i = 0; i < 5; i++) StampSpriteTo(corvette, c[i * 3], c[i * 3 + 1], c[i * 3 + 2]);
}

//Scenario used for evaluation: Large wall of water dropped in front of buildings
void WaterWorld::InitialiseTsunami()
{
	ship = LoadSprite("assets/flyingapts.vx"), corvette = LoadSprite("assets/corvette.vx");
	StampSpriteTo(ship, 20, 0, 400);
	int c[15] = { 100, 280, 400, 500, 320, 350, 600, 480, 250, 390, 350, 530, 290, 310, 30 };
	//for (int i = 0; i < 5; i++) StampSpriteTo(corvette, c[i * 3], c[i * 3 + 1], c[i * 3 + 2]);

	SetStaticBlock(00, 0, 400, 1024, 200, 1, WHITE);
	SetStaticBlock(00, 0, 700, 1024, 200, 1, WHITE);
	SetStaticBlock(400, 0, 400, 1, 200, 300, WHITE);
	//water
	SetMaterialBlock(20, 0, 400, 100, 200, 300, 1, false);
}

// -----------------------------------------------------------
// Initialize the application
// -----------------------------------------------------------
void WaterWorld::Init()
{
	ShowCursor( false );
	skyDomeLightScale = 6.0f;
	// default scene is a box; punch a hole in the ceiling
	Box( 256, 240, 256, 768, 260, 768, 0 );
	// add some objects

	InitCAPE(100);

	//A few scenario's to choose from
	//InitialiseDamHoleScenario();
	InitialiseWaterBlockDropScenario();
	//InitialiseDamBreakScenario();
	//InitialiseWaterLevelScenario();
	//InitialiseBuildingDropScenario();
	//InitialiseTsunami();
}

// -----------------------------------------------------------
// HandleInput: reusable input handling / free camera
// -----------------------------------------------------------
void WaterWorld::HandleInput( float deltaTime )
{
	// free cam controls
	float3 tmp( 0, 1, 0 ), right = normalize( cross( tmp, D ) ), up = cross( D, right );
	float speed = deltaTime * 0.1f;
	if (GetAsyncKeyState( 'W' )) O += speed * D; else if (GetAsyncKeyState( 'S' )) O -= speed * D;
	if (GetAsyncKeyState( 'A' )) O -= speed * right; else if (GetAsyncKeyState( 'D' )) O += speed * right;
	if (GetAsyncKeyState( 'R' )) O += speed * up; else if (GetAsyncKeyState( 'F' )) O -= speed * up;
	if (GetAsyncKeyState( VK_LEFT )) D = normalize( D - right * 0.025f * speed );
	if (GetAsyncKeyState( VK_RIGHT )) D = normalize( D + right * 0.025f * speed );
	if (GetAsyncKeyState( VK_UP )) D = normalize( D - up * 0.025f * speed );
	if (GetAsyncKeyState( VK_DOWN )) D = normalize( D + up * 0.025f * speed );
#if 1
	// enable to set spline path points using P key
	static bool pdown = false;
	static FILE* pf = 0;
	if (!GetAsyncKeyState( 'P' )) pdown = false; else
	{
		if (!pdown) // save a point for the spline
		{
			if (!pf) pf = fopen( "spline.bin", "wb" );
			fwrite( &O, 1, sizeof( O ), pf );
			float3 t = O + D;
			fwrite( &t, 1, sizeof( t ), pf );
		}
		pdown = true;
	}
	LookAt( O, O + D );
#else
	// playback of recorded spline path
	const size_t N = splinePath.size();
	PathPoint p0 = splinePath[(pathPt + (N - 1)) % N], p1 = splinePath[pathPt];
	PathPoint p2 = splinePath[(pathPt + 1) % N], p3 = splinePath[(pathPt + 2) % N];
	LookAt( CatmullRom( p0.O, p1.O, p2.O, p3.O ), CatmullRom( p0.D, p1.D, p2.D, p3.D ) );
	if ((t += deltaTime * 0.0005f) > 1) t -= 1, pathPt = (pathPt + 1) % N;
#endif
}

// -----------------------------------------------------------
// Main application tick function
// -----------------------------------------------------------
void WaterWorld::Tick( float deltaTime )
{
	// update camera
	HandleInput( deltaTime );

	UpdateCAPE(deltaTime);
}
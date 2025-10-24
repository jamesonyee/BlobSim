/** 
 * StarterCode for "Attack of the Blobs!" 😱
 * CS248B Fundamentals of Computer Graphics: Animation & Simulation
 * 
 * Fill in the the missing code (see TODO items) + Comment Section (Names/Features/Bugs/Sources).
 * Starter tips: 
 *    1. Try reducing MAX_BLOBS to 1 to get started. 
 *    2. Set PEGS=false for fewer obstacles.
 *    3. Set SPIKES=false for a less pokey environment. 
 *    4. Look for 👀 & TODO items below for things to adjust & do.
 * Good luck!!
 * 
 * @author Doug L. James <djames@cs.stanford.edu> 
 * @date 10/28/2022; Updated 10/16/2025 (FINER EDITION, F25)
 */

// 🌲 Stanford CS248B "Attack of the Blobs" Assignment 😱
// 📜 Startercode: Doug L. James, djames@cs.stanford.edu
// 🍂 Fall 2025
/****************************************************************************************
 🙋 Names: Jameson Yee, Basant Khalil

 ✨ Features:
  - [Explain features you implemented so we don’t miss them]

 🐛 Bugs / Issues ⚠️:
  - [Describe problems or things that don’t work]

 📚 Sources:
  - [List any external code or images you used, if any]
*****************************************************************************************/
// >>> REMOVE the following console.log() line after filling out the writeup above <<<
console.log("⚠️ Fill out the writeup comment block! ⚠️");

p5.friendlyErrorSystem = false;

const MAX_BLOBS = 1; // 👀 TODO: 100 or more to complete "Attack of the Blobs!" challenge. Use just a few for testing. 

/// TIP: TURN OFF PEGS + SPIKES UNTIL YOU ARE READY:
const PEGS = false; // 👀 Turns on pachinko pegs. Bonk!
const SPIKES = false; // 👀 Turns on the spike pit and replaces pachinko peges with spikes!! 😲 Bwaahahahahaha!!!
const DRAW_BLOB_PARTICLES = true;

/// STIFFNESS PARAMETERS TO TWEAK: 🤨
const STIFFNESS_STRETCH = -20000.0; // 👀 TODO: Set as you wish
const STIFFNESS_BEND = 1.0; //    👀 TODO: Set as you wish
const STIFFNESS_AREA = 1.0; //    👀 TODO: Set as you wish

/// WORLD PARAMETERS (DON'T CHANGE)
const CANVAS_SIZE = 1024;
const WIDTH = 1.0;
const HEIGHT = 1.0;
const PARTICLE_RADIUS = WIDTH / 400.0; // for rendering
const PARTICLE_MASS = 1.0;
const BLOB_PARTICLES = 25;//(F25)  //15 (F24) // 12 (F22) // #particles per blob
const BLOB_RADIUS = WIDTH / 21; // WIDTH / 23; (F24) // WIDTH / 25; (F22/23)

//////// IMPORTANT ARRAYS OF THINGS /////////
let particles = []; // All particles in the scene (rigid + blobs)
let edges = []; //     All edges in the scene (rigid + blobs)
let blobs = []; //     All blobs in the scene (increases over time)
let environment; //    Environment with all rigid edges available as getEdges()

let isPaused = true;
let nTimesteps = 0; // #frame-length timesteps taken
let detectedEdgeEdgeFailure = false; // Halts simulation and turns purple if true -- blobs win!

// Graph paper texture map:
let bgImage;

function preload() {
	bgImage = loadImage('graphpaper.jpg');
}

function setup() {
	createCanvas(CANVAS_SIZE, CANVAS_SIZE); // SQUARE
	background(100);
	ellipseMode(RADIUS);
	environment = new Environment();
	//print("|particles|=" + particles.length + ",  |edge|=" + edges.length + ",  |blobs|=" + blobs.length);
}

/// Timesteps (w/ substeps) and draws everything.
function draw() {

	push();
	//noStroke();
	scale(height / HEIGHT); { // DRAW IN NORMALIZED COORDS

		///// SIMULATE /////
		if (!isPaused) {
			overlapEdges = null; // clear previous error visuals once simulating again.

			/// CREATE BLOBS ///
			if (nTimesteps % 10 == 0) {
				if (blobs.length < MAX_BLOBS)
					createRandomBlob(); // tries to create one if free space available
			}

			/// TIMESTEP! ///
			{
				let dtFrame = 0.01;
				let nSubsteps = 1; // 👀 #times to split dtFrame
				for (let step = 0; step < nSubsteps; step++)
					advanceTime(dtFrame / nSubsteps);
				nTimesteps++;
			}
		}

		///// RENDER /////
		{
			push();
			background(0);
			environment.draw(); // draws edgeOverlap error, too.
			for (let blob of blobs)
				blob.draw();
			pop();
			drawMouseForce();
			drawOverlapEdges(); /// draws any edge-edge overlap + error circle
		}
	} // END OF NORMALIZED-COORDINATES DRAWING
	pop();

	/// TEXT GUI (IN PIXEL COORDINATES)
	{
		push();
		textSize(18);
		noStroke();
		fill(0);
		text("#BLOBS: " + blobs.length, 10, 20);
		text("#EDGES: " + edges.length, 10, 40);
		text("#PARTICLES: " + particles.length, 10, 60);
		pop();
	}
}

function keyPressed() {
	if (keyCode == 32) // spacebar
		isPaused = !isPaused;
}

function advanceTime(dt) {
	environment.advanceTime(dt);

	//////////////////////////////////////
	////// GATHER PARTICLE FORCES ////////
	{
		// Clear forces:
		for (let particle of particles)
			particle.f.set(0, 0);

		gatherParticleForces_Gravity();

		// 👀 Damping (springs or otherwise -- you can add some if you want): 

		// Blob springs: 
		for (let blob of blobs) {
			blob.gatherForces_Stretch();
			blob.gatherForces_Bend();
			blob.gatherForces_Area();
		}

		gatherParticleForces_Penalty();

		// Mouse force (modify if you want):
		applyMouseForce();
	}

	//////////////////////////////////////////
	// Update velocity (using mass filtering):
	for (let particle of particles)
		vacc(particle.v, dt * particle.invMass(), particle.f)

	//////////////////////////////////////////
	// Collision filter: Correct velocities //
	applyPointEdgeCollisionFilter();

	//////////////////////////////////////////
	// Update positions:
	for (let particle of particles)
		vacc(particle.p, dt, particle.v)

	// VERIFY NO OVERLAP (USEFUL FOR DEBUGGING. SPEED-UP OR TURN-OFF FOR BIG RUNS.)
	verifyNoEdgeEdgeOverlap();
}

function applyPointEdgeCollisionFilter() {
	// 👀 TEMP HACK (remove!): rigid bounce off walls so they don't fly away
	for (let blob of blobs) blob.nonrigidBounceOnWalls();

	// 👀 TODO: Process all point-edge CCD impulses 
	// 👀 FIRST: Just rigid edges.
	let edgesToCheck = environment.getEdges();
	// 👀 SECOND: All rigid + blob edges (once you get this ^^ working)
	// edgesToCheck = edges;

	// Initially just brute force all-pairs checks, later use bounding volumes or better broad phase.

}

// Computes penalty forces between all point-edge pairs
function gatherParticleForces_Penalty() {

	let warmup = true;
	if (warmup) { // First just consider rigid environment edges:
		for (let edge of environment.getEdges()) {
			// 👀 TODO (part1): Apply point-edge force (if pt not on edge!)
		}
	} else { // Consider all rigid + blob edges:
		for (let edge of edges) {// or something smarter
			// 👀 TODO (part2): Apply point-edge force (if pt not on edge!)
		}
	}
}

function gatherParticleForces_Gravity() {
	let g = vec2(0, 0.1); // weak grav accel -- they kinda float (my uncle says they are hollow and made of f*rts)
	for (let particle of particles)
		vacc(particle.f, particle.mass, g); // f += m g
}


// Creates a default particle and adds it to particles list
function createParticle(x, y) {
	let p = new Particle(vec2(x, y), 1.0, PARTICLE_RADIUS);
	particles.push(p);
	return p;
}

// Creates edge and adds to edge list
function createEdge(particle0, particle1) {
	let edge = new Edge(particle0, particle1);
	edges.push(edge);
	return edge;
}
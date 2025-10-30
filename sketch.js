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
 🙋 Names:
  - 

 ✨ Features:
  - [Explain features you implemented so we don’t miss them]

 🐛 Bugs / Issues ⚠️:
  - [Describe problems or things that don’t work]

 📚 Sources:
  - [List any external code or images you used, if any]
*****************************************************************************************/

p5.friendlyErrorSystem = false;

const MAX_BLOBS = 10; // 👀 TODO: 100 or more to complete "Attack of the Blobs!" challenge. Use just a few for testing. 

/// TIP: TURN OFF PEGS + SPIKES UNTIL YOU ARE READY:
const PEGS = true;   // Turns on pachinko pegs
const SPIKES = true; // Turns on the spike pit!!
const DRAW_BLOB_PARTICLES = true;

/// STIFFNESS PARAMETERS TO TWEAK: 🤨
const STIFFNESS_STRETCH = 50.0; // 👀 TODO: Set as you wish
const k_d = 10.0;
const STIFFNESS_BEND = 50; //    👀 TODO: Set as you wish
const STIFFNESS_AREA = 50; //    👀 TODO: Set as you wish

/// WORLD PARAMETERS (DON'T CHANGE)
const CANVAS_SIZE = 1024;
const WIDTH = 1.0;
const HEIGHT = 1.0;
const PARTICLE_RADIUS = WIDTH / 400.0; // for rendering
const PARTICLE_MASS = 1.0;
const BLOB_PARTICLES = 10;//(F25)  //15 (F24) // 12 (F22) // #particles per blob
const BLOB_RADIUS = WIDTH / 21; // WIDTH / 23; (F24) // WIDTH / 25; (F22/23)

const STIFFNESS_PENALTY = 500.0; // Spring stiffness for collisions
const PENALTY_DISTANCE = PARTICLE_RADIUS * 8.0; // d0: force field range

const COEFFICIENT_OF_RESTITUTION = 0.5; // How bouncy (0 = dead, 1 = perfect)
const MAX_IMPULSE_ITERATIONS = 10;     // How many times to re-check collisions


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
}

/// Timesteps (w/ substeps) and draws everything.
function draw() {

	push();
	scale(height / HEIGHT); { // DRAW IN NORMALIZED COORDS
		
		///// SIMULATE /////
		if (!isPaused) {
			overlapEdges = null; // clear previous error visuals once simulating again.
			detectedEdgeEdgeFailure = false; // RESET THE FAILURE FLAG!

			/// CREATE BLOBS ///
			if (nTimesteps % 10 == 0) {
				if (blobs.length < MAX_BLOBS)
					createRandomBlob(); // tries to create one if free space available
			}

			/// TIMESTEP! ///
			{
				let dtFrame = 0.01;
				let nSubsteps = 20; // 👀 #times to split dtFrame
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
		console.log("end")
	if (keyCode == ENTER) { // noLoop debugging
		redraw();
	} else if (key == 'q') {
		clear();
		lineIndex = 0.0;
	}
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


	// Update blob bounds for collision detection
	for (let blob of blobs) {
		blob.updateBound(dt);
	}

	//////////////////////////////////////////
	// Collision filter: Correct velocities //
	applyPointEdgeCollisionFilter(dt); // Pass dt();

	//////////////////////////////////////////
	// Update positions:
	for (let particle of particles)
		vacc(particle.p, dt, particle.v)

	// VERIFY NO OVERLAP (USEFUL FOR DEBUGGING. SPEED-UP OR TURN-OFF FOR BIG RUNS.)
	verifyNoEdgeEdgeOverlap();
}


function applyPointEdgeCollisionFilter(dt) { // Receive dt
	const e = COEFFICIENT_OF_RESTITUTION;
	let envEdges = environment.getEdges();
	
	for (let iter = 0; iter < MAX_IMPULSE_ITERATIONS; iter++) {

		// === PART 1: BLOB-vs-ENVIRONMENT ===
		for (let blob of blobs) {
			for (let edge of envEdges) {
				// BROAD PHASE
				let edgeAABB = getEdgeAABB(edge);
				if (!aabbOverlap(blob.aabb_min, blob.aabb_max, edgeAABB.min, edgeAABB.max)) {
					continue;
				}
				// NARROW PHASE
				for (let p of blob.BP) {
					processImpulse(p, edge, e, dt); // Pass dt				
				}
			}
		}

		// === PART 2: BLOB-vs-BLOB ===
		for (let i = 0; i < blobs.length; i++) {
			for (let j = i + 1; j < blobs.length; j++) {
				let blobA = blobs[i];
				let blobB = blobs[j];

				// BROAD PHASE
				if (!aabbOverlap(blobA.aabb_min, blobA.aabb_max, blobB.aabb_min, blobB.aabb_max)) {
					continue;
				}

				// NARROW PHASE (A vs B)
				for (let p of blobA.BP) {
					for (let edge of blobB.BE) {
						processImpulse(p, edge, e, dt); // Pass dt					
					}
				}

				// NARROW PHASE (B vs A)
				for (let p of blobB.BP) {
					for (let edge of blobA.BE) {
						processImpulse(p, edge, e, dt); // Pass dt					
					}
				}
			}
		}
	} // end iterations
}

// Creates edge and adds to edge list
function createEdge(particle0, particle1) {
	let edge = new Edge(particle0, particle1);
	edges.push(edge);
	return edge;
}

function closestPointOnSegment(p, a, b) {
	let ab = vsub(b, a);
	let ap = vsub(p, a);

	// Handle zero-length edges (e.g., a spike tip)
	let lenSq = dot2(ab);
	if (lenSq < 1e-12) {
		return { point: a.copy(), t: 0.0 };
	}

	let t = ap.dot(ab) / lenSq;

	// Clamp t to [0, 1] for the segment
	t = clamp(t, 0.0, 1.0);

	let closest = vadd(a, vmult(ab, t));
	return { point: closest, t: t };
}

/**
 * Helper function to check for overlap between two AABBs.
 */
function aabbOverlap(minA, maxA, minB, maxB) {
	if (maxA.x < minB.x || minA.x > maxB.x) return false;
	if (maxA.y < minB.y || minA.y > maxB.y) return false;
	return true; // Overlapping on both axes
}

/**
 * Helper function to get the AABB for a single edge.
 */
function getEdgeAABB(edge) {
	let min_x = min(edge.q.p.x, edge.r.p.x) - PARTICLE_RADIUS;
	let min_y = min(edge.q.p.y, edge.r.p.y) - PARTICLE_RADIUS;
	let max_x = max(edge.q.p.x, edge.r.p.x) + PARTICLE_RADIUS;
	let max_y = max(edge.q.p.y, edge.r.p.y) + PARTICLE_RADIUS;
	return {
		min: vec2(min_x, min_y),
		max: vec2(max_x, max_y)
	};
}

/**
 * NARROW PHASE (PENALTY):
 * Applies a penalty force between a single particle and a single edge.
 */
function processPenalty(p, edge, k_p, k_d, d0) {
	// Don't collide a particle with an edge it belongs to
	if (p === edge.q || p === edge.r) {
		return;
	}
	// Don't collide a particle with an edge from its OWN blob
	if (p.blob && edge.q.blob && p.blob === edge.q.blob) {
		return;
	}

	// 1. Find closest point 'c' on the edge segment
	let result = closestPointOnSegment(p.p, edge.q.p, edge.r.p);
	let c = result.point;
	let t = result.t;

	// 2. Calculate penetration
	let v_pc = vsub(p.p, c);
	let d = v_pc.mag();
	if (d < 1e-6) return; // Safety check

	let penetration = d0 - d;

	// 3. If penetrating, apply forces
	if (penetration > 0) {
		let n = vmult(v_pc, 1.0 / d); // Normal

		// 5. Spring force
		let f_spring = k_p * penetration;

		// 6. Damping force
		let v_rel = p.v;
		if (!edge.isRigid()) { // Blob-vs-blob
			let v_q = edge.q.v;
			let v_r = edge.r.v;
			let v_edge_point = vadd(vmult(v_q, 1.0 - t), vmult(v_r, t));
			v_rel = vsub(p.v, v_edge_point);
		}
		let v_n = v_rel.dot(n);
		let f_damp = -k_d * v_n;

		// 7. Total force
		let f_mag = f_spring + f_damp;
		if (f_mag < 0) f_mag = 0; // Only apply repulsive push forces
		let f_penalty = vmult(n, f_mag);

		// 9. Apply force to particle
		p.f.add(f_penalty);

		// 10. Apply equal & opposite force to edge's particles
		if (!edge.isRigid()) {
			if (!edge.q.pin) {
				vacc(edge.q.f, -(1.0 - t), f_penalty);
			}
			if (!edge.r.pin) {
				vacc(edge.r.f, -t, f_penalty);
			}
		}
	}
}

/**
 * NARROW PHASE (IMPULSE):
 * Applies a collision impulse between a single particle and a single edge.
 * This version includes Baumgarte stabilization to resolve penetration.
 */
function processImpulse(p, edge, e, dt) { 
	// Don't collide a particle with an edge it belongs to
	if (p === edge.q || p === edge.r) {
		return;
	}
	// Don't collide a particle with an edge from its OWN blob
	if (p.blob && edge.q.blob && p.blob === edge.q.blob) {
		return;
	}

	// 1. Find closest point 'c'
	let result = closestPointOnSegment(p.p, edge.q.p, edge.r.p);
	let c = result.point;
	let t = result.t;

	// 2. Check for penetration
	let v_pc = vsub(p.p, c);
	let d = v_pc.mag();
	let penetration = PARTICLE_RADIUS - d;

	if (d < 1e-9) return; // Avoid division by zero if d is tiny

	let n = vmult(v_pc, 1.0 / d); // Collision normal

	// 3. Check relative velocity
	let v_rel = p.v;
	if (!edge.isRigid()) {
		let v_q = edge.q.v;
		let v_r = edge.r.v;
		let v_edge_point = vadd(vmult(v_q, 1.0 - t), vmult(v_r, t));
		v_rel = vsub(p.v, v_edge_point);
	}
	let v_n = v_rel.dot(n);

	// 4. If penetrating, apply impulse
	if (penetration > 0) { 
		// 5. Calculate impulse magnitude (j)
		let invMass_p = p.invMass();
		let invMass_edge = 0;
		if (!edge.isRigid()) {
			invMass_edge = sq(1.0 - t) * edge.q.invMass() + sq(t) * edge.r.invMass();
		}
		let invMass_eff = invMass_p + invMass_edge;
		if (invMass_eff < 1e-9) return; // Both are effectively pinned


		// A. Calculate target velocity for restitution (bounce)
		//    This is only applied if objects are moving towards each other.
		let v_restitution = 0;
		if (v_n < 0) {
			v_restitution = -e * v_n;
		}

		// B. Calculate target velocity for penetration correction (Baumgarte)
		//    This adds a "push" velocity to correct for penetration.
		const beta = 0.2; // Stabilization factor (0.1-0.3 is typical)
		const slop = 0.001; // Allow a tiny bit of penetration to avoid jitter
		let v_correction = 0;
		if (dt > 1e-9) { // Avoid division by zero if dt is tiny
			v_correction = (beta / dt) * max(0, penetration - slop);
		}

		// C. Final target velocity is the *larger* of the two.
		//    We either bounce, or we push, whichever is stronger.
		let v_target = max(v_restitution, v_correction);

		// D. Calculate the change in velocity and the final impulse magnitude
		let delta_v = v_target - v_n;
		let j = delta_v / invMass_eff;


		let J = vmult(n, j); // Final impulse vector

		// 6. Apply impulse
		vacc(p.v, p.invMass(), J);
		if (!edge.isRigid()) {
			vacc(edge.q.v, -(1.0 - t) * edge.q.invMass(), J);
			vacc(edge.r.v, -t * edge.r.invMass(), J);
		}
	}
}
// Computes penalty forces between all point-edge pairs
function gatherParticleForces_Penalty() {
	const k_p = STIFFNESS_PENALTY;
	const k_d_penalty = k_d;
	const d0 = PENALTY_DISTANCE;
	
	let envEdges = environment.getEdges();

	// === PART 1: BLOB-vs-ENVIRONMENT ===
	for (let blob of blobs) {
		for (let edge of envEdges) {
			// BROAD PHASE: Check if blob's AABB overlaps edge's AABB
			let edgeAABB = getEdgeAABB(edge);
			if (!aabbOverlap(blob.aabb_min, blob.aabb_max, edgeAABB.min, edgeAABB.max)) {
				continue; // Skip this edge, it's too far
			}

			// NARROW PHASE: Check all particles in this blob against this edge
			for (let p of blob.BP) {
				processPenalty(p, edge, k_p, k_d_penalty, d0);
			}
		}
	}

	// === PART 2: BLOB-vs-BLOB ===
	for (let i = 0; i < blobs.length; i++) {
		for (let j = i + 1; j < blobs.length; j++) {
			let blobA = blobs[i];
			let blobB = blobs[j];

			// BROAD PHASE: Check if blob A's AABB overlaps blob B's AABB
			if (!aabbOverlap(blobA.aabb_min, blobA.aabb_max, blobB.aabb_min, blobB.aabb_max)) {
				continue; // Skip this pair, they're too far
			}

			// NARROW PHASE (A vs B):
			// Check blob A's particles against blob B's edges
			for (let p of blobA.BP) {
				for (let edge of blobB.BE) {
					processPenalty(p, edge, k_p, k_d_penalty, d0);
				}
			}

			// NARROW PHASE (B vs A):
			// Check blob B's particles against blob A's edges
			for (let p of blobB.BP) {
				for (let edge of blobA.BE) {
					processPenalty(p, edge, k_p, k_d_penalty, d0);
				}
			}
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

/*
// "press enter to continue" utility for noLoop() mode
function keyPressed() {
	if (keyCode == ENTER) {
		redraw();
	} else if (key == 'q') {
		clear();
		lineIndex = 0.0;
	}
}*/

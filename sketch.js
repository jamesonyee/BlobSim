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

const MAX_BLOBS = 3;
const SPIKES = true;
const PEGS = true;

const DRAW_BLOB_PARTICLES = true;

const STIFFNESS_STRETCH  = 180.0;   
const STIFFNESS_AREA = 30.0;
const STIFFNESS_BEND = 2.5;

const STIFFNESS_PENALTY = 2000.0;
const k_d                = 1.2;

const COEFFICIENT_OF_RESTITUTION = 0.8;

const MAX_IMPULSE_ITERATIONS = 16;

/// WORLD PARAMETERS (DON'T CHANGE)
const CANVAS_SIZE = 1024;
const WIDTH = 1.0;
const HEIGHT = 1.0;
const PARTICLE_RADIUS = WIDTH / 400.0; // for rendering
const PARTICLE_MASS = 1.0;
const BLOB_PARTICLES = 10; //(F25) //15 (F24) // 12 (F22)
const BLOB_RADIUS = WIDTH / 21; // WIDTH / 23; (F24)

const PENALTY_DISTANCE = PARTICLE_RADIUS * 2.0;   // was 8.0 → too large

// Broad-phase safety pads (decoupled from penalty)
const AABB_PAD_PARTICLE = PARTICLE_RADIUS * 8.0;   // blob bound expansion
const AABB_PAD_EDGE     = PARTICLE_RADIUS * 2.5;   // edge bound expansion

//////// IMPORTANT ARRAYS OF THINGS /////////
let particles = [];
let edges = [];
let blobs = [];
let environment;
let isPaused = true;
let nTimesteps = 0;
let detectedEdgeEdgeFailure = false;

// Graph paper texture map:
let bgImage;

function preload() {
	bgImage = loadImage('graphpaper.jpg');
}

function setup() {
	createCanvas(CANVAS_SIZE, CANVAS_SIZE);
	background(100);
	ellipseMode(RADIUS);
	environment = new Environment();
}

/// Timesteps (w/ substeps) and draws everything.
function draw() {
	push();
	scale(height / HEIGHT);
	{
		///// SIMULATE /////
		if (!isPaused) {
			overlapEdges = null;
			detectedEdgeEdgeFailure = false;

			if (nTimesteps % 10 == 0) {
				if (blobs.length < MAX_BLOBS)
					createRandomBlob();
			}

			let dtFrame = 0.004;   // smaller step
			let nSubsteps = 20;    // more sub-iterations

			for (let step = 0; step < nSubsteps; step++)
				advanceTime(dtFrame / nSubsteps);
			nTimesteps++;
		}

		///// RENDER /////
		{
			push();
			background(0);
			environment.draw();
			for (let blob of blobs) blob.draw();
			pop();
			drawMouseForce();
			drawOverlapEdges();
		}
	}
	pop();

	push();
	textSize(18);
	noStroke();
	fill(0);
	text("#BLOBS: " + blobs.length, 10, 20);
	text("#EDGES: " + edges.length, 10, 40);
	text("#PARTICLES: " + particles.length, 10, 60);
	pop();

	// --- global damping once per frame ---
	for (let p of particles)
		p.v.mult(0.99985);

}

function keyPressed() {
	if (keyCode == 32)
		isPaused = !isPaused;
	console.log("end");
	if (keyCode == ENTER) {
		redraw();
	} else if (key == 'q') {
		clear();
		lineIndex = 0.0;
	}
}

function advanceTime(dt) {
	environment.advanceTime(dt);

	for (let particle of particles)
		particle.f.set(0, 0);

	gatherParticleForces_Gravity();

	for (let blob of blobs) {
		blob.gatherForces_Stretch();
		blob.gatherForces_Bend();
		blob.gatherForces_Area();
		blob.applyShapeMemory();
	}

	gatherParticleForces_Penalty();
	applyMouseForce();

	for (let particle of particles)
		vacc(particle.v, dt * particle.invMass(), particle.f);

	for (let blob of blobs)
		blob.updateBound(dt);

	applyPointEdgeCollisionFilter(dt);

	for (let particle of particles)
		vacc(particle.p, dt, particle.v);

	verifyNoEdgeEdgeOverlap();
}

function applyPointEdgeCollisionFilter(dt) {
	const e = COEFFICIENT_OF_RESTITUTION;
	let envEdges = environment.getEdges();

	for (let iter = 0; iter < MAX_IMPULSE_ITERATIONS; iter++) {
		// === PART 1: BLOB-vs-ENVIRONMENT ===
		for (let blob of blobs) {
			for (let edge of envEdges) {
				let edgeAABB = getEdgeAABB(edge);
				if (!aabbOverlap(blob.aabb_min, blob.aabb_max, edgeAABB.min, edgeAABB.max))
					continue;
				for (let p of blob.BP)
					processImpulse(p, edge, e, dt);
			}
		}

		// === PART 2: BLOB-vs-BLOB ===
		for (let i = 0; i < blobs.length; i++) {
			for (let j = i + 1; j < blobs.length; j++) {
				let blobA = blobs[i];
				let blobB = blobs[j];
				if (!aabbOverlap(blobA.aabb_min, blobA.aabb_max, blobB.aabb_min, blobB.aabb_max))
					continue;

				for (let p of blobA.BP)
					for (let edge of blobB.BE)
						processImpulse(p, edge, e, dt);

				for (let p of blobB.BP)
					for (let edge of blobA.BE)
						processImpulse(p, edge, e, dt);
			}
		}
	}
}

function createEdge(particle0, particle1) {
	let edge = new Edge(particle0, particle1);
	edges.push(edge);
	return edge;
}

function closestPointOnSegment(p, a, b) {
	let ab = vsub(b, a);
	let ap = vsub(p, a);
	let lenSq = dot2(ab);
	if (lenSq < 1e-12) return { point: a.copy(), t: 0.0 };
	let t = ap.dot(ab) / lenSq;
	t = clamp(t, 0.0, 1.0);
	let closest = vadd(a, vmult(ab, t));
	return { point: closest, t: t };
}

function aabbOverlap(minA, maxA, minB, maxB) {
	if (maxA.x < minB.x || minA.x > maxB.x) return false;
	if (maxA.y < minB.y || minA.y > maxB.y) return false;
	return true;
}


function getEdgeAABB(edge) {
  const pad = AABB_PAD_EDGE;
  let min_x = Math.min(edge.q.p.x, edge.r.p.x) - pad;
  let min_y = Math.min(edge.q.p.y, edge.r.p.y) - pad;
  let max_x = Math.max(edge.q.p.x, edge.r.p.x) + pad;
  let max_y = Math.max(edge.q.p.y, edge.r.p.y) + pad;
  return { min: vec2(min_x, min_y), max: vec2(max_x, max_y) };
}

function processPenalty(p, edge, k_p, k_d, d0) {
	if (p === edge.q || p === edge.r) return;
	if (p.blob && edge.q.blob && p.blob === edge.q.blob && !edge.isRigid()) return;

	let result = closestPointOnSegment(p.p, edge.q.p, edge.r.p);
	let c = result.point;
	let t = result.t;
	let v_pc = vsub(p.p, c);
	let d = v_pc.mag();
	if (d < 1e-6) return;
	let penetration = d0 - d;

	if (edge.length() < 0.05 && edge.isRigid()) {
	    k_p *= 0.4; // reduce stiffness for short "spike" edges
	}

	if (penetration > 0 && d < d0 * 0.95) {
		let n = vmult(v_pc, 1.0 / d);
		let depthScale = map(penetration, 0, d0, 0, 3.0, true);
		let f_spring = k_p * penetration * (0.5 + depthScale);
		let v_rel = p.v;

		if (!edge.isRigid()) {
			let v_q = edge.q.v;
			let v_r = edge.r.v;
			let v_edge_point = vadd(vmult(v_q, 1.0 - t), vmult(v_r, t));
			v_rel = vsub(p.v, v_edge_point);
		}

		if (edge.isRigid()) {
		  stroke('red');
		  line(edge.q.p.x, edge.q.p.y, edge.r.p.x, edge.r.p.y);
		  stroke('white');
		}

		let v_n = v_rel.dot(n);
		let f_damp = -k_d * v_n;
		let f_mag = f_spring + f_damp;
		if (f_mag < 0) f_mag = 0;
		let f_penalty = vmult(n, f_mag);
		if (edge.isRigid() && penetration > d0 * 0.6) {
		    let push = vmult(n, 0.2 * penetration);
		    p.p.add(push);
		}

		p.f.add(f_penalty);
			
		// --- Tangential friction to allow deformation around spikes ---
		let tangent = vec2(-n.y, n.x);
		let v_t = v_rel.dot(tangent);
		let f_friction = vmult(tangent, -0.5 * k_d * v_t);
		p.f.add(f_friction);


		if (!edge.isRigid()) {
			if (!edge.q.pin) vacc(edge.q.f, -(1.0 - t), f_penalty);
			if (!edge.r.pin) vacc(edge.r.f, -t, f_penalty);
		}
	}
}

function processImpulse(p, edge, e, dt) {
    if (p === edge.q || p === edge.r) return;
    if (p.blob && edge.q.blob && p.blob === edge.q.blob) return;

    let result = closestPointOnSegment(p.p, edge.q.p, edge.r.p);
    let c = result.point;
    let t = result.t;
    let v_pc = vsub(p.p, c);
    let d = v_pc.mag();
    if (d < 1e-9) return;

    let n = vmult(v_pc, 1.0 / d);
    let v_rel = p.v;
    if (!edge.isRigid()) {
        let v_q = edge.q.v;
        let v_r = edge.r.v;
        let v_edge_point = vadd(vmult(v_q, 1.0 - t), vmult(v_r, t));
        v_rel = vsub(p.v, v_edge_point);
    }

    const penetration = PARTICLE_RADIUS - d;
    const v_n = v_rel.dot(n);
	if (penetration <= -PARTICLE_RADIUS * 0.3 || v_n > 0) return;

    const invMass_p = p.invMass();
    let invMass_edge = 0;
    if (!edge.isRigid()) {
        invMass_edge = sq(1.0 - t) * edge.q.invMass() + sq(t) * edge.r.invMass();
    }
    const invMass_eff = invMass_p + invMass_edge;
    if (invMass_eff < 1e-9) return;

    // --- Restitution term (elastic bounce) ---
    let v_restitution = -e * v_n;

    // --- Position correction term (very light now to prevent blowup) ---
    const beta = 0.06;  // smaller than before
    const slop = 0.0005;
    let v_correction = 0;
    if (dt > 1e-9)
        v_correction = (beta / dt) * max(0, penetration - slop);

    const v_target = v_restitution + v_correction;
    const delta_v = v_target - v_n;
    const j = delta_v / invMass_eff;
    const J = vmult(n, j);

    // --- Apply impulse ---
    vacc(p.v, invMass_p, J);
    if (!edge.isRigid()) {
        vacc(edge.q.v, -(1.0 - t) * edge.q.invMass(), J);
        vacc(edge.r.v, -t * edge.r.invMass(), J);
    }
}

function gatherParticleForces_Penalty() {
	const k_p = STIFFNESS_PENALTY;
	const k_d_penalty = k_d * 0.5; // lighter damping for environment
	let envEdges = environment.getEdges();

	// --- blob vs environment ---
	for (let blob of blobs) {
		for (let edge of envEdges) {
			const d0 = edge.isRigid() ? PARTICLE_RADIUS * 3.5 : PENALTY_DISTANCE;

			let edgeAABB = getEdgeAABB(edge);
			if (!aabbOverlap(blob.aabb_min, blob.aabb_max, edgeAABB.min, edgeAABB.max))
				continue;
			for (let p of blob.BP)
				processPenalty(p, edge, k_p, k_d_penalty, d0);
		}
	}

	// --- blob vs blob ---
	for (let i = 0; i < blobs.length; i++) {
		for (let j = i + 1; j < blobs.length; j++) {
			let blobA = blobs[i];
			let blobB = blobs[j];
			if (!aabbOverlap(blobA.aabb_min, blobA.aabb_max, blobB.aabb_min, blobB.aabb_max))
				continue;

			for (let p of blobA.BP)
				for (let edge of blobB.BE)
					processPenalty(p, edge, k_p, k_d_penalty, PENALTY_DISTANCE);

			for (let p of blobB.BP)
				for (let edge of blobA.BE)
					processPenalty(p, edge, k_p, k_d_penalty, PENALTY_DISTANCE);
		}
	}
}


function gatherParticleForces_Gravity() {
	 let g = vec2(0, 1.1);
    for (let particle of particles)
        vacc(particle.f, particle.mass, g);
}


function createParticle(x, y) {
	let p = new Particle(vec2(x, y), 1.0, PARTICLE_RADIUS);
	particles.push(p);
	return p;
}

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
}
*/

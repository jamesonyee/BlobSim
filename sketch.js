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

// 🎭 Title animation state
let titleFade = 255;        // current title opacity
let titleDrip = [];         // active blood drips
let titleEvaporating = false; // whether the title is currently fading

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

let titleFont;
let subFont; // 👻 new HUD / subtitle font (Nosifer)


function preload() {
  bgImage = loadImage('night_forest_pumpkin.png');
  titleFont = loadFont('creepster.ttf'); // 🎃 spooky title font
  subFont = loadFont('nosifer.ttf'); // 💀 dripping-blood HUD font
}


// apply image filters *after* it’s loaded (not inside preload)
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
			// 🌫️ Parallax drifting fog layer (cinematic depth)
push();
blendMode(SOFT_LIGHT);
noStroke();
for (let i = 0; i < 3; i++) {
  let speed = 0.0001 + 0.0002 * i;
  let offset = (frameCount * speed) % 1;
  let y = (i * 0.3 + offset) * HEIGHT;
  fill(255, 180, 80, 15 - i * 5);
  rect(0, y, WIDTH, HEIGHT * 0.35);
}
pop();

// 👻 ghostly orbs drifting through mist
if (!window.ghosts) {
  window.ghosts = Array.from({length:10},()=>({
    x: random(WIDTH), y: random(HEIGHT),
    r: random(0.01,0.03), s: random(0.0004,0.0008)
  }));
}
push();
noStroke();
for (let g of window.ghosts){
  g.y -= g.s;
  if(g.y < -g.r){ g.y = HEIGHT; g.x = random(WIDTH); }
fill(255,255,255,0); // make all ghosts invisible
  if (g.x < 0.3 && g.y < 0.15) continue; // skip drawing ghosts overlapping HUD
	ellipse(g.x, g.y, g.r);
}
pop();


environment.draw();

// 🦇 flying bats layer (spooky silhouettes)
if (!window.bats) {
  window.bats = Array.from({length:8},()=>({
    x: random(WIDTH), y: random(HEIGHT*0.5),
    speed: random(0.0015,0.0025),
    flap: random(TWO_PI)
  }));
}
push();
fill(0,0,0,200);
noStroke();
for (let b of window.bats){
  b.x += b.speed;
  if(b.x>WIDTH+0.05){ b.x=-0.05; b.y=random(HEIGHT*0.5); }
  b.flap += 0.3;
  let wing = 0.02 + 0.005*sin(b.flap);
  triangle(b.x, b.y, b.x+0.03, b.y+wing, b.x-0.03, b.y+wing);
}
pop();

			// 🎃 Local glow around each blob (reactive pumpkin light)
push();
blendMode(ADD);
noStroke();
for (let blob of blobs) {
  let c = blob.centerOfMass();
  let pulse = 80 + 50 * sin(frameCount * 0.05 + blob.blobIndex);
  fill(255, 120, 30, pulse);
  ellipse(c.x * width, c.y * height, width * 0.12);
}
pop();


for (let blob of blobs) blob.draw();
// 👀 Detect when blobs start moving => trigger title evaporation
if (!titleEvaporating && blobs.length > 0 && blobs[0].centerOfMass().y > 0.2) {
  titleEvaporating = true;
}

			pop();
			drawMouseForce();
			drawOverlapEdges();
		}

		// 💡 subtle screen-wide light flicker
push();
blendMode(ADD);

fill(255, 120, 40, 10 + 10 * sin(frameCount * 0.07));
ellipse(width/2, height/2, width * 1.2);
fill(180, 180, 255, 4 + 4 * sin(frameCount * 0.03));
ellipse(width/2, height/2, width * 1.4);

pop();

	}

// 💡 Neon glow enhancer (sharp + saturated pop look)
push();
// use normal blending instead of ADD → prevents haze stacking
blendMode(BLEND);
noStroke();

// a tight, subtle orange center light
fill(255, 130, 40, 25);
ellipse(width/2, height/2, width * 0.4);

// optional tiny blue accent (barely visible)
fill(100, 200, 255, 10);
ellipse(width/2, height * 0.7, width * 0.5);
pop();


// 🌈 Haunted world color pulse (slow hue drift)
push();
blendMode(OVERLAY);
let hueShift = 40 + 40 * sin(frameCount * 0.003);
fill(255, 100 + hueShift, 40, 6);
ellipse(width/2, height/2, width * 1.5);

pop();


	pop();

// 🩸 TITLE EVAPORATION + BLOOD DRIP EFFECT
push();
textAlign(CENTER);
textFont(titleFont);
textSize(90);

// gradually fade out when blobs descend
if (titleEvaporating && titleFade > 0) {
titleFade -= 8.4; // faster fade (half the time)
  // spawn new blood drips randomly
if (random() < 0.55) {
    titleDrip.push({
      x: width/2 - 200 + random(400),
      y: 150,
      speed: random(1,3),
      alpha: 255,
      size: random(4,8)
    });
  }
}

// main glowing title
for (let i = 0; i < 4; i++) {
  let glow = 100 - i * 20;
  fill(255, 130 + i * 20, 0, titleFade * (glow / 150));
  text("ATTACK OF THE BLOBS", width/2 + i * 0.5, 90 - i);
}

// flickering light shimmer per letter
let title = "ATTACK OF THE BLOBS";
for (let i = 0; i < title.length; i++) {
  let ch = title[i];
  let x = width / 2 - textWidth(title) / 2 + textWidth(title.substring(0, i));
  let y = 90 + 5 * sin(frameCount * 0.05 + i * 0.8);
  let flicker = 150 + 50 * sin(frameCount * 0.03 + i);
  fill(255, 160, 40, min(titleFade, flicker));
  stroke(255, 220, 120, titleFade * 0.8);
  strokeWeight(2.5);
  text(ch, x, y);
}

// 💀 subtitle: slowly fades too
textFont(subFont);
textSize(28);
fill(255, 40, 0, titleFade * 0.8);
noStroke();
text("Haunted Harvest Edition", width/2, 135 + 2*sin(frameCount*0.05));

// 🩸 dripping blood simulation
noStroke();
for (let i = titleDrip.length - 1; i >= 0; i--) {
  let d = titleDrip[i];
  d.y += d.speed;
  d.alpha -= 4;
  fill(255, 0, 0, d.alpha);
  ellipse(d.x, d.y, d.size, d.size * 1.3);
  if (d.alpha <= 0 || d.y > height * 0.9) titleDrip.splice(i, 1);
}

// 🔥 optional subtle smoke fade effect around text
noStroke();
for (let i = 0; i < 5; i++) {
  fill(255, 150, 80, titleFade * 0.03);
  ellipse(width/2 + random(-200, 200), 100 + random(-20, 20), random(100, 200), 30);
}
pop();

// 🎃 Haunted HUD — ghostly carved pumpkin text
push();
textFont(subFont);
textSize(26);
textAlign(LEFT);

// 👻 floating fog behind the HUD (completely transparent to prevent white circles)
push();
noStroke();
blendMode(SOFT_LIGHT);
for (let i = 0; i < 4; i++) {
  // completely invisible background layer to disable HUD glow
  fill(255, 140, 30, 0);
  ellipse(110 + i * 10, 65 + i * 5, 80, 30);
}
pop();


// 💀 eerie flicker colors (orange–amber glow)
let baseGlow = 140 + 70 * sin(frameCount * 0.07);
let emberShift = 100 + 50 * sin(frameCount * 0.15);
if (random() < 0.003) baseGlow = 255; // rare lightning flash over HUD

// carved glowing look
fill(255, baseGlow, 0, 240);
stroke(255, 180 + 40 * sin(frameCount * 0.1), 60, 200);
strokeWeight(3.2);
text(`🎃 BLOBS: ${blobs.length}`, 25, 40);

// ghostly green shimmer overlay
fill(120, 255, 180, 60 + 40 * sin(frameCount * 0.09));
noStroke();
text(`🎃 BLOBS: ${blobs.length}`, 25 + 1.5 * sin(frameCount * 0.05), 40 + 1.5 * cos(frameCount * 0.05));

// 🕸️ edges — subtle red flicker
fill(255, 100 + 80 * sin(frameCount * 0.08), 0, 200);
stroke(255, 180, 50, 150);
strokeWeight(2.5);
text(`🕸️ EDGES: ${edges.length}`, 25, 70);

// 💀 particles — eerie ghost-green decay glow
let decayPulse = 120 + 60 * sin(frameCount * 0.09);
fill(140, 255, 100, 220);
stroke(90, decayPulse, 40, 180);
strokeWeight(2.5);
text(`💀 PARTICLES: ${particles.length}`, 25, 100);
// 🕯️👻🦇 HALLOWEEN HUD EFFECTS — stacked visual layers
push();

// 🔒 restrict visuals to the HUD area only (no white circles bleeding out)
drawingContext.save();
drawingContext.beginPath();
drawingContext.rect(0, 0, 300, 130); // clip region for top-left HUD
drawingContext.clip();

blendMode(ADD);
noStroke();

// === 1️⃣ Floating Ghost Sprites ===
if (!window.hudGhosts) {
  window.hudGhosts = Array.from({length:3},()=>({
    x: random(40,180),
    y: random(30,110),
    drift: random(0.002,0.004),
    phase: random(TWO_PI)
  }));
}
for (let g of window.hudGhosts) {
  g.x += 0.15*sin(frameCount*g.drift);
  g.y += 0.1*cos(frameCount*g.drift + g.phase);
  fill(255,255,255,50 + 40*sin(frameCount*0.03 + g.phase));
  ellipse(g.x, g.y, 18 + 4*sin(frameCount*0.05 + g.phase), 14);
  fill(0,0,0,80);
  ellipse(g.x-4, g.y-1, 3,3);
  ellipse(g.x+4, g.y-1, 3,3);
}

// === 2️⃣ Tiny Flying Bats ===
if (!window.hudBats) {
  window.hudBats = Array.from({length:5},()=>({
    x: random(15,220),
    y: random(20,120),
    speed: random(0.8,1.3),
    flap: random(TWO_PI)
  }));
}
fill(0,0,0,150);
for (let b of window.hudBats){
  b.x += 0.3*b.speed;
  if(b.x>230){ b.x=10; b.y=random(30,120); }
  b.flap += 0.4*b.speed;
  let wing = 6 + 2*sin(b.flap);
  triangle(b.x, b.y, b.x+wing, b.y+2, b.x-wing, b.y+2);
}

// === 3️⃣ Dripping Candle Flames next to stats ===
for (let i=0;i<3;i++){
  let cx = 250;
  let cy = 38 + i*30;
  let flicker = 150 + 80*sin(frameCount*0.2 + i);
  fill(255, 180, 60, flicker);
  ellipse(cx, cy, 5 + 1*sin(frameCount*0.4 + i));
  fill(255, 100, 20, 180);
  rect(cx-1, cy, 2, 10, 1);
  if (random()<0.01){
    let dropY = cy + 10 + random(10);
    fill(255,150,40,180);
    ellipse(cx, dropY, 2,4);
  }
}

// === 4️⃣ Floating Eyeballs Watching Player ===
if (!window.hudEyes) {
  window.hudEyes = Array.from({length:2},()=>({
    baseX: random(70,150),
    baseY: random(40,90),
    r: 10
  }));
}
for (let e of window.hudEyes){
  let lookX = map(mouseX,0,width,-2,2);
  let lookY = map(mouseY,0,height,-1,1);
  fill(255,255,240,200);
  ellipse(e.baseX, e.baseY, e.r*2);
  fill(100,0,0,220);
  ellipse(e.baseX+lookX, e.baseY+lookY, e.r*0.8);
  fill(0);
  ellipse(e.baseX+lookX*1.2, e.baseY+lookY*1.2, e.r*0.4);
}

// === 5️⃣ Pumpkin Sparks Burst ===
for (let i=0; i<4; i++){
  let px = 60 + random(-20,200);
  let py = 40 + random(-10,90);
  fill(255, 140 + random(60), 30, 100);
  ellipse(px, py, random(2,4));
}

// 🧹 cleanup: stop bleed outside HUD
drawingContext.restore();
pop();



	// --- global damping once per frame ---
	for (let p of particles)
		p.v.mult(0.99985);

	// 🔥 Floating embers
if (!window.embers)
  window.embers = Array.from({length:50},()=>({x:random(width),y:random(height),s:random(0.5,1.5)}));
push();
blendMode(ADD);
noStroke();
for (let e of window.embers){
  e.y -= 0.5*e.s; if (e.y<0) {e.y=height; e.x=random(width);}
  fill(255,140+random(60),30,150);
  ellipse(e.x,e.y,2*e.s);
}
pop();

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

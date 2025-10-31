////////////////////////////////////////////////////
/// SUPPORTING CODE FOR ENVIRONMENT SETUP & DRAW /// 
////////////////////////////////////////////////////
// Nothing 👀 TODO here. 

// RIGID ENVIRONMENT COMPOSED OF LINE SEGMENTS (pinned Edges)
class Environment {

	constructor() {
		this.envParticles = [];
		this.envEdges = [];

		///// BOX /////
		let r = PARTICLE_RADIUS;
		this.p00 = createParticle(r, r);
		this.p01 = createParticle(r, HEIGHT - r);
		this.p11 = createParticle(WIDTH - r, HEIGHT - r);
		this.p10 = createParticle(WIDTH - r, r);
		this.p00.pin = this.p01.pin = this.p11.pin = this.p10.pin = true;
		this.envParticles.push(this.p00);
		this.envParticles.push(this.p01);
		this.envParticles.push(this.p11);
		this.envParticles.push(this.p10);
		this.envEdges.push(createEdge(this.p00, this.p01));
		this.envEdges.push(createEdge(this.p01, this.p11));
		this.envEdges.push(createEdge(this.p11, this.p10));
		this.envEdges.push(createEdge(this.p10, this.p00));

		///// OBSTACLES FOR FUN /////
		{
			// (F22) ANGLED LINES: 
			//for (let i = 0.5; i < 4; i++) this.createEnvEdge(i * width / 5, height / 2, (i + 1) * width / 5, height * 0.75);

			// (F23) PACHINKO PEGS:
			// (F23) PACHINKO SPIKES 😲: 
			if (PEGS) {
				let n = 8;
				for (let i = 1; i < n; i++) {
					for (let j = 1; j < n; j++) {
						if ((i + j) % 2 == 0) { // alternating pegs
							if (SPIKES)
								this.createSpike(WIDTH * (i / n), HEIGHT / 5 + HEIGHT / 4 * 3 * (j / n), 0.02); // spikes 😲
							else {
								//this.createPachinkoPeg(WIDTH * (i / n), HEIGHT / 5 + HEIGHT / 4 * 3 * (j / n), 0.01); // round 4-edge peg
								this.createPachinkoWedge(WIDTH * (i / n), HEIGHT / 5 + HEIGHT / 4 * 3 * (j / n), 0.01); // cheap 2-edge peg
							}
						}
					}
				}
			}

			// (F24) SPIKES AT BOTTOM, TOO, IS SO 2024 😲: 
			if (SPIKES) {
				let n = 20;
				for (let i = 1; i < n; i++) {
					this.createSpike(WIDTH * (i / n), HEIGHT, 0.02);
				}
			}
		}
	}

	// Returns all rigid-environment edges.
	getEdges() {
		return this.envEdges;
	}

	// Creates a lone rigid edge.
	createEnvEdge(x0, y0, x1, y1) {
		let p0 = createParticle(x0, y0);
		let p1 = createParticle(x1, y1);
		p0.pin = true;
		p1.pin = true;
		let e = createEdge(p0, p1);
		this.envParticles.push(p0);
		this.envParticles.push(p1);
		this.envEdges.push(e);
	}
	// Create a lone roundish peg at (x0,y0) with radius r.
	createPachinkoPeg(x, y, r) {
		let p0 = createParticle(x + r, y);
		let p1 = createParticle(x, y + r);
		let p2 = createParticle(x - r, y);
		let p3 = createParticle(x, y - r);
		p0.pin = p1.pin = p2.pin = p3.pin = true;
		let e01 = createEdge(p0, p1);
		let e12 = createEdge(p1, p2);
		let e23 = createEdge(p2, p3);
		let e30 = createEdge(p3, p0);
		this.envParticles.push(p0);
		this.envParticles.push(p1);
		this.envParticles.push(p2);
		this.envParticles.push(p3);
		this.envEdges.push(e01);
		this.envEdges.push(e12);
		this.envEdges.push(e23);
		this.envEdges.push(e30);
	}
	// Create a lighter-weight wedge-shaped peg at (x0,y0) with radius r.
	createPachinkoWedge(x, y, r) {
		let p0 = createParticle(x + r, y);
		let p1 = createParticle(x, y - r);
		let p2 = createParticle(x - r, y);
		p0.pin = p1.pin = p2.pin = true;
		let e01 = createEdge(p0, p1);
		let e12 = createEdge(p1, p2);
		this.envParticles.push(p0);
		this.envParticles.push(p1);
		this.envParticles.push(p2);
		this.envEdges.push(e01);
		this.envEdges.push(e12);
	}
	// Create a small upward spike placed at (x0,y0) with length r.
	createSpike(x, y, r) {
		let p0 = createParticle(x, y);
		let p1 = createParticle(x, y - r);
		p0.pin = p1.pin = true;
let e01 = createEdge(p0, p1);
e01.glowPhase = random(TWO_PI); // for flicker phase
		this.envParticles.push(p0);
		this.envParticles.push(p1);
		this.envEdges.push(e01);
	}

	// Updates any moveable infinite-mass rigid elements
	advanceTime(dt) {}

	// Makes popcorn <jk> no it doesn't... 
	draw() {
		push();
// 🌲 Full Halloween forest backdrop (draw first)
push();
tint(255, 255 * 1.25); // boost brightness
image(bgImage, 0, 0, WIDTH, HEIGHT);
tint(255, 200, 160); // slightly brighter
pop();

push();
blendMode(ADD);
fill(255, 140, 60, 10);
rect(0, 0, width, height);
pop();

		
// 🌌 Deep contrast enhancer (no fog, boosts color & brightness)
push();
blendMode(ADD);
noStroke();

// 🔥 warm light center for “candle” mood
// Subtle lighting overlay (adds depth without fog)
push();
blendMode(ADD);
noStroke();

// Warm candle glow (localized)
fill(255, 140, 60, 18);
ellipse(width / 2, height * 0.65, width * 0.4);

// Cool moon rim-light (faint)
fill(100, 120, 255, 10);
ellipse(width / 2, height * 0.25, width * 0.6);
pop();


pop();

// ✅ NOW draw spikes *after* the background and fog so they show on top
push();
// 🕯️ tall glowing spikes (visible candle effect)
for (let edge of this.envEdges) {
  if (edge.length() < 0.05) {
    let flicker = 180 + 75 * sin(frameCount * 0.3 + edge.glowPhase);
stroke(255, 200, 100, flicker);
strokeWeight(PARTICLE_RADIUS * 6.0);

    drawingContext.shadowBlur = 30;
drawingContext.shadowColor = frameCount % 120 < 60 ? "#ff6600" : "#00e5ff";
    line(edge.q.p.x, edge.q.p.y, edge.r.p.x, edge.r.p.y - 0.015);
    fill(255, 220, 150, 220);
    noStroke();
    ellipse(edge.r.p.x, edge.r.p.y - 0.02, 0.012 + 0.004 * sin(frameCount * 0.5));
	  // faint candle aura
drawingContext.shadowBlur = 60;
drawingContext.shadowColor = "rgba(255,160,60,0.8)";

  }
}
		
pop();


// 👻 glowing blue-green ghost pegs (keep both effects active)
for (let edge of this.envEdges) {
  if (edge.length() >= 0.05 && edge.length() < 0.2) {
    stroke(120, 255, 255, 80 + 60 * sin(frameCount * 0.07 + edge.glowPhase));
    strokeWeight(PARTICLE_RADIUS);
    drawingContext.shadowBlur = 20;
    drawingContext.shadowColor = "cyan";
    edge.draw();
  }
}
drawingContext.shadowBlur = 0;


// 🌕 Pegs = cold ghostly lanterns
for (let edge of this.envEdges) {
  if (edge.length() >= 0.05) {
    stroke(80, 200, 255, 100 + 50 * sin(frameCount * 0.05 + edge.glowPhase));
    strokeWeight(PARTICLE_RADIUS * 0.8);
    drawingContext.shadowBlur = 10;
    drawingContext.shadowColor = "cyan";
    edge.draw();
  }
}
drawingContext.shadowBlur = 0;

		noStroke();
		for (let particle of this.envParticles) {
			particle.draw();
		}
stroke("#00ffcc"); // ghostly glint
		if (floor(frameCount * 0.3) % 2 == 0) {
			for (let edge of this.envEdges) {
				if (random() > 0.95)
					if (edge.length() < 0.2) { //short edges:
						let a = edge.q.p;
						let b = edge.r.p;
						line(a.x, a.y, b.x, b.y);
					}
			}
		}
		drawingContext.shadowBlur = 0;
for (let edge of this.envEdges) {
  if (edge.length() < 0.05) { // short spikes
    let glow = 100 + 80 * sin(frameCount * 0.1 + (edge.glowPhase || 0));
    stroke(255, 180, 80, glow);
    line(edge.q.p.x, edge.q.p.y, edge.r.p.x, edge.r.p.y);
  }
}

		pop(); // wait, it does pop :/ 
	}
}

// Creates a blob centered at (x,y), and adds things to lists (blobs, edges, particles).
function createBlob(x, y) {
	let b = new Blob(vec2(x, y));
	b.blobIndex = blobs.length; // useful for avoiding self-collision checks.
	blobs.push(b);
	return b;
}

// Tries to create a new blob at the top of the screen. 
function createRandomBlob() {
	for (let attempt = 0; attempt < 5; attempt++) {
		let center = vec2(random(2 * BLOB_RADIUS, WIDTH - 2 * BLOB_RADIUS), BLOB_RADIUS * 1.3); //random horizontal spot
		// CHECK TO SEE IF NO BLOBS NEARBY:
		let tooClose = false;
		for (let blob of blobs) {
			let com = blob.centerOfMass();
			if (com.dist(center) < 3 * blob.radius) { // too close
				tooClose = true;
				break;
			}
		}
		// if we got here, then center is safe:
		if (!tooClose) {
			createBlob(center.x, center.y);
			return;
		}
	}
}

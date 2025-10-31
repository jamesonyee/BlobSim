///////////////////////////////////////////////
/// BLOB CLASS. MANY 👀 TODO ITEMS BELOW ///
//////////////////////////////////////////////
class Blob {
	constructor(centerRest) {
		//this.blobId = -1;// undefined until set. (you can give them IDs if you want)
		this.radius = BLOB_RADIUS;
		this.centerRest = centerRest; // original location
		this.blobIndex = -1;// set to 0...(nBlobs-1) when created
		this.A0 = 0;

		// CREATE PARTICLES:
		this.hit = false;
		this.BP = []; //blob particles
		this.n = BLOB_PARTICLES;
		//let v0 = vec2(1, 0)
		let v0 = vec2(random(-0.1, 0.1), random(0.20, 0.22));
		let blobAngle = TWO_PI * random();
		for (let i = 0; i < this.n; i++) {
			let bump = 1.0 - 0.1 * abs(sin(2.5 * TWO_PI * i / this.n)); //nonspherical
			let ri = this.radius * bump;
			let xi = ri * cos(i / this.n * TWO_PI);
			let yi = ri * sin(i / this.n * TWO_PI);
			let v = vec2(xi, yi).rotate(blobAngle);
			v.add(centerRest);
			let particle = createParticle(v.x, v.y);
			particle.v.set(v0);
			particle.blob = this; // set blob ref
			this.BP.push(particle);
		}

		// CREATE EDGES FOR STRETCH SPRINGS + COLLISIONS:
		this.BE = []; // blob edges
		for (let i = 0; i < this.n; i++) {
			let p0 = this.BP[i];
			let p1 = this.BP[(i + 1) % this.n];
			this.BE.push(createEdge(p0, p1));
		}

		// SETUP YOUR APPEARANCE/FACIAL ELEMENTS:
		// 👀 TODO
		// 🎃 Halloween pumpkin orange tones
		let dc = 25;
		// 🎃 Randomly choose pumpkin, ghost, or slime tone
		const palette = [
		    color(255, 130 + random(-10, 10), 0),      // pumpkin
		    color(200 + random(-20, 20), 255, 240),   // ghostly white
		    color(100 + random(-20, 10), 255, 100)    // eerie green slime
		];
		this.fillColor = random(palette);
		

		// Rest Area
		let sum = 0;
		for (let i = 0; i < this.n; i++) {
			console.log("index: " + i)
			let particle1 = this.BP[i];
			let particle2 = this.BP[(i + 1) % this.n];
			let particle0;
			if (i - 1 < 0) {
				particle0 = this.BP[this.n - 1];
			} else {
				particle0 = this.BP[i - 1];
			}
			let p0 = particle0.p;
			let p1 = particle1.p;
			let p2 = particle2.p;
			console.log("particle0: " + p0)
			console.log("particle1: " + p1)
			console.log("particle2: " + p2)
			let b = vsub(p0, p2);
			sum += (b.x * p1.y) - (b.y * p1.x);
			console.log("sum: " + sum)
		}
		this.A0 = 0.5 * sum;

		this.aabb_min = vec2(0, 0);
		this.aabb_max = vec2(0, 0);

		this.eyeSize = 4 * PARTICLE_RADIUS;
		this.pupilSize = 1.5 * PARTICLE_RADIUS;
		this.eyeSpacing = 4 * PARTICLE_RADIUS;
		this.blinkTimer = 0;
		this.blinkDuration = 10; // in frames
	}

	getBlobIndex() { return this.blobIndex; }

	blobParticles() { return this.BP; }

	gatherForces_Stretch() {
		let k_s = STIFFNESS_STRETCH;
		let k_d_spring = k_d;
		for (let edge of this.BE) {
			let p0 = edge.q;
			let p1 = edge.r;
			let L0 = edge.restLength;
			let v = vsub(p1.p, p0.p);
			let L = v.mag();
			if (L < 1e-9) continue;
			let n = vmult(v, 1.0 / L);
			let f_spring_mag = k_s * (L - L0);
			let v_rel = vsub(p1.v, p0.v);
			let v_n = dot(v_rel, n);
			let f_damp_mag = k_d_spring * v_n;
			let f_total_mag = f_spring_mag + f_damp_mag;
			let f_total = vmult(n, f_total_mag);
			p0.f.add(f_total);
			p1.f.sub(f_total);
		}
	}

	gatherForces_Bend() {
		let k = STIFFNESS_BEND;
		for (let i = 0; i < this.n; i++) {
			this.hit = true;
			let p1 = this.BP[i];
			let p2 = this.BP[(i + 1) % this.n];
			let p0;
			if (i - 1 < 0) { p0 = this.BP[this.n - 1]; }
			else { p0 = this.BP[i - 1]; }

			let vec_a = vec2(p1.p.x - p0.p.x, p1.p.y - p0.p.y)
			let length_a = length(vec_a);
			let a_hat = vec_a.div(length_a);
			let vec_b = vec2(p2.p.x - p1.p.x, p2.p.y - p1.p.y);
			let length_b = length(vec_b)
			let b_hat = vec_b.div(length_b);

			let f0 = vax(-k / (2 * length_a), (vsub(b_hat, vax(dot(a_hat, b_hat), a_hat))))
			let f2 = vax(k / (2 * length_b), (vsub(a_hat, vax(dot(a_hat, b_hat), b_hat))));
			let neg_f0 = vax(-1, f0)
			let f1 = vsub(neg_f0, f2);

			p0.f = vadd(p0.f, f0);
			p1.f = vadd(p1.f, f1);
			p2.f = vadd(p2.f, f2);
		}
	}
	
	gatherForces_Area() {
	    const k = STIFFNESS_AREA;
	    // --- compute current signed area ---
	    let A = 0;
	    for (let i = 0; i < this.n; i++) {
	        const p0 = this.BP[i].p;
	        const p1 = this.BP[(i + 1) % this.n].p;
	        A += (p0.x * p1.y - p1.x * p0.y);
	    }
	    A *= 0.5;
	
	    const areaDiff = (A - this.A0);
	
	    // --- apply corrective pressure to each vertex ---
	    for (let i = 0; i < this.n; i++) {
	        const pPrev = this.BP[(i - 1 + this.n) % this.n];
	        const pNext = this.BP[(i + 1) % this.n];
	        const edge = vsub(pNext.p, pPrev.p);
	        const n = vec2(-edge.y, edge.x);        // outward normal
	        const f = vmult(n, -k * areaDiff / this.n); // distribute evenly
	        this.BP[i].f.add(f);
	    }
	}
	
	updateBound(dt) {
	  if (this.BP.length === 0) return;
	  let min_x =  Infinity, min_y =  Infinity;
	  let max_x = -Infinity, max_y = -Infinity;
	
	  for (let p of this.BP) {
	    min_x = Math.min(min_x, p.p.x);
	    min_y = Math.min(min_y, p.p.y);
	    max_x = Math.max(max_x, p.p.x);
	    max_y = Math.max(max_y, p.p.y);
	
	    // include motion over the substep
	    const nx = p.p.x + dt * p.v.x;
	    const ny = p.p.y + dt * p.v.y;
	    min_x = Math.min(min_x, nx);
	    min_y = Math.min(min_y, ny);
	    max_x = Math.max(max_x, nx);
	    max_y = Math.max(max_y, ny);
	  }
	
	  const padding = AABB_PAD_PARTICLE; // <- decoupled from PENALTY_DISTANCE
	  min_x -= padding; min_y -= padding;
	  max_x += padding; max_y += padding;
	
	  this.aabb_min.set(min_x, min_y);
	  this.aabb_max.set(max_x, max_y);
	}
	

	drawBound() {
		push();
		noFill();
		stroke("magenta");
		strokeWeight(0.002);
		rectMode(CORNERS);
		rect(this.aabb_min.x, this.aabb_min.y, this.aabb_max.x, this.aabb_max.y);
		pop();
	}

	centerOfMass() {
		let com = vec2(0, 0);
		for (let particle of this.BP) {
			vacc(com, 1 / this.BP.length, particle.p);
		}
		return com;
	}

	centerOfVelocity() {
		let cov = vec2(0, 0);
		for (let particle of this.BP) vacc(cov, 1 / this.BP.length, particle.v);
		return cov;
	}

	rigidBounceOnWalls() {
		let pos = this.centerOfMass();
		let vel = this.centerOfVelocity();
		let R = BLOB_RADIUS + PARTICLE_RADIUS;
		if ((pos.x < R && vel.x < 0) || (pos.x > WIDTH - R && vel.x > 0)) {
			for (let particle of this.BP) particle.v.x *= -0.4;
		}
		if ((pos.y < R && vel.y < 0) || (pos.y > HEIGHT - R && vel.y > 0)) {
			for (let particle of this.BP) particle.v.y *= -0.4;
		}
	}

	nonrigidBounceOnWalls() {
		let R = PARTICLE_RADIUS;
		for (let particle of this.BP) {
			let pos = particle.p;
			let vel = particle.v;
			if ((pos.x < R && vel.x < 0) || (pos.x > WIDTH - R && vel.x > 0)) {
				vel.x *= -0.4;
			}
			if ((pos.y < R && vel.y < 0) || (pos.y > HEIGHT - R && vel.y > 0)) {
				vel.y *= -0.4;
			}
		}
	}

	draw() {
	  push();
	
	  strokeWeight(PARTICLE_RADIUS);
	  stroke(50, 0, 80, 80);
	
	  // keep particles visible + glowing
	  let alpha = this.hit ? 180 : 230;
	  let flicker = 0.6 + 0.4 * sin(frameCount * 0.15 + this.blobIndex);
	  let theme = this.theme || ["ghost","pumpkin","skull"][floor(random(3))];
	  this.theme = theme; // freeze per blob
	
	  // === choose palette per type ===
	  if (theme === "ghost") this.fillColor = color(200, 240, 255, alpha * 0.7);
	  if (theme === "pumpkin") this.fillColor = color(255, 120 + 30*sin(frameCount*0.2), 0, alpha);
	  if (theme === "skull") this.fillColor = color(245, 245, 230, alpha * 0.8);
	
	  // === draw outer shape ===
	  drawingContext.shadowBlur = 30;
	  drawingContext.shadowColor =
	    theme === "ghost" ? "rgba(180,220,255,0.9)" :
	    theme === "pumpkin" ? "rgba(255,130,0,0.9)" :
	    "rgba(255,255,255,0.6)";
	
	  beginShape();
	  for (let p of this.BP) vertex(p.p.x, p.p.y);
	  endShape(CLOSE);
	
	  // --- visible particles integrated into body ---
	  for (let p of this.BP) {
	    noStroke();
	    if (theme === "ghost") fill(220, 255, 255, 80 + 40*sin(frameCount*0.3));
	    if (theme === "pumpkin") fill(255, 100, 0, 120 + 60*sin(frameCount*0.3));
	    if (theme === "skull") fill(255, 255, 255, 100 + 40*sin(frameCount*0.3));
	    circle(p.p.x, p.p.y, PARTICLE_RADIUS * 1.5);
	  }
	
	  // === per-theme facial design ===
	  if (theme === "ghost") this.drawGhostFace();
	  else if (theme === "pumpkin") this.drawPumpkinFace();
	  else this.drawSkullFace();
	
	  this.drawBound();
	  drawingContext.shadowBlur = 0;
	  pop();
	}
	
	
	//------------------------------------------------------
	// 👻 GHOST FACE — floating soul
	//------------------------------------------------------
	drawGhostFace() {
	  push();
	  let c = this.centerOfMass();
	  noStroke();
	
	  // hollow eyes
	  fill(0, 0, 0, 220);
	  ellipse(c.x - 0.015, c.y - 0.01, 0.018, 0.024);
	  ellipse(c.x + 0.015, c.y - 0.01, 0.018, 0.024);
	
	  // soft mouth
	  fill(0, 0, 0, 160);
	  ellipse(c.x, c.y + 0.018, 0.022, 0.012);
	
	  // faint cyan glow around eyes
	  drawingContext.shadowBlur = 25;
	  drawingContext.shadowColor = "rgba(180,255,255,0.8)";
	  noFill();
	  stroke(180, 255, 255, 80);
	  strokeWeight(0.003);
	  ellipse(c.x - 0.015, c.y - 0.01, 0.025);
	  ellipse(c.x + 0.015, c.y - 0.01, 0.025);
	  pop();
	}
	
	
	//------------------------------------------------------
	// 🎃 PUMPKIN FACE — jack-o’-lantern
	//------------------------------------------------------
	drawPumpkinFace() {
	  push();
	  let c = this.centerOfMass();
	  drawingContext.shadowBlur = 30;
	  drawingContext.shadowColor = "rgba(255,120,0,0.9)";
	  stroke(0);
	  strokeWeight(0.001);
	  fill(0);
	
	  // triangle eyes
	  beginShape();
	  vertex(c.x - 0.02, c.y - 0.01);
	  vertex(c.x - 0.01, c.y - 0.03);
	  vertex(c.x, c.y - 0.01);
	  endShape(CLOSE);
	  beginShape();
	  vertex(c.x + 0.02, c.y - 0.01);
	  vertex(c.x + 0.01, c.y - 0.03);
	  vertex(c.x, c.y - 0.01);
	  endShape(CLOSE);
	
	  // jagged grin
	  noFill();
	  stroke(255, 120 + 60*sin(frameCount*0.3), 0);
	  strokeWeight(0.004);
	  beginShape();
	  for (let i = -3; i <= 3; i++) {
	    let x = c.x + i * 0.006;
	    let y = c.y + 0.02 + 0.004 * ((i % 2 == 0) ? 1 : -1);
	    vertex(x, y);
	  }
	  endShape();
	  pop();
	}
	
	
	//------------------------------------------------------
	// 💀 SKULL FACE — eerie bone wisp
	//------------------------------------------------------
	drawSkullFace() {
	  push();
	  let c = this.centerOfMass();
	  drawingContext.shadowBlur = 20;
	  drawingContext.shadowColor = "rgba(255,255,255,0.7)";
	  noStroke();
	  fill(0);
	  ellipse(c.x - 0.015, c.y - 0.01, 0.018);
	  ellipse(c.x + 0.015, c.y - 0.01, 0.018);
	
	  // nose hole
	  triangle(c.x - 0.004, c.y + 0.002, c.x + 0.004, c.y + 0.002, c.x, c.y + 0.008);
	
	  // teeth lines
	  stroke(0);
	  strokeWeight(0.002);
	  for (let i = -3; i <= 3; i++) {
	    line(c.x + i*0.004, c.y + 0.018, c.x + i*0.004, c.y + 0.023);
	  }
	  pop();
	}


	drawBlobFace() {
		push();
		// Candle glow pulsation
		let flicker = 0.8 + 0.2 * sin(frameCount * 0.3 + random(0, 1));
		drawingContext.shadowBlur = 25 * flicker;
		drawingContext.shadowColor = "orange";

		let com = this.centerOfMass();
		let cov = this.centerOfVelocity();
		let up = vec2(0, 0);
		for (let p of this.BP) { up.add(vsub(p.p, com)); }
		up.normalize();
		if (up.magSq() < 1e-6) { up.set(0, -1); }
		let right = vec2(-up.y, up.x);

		this.blinkTimer--;
		if (this.blinkTimer < 0 && random() < 0.01) {
			this.blinkTimer = this.blinkDuration;
		}
		let isBlinking = this.blinkTimer > 0;

		let leftEyePos = vadd(com, vmult(right, -this.eyeSpacing));
		let rightEyePos = vadd(com, vmult(right, this.eyeSpacing));
		let lookDir = cov.copy();
		if (lookDir.magSq() < 0.0001) {
			lookDir.set(0, 0);
		} else {
			lookDir.normalize().mult(this.eyeSize * 0.3);
		}
		// Halloween glow
drawingContext.shadowBlur = 20;
drawingContext.shadowColor = "rgba(255, 140, 0, 0.8)";

		stroke(0);
		strokeWeight(0.001);
fill(255, 140, 0); // orange glow eyes
		circle(leftEyePos.x, leftEyePos.y, this.eyeSize);
		if (isBlinking) {
			strokeWeight(0.002);
			line(leftEyePos.x - this.eyeSize, leftEyePos.y, leftEyePos.x + this.eyeSize, leftEyePos.y);
		} else {
			fill(0);
			circle(leftEyePos.x + lookDir.x, leftEyePos.y + lookDir.y, this.pupilSize);
		}
fill(255, 120 + 40 * flicker, 0);
		circle(rightEyePos.x, rightEyePos.y, this.eyeSize);
		if (isBlinking) {
			strokeWeight(0.002);
			line(rightEyePos.x - this.eyeSize, rightEyePos.y, rightEyePos.x + this.eyeSize, rightEyePos.y);
		} else {
			fill(0);
			circle(rightEyePos.x + lookDir.x, rightEyePos.y + lookDir.y, this.pupilSize);
		}
// 🎃 spooky grin when compressed
let mouthCenter = this.centerOfMass(); // renamed variable
let mouthWidth = this.eyeSpacing * 2;
let mouthHeight = this.eyeSize * 0.6;
noFill();
stroke(255, 100 + 60 * sin(frameCount * 0.3), 0);
strokeWeight(0.004);
arc(mouthCenter.x, mouthCenter.y + this.eyeSpacing * 0.8, mouthWidth, mouthHeight, 0, PI);

		pop();
	}

	applyShapeMemory() {
		 const k_mem = 4.5; // stronger “spring” back toward rest
	    const com = this.centerOfMass();
	    for (let i = 0; i < this.n; i++) {
	        const restVec = vsub(this.BP[i].p, com);
	        const target = vadd(this.centerRest, restVec);
	        const correction = vsub(target, this.BP[i].p);
	        vacc(this.BP[i].f, k_mem, correction);
	    }
	}


}

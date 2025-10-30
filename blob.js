///////////////////////////////////////////////
/// BLOB CLASS.   MANY 👀 TODO ITEMS BELOW ///
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
			let bump = 1.0 - 0.1 * abs(sin(2.5 * TWO_PI * i / this.n)); //nonspherical bump (-0.2 (F24))
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
		let dc = 26;
		this.fillColor = color([221 + random(-dc, dc), 160 + random(-dc, dc), 221 + random(-dc, dc), 255]); // ("Plum"); // 221, 160, 221

		// Rest Area
		let sum = 0;
		for (let i = 0; i < this.n; i++) {
			console.log("index: " + i)
			let particle1 = this.BP[i];
			let particle2 = this.BP[(i + 1) % this.n];
			let particle0;
			if (i -1 < 0) {
				particle0 = this.BP[this.n - 1];
			}
			else {
				particle0 = this.BP[i-1];
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
	
	getBlobIndex() {
		return this.blobIndex;
	}

	blobParticles() {
		return this.BP;
	}

	// Loops over blob edges and accumulates stretch forces (Particle.f += ...)
	gatherForces_Stretch() {
		let k_s = STIFFNESS_STRETCH; // Spring stiffness
		let k_d_spring = k_d;       // Damping stiffness

		for (let edge of this.BE) {
			let p0 = edge.q;
			let p1 = edge.r;
			let L0 = edge.restLength; // Rest length

			// 1. Get current vector and length
			let v = vsub(p1.p, p0.p); // Vector from p0 to p1
			let L = v.mag();         // Current length

			if (L < 1e-9) {
				continue; // Avoid division by zero
			}

			// 2. Calculate force direction
			let n = vmult(v, 1.0 / L); // Normalized direction vector

			// 3. Calculate spring force (Hooke's Law)
			let f_spring_mag = k_s * (L - L0);

			// 4. Calculate damping force
			let v_rel = vsub(p1.v, p0.v); // Relative velocity
			let v_n = dot(v_rel, n);      // Relative velocity along the spring
			let f_damp_mag = k_d_spring * v_n;

			// 5. Total force magnitude
			let f_total_mag = f_spring_mag + f_damp_mag;
			let f_total = vmult(n, f_total_mag);

			// 6. Apply forces to particles
			p0.f.add(f_total);
			p1.f.sub(f_total);
		}
	}
	
	// Loops over blob particles and accumulates bending forces (Particle.f += ...)
	gatherForces_Bend() {
		let k = STIFFNESS_BEND;
		for (let i = 0; i < this.n; i++) {
			this.hit = true;
			let p1 = this.BP[i];
			let p2 = this.BP[(i + 1) % this.n];
			let p0;
			if (i - 1 < 0) {
				p0 = this.BP[this.n -1];
			}
			else {
				p0 = this.BP[i-1];
			} 
			let vec_a = vec2(p1.p.x - p0.p.x, p1.p.y - p0.p.y) // vector a = p1 - p0
			let length_a = length(vec_a);
			let a_hat = vec_a.div(length_a);

			let vec_b = vec2(p2.p.x - p1.p.x, p2.p.y - p1.p.y); // vector b = p2 - p1
			let length_b = length(vec_b)
			let b_hat = vec_b.div(length_b);

			// Calculate forces
			let f0 = vax(-k / (2*length_a),(vsub(b_hat, vax(dot(a_hat, b_hat), a_hat))))
			let f2 = vax(k / (2*length_b),(vsub(a_hat, vax(dot(a_hat, b_hat), b_hat))));
			let neg_f0 = vax(-1, f0)
			let f1 = vsub(neg_f0, f2);

			// add bending forces to particles
			p0.f = vadd(p0.f, f0);
			p1.f = vadd(p1.f, f1);
			p2.f = vadd(p2.f, f2);

		}
	}

	// Loops over blob particles and gathers area compression forces (Particle.f += ...)
	gatherForces_Area() {
		let k = STIFFNESS_AREA;

		// Deformed Area
		let sum = 0;
		for (let i = 0; i < this.n; i++) {
			let particle1 = this.BP[i];
			let particle2 = this.BP[(i + 1) % this.n];
			let particle0;
			if (i -1 < 0) {
				particle0 = this.BP[this.n - 1];
			}
			else {
				particle0 = this.BP[i-1];
			} 

			let p0 = particle0.p;
			let p1 = particle1.p;
			let p2 = particle2.p;
			let b = vsub(p0, p2);
			sum += (b.x * p1.y) - (b.y * p1.x);
		}
		let A = 0.5 * sum;

		// partial derivative of A with respect to p_i
		for (let i = 0; i < this.n; i++) {
			let particle1 = this.BP[i];
			let particle2 = this.BP[(i + 1) % this.n];
			let particle0;
			if (i -1 < 0) {
				particle0 = this.BP[this.n - 1];
			}
			else {
				particle0 = this.BP[i-1];
			} 

			let p0 = particle0.p;
			let p1 = particle1.p;
			let p2 = particle2.p;

			let b = vsub(p0, p2);
			let b_perp = vec2(-1 * b.y, b.x);
			let di_A = b_perp.div(2)
			let f_i = vax(-k * (A - this.A0), di_A);
			particle1.f = vadd(particle1.f, f_i);
		}
	}

	/**
	 * Updates the blob's Axis-Aligned Bounding Box (AABB).
	 * This bound encloses the swept motion for the *next* timestep (dt)
	 * and is inflated by the penalty distance.
	 */
	updateBound(dt) {
		if (this.BP.length === 0) return;

		// 1. Find the min/max of all current particle positions (p)
		//    and predicted next positions (p_next = p + dt*v)
		let min_x = Infinity, min_y = Infinity;
		let max_x = -Infinity, max_y = -Infinity;

		for (let p of this.BP) {
			// Check current position
			min_x = min(min_x, p.p.x);
			min_y = min(min_y, p.p.y);
			max_x = max(max_x, p.p.x);
			max_y = max(max_y, p.p.y);

			// Calculate predicted next position
			let p_next_x = p.p.x + dt * p.v.x;
			let p_next_y = p.p.y + dt * p.v.y;

			// Check predicted next position
			min_x = min(min_x, p_next_x);
			min_y = min(min_y, p_next_y);
			max_x = max(max_x, p_next_x);
			max_y = max(max_y, p_next_y);
		}
		
		// 2. Inflate the bound by the penalty distance (d0)
		//    (and a little extra just to be safe)
		let padding = PENALTY_DISTANCE + PARTICLE_RADIUS * 2;
		min_x -= padding;
		min_y -= padding;
		max_x += padding;
		max_y += padding;

		// 3. Store the final AABB
		this.aabb_min.set(min_x, min_y);
		this.aabb_max.set(max_x, max_y);
	}

	/**
	 * (For Debugging) Draws the blob's AABB.
	 */
	drawBound() {
		push();
		noFill();
		stroke("magenta"); // A bright color so we can see it
		strokeWeight(0.002);
		rectMode(CORNERS);
		rect(this.aabb_min.x, this.aabb_min.y, this.aabb_max.x, this.aabb_max.y);
		pop();
	}

	// Center of mass of all blob particles
	centerOfMass() {
		let com = vec2(0, 0);
		for (let particle of this.BP) {
			vacc(com, 1 / this.BP.length, particle.p); // assumes equal mass
		}
		return com;
	}

	// Center of velocity of all blob particles
	centerOfVelocity() {
		let cov = vec2(0, 0);
		for (let particle of this.BP)
			vacc(cov, 1 / this.BP.length, particle.v); // assumes equal mass
		return cov;
	}

	// Something simple to keep rigid blobs inside the box:
	rigidBounceOnWalls() {
		let pos = this.centerOfMass();
		let vel = this.centerOfVelocity();

		let R = BLOB_RADIUS + PARTICLE_RADIUS;

		// Boundary reflection (only if outside domain AND still moving outward):
		if ((pos.x < R && vel.x < 0) ||
			(pos.x > WIDTH - R && vel.x > 0)) {
			for (let particle of this.BP)
				particle.v.x *= -0.4;
		}
		if ((pos.y < R && vel.y < 0) ||
			(pos.y > HEIGHT - R && vel.y > 0)) {
			for (let particle of this.BP)
				particle.v.y *= -0.4;
		}
	}

	// Something simple to keep nonrigid blob particles inside the box:
	nonrigidBounceOnWalls() {
		let R = PARTICLE_RADIUS;
		for (let particle of this.BP) {
			let pos = particle.p;
			let vel = particle.v;

			// Boundary reflection (only if outside domain AND still moving outward):
			if ((pos.x < R && vel.x < 0) ||
				(pos.x > WIDTH - R && vel.x > 0)) {
				vel.x *= -0.4;
			}
			if ((pos.y < R && vel.y < 0) ||
				(pos.y > HEIGHT - R && vel.y > 0)) {
				vel.y *= -0.4;
			}
		}
	}

	draw() {
		push();
		strokeWeight(PARTICLE_RADIUS);
		stroke("DarkOrchid"); 
		fill(this.fillColor); { // draw blob
			beginShape(TESS);
			for (let particle of this.BP)
				vertex(particle.p.x, particle.p.y);
			endShape(CLOSE);
		}

		if (DRAW_BLOB_PARTICLES) {
			fill("DarkOrchid");
			if (this.hit == true) {
				stroke("yellow")
				fill("yellow"); // yellow
			}
			for (let particle of this.BP)
				circle(particle.p.x, particle.p.y, PARTICLE_RADIUS);
		}

		this.drawBlobFace();
		this.drawBound(); // (For debugging)
		pop();
	}


	drawBlobFace() {
		push();

		// --- 1. Find a stable center and orientation ---
		let com = this.centerOfMass();
		let cov = this.centerOfVelocity();
		let r = BLOB_RADIUS * 0.2; // Scale for face features

		// Calculate an "up" vector. We can average the vectors from the
		// center of mass to each particle. This gives a stable orientation.
		let up = vec2(0, 0);
		for (let p of this.BP) {
			up.add(vsub(p.p, com));
		}
		up.normalize();

		// If the blob is perfectly symmetrical, up might be (0,0).
		// In that case, just default to pointing straight up.
		if (up.magSq() < 1e-6) {
			up.set(0, -1);
		}

		// Get a "right" vector using the 2D cross product (perp)
		let right = vec2(-up.y, up.x);

		// --- 2. Handle Blinking ---
		this.blinkTimer--;
		if (this.blinkTimer < 0 && random() < 0.01) {
			// Start a blink
			this.blinkTimer = this.blinkDuration;
		}
		let isBlinking = this.blinkTimer > 0;

		// --- 3. Calculate Eye Positions ---
		let leftEyePos = vadd(com, vmult(right, -this.eyeSpacing));
		let rightEyePos = vadd(com, vmult(right, this.eyeSpacing));

		// --- 4. Calculate Pupil Position ---
		// Pupils look in the direction of velocity
		let lookDir = cov.copy();
		if (lookDir.magSq() < 0.0001) {
			lookDir.set(0, 0); // Not moving, look forward
		} else {
			lookDir.normalize().mult(this.eyeSize * 0.3);
		}
		
		// --- 5. Draw Eyes ---
		stroke(0);
		strokeWeight(0.001);
		
		// Draw Left Eye
		fill(255); // White part
		circle(leftEyePos.x, leftEyePos.y, this.eyeSize);
		if (isBlinking) {
			strokeWeight(0.002);
			line(leftEyePos.x - this.eyeSize, leftEyePos.y, leftEyePos.x + this.eyeSize, leftEyePos.y);
		} else {
			fill(0); // Pupil
			circle(leftEyePos.x + lookDir.x, leftEyePos.y + lookDir.y, this.pupilSize);
		}

		// Draw Right Eye
		fill(255); // White part
		circle(rightEyePos.x, rightEyePos.y, this.eyeSize);
		if (isBlinking) {
			strokeWeight(0.002);
			line(rightEyePos.x - this.eyeSize, rightEyePos.y, rightEyePos.x + this.eyeSize, rightEyePos.y);
		} else {
			fill(0); // Pupil
			circle(rightEyePos.x + lookDir.x, rightEyePos.y + lookDir.y, this.pupilSize);
		}
		
		// ---- MOUTH ----
		// Position mouth relative to blob orientation (like the eyes)
		let mouthOffset = vmult(up, 1.5 * r); // Move down from center in blob's "up" direction
		let mouthPos = vadd(com, mouthOffset);

		// Calculate the angle of the blob's orientation
		let angle = atan2(up.y, up.x) - PI/2; // up vector angle + 90° to make it vertical

		// Save current transformation
		push();
		// Move to mouth position and rotate
		translate(mouthPos.x, mouthPos.y);
		rotate(angle);

		// Outer mouth (black)
		fill(0);
		let mouthW = 1.5 * r;
		let mouthH = 1.5 * r;
		arc(0, 0, mouthW, mouthH, 0, PI, CHORD); // Draw at (0,0) since we translated

		// Inner mouth (pink)
		fill(255, 192, 203); // Pink color
		let innerMouthW = mouthW * 0.7; // Slightly smaller
		let innerMouthH = mouthH * 0.6;
		arc(0, 0, innerMouthW, innerMouthH, 0, PI, CHORD); // Draw at (0,0)

		pop();
	}

}

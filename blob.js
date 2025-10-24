///////////////////////////////////////////////
/// BLOB CLASS.   MANY 👀 TODO ITEMS BELOW ///
//////////////////////////////////////////////
class Blob {
	constructor(centerRest) {
		//this.blobId = -1;// undefined until set. (you can give them IDs if you want)
		this.radius = BLOB_RADIUS;
		this.centerRest = centerRest; // original location
		this.blobIndex = -1;// set to 0...(nBlobs-1) when created

		// CREATE PARTICLES:
		this.BP = []; //blob particles
		this.n = BLOB_PARTICLES;
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
	}
	
	getBlobIndex() {
		return this.blobIndex;
	}

	blobParticles() {
		return this.BP;
	}

	// Loops over blob edges and accumulates stretch forces (Particle.f += ...)                      
	gatherForces_Stretch() {
	let k = STIFFNESS_STRETCH;  
	for (let edge of this.BE) {
		let p0 = edge.q;
		let p1 = edge.r;

		let delta = vsub(p1.p, p0.p);
		let L = delta.mag();
		let L0 = edge.restLength;
		if (L === 0) continue;
		
			let dir = vmult(delta, 1 / L);

			let stretch = L - L0;
			let force = vmult(dir, k * stretch);

			vacc(p0.f, 1, force);
			vacc(p1.f, -1, force);
		}
	}

	// Loops over blob particles and accumulates bending forces (Particle.f += ...)
	gatherForces_Bend() {
		let k = STIFFNESS_BEND;
		for (let i = 0; i < this.n; i++) {
			// FIND particles before (p0) and after (p2) particle i (p1): 
			//let p0 = ????
			let p1 = this.BP[i];
			let p2 = this.BP[(i + 1) % this.n];

			// 👀 TODO

		}
	}
	// Loops over blob particles and gathers area compression forces (Particle.f += ...)
	gatherForces_Area() {
		let k = STIFFNESS_AREA;
		// TODO 

	}

	// Center of mass of all blob particles
	centerOfMass() {
		let com = vec2(0, 0);
		for (let particle of this.BP)
			vacc(com, 1 / this.BP.length, particle.p); // assumes equal mass
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
		stroke("DarkOrchid"); //BlueViolet");
		fill(this.fillColor); { // draw blob
			beginShape(TESS);
			for (let particle of this.BP)
				vertex(particle.p.x, particle.p.y);
			endShape(CLOSE);
		}

		if (DRAW_BLOB_PARTICLES) {
			fill("DarkOrchid");
			for (let particle of this.BP)
				circle(particle.p.x, particle.p.y, PARTICLE_RADIUS);
		}

		this.drawBlobFace();
		pop();
	}

	drawBlobFace() {
		push();
		// 👀 TODO: Draw your face here! :D 

		// CENTER OF MASS eyeball for now :/
		let com = this.centerOfMass();
		let cov = this.centerOfVelocity();
		stroke(0);
		fill(255);
		circle(com.x, com.y, 5 * PARTICLE_RADIUS);
		fill(0);
		circle(com.x + 0.01 * cov.x + 0.005 * sin(nTimesteps / 3), com.y + 0.01 * cov.y + 0.001 * random(-1, 1), PARTICLE_RADIUS);
		pop();
	}

}
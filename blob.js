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
	
	  // --- flicker & theme selection ---
	  let alpha = this.hit ? 180 : 230;
	  let flicker = 0.6 + 0.4 * sin(frameCount * 0.1 + this.blobIndex);
	  let theme = this.theme || ["ghost","pumpkin","skull"][floor(random(3))];
	  this.theme = theme;
	
	  // --- body contour glow ---
	  drawingContext.shadowBlur = 40;
	  drawingContext.shadowColor =
	    theme === "ghost" ? "rgba(170,240,255,0.9)" :
	    theme === "pumpkin" ? "rgba(255,100,0,0.9)" :
	    "rgba(255,240,180,0.8)";
	  noStroke();
	
	  // --- dynamic fill flicker ---
	  if (theme === "ghost") fill(210, 250, 255, 80 + 50*flicker);
	  if (theme === "pumpkin") fill(255, 120 + 40*sin(frameCount*0.2), 0, 100 + 50*flicker);
	  if (theme === "skull") fill(255, 255, 245, 90 + 40*flicker);
	
	  beginShape(); for (let p of this.BP) vertex(p.p.x, p.p.y); endShape(CLOSE);
	
	  // --- visible particles (integrated) ---
	  for (let p of this.BP) {
	    noStroke();
	    if (theme === "ghost") fill(160, 255, 255, 70 + 40*sin(frameCount*0.3));
	    if (theme === "pumpkin") fill(255, 120, 0, 90 + 40*sin(frameCount*0.3));
	    if (theme === "skull") fill(255, 255, 255, 80 + 40*sin(frameCount*0.3));
	    circle(p.p.x, p.p.y, PARTICLE_RADIUS * 1.6);
	  }
	
	  // --- spectral interior (new core) ---
	  let c = this.centerOfMass();
	  blendMode(ADD);
	  for (let i = 0; i < 3; i++) {
	    if (theme === "ghost") fill(120, 200, 255, 30 + 20*i);
	    if (theme === "pumpkin") fill(255, 80 + 60*i, 0, 25 + 15*i);
	    if (theme === "skull") fill(255, 220, 180, 25 + 15*i);
	    ellipse(c.x, c.y, this.radius * (0.6 + 0.1 * i) + 0.002*sin(frameCount*0.5+i));
	  }
	
	  // --- pulsating veins (motion energy) ---
	  blendMode(BLEND);
	  stroke(255, 200, 60, 35);
	  strokeWeight(0.0009);
	  for (let i=0;i<this.BP.length;i+=2){
	    let a=this.BP[i].p,b=this.BP[(i+3)%this.BP.length].p;
	    line(a.x,a.y,b.x,b.y);
	  }
	
	  // --- face ---
	  if (theme==="ghost") this.drawGhostFace(c,flicker);
	  else if (theme==="pumpkin") this.drawPumpkinFace(c,flicker);
	  else this.drawSkullFace(c,flicker);
	
	  drawingContext.shadowBlur=0;
	  pop();
	}
	
	// 👻 GHOST FACE — hollow eyes + eerie light mouth
	drawGhostFace(c, flicker){
		push();
		if (this.hit) {
		  drawingContext.shadowColor = "rgba(255,255,255,1)";
		  drawingContext.shadowBlur = 80;
		}
		
		drawingContext.shadowBlur = 50;
		drawingContext.shadowColor = "rgba(160,255,255,0.9)";
		noStroke();
	
		// 👁 Animated hollow eyes (haunting version)
		let look = 0.002 * sin(frameCount * 0.5); // subtle eye movement
		let pulse = 0.002 * sin(frameCount * 0.3); // breathing size change
		let glowPulse = 60 + 40 * sin(frameCount * 0.4); // glow intensity
			
		// moving hollow eyes
		fill(0, 0, 0, 220);
		ellipse(c.x - 0.015 + look, c.y - 0.01, 0.018 + pulse, 0.024 + pulse);
		ellipse(c.x + 0.015 + look, c.y - 0.01, 0.018 + pulse, 0.024 + pulse);
		
		// flickering cyan glow
		drawingContext.shadowBlur = 25 + 10 * sin(frameCount * 0.2);
		drawingContext.shadowColor = `rgba(180,255,255,${0.5 + 0.3 * sin(frameCount * 0.3)})`;
		noFill();
		stroke(180, 255, 255, glowPulse);
		strokeWeight(0.003 + 0.001 * sin(frameCount * 0.4));
		ellipse(c.x - 0.015 + look, c.y - 0.01, 0.025 + pulse);
		ellipse(c.x + 0.015 + look, c.y - 0.01, 0.025 + pulse);

		// Wailing mouth with ghostly teeth glow
	   fill(0,0,0,230); ellipse(c.x,c.y+0.018,0.034,0.018);
	   stroke(180,255,255,120); strokeWeight(0.0015);
	   for(let i=-3;i<=3;i++){
		   let x=c.x+i*0.004, y=c.y+0.018;
		   line(x,y,x,y+0.004*sin(frameCount*0.5+i));
	   }
		pop();
	}

	// 🎃 PUMPKIN FACE — infernal eyes + sawtooth grin
	drawPumpkinFace(c,flicker){
	  push();
	  drawingContext.shadowBlur = 45;
	  drawingContext.shadowColor = "rgba(255,100,0,0.9)";
	  noStroke();
	
	  // Eyes
	  for (let s of [-1,1]){
	    let ex=c.x+s*0.018, ey=c.y-0.015;
	    fill(255,90,0,100); ellipse(ex,ey,0.022);
	    fill(0); ellipse(ex,ey,0.018);
	    fill(255,0,0,220);
	    ellipse(ex+0.001*sin(frameCount*0.6+s),ey,0.006,0.012);
	  }
	
	  // Old-style jagged glowing teeth (restored)
	  noFill();
	  stroke(255,120+80*flicker,0);
	  strokeWeight(0.004);
	  beginShape();
	  for(let i=-3;i<=3;i++){
	    let x=c.x+i*0.006;
	    let y=c.y+0.02+0.005*((i%2==0)?1:-1);
	    vertex(x,y);
	  }
	  endShape();
	
	  // Inner ember flicker behind teeth
	  noStroke(); fill(255,100,0,50+50*flicker);
	  ellipse(c.x,c.y+0.02,0.04,0.012);
	  pop();
	}
	
	// 💀 SKULL — bone shell, fire eyes, black skeletal fangs
	drawSkullFace(c,flicker){
	  push();
	  drawingContext.shadowBlur = 55;
	  drawingContext.shadowColor = "rgba(255,180,80,0.8)";
	  noStroke();
	
	  // Eyes – hellfire cores
	  for (let s of [-1,1]){
	    let ex=c.x+s*0.018, ey=c.y-0.012;
	    fill(255,140,0,130+80*flicker); ellipse(ex,ey,0.026);
	    fill(0); ellipse(ex,ey,0.018);
	    fill(255,30,0,240);
	    ellipse(ex,ey,0.007+0.003*sin(frameCount*0.5+s),0.014);
	  }
	
	  // Nose cavity
	  fill(0); triangle(c.x-0.005,c.y,c.x+0.005,c.y,c.x,c.y+0.012);
	
	  // Black teeth etched into bone
	  stroke(0,0,0,230);
	  strokeWeight(0.0025);
	  for(let i=-4;i<=4;i++){
	    let jitter=0.001*sin(frameCount*0.4+i);
	    line(c.x+i*0.004,c.y+0.02,c.x+i*0.004+jitter,c.y+0.026);
	  }
	
	  // Subtle bone glow around jaw
	  noStroke(); fill(255,200,120,40);
	  ellipse(c.x,c.y+0.024,0.04,0.014);
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

// ONE LOUSY PARTICLE. Modify as needed.
// Nothing 👀 TODO here. 
class Particle {
	constructor(pRest, mass, radius) {
		this.pRest = pRest.copy();
		this.p = pRest.copy();
		this.v = vec2(0, 0);
		this.pin = false; // true if doesn't respond to forces
		this.mass = mass;
		this.radius = radius;
		this.f = vec2(0, 0);
		this.blob = null; // blob reference (if blob particle) -- useful for avoiding self-collision checks.
	}
	invMass() {
		return (this.pin ? 0.0 : 1.0 / this.mass);
	}
	// Emits a circle
	draw() {
		// nobody = (this.pin ? fill("red") : fill(0)); // default colors (red if pinned)
		circle(this.p.x, this.p.y, this.radius); //ellipseMode(RADIUS);
	}
}
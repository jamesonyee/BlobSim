// EDGE CLASS.
// Nothing 👀 TODO here. 
class Edge {
	// Creates edge spring of default stiffness, STIFFNESS_STRETCH
	constructor(particle0, particle1) {
		this.q = particle0;
		this.r = particle1;
		this.restLength = this.q.pRest.dist(this.r.pRest);
		this.stiffness = STIFFNESS_STRETCH;
	}
	// True if both particles are pinned
	isRigid() {
		return (this.q.pin && this.r.pin);
	}
	// Current length of edge spring
	length() {
		return this.q.p.dist(this.r.p);
	}
	midpoint() {
		return vadd(this.q.p, this.r.p).mult(0.5);
	}
	// Rest length of edge spring
	lengthRest() {
		return this.restLength;
	}
	// Draws the unstylized line 
	draw() {
		let a = this.q.p;
		let b = this.r.p;
		line(a.x, a.y, b.x, b.y);
	}
}

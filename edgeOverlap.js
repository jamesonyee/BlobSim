////////////////////////////////////////////////////
/// SUPPORTING CODE FOR EDGE-EDGE OVERLAP CHECKS /// 
////////////////////////////////////////////////////
// Nothing 👀 TODO here. 

// Efficiently checks that no pair of edges overlap, where the pairs do not share a particle in common.
function verifyNoEdgeEdgeOverlap() {
	if (detectedEdgeEdgeFailure) {
		return; // already done
	}

	// TODO: Optional: Make faster with broad phase
	// SIMPLE: Brute force check on edges i<j:
	for (let i = 0; i < edges.length - 1; i++) {
		let ei = edges[i];
		for (let j = i + 1; j < edges.length; j++) {
			let ej = edges[j];
			if (checkEdgeEdgeOverlap(ei, ej)) { // HALT!
				detectedEdgeEdgeFailure = true;
				isPaused = true;
				return;
			}
		}
	}
}

let overlapEdges = null;
/// IMPLEMENTATION (yeinj (F23))
function checkEdgeEdgeOverlap(ei, ej) {

	if (ei.isRigid() && ej.isRigid()) {
		return false;
	}

	// if (ei.blobId === ej.blobId) {
	// 	return false;
	// }

	let r = vsub(ei.r.p, ei.q.p);
	let s = vsub(ej.r.p, ej.q.p);
	let p = ei.q.p;
	let q = ej.q.p;

	let denominator = cross(r, s).z;
	let numerator1 = cross(vsub(q, p), s).z;
	let numerator2 = cross(vsub(q, p), r).z;

	if (denominator == 0) {
		let collinear = (numerator1 == 0 && numerator2 == 0);
		if (collinear) {
			let rdotr = dot(r, r);
			let t0 = dot(vsub(q, p), r) / rdotr;
			let t1 = t0 + dot(s, r) / rdotr;
			// are they in the same line seg?
			if (0 <= t0 && t0 <= 1 && 0 <= t1 && t1 <= 1) {
				overlapEdges = [ei, ej];
				return true;
			}
		} else {
			// colinear but disjoint.
			return false;
		}
	}

	// non-colinear case
	let t = numerator1 / denominator;
	let u = numerator2 / denominator;

	let eps = 0; //1e-4;

	if (eps < t && t < 1 - eps && eps < u && u < 1 - eps) {
		print("OVERLAP: ");
		print("  ei.r.p: " + ei.r.p + ", ei.q.p: " + ei.q.p);
		print("  ej.r.p: " + ej.r.p + ", ej.q.p: " + ej.q.p);
		print("  t: " + t + ", u: " + u);

		overlapEdges = [ei, ej];
		return true;
	}

	return false;
}

// Draws the overlapping edge pair (if any) found by verifyNoEdgeEdgeOverlap().
function drawOverlapEdges() {
	if (overlapEdges) {
		push();
		noFill();
		stroke("red");
		strokeWeight(0.01);
		let edgeShortest = (overlapEdges[0].length() < overlapEdges[1].length() ? overlapEdges[0] : overlapEdges[1]);
		let circleCenter = edgeShortest.midpoint();
		circle(circleCenter.x, circleCenter.y, 0.05); // Circle at midpoint of shortest edge (close enough)
		strokeWeight(0.0025);
		stroke("yellow");
		if (floor(frameCount / 4) % 2 == 0) /// Draw flashing yellow edges:
			for (let e of overlapEdges) e.draw();
		pop();
	}
}
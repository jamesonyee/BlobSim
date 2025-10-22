////////////////////////////////////////////////////
/// SUPPORTING CODE FOR MOUSE FORCE INTERACTIONS /// 
////////////////////////////////////////////////////
// Nothing 👀 TODO here. (but feel free to use the mouse for debugging/etc.)

// Blob currently being dragged by mouse forces, or undefined.
let mouseBlob;
// Selects closest blob for mouse forcing (mask out if using a GUI)
function mousePressed() {
	if (blobs.length == 0 || isPaused) return;

	// Set mouseBlob to blob closest to the mouse:
	let m = vec2(mouseX / CANVAS_SIZE, mouseY / CANVAS_SIZE);
	let minDist = 1000000000;
	let minCOM;
	let minBlob;
	for (let blob of blobs) {
		let com = blob.centerOfMass();
		if (com.dist(m) < minDist) {
			minDist = com.dist(m);
			minBlob = blob;
			minCOM = com;
		}
	}
	mouseBlob = minBlob;
}

function mouseReleased() {
	mouseBlob = undefined;
}

// Applies spring + damping force to all mouseBlob particles
function applyMouseForce() {
	if (mouseIsPressed && mouseBlob) {
		if (blobs.length < 1) return;
		let m = vec2(mouseX / CANVAS_SIZE, mouseY / CANVAS_SIZE);
		let blobCOM = mouseBlob.centerOfMass();
		let blobDist = blobCOM.dist(m);
		let mforce = vsub(m, blobCOM).normalize().mult(100.0 * clamp(blobDist, 0, 0.05));

		// Apply force to blob particles:
		let P = mouseBlob.blobParticles();
		for (let part of P) {
			part.f.add(mforce);
			vacc(part.f, -0.9, part.v); //some damping
		}
	}
}

// Draws line from the mouse to any forced mouseBlob
function drawMouseForce() {
	if (mouseIsPressed && mouseBlob) {
		if (blobs.length < 1) return;
		let m = vec2(mouseX / CANVAS_SIZE, mouseY / CANVAS_SIZE);
		let blobCOM = mouseBlob.centerOfMass();
		push();
		stroke(0);
		strokeWeight(0.003);
		line(m.x, m.y, blobCOM.x, blobCOM.y);
		pop();
	}
}
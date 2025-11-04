// A simple uniform 2D subdivision for interval-overlap tests on [0,width]x[0,height].
// Given an object and rectangular interval, [xmin,xmax]x[ymin,ymax], the object is stored in cells it overlaps.
// A test rectangle can be used to find objects with overlapping intervals.
// Doug L. James, Sept 2020, http://graphics.stanford.edu/~djames
class GridCheck2D {
	
	constructor(nXCells, width, nYCells, height) {
		this.nXCells       = nXCells;
		this.nYCells       = nYCells;
		this.nXCellsMinus1 = nXCells-1; // used often
		this.nYCellsMinus1 = nYCells-1; // used often
		this.dx     = width /this.nXCells;
		this.dy     = height/this.nYCells;
		this.invDX  = 1/this.dx; // used often
		this.invDY  = 1/this.dy; // used often
		this.width  = width;
		this.height = height;

		if(nXCells<1 || nYCells<1) 
			print("ERROR: GridCheck2D: nXCells & nYCells must be positive but were "+nXCells+" & "+nYCells);

		// Build cells:
		this.cells = [];// [i*nXCells + j] row-major form
		for(let i=0; i<nXCells; i++) {
			let xmin = i*this.dx;
			let xmax = xmin+this.dx;
			for(let j=0; j<nYCells; j++) {
				let ymin = j*this.dy;
				let ymax = ymin+this.dy;
				this.cells.push(new Cell2D());
			}
		}
	}
	
	// Clears all cells.
	clear() {
		for(let cell of this.cells) cell.clear();	
	}
	
	// Maps x to a cell index that is clamped to 0 to nXCells-1.
	getCellIndexX(x) {
		return Math.min(this.nXCellsMinus1, Math.max(0, Math.floor(x*this.invDX)));
	}
	// Maps y to a cell index that is clamped to 0 to nYCells-1.
	getCellIndexY(y) {
		return Math.min(this.nYCellsMinus1, Math.max(0, Math.floor(y*this.invDY)));
	}
	
	// Primary method to add an object on [xmin,xmax]x[ymin,ymax] to all overlapping-cell sets.
	addEntry(xmin, xmax, ymin, ymax, object) {
		let imin = this.getCellIndexX(xmin); 
		let imax = this.getCellIndexX(xmax);
		let jmin = this.getCellIndexY(ymin); 
		let jmax = this.getCellIndexY(ymax);
		for(let i=imin; i<=imax; i++) {
			for(let j=jmin; j<=jmax; j++) {
				let cell = this.cells[i*this.nYCells+j];
				cell.addEntry(object);
			}
		}
	}
	
	// Set of all objects in cells that overlap this interval.
	// NOTE: Contains no duplicates.
	getOverlappingCellObjects(xmin, xmax, ymin, ymax) {
		let imin = this.getCellIndexX(xmin); 
		let imax = this.getCellIndexX(xmax);
		let jmin = this.getCellIndexY(ymin); 
		let jmax = this.getCellIndexY(ymax);
		let result = new Set();
		for(let i=imin; i<=imax; i++) {
			for(let j=jmin; j<=jmax; j++) {
				let objs = this.cells[i*this.nYCells+j].getObjects();
				for(let obj of objs) 
					result.add(obj);
			}
		}
		return result;
	}
	
	// Returns a Set of objects in overlapping cells that return true for the testFunction(obj).
	// NOTE: Contains no duplicates.
	overlapsTest(xmin, xmax, ymin, ymax, testFunction) {
		let imin = this.getCellIndexX(xmin); 
		let imax = this.getCellIndexX(xmax);
		let jmin = this.getCellIndexY(ymin); 
		let jmax = this.getCellIndexY(ymax);
		let result = new Set();// ensures no duplicates
		// Apply test function to cells objects. 
		for(let i=imin; i<=imax; i++) {
			for(let j=jmin; j<=jmax; j++) { 
				let objs = this.cells[i*this.nYCells+j].overlapTest(testFunction);
				for(let obj of objs) 
					result.add(obj);
			}
		}
		return result;
	}

    populateGrid(blobs, environment) {
        // Clear the grid each frame
        this.clear();
        
        // Add all blob particles to the grid
        for (let blob of blobs) {
            for (let p of blob.BP) {
                // Calculate the AABB for this particle with padding
                const xmin = p.p.x - PARTICLE_RADIUS - AABB_PAD_PARTICLE;
                const xmax = p.p.x + PARTICLE_RADIUS + AABB_PAD_PARTICLE;
                const ymin = p.p.y - PARTICLE_RADIUS - AABB_PAD_PARTICLE;
                const ymax = p.p.y + PARTICLE_RADIUS + AABB_PAD_PARTICLE;
                
                this.addEntry(xmin, xmax, ymin, ymax, p);
            }
            
            // Add all blob edges to the grid
            for (let edge of blob.BE) {
                const edgeAABB = getEdgeAABB(edge);
                this.addEntry(
                    edgeAABB.min.x - AABB_PAD_EDGE, 
                    edgeAABB.max.x + AABB_PAD_EDGE,
                    edgeAABB.min.y - AABB_PAD_EDGE,
                    edgeAABB.max.y + AABB_PAD_EDGE,
                    edge
                );
            }
        }
        
        // Add environment edges to the grid
        let envEdges = environment.getEdges();
        for (let edge of envEdges) {
            const edgeAABB = getEdgeAABB(edge);
            this.addEntry(
                edgeAABB.min.x - AABB_PAD_EDGE,
                edgeAABB.max.x + AABB_PAD_EDGE,
                edgeAABB.min.y - AABB_PAD_EDGE,
                edgeAABB.max.y + AABB_PAD_EDGE,
                edge
            );
        }
    }
    
    // Optional: Helper method to get objects near a particle
    getObjectsNearParticle(particle, padding = AABB_PAD_PARTICLE) {
        const xmin = particle.p.x - PARTICLE_RADIUS - padding;
        const xmax = particle.p.x + PARTICLE_RADIUS + padding;
        const ymin = particle.p.y - PARTICLE_RADIUS - padding;
        const ymax = particle.p.y + PARTICLE_RADIUS + padding;
        
        return this.getOverlappingCellObjects(xmin, xmax, ymin, ymax);
    }
    
    // Optional: Helper method to get objects near an edge
    getObjectsNearEdge(edge, padding = AABB_PAD_EDGE) {
        const edgeAABB = getEdgeAABB(edge);
        const xmin = edgeAABB.min.x - padding;
        const xmax = edgeAABB.max.x + padding;
        const ymin = edgeAABB.min.y - padding;
        const ymax = edgeAABB.max.y + padding;
        
        return this.getOverlappingCellObjects(xmin, xmax, ymin, ymax);
    }

	// Draws interior cell boundaries
	draw() {
		/// DRAW CELL BOUNDARIES AS LINES:
		push();
		stroke(100,0,0);
        strokeWeight(0.001);
		for(let x=0; x<=this.width; x+=this.dx) {
			line(x,0, x,this.height);	
		}
		for(let y=0; y<=this.height; y+=this.dy) {
			line(0,y, this.width,y);
		}
		pop();
	}

	// Prints debugging info.
	debug() {
		let info = "";
		for(let k=0; k<this.cells.length; k++) {
			info = info + " " + this.cells[k].objects.size;
		}
		print("GridCheck2D: "+this.nXCells+"x"+this.nYCells+": cells: "+info);
	}

    drawDebug() {
        push();
        noStroke();
        
        for(let i=0; i<this.nXCells; i++) {
            for(let j=0; j<this.nYCells; j++) {
                let cell = this.cells[i*this.nYCells+j];
                let objectCount = cell.objects.size;
                
                if(objectCount > 0) {
                    // Red color with opacity based on object count
                    let alpha = map(objectCount, 0, 10, 50, 200, true);
                    fill(255, 0, 0, alpha);
                    rect(i*this.dx, j*this.dy, this.dx, this.dy);
                }
            }
        }
        pop();
    }
}

// Set implementation of objects in a 2D cell.
class Cell2D {
	constructor() {
		this.objects = new Set();
	}

	// Adds entry (assumed on this cell interval):
	addEntry(object){
		this.objects.add(object);	
	}

	// Returns reference to Set of all objects in this cell.	
	getObjects() {
		return this.objects; 
	}
	
	// Clears cell contents.
	clear() { this.objects.clear(); }
	
	// Returns list of all objects that return true for the testFunction(obj).
	overlapTest(testFunction) {
		let result = [];
		for(let obj of this.objects) {
			if(testFunction(obj)) // e.g., circle-circle overlap
				result.push(obj);
		}
		return result;
	}
}
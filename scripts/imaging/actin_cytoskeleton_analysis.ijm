#@ File (label="Select a directory", style="directory") myDir

// === GLOBAL VARIABLES ===
outputDirPath = "";
originalName = "";


// === MAIN EXECUTION ===
chosenDirectory = myDir;
if (!endsWith(chosenDirectory, File.separator))
    chosenDirectory = chosenDirectory + File.separator;

processFiles(chosenDirectory);


// === FUNCTION DEFINITIONS ===

/*
 * Recursively process all .tif images in the selected directory and its subfolders.
 * Creates an Analysis_<n> subfolder for output.
 */
function processFiles(directory) {

	fileList = getFileList(directory);
	
	outputDirName = directory + "Analysis";
	folderCount = 1;
	while (File.exists(outputDirName)) {
		outputDirName = directory + "Analysis_" + folderCount;
		folderCount++;
	}
	
	outputDirPath = outputDirName + File.separator;
	File.makeDirectory(outputDirName);

	for (file = 0; file < fileList.length; file++) {
		if (endsWith(fileList[file], ".tif")) {
			open(directory + fileList[file]);
			analyzeImage();	
		} else if (endsWith(fileList[file], "/") && !matches(fileList[file], ".*Analysis.*")) {
			processFiles(directory + fileList[file]);
		}
	}
}


/*
 * Main image analysis routine: nuclei, skeleton, tables, cleanup
 */
function analyzeImage() {
	originalName = getTitle();
	subFolder = outputDirPath + originalName + File.separator;
	File.makeDirectory(subFolder);
	
	Nuclei(subFolder);
	Skeleton(subFolder);
	FinalTable(subFolder);
	cleanUp();
}


/*
 * Segment and count nuclei
 */
function Nuclei(subFolder){
	run("Duplicate...", "duplicate channels=1");
	rename("Nuclei");
	run("Subtract Background...", "rolling=200 stack");
	run("Gaussian Blur...", "sigma=5 stack");
	run("Convert to Mask", "method=Otsu background=Dark calculate black");
	run("Fill Holes", "stack");
	run("Watershed", "stack");
	run("Analyze Particles...", "size=200-Infinity pixel show=Masks summarize stack");

	saveAs("Tiff", subFolder + "Mask_of_Nuclei_" + originalName + ".tif");

	if (isOpen("Summary of Nuclei")) {
		selectWindow("Summary of Nuclei"); 
		saveAs("Results", subFolder + "Nuclei_" + originalName + ".tsv");
		run("Close");
	}
}


/*
 * Segment phalloidin, skeletonize, measure branch data
 */
function Skeleton(subFolder) {
	selectWindow(originalName); 
	run("Duplicate...", "duplicate channels=2");
	rename("Phalloidin");
	run("Subtract Background...", "rolling=10 stack");
	setAutoThreshold("Otsu dark 16-bit stack");
	run("Convert to Mask", "method=Otsu background=Dark calculate black");
	run("Analyze Particles...", "size=50-Infinity pixel circularity=0.00-0.50 show=Masks stack");
	run("Gaussian Blur...", "sigma=1 stack");
	run("Convert to Mask", "method=Otsu background=Dark black");
	run("Skeletonize", "stack");

	saveAs("Tiff", subFolder + "Skeleton_" + originalName + ".tif");

	stackImg = getImageID();
	getDimensions(width, height, channels, slices, frames);

	BranchTot = newArray(slices);
	SkeletonLengths = newArray(slices);

	for (i = 1; i <= slices; i++) {
	    selectImage(stackImg);
	    setSlice(i);
	    run("Duplicate...", "");
	    sliceImg = getImageID();
	    run("Analyze Skeleton (2D/3D)", "prune=none");
	    saveAs("Results", subFolder + originalName + "_SkeletonLength_frame-" + i + ".tsv");
	    
	    if (isOpen("Results")) {
	    	selectWindow("Results");
	    	Table.rename("Results", "Results_" + i);
	    }

	    if (isOpen("Tagged skeleton")) {
	        selectWindow("Tagged skeleton");
	        rename("Tagged_skeleton_" + i);
	    }

	    selectImage(sliceImg);
	    close();
	}

	// Summarize branch and skeleton length per slice
	for (i = 1; i <= slices; i++) {
		Table.rename("Results_" + i, "Results");
		Skel = newArray(nResults);
		Branch = newArray(nResults);
		for(h = 0; h < nResults; h++) {
			Skel[h] = getResult("# Branches", h) * getResult("Average Branch Length", h);
			Branch[h] = getResult("# Branches", h);
		}
		skeletonSum = 0;
		branchSum = 0;
		for (j = 0; j < Skel.length; j++) {
			skeletonSum += Skel[j];
			branchSum += Branch[j];
		}
		SkeletonLengths[i-1] = skeletonSum;
		BranchTot[i-1] = branchSum;
		selectWindow("Results");
		run("Close");
	}

	Table.create("Lengths");
	Table.setColumn("Skeleton Lengths", SkeletonLengths);
	Table.setColumn("Branches", BranchTot);

	run("Merge Channels...", "c3=[Mask_of_Nuclei_" + originalName + "] c2=[Phalloidin] c1=[" + originalName + "] create ignore");
	saveAs("Tiff", subFolder + "Composite_" + originalName + ".tif");
}


/*
 * Combine nuclei counts and skeleton data into final summary table
 */
function FinalTable(subFolder) { 	
	Table.rename("Lengths", "Results");
	Lengths = newArray(nResults);
	BranchTotal = newArray(nResults);
	ABL = newArray(nResults);

	for(i=0; i<nResults; i++) {
		Lengths[i] = getResult("Skeleton Lengths", i);
		BranchTotal[i] = getResult("Branches", i);
		ABL[i] = getResult("Skeleton Lengths", i) / getResult("Branches", i);
	}
	selectWindow("Results");
	run("Close");

	// Load nuclei counts
	open(subFolder + "Nuclei_" + originalName + ".tsv");
	Table.rename("Results", "NucleiResults");

	Cells = newArray(nResults);
	for(i=0; i<nResults; i++) {
		Cells[i] = getResult("Count", i);
	}
	selectWindow("NucleiResults");
	run("Close");

	Table.create("FinalTable");
	Table.setColumn("Nuclei", Cells);
	Table.setColumn("Skeleton Lengths", Lengths);
	Table.setColumn("# Branches", BranchTotal);
	Table.setColumn("Average Branch Length", ABL);
	saveAs("Results", subFolder + "FinalTable_" + originalName + ".tsv");
}


/*
 * Close all open windows and tables
 */
function cleanUp() {
    if (isOpen("Results")) {
         selectWindow("Results"); 
         run("Close");
    }
    if (isOpen("Log")) {
         selectWindow("Log");
         run("Close");
    }
    while (nImages() > 0) {
          selectImage(nImages());  
          run("Close");
    }
    return;
}

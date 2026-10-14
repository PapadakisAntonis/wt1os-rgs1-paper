#@ Integer (label="Threshold channel 1 (Nucleus)", value="25") T1
#@ Integer (label="Prominence dots (lncRNA)", value="10") P1
#@ File (label="Select a directory", style="directory") myDir


//parameter SetUp
run("Set Measurements...", "area mean min perimeter integrated median display redirect=None decimal=3");


chosenDirectory = myDir + File.separator;

//Open the directory
processFiles(chosenDirectory);


/*
 * this function will create a directory named "analysis" inside the folder 
 * with the pictures and create the path to find that directory  
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
	saveAs("Measurements", outputDirPath +"Results.csv");
	run("Close");
}


//Change name of the figure with VSGV and Golgi
function normaliseNames(){
  run("Split Channels");
  selectWindow ("C2-"+originalName);
  rename (originalName + "- lncRNA");
  selectWindow ("C1-"+originalName);
  rename (originalName + "- Nucleus");
}

 



// Create a mask using the chosen threshold for the Golgi channel
function AreaIsolation(low, image){
	selectWindow(image);
	run("Gaussian Blur...", "sigma=5 slice");
    setThreshold(low, 255);
    setOption("BlackBackground", true);
    run("Convert to Mask");
}


//Create a ROI using the thrasholded Golgi channel 
function ROICalc(image){
    selectWindow(image);
    run("Create Selection");
    roiManager("add");
    roiManager("Select", 0);
    roiManager("Rename", " Nuclear");
    selectWindow(image);
    close();
	
}

//Dots calculator
function DotCalc(image, Prom) {
	selectWindow(image);
	run("Find Maxima...", "prominence=Prom output=Count");
	run("Find Maxima...", "prominence=Prom output=[Point Selection]");
	roiManager("Add");
	roiManager("Select", 1);
	roiManager("Rename", " Dots");
	roiManager("Select", 0);
    run("Find Maxima...", "prominence=Prom output=Count");
    roiManager("save", outputDirPath + originalName + ".zip");
    selectWindow(image);
    close();
    roiManager("reset");
}





function analyzeImage() {

originalName = getTitle ();
  

// Call first function
normaliseNames();

// Call second function
AreaIsolation(T1, originalName + "- Nucleus");

// Call third function and fourth
ROICalc(originalName + "- Nucleus");
DotCalc(originalName + "- lncRNA", P1);

}
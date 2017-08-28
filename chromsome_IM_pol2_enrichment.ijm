/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////  MACRO TO IDENTIFY AND CHARACTERISE INTERMINGLING REGIONS BETWEEN CHROMOSOMES                											                    /////////////
///////  WRITTEN BY: SARADHA VENKATACHALAPATHY                                                                                                                      /////////////
///////  DATE LAST EDITED: April 5th, 2017                                                                                                                          /////////////
///////  ASSUMPTIONS: The input image is a confocal zstack (slice height=0.5microns) with 4 channels, the channel One contains the nucleus and then there are two 
///////               channels for the chromosomes and one channel for RNA polymerase2. The channel numbers for the chromosomes and Pol2 are specified.             /////////////
///////  DESCRIPTION: Two dialog box opens where the user inputs the source (where the images are strored) the results(where results need to be saved) directory
///////               The program opens nucleus channel, uses the gaussian filter (sigma=0.15, scaled) to the stack and the thresholds the image and asks the user to
///////               check the thresholded image(optional). The user can then adjust the threshold for prescision. The program then converts the thresholded image to 
///////               a binary image and these images are stored in the results folder. For the chromsome channels, it identifies nuclear region and then it performs
///////               the aforementioned processes. The threshold values are stored in the datafiles directory. The binary chromsome images are opened in combination  
///////               (2-3) and passed through the AND filter to identify the intermingling regions (IMR) and the images saved in the results directory. Next, the 
///////				  volumes of the binary nucleus, chromosomes and IMRs are measured and saved in the datafiles directory. Then, the amount of pol2 in the total image, 
///////				  nucleus, and the IMR are calculated and stored. 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

nucleus=1;
chromosomeA=3;
chromosomeB=4;


dir=getDirectory("Please choose the source directory");
dir2=getDirectory("Please choose the results directory");
newDir1 = dir2 + "ChromosomeA_gaussian_blur" + File.separator;
newDir11 = dir2 + "ChromosomeA_binary" + File.separator;
newDir13 = dir2 + "ChromosomeA_raw" + File.separator;
newDir2 = dir2 + "ChromosomeB_gaussian_blur" + File.separator;
newDir21 = dir2 + "ChromosomeB_binary" + File.separator;
newDir23 = dir2 + "ChromosomeB_raw" + File.separator;
newDir9 = dir2 + "Ch1_nucleus_gaussian_blur" + File.separator;
newDir91 = dir2 + "Ch1_nucleus_binary" + File.separator;
newDir93 = dir2 + "Ch1_nucleus_raw" + File.separator;
newDir4 = dir2 + "ChrA_AND_ChrB_binary" + File.separator;
newDir7 = dir2 + "datafiles" + File.separator;

File.makeDirectory(newDir1); 
File.makeDirectory(newDir11); 
File.makeDirectory(newDir13); 
File.makeDirectory(newDir2); 
File.makeDirectory(newDir21); 
File.makeDirectory(newDir23); 
File.makeDirectory(newDir9); 
File.makeDirectory(newDir91); 
File.makeDirectory(newDir93); 
File.makeDirectory(newDir4); 
File.makeDirectory(newDir7);
setBatchMode(true);

filenames=getFileList(dir);
lower_ch1=newArray(filenames.length);
upper_ch1=newArray(filenames.length);
lower_ch2=newArray(filenames.length);
upper_ch2=newArray(filenames.length);
lower_ch3=newArray(filenames.length);
upper_ch3=newArray(filenames.length);
title=newArray(filenames.length);

for(i=0; i<filenames.length; i++){
	path=dir+filenames[i];

	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=nucleus c_end=nucleus c_step=1");
	imgName=getTitle(); 
	title[i]=getTitle(); 
	baseNameEnd=indexOf(imgName, ".nd2"); 
	baseName=substring(imgName, 0, baseNameEnd); 
	saveAs("tiff", newDir93+baseName);
	run("Gaussian Blur...", "sigma=0.15 scaled stack");
	saveAs("tiff", newDir9+baseName);
	setAutoThreshold("Otsu dark stack");
	setSlice(nSlices/2);
	//waitForUser("check please");
	getThreshold(lower_ch1[i], upper_ch1[i]);
	run("Convert to Mask", "method=Otsu background=Dark black");
	saveAs("tiff", newDir91+baseName);	
	run("Close All");

}

//NUCLEAR MASKS//
filenames=getFileList(dir);
filenames1=getFileList(newDir91);
for(i=0; i<filenames.length; i++){
	
	path1=newDir91+filenames1[i];
	open(path1);
	imgName1=getTitle(); 
	
	path=dir+filenames[i];
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=chromosomeA c_end=chromosomeA c_step=1");
	imgName2=getTitle(); 
	baseNameEnd=indexOf(imgName2,".nd2"); 
	baseName=substring(imgName2, 0, baseNameEnd); 
	saveAs("tiff", newDir13+baseName+"CH2");
	imgName2=getTitle(); 
	imageCalculator("AND create stack", imgName1,imgName2);
	saveAs("tiff", newDir13+baseName+"nuclear_mask");
	run("Gaussian Blur...", "sigma=0.15 scaled stack");
	saveAs("tiff", newDir1+baseName);
	setAutoThreshold("RenyiEntropy dark stack");
	setSlice(nSlices/2);
	//waitForUser("check please chromosomeA");
	getThreshold(lower_ch2[i], upper_ch2[i]);
	run("Convert to Mask", "method=RenyiEntropy background=Dark black");
	saveAs("tiff", newDir11+baseName);	
	run("Close All");


	path1=newDir91+filenames1[i];
	open(path1);
	imgName1=getTitle();
	
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=chromosomeB c_end=chromosomeB c_step=1");
	imgName2=getTitle(); 
	baseNameEnd=indexOf(imgName2, ".nd2"); 
	baseName=substring(imgName2, 0, baseNameEnd); 
	saveAs("tiff", newDir23+baseName+"CH3");
	imgName2=getTitle(); 
	imageCalculator("AND create stack", imgName1,imgName2);
	saveAs("tiff", newDir23+baseName+"nuclear_mask");
	run("Gaussian Blur...", "sigma=0.15 scaled stack");
	saveAs("tiff", newDir2+baseName);
	setAutoThreshold("RenyiEntropy dark stack");
	setSlice(nSlices/2);
	//waitForUser("check please chromsome B");
	getThreshold(lower_ch3[i], upper_ch3[i]);
	run("Convert to Mask", "method=RenyiEntropy background=Dark black");
	saveAs("tiff", newDir21+baseName);	
	run("Close All");

}


// INTERMINGLING //

//calculate intermingling for the binary channel(2-3)//
filenames1=getFileList(newDir11);
filenames2=getFileList(newDir21);

for(i=0; i<filenames1.length; i++){
	path1=newDir11+filenames1[i];
	open(path1);
	imgName1=getTitle(); 
	baseNameEnd1=indexOf(imgName1, ".tif"); 
	baseName1=substring(imgName1, 0, baseNameEnd1); 
	
	path2=newDir21+filenames2[i];
	open(path2);
	imgName2=getTitle(); 
	
	imageCalculator("AND create stack", imgName2,imgName1);
	waitForUser("ok?");
	saveAs("tiff", newDir4+baseName1);
	
	run("Z Project...", "projection=[Sum Slices]");
	setThreshold(1, 65025);
	run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis area_fraction stack  limit display redirect=None decimal=4");
	run("Measure");
	run("Close All");

}
saveAs("Results", newDir7+"intermingling_chr20_chr3.csv");
run("Close All");
run("Clear Results");

for(i=0; i<filenames1.length; i++){
	setResult("Title",i,title[i]);
	setResult("Nuclear_lower_threshold",i,lower_ch1[i]);
	setResult("Nuclear_upper_threshold",i,upper_ch1[i]);
	setResult("Chromosome1_lower_threshold",i,lower_ch2[i]);
	setResult("Chromosome1_upper_threshold",i,upper_ch2[i]);
	setResult("Chromosome2_lower_threshold",i,lower_ch3[i]);
	setResult("Chromosome2_upper_threshold",i,upper_ch3[i]);
	

}
saveAs("Results", newDir7+"thresholds.csv");


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Calculating the volumes of the nucleus, the chromsomes and the imtermingling regions


//Binary nucleus
list1 = getFileList(newDir91);
n= list1.length;
volume=newArray(n);
label=newArray(n);
for (i=0; i<list1.length; i++) {
		path = newDir91+list1[i];
		open(path);
		run("Clear Results");
		//run("Threshold...");
		setThreshold(1, 65025);
		run("Set Measurements...", "area limit display redirect=None decimal=4");
		run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing display stack");
		label[i]=File.getName(path);
		volume[i]=0;
		for (j=0; j<nResults; j++){
			volume[i]=volume[i]+getResult("Area",j)*0.5;
		}
		close();
}
run("Clear Results");
for (k=0; k<list1.length; k++) {
	setResult("Label",k ,label[k]);
	setResult("Volume",k, volume[k]);
}
saveAs("Results", newDir7+"binary_nucleus.csv");

//Binary chromsome 2
list1 = getFileList(newDir11);
n= list1.length;
volume=newArray(n);
label=newArray(n);
for (i=0; i<list1.length; i++) {
		path = newDir11+list1[i];
		open(path);
		run("Clear Results");
		//run("Threshold...");
		setThreshold(1, 65025);
		run("Set Measurements...", "area limit display redirect=None decimal=4");
		run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing display stack");
		label[i]=File.getName(path);
		volume[i]=0;
		for (j=0; j<nResults; j++){
			volume[i]=volume[i]+getResult("Area",j)*0.5;
		}
		close();
}
run("Clear Results");
for (k=0; k<list1.length; k++) {
	setResult("Label",k ,label[k]);
	setResult("Volume",k, volume[k]);
}
saveAs("Results", newDir7+"binary_chromsome_2.csv");
//Binary chromsome 3
list1 = getFileList(newDir21);
n= list1.length;
volume=newArray(n);
label=newArray(n);
for (i=0; i<list1.length; i++) {
		path = newDir21+list1[i];
		open(path);
		run("Clear Results");
		//run("Threshold...");
		setThreshold(1, 65025);
		run("Set Measurements...", "area limit display redirect=None decimal=4");
		run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing display stack");
		label[i]=File.getName(path);
		volume[i]=0;
		for (j=0; j<nResults; j++){
			volume[i]=volume[i]+getResult("Area",j)*0.5;
		}
		close();
}
run("Clear Results");
for (k=0; k<list1.length; k++) {
	setResult("Label",k ,label[k]);
	setResult("Volume",k, volume[k]);
}
saveAs("Results", newDir7+"binary_chromsome_3.csv");
//Binary IMR
list1 = getFileList(newDir4);
n= list1.length;
volume=newArray(n);
label=newArray(n);
for (i=0; i<list1.length; i++) {
		path = newDir4+list1[i];
		open(path);
		run("Clear Results");
		//run("Threshold...");
		setThreshold(1, 65025);
		run("Set Measurements...", "area limit display redirect=None decimal=4");
		run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing display stack");
		label[i]=File.getName(path);
		volume[i]=0;
		for (j=0; j<nResults; j++){
			volume[i]=volume[i]+getResult("Area",j)*0.5;
		}
		close();
}
run("Clear Results");
for (k=0; k<list1.length; k++) {
	setResult("Label",k ,label[k]);
	setResult("Volume",k, volume[k]);
}
saveAs("Results", newDir7+"binary_IMR.csv");


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Calculating the amount of pol2 in the nucleus and the IMR


pol2=2;

newDir1=dir2 + "pol2_tiff_raw" + File.separator;
newDir2=dir2 + "pol2_tiff_nuclear_mask" + File.separator;
newDir21=dir2 + "pol2_tiff_IMR_mask_binary" + File.separator;

File.makeDirectory(newDir1); 
File.makeDirectory(newDir2);
File.makeDirectory(newDir21);

//total_raw_pol2
filenames=getFileList(dir);
for(i=0; i<filenames.length; i++){
	path=dir+filenames[i];
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=pol2 c_end=pol2 c_step=1");
	imgName=getTitle(); 
	baseNameEnd=indexOf(imgName, ".nd2"); 
	baseName=substring(imgName, 0, baseNameEnd); 
	saveAs("tiff", newDir1+baseName);

	run("Z Project...", "projection=[Sum Slices]");
	run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis area_fraction stack display redirect=None decimal=4");
	run("Measure");
	run("Close All");	
}
saveAs("Results", newDir7+"pol2_raw.csv");
run("Close All");
run("Clear Results");

//total_pol2_in_nucleus_binary
filenames1=getFileList(newDir91);
filenames2=getFileList(newDir1);

for(i=0; i<filenames1.length; i++){
	path1=newDir91+filenames1[i];
	open(path1);
	imgName1=getTitle(); 
	baseNameEnd1=indexOf(imgName1, ".tif"); 
	baseName1=substring(imgName1, 0, baseNameEnd1); 
	
	path2=newDir1+filenames2[i];
	open(path2);
	imgName2=getTitle(); 
	
	imageCalculator("AND create stack", imgName2,imgName1);
	saveAs("tiff", newDir2+baseName1);
	
	run("Z Project...", "projection=[Sum Slices]");
	run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis area_fraction stack display redirect=None decimal=4");
	run("Measure");
	run("Close All");

}
saveAs("Results", newDir7+"pol2_nucelus_binary.csv");
run("Close All");
run("Clear Results");

//total_pol2_in_IMR_binary
filenames1=getFileList(newDir4);
filenames2=getFileList(newDir1);

for(i=0; i<filenames1.length; i++){
	path1=newDir4+filenames1[i];
	open(path1);
	imgName1=getTitle(); 
	baseNameEnd1=indexOf(imgName1, ".tif"); 
	baseName1=substring(imgName1, 0, baseNameEnd1); 
	
	path2=newDir1+filenames2[i];
	open(path2);
	imgName2=getTitle(); 
	
	imageCalculator("AND create stack", imgName2,imgName1);
	saveAs("tiff", newDir21+baseName1);
	
	run("Z Project...", "projection=[Sum Slices]");
	run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis area_fraction stack display redirect=None decimal=4");
	run("Measure");
	run("Close All");

}
saveAs("Results", newDir7+"pol2_IMR_binary.csv");
run("Close All");
run("Clear Results");



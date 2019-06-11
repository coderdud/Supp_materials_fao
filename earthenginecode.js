// define country boarders
var lesotho = ee.FeatureCollection('ft:1tdSwUL7MVpOauSgRzqVTOwdfy17KDbw-1d9omPw')
  .filter(ee.Filter.eq('Country', 'Lesotho'));
var geometry = lesotho.geometry();
var clip2country = function(inputIm){
  return inputIm.clipToCollection(lesotho);
};

// Load all Sentinel 2 images as ImageCollection
// limmit the images to Lesotho boarders
var selected_bands = ['B2','B3','B4','B5','B6','B7','B8','B8A',/*'B9',*/'B10','B11','B12','QA60'];
var LesothoImages = ee.ImageCollection('COPERNICUS/S2').select(selected_bands)
                      .filterBounds(lesotho.geometry())
                      .filterDate('2015-10-01','2018-01-28')       //fixed dates until end of my work period in FAO
                      .map(clip2country).sort('system:time_start');

// sort and select images with less than 20% cloud cover
var lesotho_less_cloud = LesothoImages
                        .filterMetadata('CLOUDY_PIXEL_PERCENTAGE','less_than',0.1);

//reduce image collection to one image for PCA analysis and show on map.
var lesotho_reduce_mean = lesotho_less_cloud.reduce(ee.Reducer.mean());



//define visualization parameters
var vizParams = {
  bands: ['B4', 'B3', 'B2'],
  min: 0,
  max: 3000,
  gamma: [1.6]
};


// visualize image tiles on map
Map.setCenter(27.6649, -29.4178,10);
/*Map.addLayer(lesotho_reduce_mean,{bands:['B4_mean','B3_mean','B2_mean'],
             min:0, max:3000,gamma:[1]},'RGB image National level');*/


//function to calculate EVI on an image
var eviOnCollection = function (inputImage){
  return inputImage.expression (
       '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
      'NIR': inputImage.select('B8'),
      'RED': inputImage.select('B4'),
      'BLUE': inputImage.select('B2')
       });
};


// function to calculate NDVI---------------------------------------------------
var ndviOnCollection = function (inputImage){
  return inputImage.normalizedDifference(['B8','B4']);
};


//PCA ANALYSIS SECTION-----------------------------------------------------

// initilaize image for pca analysis:
// Mean center the data to enable a faster covariance reducer
// and an SD stretch of the principal components.
var bandNames = lesotho_reduce_mean
    .select(['B8_mean','B4_mean','B3_mean','B2_mean']) ///'B5_mean','B6_mean','B7_mean','B8A_mean',,'B12_mean' /// specify bands which will be use for PCA analysis
    .bandNames();
var scale = 10;
var region = geometry;

//=------------------------
//RUN WHEN APPLY TO NEW country

//following code return the mean value for each band in bandNames for whole country
//it is time consuming and causing "computation time out" error in some cases
//therefore after successfully calculating we use the values directly instead of runing
// the code again

/*var meanDict = lesotho_reduce_mean.select(bandNames).reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region,
    scale: scale,
    maxPixels: 1e9
});//.aside(print);*/
//=--------------------------

// same values from last part. in case you apply to new regio on earth please run the 
// previous part instead.
var meanDict = ee.Dictionary({
  "B2_mean": 1010.4948406068576,
  "B3_mean": 994.271157638902,
  "B4_mean": 1151.1745636536798,
  "B8_mean": 1843.391798410535
});

// mean of each band in image
var means = ee.Image.constant(meanDict.values(bandNames)).clipToCollection(lesotho);


//subtract the mean value of each band from pixel values to mean sentering the image
var centered = lesotho_reduce_mean.select(bandNames).subtract(means).clipToCollection(lesotho);


// this funcion Mean center the data to enable a faster covariance reducer
// and an SD stretch of the principal components.
var getNewBandNames = function(prefix) {
  var seq = ee.List.sequence(1, bandNames.length());
  return seq.map(function(b) {
    return ee.String(prefix).cat(ee.Number(b).int());
  });
};
//*******************begin: FUNCTION Decleration************************
// This function accepts mean centered imagery, a scale and
// a region in which to perform the analysis.  It returns the
// Principal Components (PC) in the region as a new image.
var getPrincipalComponents = function(centered, scale, region) {
  // Collapse the bands of the image into a 1D array per pixel.
  var arrays = centered.toArray();
  
  
//=------------------------------ 
  //RUN WHEN APPLY NEW country
  
  //following code return the covariance value for each band in bandNames for whole country
  //it is time consuming and causing "computation time out" error in some cases
  //therefore after successfully calculating we use the values directly instead of runing
  //the code again
  
  // Compute the covariance of the bands within the region.
  /*var covar = arrays.reduceRegion({
    reducer: ee.Reducer.centeredCovariance(),
    geometry: region,
    scale: scale,
    maxPixels: 1e9
  });*/
//=-------------------------------  
  
  
  // FOR LESOTHO: same values from last part. in case you apply to new regio on earth please run the 
  // previous part instead.
  var covar = ee.Dictionary({
  "array": [
    [
      318408.4475468951,
      211104.85235648695,
      137335.41840500684,
      94909.8942295898
    ],
    [
      211104.85235648695,
      173019.05147897548,
      109767.29567174338,
      77495.3009322425
    ],
    [
      137335.41840500684,
      109767.29567174338,
      73366.56275666748,
      53144.75064914523
    ],
    [
      94909.8942295898,
      77495.3009322425,
      53144.75064914523,
      39933.70570829929
    ]
  ]
  });
  


  // Get the 'array' covariance result and cast to an array.
  // This represents the band-to-band covariance within the region.
  var covarArray = ee.Array(covar.get('array'));

  // Perform an eigen analysis and slice apart the values and vectors.
  var eigens = covarArray.eigen();

  // This is a P-length vector of Eigenvalues.
  var eigenValues = eigens.slice(1, 0, 1);
  // This is a PxP matrix with eigenvectors in rows.
  var eigenVectors = eigens.slice(1, 1);

  // Convert the array image to 2D arrays for matrix computations.
  var arrayImage = arrays.toArray(1);

  // Left multiply the image array by the matrix of eigenvectors.
  var principalComponents = ee.Image(eigenVectors).matrixMultiply(arrayImage);

  // Turn the square roots of the Eigenvalues into a P-band image.
  var sdImage = ee.Image(eigenValues.sqrt())
    .arrayProject([0]).arrayFlatten([getNewBandNames('sd')]);

  // Turn the PCs into a P-band image, normalized by SD.
  return principalComponents
    // Throw out an an unneeded dimension, [[]] -> [].
    .arrayProject([0])
    // Make the one band array image a multi-band image, [] -> image.
    .arrayFlatten([getNewBandNames('pc')])
    // Normalize the PCs by their SDs.
    .divide(sdImage);
};
//*******************END: FUNCTION Decleration************************

// =====--------------------------------------------------------------------------------------------------->
//calculate PCA component for images

// producws the PCA for country's image
var pca_oneImage = getPrincipalComponents(centered, 10 ,region)    //the region var is all the country's region
                   .clipToCollection(lesotho);


//print(pca_oneImage,'pca');
//Map.addLayer(pca_oneImage.select(['pc1']),{min: -2.5, max: 2.5, palette: ['ff0000','00ff00','0000ff']});// 2018-1-17,16:45 it worked perfectly!! //modified on 2018-09-23





//VISUALIZING SHARPEN IMAGE
// sharpening the image by pca components 
                                                              //coefficiens can be modify 
var sharp1 = lesotho_reduce_mean.add(pca_oneImage.select(['pc1']).multiply(200));
Map.addLayer(sharp1.select(['B4_mean','B3_mean','B2_mean']),{min:0,max:3000}, 'sharp image - PC1');

//var sharp2 = lesotho_reduce_mean.add(pca_oneImage.select(['pc2']).multiply(400));
//Map.addLayer(sharp2.select(['B4_mean','B3_mean','B2_mean']),{min:0,max:3000});

//var sharp3 = lesotho_reduce_mean.add(pca_oneImage.select(['pc3']).multiply(100));
//Map.addLayer(sharp3.select(['B4_mean','B3_mean','B2_mean']),{min:0,max:3000});

//var sharp4 = lesotho_reduce_mean.add(pca_oneImage.select(['pc4']).multiply(150));
//Map.addLayer(sharp4.select(['B4_mean','B3_mean','B2_mean']),{min:0,max:3000});

//------------------------------------------------------------------------------------------------------------------
// Compute the gray-level co-occurrence matrix (GLCM), 

//claculate GLCM of PCA in 300 times 
//var glcm_pca = pca_oneImage.multiply(300).toInt().glcmTexture({size: 8});  // return glcm 
//print(glcm_pca,'glcm of pca');
/*var contrast = glcm_pca.select('pc2_contrast');
Map.addLayer(contrast,
             {min: 0, max: 1500, palette: ['0000CC', 'CC0000']},
             'contrast');
*/


// GLCM of sharpen images
var onBands = ['B2_mean','B3_mean','B4_mean',/*'B5_mean',         // specify bands which should be use for GLCM analysis
               'B6_mean','B7_mean',*/'B8_mean',/*'B8A_mean',
               'B10_mean','B11_mean',*/'B12_mean'];



// extract glcm texture properties for each band of pca sharpen image 
var glcm_sharp_1 = sharp1.select(onBands).toInt().glcmTexture({size:8});
//var glcm_sharp_2 = sharp2.select(onBands).toInt().glcmTexture({size:8}).aside(print);
/*var glcm_sharp_3 = sharp3.select(onBands).toInt().glcmTexture({size:8});      deactivated on 2018-09-23
var glcm_sharp_4 = sharp4.select(onBands).toInt().glcmTexture({size:8});*/

// visualize contrast of one band glcm for reference!
//var contrast = glcm_sharp_1.select('B8_mean_contrast');
/*Map.addLayer(contrast,
             {min: 0, max: 1500, palette: ['0000CC', 'CC0000']},
             'contrast');
*/

 

//<__________________________________spectral features _________________________________________________________>>

// define country boarders
//var geometry = lesotho.geometry();  
// Load all Sentinel 2 images as ImageCollection and sor by time
// filter more than 35% cloudy images


// cloud function to remove clouds
var cloudfunction_ST2 = function(image){
  //use add the cloud likelihood band to the image
  var quality = image.select("QA60").unmask();
  //get pixels above the threshold
  var cloud01 = quality.gt(3276);       //1000 out of 65535 means 1.5% (0.0152590219) cloudy pixels     3276=5%
  //create a mask from high likelihood pixels
  var cloudmask = image.mask().and(cloud01.not());
  //mask those pixels from the image
  return image.updateMask(cloudmask);
};


//function to calculate EVI on an image
var eviOnImage = function (inputImage){
  return inputImage.expression (
       '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
      'NIR': inputImage.select('B8_mean'),
      'RED': inputImage.select('B4_mean'),
      'BLUE': inputImage.select('B2_mean')
})
};

// following process repeat for all 4 season ------------------------------------------------------------------------------------>


// October to end of December-------------------------------------------------------------------------------------------------------------
//NDVI
var visual_param = {
                    bands: 'nd',
                    min: -0.25,
                    max: 0.4, 
                    palette: ['FF0000','000000', '00FF00']
};

var visual_param_evi = {
                    bands: 'constant',
                    min: -5,
                    max: 5, 
                    palette: ['FF0000','000000', '00FF00']
};


var OD2015ndvi = LesothoImages.filterDate('2015-10-01','2015-12-31')  //filter collection to select images of season
    .map(cloudfunction_ST2)      //get ride of cloudy pixels
    .reduce(ee.Reducer.mean())   // reduce to one image
    .normalizedDifference(['B8_mean','B4_mean']);   //calculate NDVI for season
var OD2016ndvi = LesothoImages.filterDate('2016-10-01','2016-12-31')
    .map(cloudfunction_ST2)
    .reduce(ee.Reducer.mean())
    .normalizedDifference(['B8_mean','B4_mean']);
var OD2017ndvi = LesothoImages.filterDate('2017-10-01','2017-12-31')
    .map(cloudfunction_ST2)
    .reduce(ee.Reducer.mean())
    .normalizedDifference(['B8_mean','B4_mean']);
//var OD2017ndvi = LesothoImages.filterDate('2017-10-01','2017-12-31')
//    .reduce(ee.Reducer.mean())
//    .normalizedDifference(['B12_mean','B8_mean']);       //NDVI in case of strong atmospheric disturbance based on: 
                                                           //https://www.indexdatabase.de/db/si-single.php?rsindx_d=59=&sensor_id=96



var ODndvi = ee.Image(OD2015ndvi)
                .addBands(OD2016ndvi)   // add NDVI of different years in same season as bands to one image 
                .addBands(OD2017ndvi)
                .reduce(ee.Reducer.mean());  // calculate mean of NDVI for season */

//EVI
var OD2015evi = eviOnImage(LesothoImages.filterDate('2015-10-01','2015-12-31')   // select images of season
                .map(cloudfunction_ST2)   // get ride of cloudy pixels 
                .reduce(ee.Reducer.mean()));  // reduce to one image and calculate EVI
var OD2016evi = eviOnImage(LesothoImages.filterDate('2016-10-01','2016-12-31')
                .map(cloudfunction_ST2)
                .reduce(ee.Reducer.mean()));
var OD2017evi = eviOnImage(LesothoImages.filterDate('2017-10-01','2017-12-31')
                .map(cloudfunction_ST2)
                .reduce(ee.Reducer.mean()));

var ODevi = ee.Image(OD2015evi)
              .addBands(OD2016evi)
              .addBands(OD2017evi)
              .reduce(ee.Reducer.mean());

//January to end of March ---------------------------------------------------------------------------------------------------
//NDVI
var JM2016ndvi = LesothoImages.filterDate('2016-01-01','2016-03-31')
                      .map(cloudfunction_ST2)
                      .reduce(ee.Reducer.mean())
                      .normalizedDifference(['B8_mean','B4_mean']);
var JM2017ndvi = LesothoImages.filterDate('2017-01-01','2017-03-31')
                      .map(cloudfunction_ST2)
                      .reduce(ee.Reducer.mean())
                      .normalizedDifference(['B8_mean','B4_mean']);

var JMndvi = ee.Image(JM2016ndvi)
               .addBands(JM2017ndvi)
               .reduce(ee.Reducer.mean());

//EVI
var JM2016evi = eviOnImage(LesothoImages.filterDate('2016-01-01','2016-03-31')
                      .map(cloudfunction_ST2)
                      .reduce(ee.Reducer.mean()));
var JM2017evi = eviOnImage(LesothoImages.filterDate('2017-01-01','2017-03-31')
                      .map(cloudfunction_ST2)
                      .reduce(ee.Reducer.mean()));

var JMevi = ee.Image(JM2016evi)
              .addBands(JM2017evi)
              .reduce(ee.Reducer.mean());

// April to end of June-----------------------------------------------------------------------------------------------------------------
//NDVI
var AJ2016ndvi = LesothoImages.filterDate('2016-04-01','2016-06-30')
                      .map(cloudfunction_ST2)
                      .reduce(ee.Reducer.mean())
                      .normalizedDifference(['B8_mean','B4_mean']);
var AJ2017ndvi = LesothoImages.filterDate('2017-04-01','2017-06-30')
                      .map(cloudfunction_ST2)
                      .reduce(ee.Reducer.mean())
                      .normalizedDifference(['B8_mean','B4_mean']);

var AJndvi = ee.Image(AJ2016ndvi)
               .addBands(AJ2017ndvi)
               .reduce(ee.Reducer.mean());

//EVI
var AJ2016evi = eviOnImage(LesothoImages.filterDate('2016-04-01','2016-06-30')
                      .map(cloudfunction_ST2)
                      .reduce(ee.Reducer.mean()));
var AJ2017evi = eviOnImage(LesothoImages.filterDate('2017-04-01','2017-06-30')
                      .map(cloudfunction_ST2)
                      .reduce(ee.Reducer.mean()));

var AJevi = ee.Image(AJ2016evi)
              .addBands(AJ2017evi)
              .reduce(ee.Reducer.mean());



// July to end of September--------------------------------------------------------------------------------------------------
//NDVI
var JS2016ndvi = LesothoImages.filterDate('2016-07-01','2016-09-30')
                      .map(cloudfunction_ST2)
                      .reduce(ee.Reducer.mean())
                      .normalizedDifference(['B8_mean','B4_mean']);
var JS2017ndvi = LesothoImages.filterDate('2017-07-01','2017-09-30')
                      .map(cloudfunction_ST2)
                      .reduce(ee.Reducer.mean())
                      .normalizedDifference(['B8_mean','B4_mean']);

var JSndvi = ee.Image(JS2016ndvi)
               .addBands(JS2017ndvi)
               .reduce(ee.Reducer.mean());

//EVI
var JS2016evi = eviOnImage(LesothoImages.filterDate('2016-07-01','2016-09-30')
                      .map(cloudfunction_ST2)
                      .reduce(ee.Reducer.mean()));
var JS2017evi = eviOnImage(LesothoImages.filterDate('2017-07-01','2017-09-30')
                      .map(cloudfunction_ST2)
                      .reduce(ee.Reducer.mean()));

var JSevi = ee.Image(JS2016evi)
              .addBands(JS2017evi)
              .reduce(ee.Reducer.mean());

//===------------------------------------------------------------------------------------------------------------------

var visual_param = {
                    bands: 'mean',
                    min: -0.25,
                    max: 0.4, 
                    palette: ['FF0000','000000', '00FF00']
};

var visual_param_evi = {
                    bands: 'mean',
                    min: -5,
                    max: 5, 
                    palette: ['FF0000','000000', '00FF00']
};



/*Map.addLayer(LesothoImages.filterDate('2017-10-01','2017-12-31').map(cloudfunction_ST2).reduce(ee.Reducer.mean())
.select('B4_mean','B3_mean','B2_mean'),{min:0,max:3000,gamma:[1.6]},'OD');*/

//Map.addLayer(ODndvi,visual_param ,'ndvi OD');
//Map.addLayer(ODevi,visual_param_evi,'evi OD')

//==----------------
//Map.addLayer(LesothoImages.filterDate('2017-01-01','2017-03-31').map(cloudfunction_ST2).reduce(ee.Reducer.mean())
//.select('B4_mean','B3_mean','B2_mean'),{min:0,max:3000,gamma:[1.6]},'JM'); 

//Map.addLayer(JMndvi,visual_param,'ndvi JM')
//Map.addLayer(JMevi,visual_param_evi,'evi JM')


//==----------------
//Map.addLayer(LesothoImages.filterDate('2017-04-01','2017-06-30').map(cloudfunction_ST2).reduce(ee.Reducer.mean())
//.select('B4_mean','B3_mean','B2_mean'),{min:0,max:3000,gamma:[1.6]},'AJ'); 

//Map.addLayer(AJndvi,visual_param,'ndvi AJ')
//Map.addLayer(AJevi,visual_param_evi,'evi AJ')

//==----------------
//Map.addLayer(LesothoImages.filterDate('2017-07-01','2017-09-30').map(cloudfunction_ST2).reduce(ee.Reducer.mean())
//.select('B4_mean','B3_mean','B2_mean'),{min:0,max:3000,gamma:[1.6]},'JS'); 

//Map.addLayer(JSndvi,visual_param,'ndvi JS')
//Map.addLayer(JSevi,visual_param_evi,'evi JS')  */



// =====------------------------------------------------------------->
// make class band from impoerted raster data

var palette = [
  'FFFFFF', //          white
  '0001FB', // BUILT-UP           =1
  '0081FF', // AGRICULTURE        =2
  '04FDFF', // TREES              =3
  '7CFE88', // HYDROLOGY,wetlands =4
  'FFFF01', // SHRUBLAND          =5
  'FF8000', // GRASSLAND          =6
  'FE0000', // BARREN LAND        =7
];

var klass = ee.ImageCollection([east,north,south,west]).mosaic().rename(['classBand']);
Map.addLayer(klass,{min:0,max:7,palette: palette},'FAO LC');


//make one imege include all features as bands



var features = ee.Image(ODndvi).addBands(JMndvi).addBands(AJndvi).addBands(JSndvi)    // NDVIs
                 .addBands(ODevi).addBands(JMevi).addBands(AJevi).addBands(JSevi)     // EVIs
                 /*.addBands(glcm_pca)*/       // GLCM of PCA.select(1,19)
                 .addBands(glcm_sharp_1)
                 //.addBands(glcm_sharp_2)     //GLCM of sharpened Image
                 //.addBands(glcm_sharp_3).addBands(glcm_sharp_4)
                 .addBands(klass).toFloat().aside(print);






///  =====__________________________________________________Export data__________________________________________________________________>
// =======-------------------------------------------------------------------------------------------------->
// export feature vector to asset in order to process in other script file 
//to privent run out of memory or timeout problems

Export.image.toDrive({
  image: features,
  description: 'exported_feature_vector',  //give the file name
  scale: 10,
  region: forExport,   // enter the polygon name te be downloaded
  maxPixels: 1e10
});




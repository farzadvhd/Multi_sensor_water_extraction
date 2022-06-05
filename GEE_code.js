var ID = 436; // Wetland ID 
var yer = 2020;
var month_ = 6;
// day1 and day2 specifies the 6-day period
var day1 = 18;
var day2 = 24;
var bnd = 'VV_filt'
var orbit = 'ASCENDING'
Map.centerObject(aoi,10)

// LC is 10m high-resolution national land cover of Sweden
//relevant classess has been extracted for further analysis
var LC = ee.Image(0).where(
  LC.eq(2)
    .or(LC.eq(121))
     .or(LC.eq(122))
      .or(LC.eq(123))
       .or(LC.eq(124))
        .or(LC.eq(125))
         .or(LC.eq(126))
          .or(LC.eq(127))
           .or(LC.eq(128))
             .or(LC.eq(61))
              .or(LC.eq(62))
    , 1)

// PeronaMalikFilter (https://mygeoblog.com/2021/01/22/perona-malik-filter/)
function PeronaMalikFilter(img) {
  
  var K = 3.5;
  var iter = 10;
  var method = 2;
  
  var dxW = ee.Kernel.fixed(3, 3,
                           [[ 0,  0,  0],
                            [ 1, -1,  0],
                            [ 0,  0,  0]]);
  
  var dxE = ee.Kernel.fixed(3, 3,
                           [[ 0,  0,  0],
                            [ 0, -1,  1],
                            [ 0,  0,  0]]);
  
  var dyN = ee.Kernel.fixed(3, 3,
                           [[ 0,  1,  0],
                            [ 0, -1,  0],
                            [ 0,  0,  0]]);
  
  var dyS = ee.Kernel.fixed(3, 3,
                           [[ 0,  0,  0],
                            [ 0, -1,  0],
                            [ 0,  1,  0]]);

  var lambda = 0.2;

  var k1 = ee.Image(-1.0/K);
  var k2 = ee.Image(K).multiply(ee.Image(K));

  for(var i = 0; i < iter; i++) {
    var dI_W = img.convolve(dxW)
    var dI_E = img.convolve(dxE)
    var dI_N = img.convolve(dyN)
    var dI_S = img.convolve(dyS)

    switch(method) {
      case 1:
        var cW = dI_W.multiply(dI_W).multiply(k1).exp();
        var cE = dI_E.multiply(dI_E).multiply(k1).exp();
        var cN = dI_N.multiply(dI_N).multiply(k1).exp();
        var cS = dI_S.multiply(dI_S).multiply(k1).exp();
    
        img = img.add(ee.Image(lambda).multiply(cN.multiply(dI_N).add(cS.multiply(dI_S)).add(cE.multiply(dI_E)).add(cW.multiply(dI_W))))
        break;
      case 2:
        var cW = ee.Image(1.0).divide(ee.Image(1.0).add(dI_W.multiply(dI_W).divide(k2)));
        var cE = ee.Image(1.0).divide(ee.Image(1.0).add(dI_E.multiply(dI_E).divide(k2)));
        var cN = ee.Image(1.0).divide(ee.Image(1.0).add(dI_N.multiply(dI_N).divide(k2)));
        var cS = ee.Image(1.0).divide(ee.Image(1.0).add(dI_S.multiply(dI_S).divide(k2)));
    
        img = img.add(ee.Image(lambda).multiply(cN.multiply(dI_N).add(cS.multiply(dI_S)).add(cE.multiply(dI_E)).add(cW.multiply(dI_W))))
        break;
    }
  }

  return img;
}

// Apply angle correction (for VV and VH)
function toGammaVV(image) {
      return image.addBands(image.select('VV').subtract(image.select('angle')
        .multiply(Math.PI/180.0).cos().log10().multiply(10.0)).rename('VV_corr'));
    }
function toGammaVH(image) {
      return image.addBands(image.select('VH').subtract(image.select('angle')
        .multiply(Math.PI/180.0).cos().log10().multiply(10.0)).rename('VH_corr'));
    }
// table contains wetland shape files
var RMS = table.filterBounds(aoi).geometry();

//Loading Sentinel1 images
var S1 = ee.ImageCollection('COPERNICUS/S1_GRD')
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
    .filter(ee.Filter.calendarRange(yer,yer,'year'))
    .filter(ee.Filter.calendarRange(month_,month_,'month'))
    .filter(ee.Filter.calendarRange(day1,day2,'day_of_month'))
//    .filter(ee.Filter.calendarRange(152,158,'day_of_year'))
    .filter(ee.Filter.eq('orbitProperties_pass', orbit)) //Filter to get images from different look angles
    .filter(ee.Filter.contains({leftField: ".geo", rightValue: RMS}))
    .filterBounds(aoi)  
print(S1,'S1')

//Applying filters 
var S1 = S1.map(toGammaVV)
    .map(toGammaVH)
    .map(PeronaMalikFilter)
    .mean()
    .select(['VV_corr','VH_corr']).rename(['VV_filt','VH_filt'])

Map.addLayer(S1.select('VV_filt').clip(aoi),{min: -30, max: 5},'S1_VV');

//creating vv-vh
var vv_vh = S1.select('VV_filt').subtract(S1.select('VH_filt')).clip(aoi)
Map.addLayer(vv_vh,{min: 5, max: 15},'vv-vh');

// CC is Interferometric coherence layer
var CC = CC.double()
Map.addLayer(CC,{min: 0.5, max: 1}, 'coherence')
print(CC,'coherence')

// sentinel2 image
function maskS2clouds(image) {
  var qa = image.select('QA60');
  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return image.updateMask(mask).divide(10000);
}

var s2 = ee.ImageCollection("COPERNICUS/S2")
          .filterBounds(aoi)
          .filter(ee.Filter.calendarRange(yer,yer,'year'))
          .filter(ee.Filter.calendarRange(month_,month_,'month'))
          .filter(ee.Filter.calendarRange(day1,day2,'day_of_month'))
//          .filter(ee.Filter.calendarRange(152,158,'day_of_year'))
          .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',5))
          .map(maskS2clouds)
          
var s2_lay = s2.mean()
          .clip(aoi)
var rgbVis = {min: 0.0, max: 0.3 , bands: ['B4', 'B3', 'B2']};
print(s2)
Map.addLayer(s2_lay, rgbVis, 's2');
// ndvi layer is created to help vegetation detection
var ndvi = s2.map(function(img){
            return img.normalizedDifference(['B8','B4']).rename('ndvi');
          })
          .mean()
          .clip(aoi)
        
var ndviViz = {min: -1, max: 1, palette: ['blue', 'white', 'green']};
var ndvi = ndvi.toDouble()
Map.addLayer(ndvi, ndviViz, 'NDVI');

var blank = ee.Image.constant(0)
print(blank)

// classification
var input = S1.select('VV_filt').addBands(ndvi).addBands(vv_vh).addBands(CC)
print(input,'input')

// Define the visualization parameters.
var vizParams = {
  bands: ['VV_filt', 'ndvi','VV_filt_1'],
    min: [-30,-1,10],
  max: [5,1,20],
  //gamma: [1.1, 1,1.2]
};

Map.addLayer(input, vizParams, 'composite');

// Make the training dataset.
var training = input.sample({
  region: aoi,
  scale: 10,
  numPixels: 10000
});

// Instantiate the clusterer and train it.
var clusterer = ee.Clusterer.wekaXMeans({
  minClusters: 15,
  maxClusters: 30
}).train(training);

// Cluster the input using the trained clusterer.
var result = input.cluster(clusterer).clip(aoi);

// Display the clusters with random colors.
Map.addLayer(result.randomVisualizer(), {}, 'clusters_cc');

//manually grouping pixels into four classes
var reclassified = ee.Image(0)
.where(result.eq(13), 1)
 .where(result.eq(10), 2)
    .where(result.eq(1), 3)
      .where(result.eq(0), 4)
    .clip(aoi);
    
var final = reclassified.mask(LC).unmask(0)

Map.addLayer(final.clip(RMS).randomVisualizer(), {}, 'final');
Map.addLayer(table, {color: 'red', strokeWidth: 10}, 'ramsar_inland_sites');

Export.image.toDrive({
  image: input.clip(RMS),
  scale: 10,
  region: RMS,
  fileNamePrefix: ID+'_composit_'+yer+'_'+month_+'_'+day1+'_'+day2
  
})

Export.image.toDrive({
  image: final.clip(RMS),
  scale: 10,
  region: RMS,
  fileNamePrefix: ID+'_final_'+yer+'_'+month_+'_'+day1+'_'+day2
  
})

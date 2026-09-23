// ============================================================================
// INPUT
// ============================================================================

var primarySchools = ee.FeatureCollection("projects/ee-juvemec/assets/EducationP");

var bangladesh = ee.Geometry.Rectangle([88.0, 20.5, 92.7, 26.7], null, false);

var school = ee.FeatureCollection("projects/ee-juvemec/assets/EducationP")
  .filter(ee.Filter.or(
    ee.Filter.eq("F_Type", "Primary School"),
    ee.Filter.eq("F_Type", "Pri mary School")
  ))
  .filterBounds(bangladesh);

print("Valid Primary Schools:", school.size());

print("Total Primary Schools:", school.size());
print("First 10 Primary Schools:", school.limit(10));

var startDate = "2026-01-01";
var endDate = "2027-01-01";

var bufferRadiusMeters = 100;
var chunkSize = 2000;


// ============================================================================
// BANGLADESH PROCESSING EXTENT
// ============================================================================


Map.addLayer(school, {color: "red"}, "Primary Schools");
Map.centerObject(school, 7);


// ============================================================================
// SENTINEL-2
// ============================================================================

function maskS2clouds(image) {
  image = ee.Image(image);
  var qa = image.select("QA60");
  var mask = qa.bitwiseAnd(1 << 10).eq(0).and(qa.bitwiseAnd(1 << 11).eq(0));
  return image.updateMask(mask).divide(10000).copyProperties(image, ["system:time_start"]);
}

var s2Col = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
  .filterDate(startDate, endDate)
  .filterBounds(bangladesh)
  .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE",5))
  .map(maskS2clouds);

print("Sentinel-2 Image Count:", s2Col.size());

var s2Median = s2Col.median();

var blue = s2Median.select("B2");
var green = s2Median.select("B3");
var red = s2Median.select("B4");
var nir = s2Median.select("B8");
var swir1 = s2Median.select("B11");

var ndvi = nir.subtract(red).divide(nir.add(red)).rename("NDVI");
var evi = nir.subtract(red).multiply(2.5).divide(nir.add(red.multiply(6)).subtract(blue.multiply(7.5)).add(1)).rename("EVI");
var ndmi = nir.subtract(swir1).divide(nir.add(swir1)).rename("NDMI");
var ndwi = green.subtract(nir).divide(green.add(nir)).rename("NDWI");
var mndwi = green.subtract(swir1).divide(green.add(swir1)).rename("MNDWI");
var ndbi = swir1.subtract(nir).divide(swir1.add(nir)).rename("NDBI");
var bsi = swir1.add(red).subtract(nir.add(blue)).divide(swir1.add(red).add(nir).add(blue)).rename("BSI");


/// ==================== LANDSAT 8 LST ====================
var l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
  .filterDate(startDate, endDate)
  .filterBounds(bangladesh)
  .filter(ee.Filter.lt("CLOUD_COVER", 20));

var lstCol = l8.map(function(img) {
  var image = ee.Image(img);
  var qa = image.select("QA_PIXEL");

  // Cloud, dilated cloud, cirrus, cloud shadow and snow mask
  var mask = qa.bitwiseAnd(1 << 1).eq(0)
    .and(qa.bitwiseAnd(1 << 2).eq(0))
    .and(qa.bitwiseAnd(1 << 3).eq(0))
    .and(qa.bitwiseAnd(1 << 4).eq(0))
    .and(qa.bitwiseAnd(1 << 5).eq(0));

  // Landsat Collection 2 Level-2 Surface Temperature
  var lst = image.select("ST_B10")
    .multiply(0.00341802)
    .add(149.0)
    .subtract(273.15)
    .rename("LST");

  return lst
    .updateMask(mask)
    .copyProperties(image, ["system:time_start"]);
});

// Annual mean LST
var lstMean = lstCol.mean().rename("LST_Mean");

// 90th percentile LST
var lstP90 = lstCol
  .reduce(ee.Reducer.percentile([90]))
  .rename("LST_P90");

print("Landsat 8 images:", l8.size());
print("LST images:", lstCol.size());
print("LST bands:", lstCol.first().bandNames());


// ============================================================================
// ERA5 EXTREME HEAT
// ============================================================================

var era5Daily = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")
  .filterDate(startDate, endDate)
  .filterBounds(bangladesh)
  .select("temperature_2m_max");

var extremeHeatDays = era5Daily.map(function(image) {
  return ee.Image(image).subtract(273.15).gte(35).rename("heat_day");
}).sum().rename("Extreme_Heat_Days");


// ============================================================================
// SENTINEL-5P NO2
// ============================================================================

var no2 = ee.ImageCollection("COPERNICUS/S5P/OFFL/L3_NO2")
  .filterDate(startDate, endDate)
  .filterBounds(bangladesh)
  .select("NO2_column_number_density")
  .mean()
  .multiply(1000000)
  .rename("NO2_umol_m2");


// ============================================================================
// SENTINEL-5P O3
// ============================================================================

var o3 = ee.ImageCollection("COPERNICUS/S5P/OFFL/L3_O3")
  .filterDate(startDate, endDate)
  .filterBounds(bangladesh)
  .select("O3_column_number_density")
  .mean()
  .multiply(2241.15)
  .rename("O3_DU");


// ============================================================================
// MODIS AOD
// ============================================================================

var aod = ee.ImageCollection("MODIS/061/MCD19A2_GRANULES")
  .filterDate(startDate, endDate)
  .filterBounds(bangladesh)
  .select("Optical_Depth_047")
  .mean()
  .multiply(0.001)
  .rename("AOD");


// ============================================================================
// DYNAMIC WORLD
// ============================================================================

var dw = ee.ImageCollection("GOOGLE/DYNAMICWORLD/V1")
  .filterDate(startDate, endDate)
  .filterBounds(bangladesh)
  .select(["trees", "water", "built"])
  .mean();

var treeFrac = dw.select("trees").rename("Tree_Fraction");
var blueFrac = dw.select("water").rename("Blue_Fraction");
var impFrac = dw.select("built").rename("Impervious_Fraction");


// ============================================================================
// META CANOPY HEIGHT
// ============================================================================

var canopyHeight = ee.ImageCollection("projects/sat-io/open-datasets/facebook/meta-canopy-height")
  .filterBounds(bangladesh)
  .mosaic();

var treeCanopyH = canopyHeight
  .updateMask(canopyHeight.gte(1))
  .rename("Tree_Canopy_Height_m");

var treeCanopyFrac = canopyHeight
  .gte(1)
  .rename("Tree_Canopy_Fraction");


// ============================================================================
// VIIRS NIGHTTIME LIGHT
// ============================================================================

var alan = ee.ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG")
  .filterDate(startDate, endDate)
  .filterBounds(bangladesh)
  .select("avg_rad")
  .mean()
  .rename("ALAN_Radiance");


// ============================================================================
// GOOGLE OPEN BUILDINGS HEIGHT
// ============================================================================

var openBuildingsTemporal = ee.ImageCollection("GOOGLE/Research/open-buildings-temporal/v1")
  .filterDate(startDate, endDate)
  .filterBounds(bangladesh)
  .select("building_height")
  .mean()
  .rename("Building_Height_m");


// ============================================================================
// MASTER IMAGE
// ============================================================================

var masterComposite = ee.Image.cat([
  ndvi,
  evi,
  ndmi,
  ndwi,
  mndwi,
  ndbi,
  bsi,
  lstMean,
  lstP90,
  extremeHeatDays,
  no2,
  o3,
  aod,
  treeFrac,
  blueFrac,
  impFrac,
  treeCanopyH,
  treeCanopyFrac,
  alan,
  openBuildingsTemporal
]);

print("Master Bands:", masterComposite.bandNames());


// ============================================================================
// OUTPUT COLUMNS
// ============================================================================

var keepColumns = [
  "OBJECTID",
  "NDVI",
  "EVI",
  "NDMI",
  "NDWI",
  "MNDWI",
  "NDBI",
  "BSI",
  "LST_Mean",
  "LST_P90",
  "Extreme_Heat_Days",
  "NO2_umol_m2",
  "O3_DU",
  "AOD",
  "Tree_Fraction",
  "Tree_Cover_Pct",
  "Blue_Fraction",
  "Blue_Space_Pct",
  "Impervious_Fraction",
  "Impervious_Pct",
  "Tree_Canopy_Height_m",
  "Tree_Canopy_Fraction",
  "Tree_Canopy_Coverage_Pct",
  "ALAN_Radiance",
  "Building_Height_m"
];


// ============================================================================
// SCHOOL EXPORT
// ============================================================================


var totalSchools = 56880;
var schoolList = school.toList(totalSchools);

for (var i = 0; i < totalSchools; i += chunkSize) {

  var endIndex = Math.min(i + chunkSize, totalSchools);

  var subCol = ee.FeatureCollection(schoolList.slice(i, endIndex));

  var bufferedSchools = subCol.map(function(feature) {
    return feature.buffer(bufferRadiusMeters);
  });

  var reduced = masterComposite.reduceRegions({
    collection: bufferedSchools,
    reducer: ee.Reducer.mean(),
    scale: 10,
    tileScale: 8,
    maxPixelsPerRegion: 100000
  });

  var formatted = reduced.map(function(feature) {

    var tree = ee.Number(feature.get("Tree_Fraction"));
    var blue = ee.Number(feature.get("Blue_Fraction"));
    var impervious = ee.Number(feature.get("Impervious_Fraction"));
    var canopy = ee.Number(feature.get("Tree_Canopy_Fraction"));

    return feature.set({
      "Tree_Cover_Pct": tree.multiply(100),
      "Blue_Space_Pct": blue.multiply(100),
      "Impervious_Pct": impervious.multiply(100),
      "Tree_Canopy_Coverage_Pct": canopy.multiply(100)
    });

  });

  var chunkNumber = Math.floor(i / chunkSize) + 1;

  Export.table.toDrive({
    collection: formatted,
    description: "PrimarySchools_2026_Part_" + chunkNumber,
    folder: "Green School",
    fileFormat: "CSV",
    selectors: keepColumns
  });

  print("Export Task:", chunkNumber, "Schools:", i, "-", endIndex);
}

print("All export tasks created.");
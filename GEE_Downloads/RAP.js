
/*
Band 1 - annual forb and grass
Band 2 - bare ground
Band 3 - litter
Band 4 - perennial forb and grass
Band 5 - shrub
Band 6 - tree
*/

var prism_ic = ee.ImageCollection('projects/sat-io/open-datasets/OREGONSTATE/PRISM_800_MONTHLY');
var first_im = prism_ic.first().select('ppt');
var scale = first_im.projection().nominalScale().getInfo();

var rap_ic = ee.ImageCollection('projects/rap-data-365417/assets/vegetation-cover-v3');
print(rap_ic.first());
var years_list = ee.List([1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005]);
var index_list = ee.List([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]);
var ndays_months = ee.List([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]);
var order_months = ee.List([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]);

var area_shp = ee.Geometry.Rectangle([-121, 30, -102, 43], 'EPSG:4326', false);

var n = years_list.size();
var dof = n.subtract(4);

var im = rap_ic.first();
var proj = im.projection().getInfo();

var transform = [
  proj['transform'][0],
  proj['transform'][1],
  proj['transform'][2],
  proj['transform'][3],
  proj['transform'][4],
  proj['transform'][5],
];

var transform_new = [
  0.0083333333,
  proj['transform'][1],
  proj['transform'][2],
  proj['transform'][3],
  0.0083333333,
  proj['transform'][5],
]

var proj = im.projection();

//PRISM STUFF

function yr_fn(yrobj){
  var year = ee.Number(yrobj);
  var start = ee.Date.fromYMD(year, 1, 1);
  var end = ee.Date.fromYMD(year.add(1), 1, 1);
  var precip_ic = prism_ic.filterDate(start, end).select('ppt');
  var precip_im = precip_ic.reduce(ee.Reducer.sum()).clip(area_shp);
  var tmax_ic = prism_ic.filterDate(start, end).select('tmax');
  var tmax_im = tmax_ic.reduce(ee.Reducer.mean()).clip(area_shp);
  var tmin_ic = prism_ic.filterDate(start, end).select('tmin');
  var tmin_im = tmin_ic.reduce(ee.Reducer.mean()).clip(area_shp);
  var year_im = precip_im.addBands(tmax_im).addBands(tmin_im);
  return year_im;
}

var clima_ic = ee.ImageCollection(years_list.map(yr_fn));



//RAP STUFF

function clip_fn(yrobj){
  var year = ee.Number(yrobj);
  var start = ee.Date.fromYMD(year, 1, 1);
  var end = ee.Date.fromYMD(year.add(1), 1, 1);
  //var yr_im = rap_ic.filterDate(start, end).select(['AFG', 'PFG']).first();
  var yr_im = rap_ic.filterDate(start, end).select(['BGR']).first();
  var im = yr_im.setDefaultProjection('EPSG:4326', transform_new);
  var im = im.reproject({crs:proj.crs(), crsTransform:transform_new});
  var clip_im = ee.Image(im).clip(area_shp);
  return clip_im;
}

var cover_ic = ee.ImageCollection(years_list.map(clip_fn));



//MERGE STUFF

var clima_ic_list = clima_ic.toList(999);
var cover_ic_list = cover_ic.toList(999);

function merge_bands_fn(iobj){
  var i = ee.Number(iobj);
  var p_im = ee.Image(clima_ic_list.get(i));
  var c_im = ee.Image(cover_ic_list.get(i));
  var merge_im = p_im.addBands(c_im);
  //var add_im = ee.Image(c_im.select('AFG')).add(ee.Image(c_im.select('PFG')));
  //var merge_im = merge_im.addBands(add_im);
  return merge_im;
}

var merge_ic = ee.ImageCollection(ee.List(index_list.map(merge_bands_fn)));
print('MERGE IC');
print(merge_ic);



//CORRELATION STUFF

function createConstantBand_fn(image){
  return ee.Image(1).addBands(image);
}

var regr_ic = merge_ic.map(createConstantBand_fn);
//var regr_ic = regr_ic.select(['constant', 'ppt_sum', 'tmax_mean', 'tmin_mean', 'AFG_1']);
var regr_ic = regr_ic.select(['constant', 'ppt_sum', 'tmax_mean', 'tmin_mean', 'BGR']);
var regr_im = regr_ic.reduce(ee.Reducer.linearRegression({numX: 4, numY: 1}));

var rmsr_im = regr_im.select('residuals').arrayProject([0]).arrayFlatten([['rmsr']]);
var rss_im = rmsr_im.pow(2).multiply(n);
var sSquared_im = rss_im.divide(dof);
//var yVariance_im = merge_ic.select('AFG_1').reduce(ee.Reducer.sampleVariance());
var yVariance_im = merge_ic.select('BGR').reduce(ee.Reducer.sampleVariance());
var rSquareAdj_im = ee.Image(1).subtract(sSquared_im.divide(yVariance_im));
var coeff_im = regr_im.select('coefficients').arrayProject([0]).arrayFlatten([['c1', 'c2', 'c3', 'c4']]);
var mlr_im = rmsr_im.addBands(rSquareAdj_im).addBands(coeff_im);

Map.addLayer(rSquareAdj_im, {min:0, max:1}, 'rsqr');

//PREDICTION INPUT IC STUFF

var input_ic = clima_ic;



//PREDICTION OUTPUT STUFF
//grass

function prediction_fn(imobj){
  var cli_im = ee.Image(imobj);
  var c1 = mlr_im.select('c1');
  var c2 = cli_im.select('ppt_sum').multiply(mlr_im.select('c2'));
  var c3 = cli_im.select('tmax_mean').multiply(mlr_im.select('c3'));
  var c4 = cli_im.select('tmin_mean').multiply(mlr_im.select('c4'));
  var pred_im = c1.add(c2).add(c3).add(c4);
  return pred_im;
}

var grass_prediction_ic = ee.ImageCollection(input_ic.map(prediction_fn));
Map.addLayer(grass_prediction_ic.first(), {min:0, max:1}, 'Grass Prediction');
//Map.addLayer(merge_ic.first().select('AFG_1'), {min:0, max:1}, 'merge_ic.first()');
Map.addLayer(merge_ic.first().select('BGR'), {min:0, max:1}, 'merge_ic.first()');


//POINT STUFF

// print('POINT STUFF');
// var pt = ee.Geometry.Point(-110, 35);
// Map.addLayer(pt, null, 'point');


// var merge_props = merge_ic.getRegion(pt, 9);
// print('MERGE PROPS');
// print(merge_props);
// var merge_props = merge_props.slice(1);
// var merge_props = merge_props.getInfo();

// var point_result = [];
// for (var i = 0; i < merge_props.length; i++) {
//   point_result.push(merge_props[i][9]);
// }

// print(point_result);


/*
Band 1 - annual forb and grass
Band 2 - bare ground
Band 3 - litter
Band 4 - perennial forb and grass
Band 5 - shrub
Band 6 - tree
*/

var prism_ic = ee.ImageCollection('OREGONSTATE/PRISM/AN81m');
var first_im = prism_ic.first().select('ppt');
var scale = first_im.projection().nominalScale().getInfo();

var rap_ic = ee.ImageCollection('projects/rap-data-365417/assets/vegetation-cover-v3');

var years_list= ee.List([1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005]);
var index_list = ee.List([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]);
var ndays_months = ee.List([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]);
var order_months = ee.List([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]);

var area_shp = ee.Geometry.Rectangle([-121, 30, -102, 43], 'EPSG:4326', false);


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
  0.041666666666,
  proj['transform'][1],
  proj['transform'][2],
  proj['transform'][3],
  0.041666666666,
  proj['transform'][5],
];

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
  var yr_im = rap_ic.filterDate(start, end).select(['AFG', 'PFG']).first();
  var im = yr_im.setDefaultProjection('EPSG:4326', transform);
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
  var add_im = ee.Image(c_im.select('AFG')).add(ee.Image(c_im.select('PFG')));
  var merge_im = merge_im.addBands(add_im);
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
var regr_ic = regr_ic.select(['constant', 'ppt_sum', 'tmax_mean', 'tmin_mean', 'AFG_1']);
var regr_im = regr_ic.reduce(ee.Reducer.linearRegression({numX: 4, numY: 1}));

Map.addLayer(regr_im);

var n = merge_ic.size();
var dof = n.subtract(4);

var rmsr = regr_im.select('residuals').arrayProject([0]).arrayFlatten([['rmsr']]);
var rss = rmsr.pow(2).multiply(n);
var sSquared = rss.divide(dof);

var yVariance = merge_ic.select('AFG_1').reduce(ee.Reducer.sampleVariance());
var rSquareAdj = ee.Image(1).subtract(sSquared.divide(yVariance));

Map.addLayer(rSquareAdj);



//POINT STUFF

print('POINT STUFF');
var pt = ee.Geometry.Point(-110, 35);
Map.addLayer(pt);


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

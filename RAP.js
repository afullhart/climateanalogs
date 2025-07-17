
var prism_ic = ee.ImageCollection('projects/sat-io/open-datasets/OREGONSTATE/PRISM_800_MONTHLY');
var first_im = prism_ic.first().select('ppt');
var scale = first_im.projection().nominalScale().getInfo();
var area_shp = ee.Geometry.Rectangle([-121, 30, -102, 43], 'EPSG:4326', false);

var years_list = ee.List.sequence(1986, 2024);
var yr_idx_list = ee.List.sequence(0, years_list.length().subtract(1));
var ndays_months = ee.List([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]);
var order_months = ee.List([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]);
var n = years_list.size();
var dof = n.subtract(4);

var years_str_list = ['1986', '1987', '1988', '1989', '1990', '1991', '1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999', '2000', '2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018', '2019', '2020', '2021', '2022', '2023', '2024'];
var cover_bands = ['AFG', 'BGR', 'LTR', 'PFG', 'SHR', 'TRE'];
var prod_bands = ['afgNPP', 'pfgNPP', 'shrNPP', 'treNPP'];
var cover_ic = ee.ImageCollection('projects/rap-data-365417/assets/vegetation-cover-v3');
var prod_ic = ee.ImageCollection('projects/rap-data-365417/assets/npp-partitioned-v3');
var metric_strs = ['rmsr', 'rsqr'];


var im = cover_ic.first();
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
  var tdew_ic = prism_ic.filterDate(start, end).select('tdmean');
  var tdew_im = tdew_ic.reduce(ee.Reducer.mean()).clip(area_shp);
  var vmin_ic = prism_ic.filterDate(start, end).select('vpdmin');
  var vmin_im = vmin_ic.reduce(ee.Reducer.mean()).clip(area_shp);
  var vmax_ic = prism_ic.filterDate(start, end).select('vpdmax');
  var vmax_im = vmax_ic.reduce(ee.Reducer.mean()).clip(area_shp);
  var year_im = precip_im.addBands(tmax_im).addBands(tmin_im).addBands(tdew_im).addBands(vmin_im).addBands(vmax_im);
  return year_im;
}

var clima_ic = ee.ImageCollection(years_list.map(yr_fn));



function main_fn(band, rap_ic, out_im_type){

  function clip_fn(yrobj){
    var year = ee.Number(yrobj);
    var start = ee.Date.fromYMD(year, 1, 1);
    var end = ee.Date.fromYMD(year.add(1), 1, 1);
    var yr_im = rap_ic.filterDate(start, end).select([band]).first();
    var im = yr_im.setDefaultProjection('EPSG:4326', transform_new);
    var im = im.reproject({crs:proj.crs(), crsTransform:transform_new});
    var clip_im = ee.Image(im).clip(area_shp);
    return clip_im;
  }
  
  var rap_ic = ee.ImageCollection(years_list.map(clip_fn));

  var clima_ic_list = clima_ic.toList(999);
  var rap_ic_list = rap_ic.toList(999);
  
  function merge_bands_fn(iobj){
    var i = ee.Number(iobj);
    var p_im = ee.Image(clima_ic_list.get(i));
    var c_im = ee.Image(rap_ic_list.get(i));
    var merge_im = p_im.addBands(c_im);
    return merge_im;
  }
  
  var merge_ic = ee.ImageCollection(ee.List(yr_idx_list.map(merge_bands_fn)));

  function createConstantBand_fn(image){
    return ee.Image(1).addBands(image);
  }
  
  var regr_ic = merge_ic.map(createConstantBand_fn);
  var regr_ic = regr_ic.select(['constant', 'ppt_sum', 'tmax_mean', 'tmin_mean', 'tdmean_mean', 'vpdmin_mean', 'vpdmax_mean', band]);
  var regr_im = regr_ic.reduce(ee.Reducer.linearRegression({numX: 7, numY: 1}));
  
  var rmsr_im = regr_im.select('residuals').arrayProject([0]).arrayFlatten([['rmsr']]);
  var rss_im = rmsr_im.pow(2).multiply(n);
  var sSquared_im = rss_im.divide(dof);
  var yVariance_im = merge_ic.select(band).reduce(ee.Reducer.sampleVariance());
  var rSquareAdj_im = ee.Image(1).subtract(sSquared_im.divide(yVariance_im));
  var coeff_im = regr_im.select('coefficients').arrayProject([0]).arrayFlatten([['c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7']]);
  var mlr_im = rmsr_im.addBands(rSquareAdj_im).addBands(coeff_im);

  function prediction_fn(imobj){
    var cli_im = ee.Image(imobj);
    var c1 = mlr_im.select('c1');
    var c2 = cli_im.select('ppt_sum').multiply(mlr_im.select('c2'));
    var c3 = cli_im.select('tmax_mean').multiply(mlr_im.select('c3'));
    var c4 = cli_im.select('tmin_mean').multiply(mlr_im.select('c4'));
    var c5 = cli_im.select('tdmean_mean').multiply(mlr_im.select('c5'));
    var c6 = cli_im.select('vpdmin_mean').multiply(mlr_im.select('c6'));
    var c7 = cli_im.select('vpdmax_mean').multiply(mlr_im.select('c7'));
    var pred_im = c1.add(c2).add(c3).add(c4).add(c5).add(c6).add(c7);
    return pred_im;
  }

  if (out_im_type == 'rmsr'){
    return rmsr_im;
  }else if (out_im_type == 'rsqr'){
    return rSquareAdj_im;
  }else if (out_im_type.slice(0, 3) == 'cov' || out_im_type.slice(0, 3) == 'pro'){
    var year = ee.Number.parse(out_im_type.slice(3));
    var year_idx = years_list.indexOf(year);
    var rap_im = ee.Image(rap_ic_list.get(year_idx));
    return rap_im;
  }else if (out_im_type == 'OnetoOne'){
    var prediction_ic = ee.ImageCollection(merge_ic.map(prediction_fn));
    var oneone_ic = rap_ic.merge(prediction_ic);
    return oneone_ic;
  }
}


///////////////////////
///////////////////////
//USER INTERFACE
///////////////////////
///////////////////////


//https://github.com/gee-community/ee-palettes
var palettes = require('users/gena/packages:palettes');

var rmsrVis = {
  min:0,
  max:10,
  palette:palettes.misc.parula[7]
};

var rsqrVis = {
  min:0,
  max:0.67,
  palette:palettes.misc.warmcool[7]
};

var covVis = {
  min:0,
  max:50,
  palette:palettes.niccoli.cubicl[7]
};

var proVis = {
  min:0,
  max:100,
  palette:palettes.niccoli.cubicl[7]
};

var widgetStyle = {
  position:'bottom-center'
};

var checkStyle = {
  position:'bottom-center'
};

var timelinePanelStyle = {
  position:'top-center', 
  stretch:'vertical',
  height:'400px',
  width:'400px',
  margin:'10px 10px'};

/////////////////////
//Global widget vars

var band_selection = 'AFG';
var ic_selection = cover_ic;
var type_selection = 'rmsr';
var bandVis = palettes.colorbrewer.Paired;

var im_to_show = main_fn(band_selection, ic_selection, 'OnetoOne');
print(im_to_show);

/////////////////////
//Change RAP Variable

function renderVariable(var_selection){
  Map.layers().reset();
  band_selection = var_selection;
  if (cover_bands.indexOf(var_selection) >= 0){
    ic_selection = cover_ic;
  } else {
    ic_selection = prod_ic;
  }
  var im_to_show = main_fn(band_selection, ic_selection, type_selection);
  if (type_selection == 'rmsr'){
    var bandVis = rmsrVis;
  }else if (type_selection == 'rsqr'){
    var bandVis = rsqrVis;
  }
  Map.addLayer(im_to_show, bandVis);
}

var variable_dropdown = ui.Select({
  items:cover_bands.concat(prod_bands), 
  placeholder:'Select Variable', 
  onChange:renderVariable,
  style:widgetStyle
});

ui.root.add(variable_dropdown);

/////////////////////
//Change Metric Layer

function renderMetric(metric_selection){
  Map.layers().reset();
  type_selection = metric_selection;
  var im_to_show = main_fn(band_selection, ic_selection, type_selection);
  if (type_selection == 'rmsr'){
    var bandVis = rmsrVis;
  }else if (type_selection == 'rsqr'){
    var bandVis = rsqrVis;
  }
  Map.addLayer(im_to_show, bandVis);
}

var metric_dropdown = ui.Select({
  items:metric_strs, 
  placeholder:'Select Metric', 
  onChange:renderMetric,
  style:widgetStyle
});

ui.root.add(metric_dropdown);

/////////////////////
//Render RAP Layer

function renderYear(year_selection){
  Map.layers().reset();
  if (cover_bands.indexOf(band_selection) != -1){
    type_selection = 'cov'.concat(year_selection);
  }else if (prod_bands.indexOf(band_selection) != -1){
    type_selection = 'pro'.concat(year_selection);
  }
  var im_to_show = main_fn(band_selection, ic_selection, type_selection);
  if (type_selection.slice(0, 3) == 'cov'){
    var bandVis = covVis;
  }else if (type_selection.slice(0, 3) == 'pro'){
    var bandVis = proVis;
  }
  Map.addLayer(im_to_show, bandVis);
}

var year_dropdown = ui.Select({
  items:years_str_list, 
  placeholder:'Select RAP Year', 
  onChange:renderYear,
  style:widgetStyle
});

ui.root.add(year_dropdown);

/////////////////////////
//Render One-to-One Chart

var timelinePanel = ui.Panel({style:timelinePanelStyle});

function makePlot(){
  
}

function clickCallback(clickInfo_obj){
  timelinePanel = ui.Panel({widgets:[], style:timelinePanelStyle});
  var lat = clickInfo_obj.lat;
  var lon = clickInfo_obj.lon;
  var pt = ee.Geometry.Point([lon, lat]);
  var ic_to_show = main_fn(band_selection, ic_selection, 'OnetoOne');
  //var 
  var datelabel = ui.Label({
    value:String(lat).concat(String(lon)),
    style:{
      padding:'1px',
      margin:'0px',
      position:'middle-left',
      fontSize:'12px'}
  });
  var row = ui.Panel({
    widgets:datelabel,
    layout:ui.Panel.Layout.Flow('horizontal')
  });
  timelinePanel.add(row);
  ui.root.add(timelinePanel);
}

function renderOnetoOne(checkbox_bool){
  if (checkbox_bool === true){
    Map.onClick(clickCallback);
  }
  else{
    ui.root.remove(timelinePanel);
  }
}

var OnetoOne_checkbox = ui.Checkbox({
  label:'One-to-One Plot on Click',
  onChange:renderOnetoOne,
  style:checkStyle
});

ui.root.add(OnetoOne_checkbox);







//PREDICTION INPUT IC STUFF

//var input_ic = clima_ic;

// function prediction_fn(imobj){
//   var cli_im = ee.Image(imobj);
//   var c1 = mlr_im.select('c1');
//   var c2 = cli_im.select('ppt_sum').multiply(mlr_im.select('c2'));
//   var c3 = cli_im.select('tmax_mean').multiply(mlr_im.select('c3'));
//   var c4 = cli_im.select('tmin_mean').multiply(mlr_im.select('c4'));
//   var c5 = cli_im.select('tdmean_mean').multiply(mlr_im.select('c5'));
//   var c6 = cli_im.select('vpdmin_mean').multiply(mlr_im.select('c6'));
//   var c7 = cli_im.select('vpdmax_mean').multiply(mlr_im.select('c7'));
//   var pred_im = c1.add(c2).add(c3).add(c4).add(c5).add(c6).add(c7);
//   return pred_im;
// }

// var grass_prediction_ic = ee.ImageCollection(input_ic.map(prediction_fn));
// Map.addLayer(grass_prediction_ic.first(), {min:0, max:70}, 'Grass Prediction');
// Map.addLayer(merge_ic.first().select(band_selection), {min:0, max:70}, 'merge_ic.first()');


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

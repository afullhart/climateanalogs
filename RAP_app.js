
var prism_ic = ee.ImageCollection('projects/sat-io/open-datasets/OREGONSTATE/PRISM_800_MONTHLY');
var first_im = prism_ic.first().select('ppt');
var scale = first_im.projection().nominalScale().getInfo();
var area_shp = ee.Geometry.Rectangle([-121, 30, -102, 43], 'EPSG:4326', false);

var years_list = ee.List.sequence(1986, 2024);
var yr_idx_list = ee.List.sequence(0, years_list.length().subtract(1));
var ndays_months = ee.List([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]);
var order_months = ee.List([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]);
function make_date_fn(year_num){
  var date = ee.Date.fromYMD(year_num, 1, 1);
  return date;
}
var dates_list = ee.List(years_list.map(make_date_fn));
var years_str_list = ['1986', '1987', '1988', '1989', '1990', '1991', '1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999', '2000', '2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018', '2019', '2020', '2021', '2022', '2023', '2024'];
var cover_bands = ['AFG', 'BGR', 'LTR', 'PFG', 'SHR', 'TRE'];
var prod_bands = ['afgNPP', 'pfgNPP', 'shrNPP', 'treNPP'];
var cover_ic = ee.ImageCollection('projects/rap-data-365417/assets/vegetation-cover-v3');
var prod_ic = ee.ImageCollection('projects/rap-data-365417/assets/npp-partitioned-v3');
var metric_strs = ['rmse', 'rsqr', 'rsqrA', 'Fconf', 'coeff', 'Tconf'];

var n = years_list.size();
var k = ee.Number(6);
var dof = n.subtract(k).subtract(1);
//for F-statistic look-up table
//n - k = 39 - 6 = 33
//k - 1 = 6 - 1 = 5
//1% = 3.630
//5% = 2.503
//10% = 2.030
//for T-statistic look-up
//df = n - 1 = 39 - 1 = 38
//1% = 2.712
//5% = 2.024
//10% = 1.686

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
  
  //Climate model statistics
  var regr_ic = merge_ic.map(createConstantBand_fn);
  var regr_ic = regr_ic.select(['constant', 'ppt_sum', 'tmax_mean', 'tmin_mean', 'tdmean_mean', 'vpdmin_mean', 'vpdmax_mean', band]);
  var regr_im = regr_ic.reduce(ee.Reducer.linearRegression({numX:7, numY:1}));
  var rmsr_im = regr_im.select('residuals').arrayProject([0]).arrayFlatten([['rmsr']]);
  var rss_im = rmsr_im.pow(2).multiply(n);
  var sSquare_im = rss_im.divide(dof);
  var yVariance_im = merge_ic.select(band).reduce(ee.Reducer.sampleVariance());
  var rSquare_im = ee.Image(1).subtract(ee.Image(rss_im.divide(yVariance_im.multiply(n.subtract(1)))));
  var rSquareAdj_im = ee.Image(1).subtract(sSquare_im.divide(yVariance_im));
  var coeff_im = regr_im.select('coefficients').arrayProject([0]).arrayFlatten([['c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7']]);
  var top_im = rSquare_im.divide(k);
  var bot_im = ee.Image(ee.Image(1).subtract(rSquare_im)).divide(n.subtract(k).subtract(1));
  var f_im = top_im.divide(bot_im);
  var zero_im = rmsr_im.lt(0.0);
  var ninenine_im = zero_im.where(f_im.gte(3.630), 1);
  var ninefive_im = zero_im.where(f_im.gte(2.503), 1);
  var ninezero_im = zero_im.where(f_im.gte(2.030), 1);
  var conf_im = zero_im.add(ninenine_im).add(ninefive_im).add(ninezero_im);
  var mlr_im = rmsr_im.addBands(rSquareAdj_im).addBands(coeff_im);

  //Trend model statistics
  function merg_bands_fn(iobj){
    var i = ee.Number(iobj);
    var p_im = ee.Image.constant(yr_idx_list.get(i)).toFloat();
    var c_im = ee.Image(rap_ic_list.get(i));
    var merge_im = c_im.addBands(p_im);
    return merge_im;
  }
  var merg_ic = ee.ImageCollection(ee.List(yr_idx_list.map(merg_bands_fn)));
  var reg_ic = merg_ic.map(createConstantBand_fn);
  var reg_ic = reg_ic.select(['constant', 'constant_1', band]);
  var reg_im = reg_ic.reduce(ee.Reducer.linearRegression({numX:2, numY:1}));
  var rmse_im = reg_im.select('residuals').arrayProject([0]).arrayFlatten([['rmsr']]);
  var top_im = ee.Image(rmse_im.pow(2).multiply(n).divide(n.subtract(2))).pow(0.5);
  var bot_im = reg_ic.select('constant_1').reduce(ee.Reducer.sampleVariance()).multiply(n.subtract(1)).pow(0.5);
  var sslp_im = top_im.divide(bot_im);
  var coef_im = reg_im.select('coefficients').arrayProject([0]).arrayFlatten([['c1', 'c2']]);
  var t_im = coef_im.select('c2').divide(sslp_im);
  var zer_im = rmse_im.lt(0.0);
  var ninenin_im = zer_im.where(t_im.gte(2.712), 1);
  var ninefiv_im = zer_im.where(t_im.gte(2.024), 1);
  var ninezer_im = zer_im.where(t_im.gte(1.686), 1);
  var con_im = zer_im.add(ninenin_im).add(ninefiv_im).add(ninezer_im);

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

  if (out_im_type == 'rmse'){
    return rmsr_im;
  }else if (out_im_type == 'rsqr'){
    return rSquare_im;
  }else if (out_im_type == 'rsqrA'){
    return rSquareAdj_im;
  }else if (out_im_type == 'Fconf'){
    return conf_im;
  }else if (out_im_type == 'coeff'){
    return coef_im.select('c2');
  }else if (out_im_type == 'Tconf'){
    return con_im;
  }else if (out_im_type.slice(0, 3) == 'cov' || out_im_type.slice(0, 3) == 'pro'){
    var year = ee.Number.parse(out_im_type.slice(3));
    var year_idx = years_list.indexOf(year);
    var rap_im = ee.Image(rap_ic_list.get(year_idx));
    return rap_im;
  }else if (out_im_type.slice(0, 3) == 'covpred' || out_im_type.slice(0, 3) == 'propred'){
    var prediction_ic = ee.ImageCollection(merge_ic.map(prediction_fn));
    var pred_ic_list = prediction_ic.toList(999);
    var year = ee.Number.parse(out_im_type.slice(7));
    var year_idx = years_list.indexOf(year);
    var pred_im = ee.Image(pred_ic_list.get(year_idx));
    return pred_im;
  }else if (out_im_type == 'OnetoOne'){
    var prediction_ic = ee.ImageCollection(merge_ic.map(prediction_fn));
    var oneone_ic = rap_ic.merge(prediction_ic);
    return oneone_ic;
  }else if (out_im_type == 'Trend'){
    return rap_ic;
  }else if (out_im_type == 'Debug'){
    var prediction_ic = ee.ImageCollection(merge_ic.map(prediction_fn));
    return ee.List([t_im, prediction_ic]);
  }
}



//START TesTING TEsTING tEStIng tESTiNG TEsTiNG
//var t_im = ee.Image(main_fn('PFG', cover_ic, 'Debug').get(0));
//print(t_im);
//Map.addLayer(t_im);

// var merge_props = merge_ic.getRegion(geometry, scale);
// var pred_props = pred_ic.getRegion(geometry, scale);
// var merge_props = merge_props.slice(1);
// var pred_props = pred_props.slice(1);

// function make_prop_feats_fnA(p_list){
//   var p_list = ee.List(p_list);
//   var a = p_list.get(4);
//   var b = p_list.get(5);
//   var c = p_list.get(6);
//   var d = p_list.get(7);
//   var e = p_list.get(8);
//   var f = p_list.get(9);
//   var g = p_list.get(10);
//   var ft = ee.Feature(null, {a:a, b:b, c:c, d:d, e:e, f:f, g:g});
//   return ft;
// }

// var out_fc = ee.FeatureCollection(merge_props.map(make_prop_feats_fnA));

// function make_prop_feats_fnB(p_list){
//   var p_list = ee.List(p_list);
//   var a = p_list.get(4);
//   var ft = ee.Feature(null, {a:a});
//   return ft;
// }

// var pred_fc = ee.FeatureCollection(pred_props.map(make_prop_feats_fnB));

// Export.table.toDrive({
//   collection:out_fc,
//   description:'vectorsToDriveExample'
// });

// Export.table.toDrive({
//   collection:pred_fc,
//   description:'vectorsPred'
// });
//END TesTING TEsTING tEStIng tESTiNG TEsTiNG



///////////////////////
///////////////////////
//USER INTERFACE
///////////////////////
///////////////////////


//https://github.com/gee-community/ee-palettes
var palettes = require('users/gena/packages:palettes');
var rmseCovVis = {
  min:0,
  max:5,
  palette:palettes.kovesi.diverging_linear_bjr_30_55_c53[7]
};

var palettes = require('users/gena/packages:palettes');
var rmseProVis = {
  min:0,
  max:100,
  palette:palettes.kovesi.diverging_linear_bjr_30_55_c53[7]
};

var palettes = require('users/gena/packages:palettes');
var rsqrVis = {
  min:0,
  max:0.6,
  palette:palettes.kovesi.diverging_linear_bjr_30_55_c53[7].reverse()
};

var palettes = require('users/gena/packages:palettes');
var rsqrAdjVis = {
  min:0,
  max:0.4,
  palette:palettes.kovesi.diverging_linear_bjr_30_55_c53[7].reverse()
};

var palettes = require('users/gena/packages:palettes');
var covVis = {
  min:0,
  max:50,
  palette:palettes.niccoli.cubicl[7]
};

var palettes = require('users/gena/packages:palettes');
var proVis = {
  min:0,
  max:800,
  palette:palettes.niccoli.cubicl[7]
};

var confVis = {
  min:0,
  max:3,
  palette:['#FFFFFF', '#d7481d', '#59f720', '#800080']
};

var palettes = require('users/gena/packages:palettes');
var covcoeffVis = {
  min:-0.5,
  max:0.5,
  palette:palettes.kovesi.diverging_gkr_60_10_c40[7].reverse()
};

var palettes = require('users/gena/packages:palettes');
var procoeffVis = {
  min:-10,
  max:10,
  palette:palettes.kovesi.diverging_gkr_60_10_c40[7]
};

var widgetStyle = {
  position:'bottom-center'
};

var checkStyle = {
  position:'bottom-center'
};

var chartPanelStyle = {
  position:'bottom-center', 
  stretch:'vertical',
  height:'400px',
  width:'400px',
  margin:'10px 10px'};

//ui.root.setLayout(ui.Panel.Layout.absolute());
var main_panel = ui.Panel({
  layout:ui.Panel.Layout.flow('vertical'),
  style: {width: '300px'}
});

/////////////////////
//Global widget vars

var band_selection = 'PFG';
var ic_selection = cover_ic;
var type_selection = 'cov2000';
var im_to_show = main_fn(band_selection, ic_selection, type_selection);
var bandVis = covVis;
var chart_panelA = ui.Panel({style:chartPanelStyle});
var chart_panelB = ui.Panel({style:chartPanelStyle});

Map.addLayer(im_to_show, bandVis);

/////////////////////
//Change RAP Variable

function renderVariable(var_selection){
  Map.layers().reset();
  band_selection = var_selection;
  if (cover_bands.indexOf(band_selection) >= 0){
    ic_selection = cover_ic;
    var bandVis = covVis;
    if (type_selection.slice(0, 3) == 'pro'){
      var year = type_selection.slice(-4);
      type_selection = type_selection.replace('pro', 'cov').slice(0, -4).concat(year);
    }    
  } else {
    ic_selection = prod_ic;
    var bandVis = proVis;
    if (type_selection.slice(0, 3) == 'cov'){
      var year = type_selection.slice(-4);
      type_selection = type_selection.replace('cov', 'pro').slice(0, -4).concat(year);
    }
  }
  var im_to_show = main_fn(band_selection, ic_selection, type_selection);
  if (type_selection == 'rmse' && cover_bands.indexOf(band_selection) >= 0){
    var bandVis = rmseCovVis;
  }else if (type_selection == 'rmse' && prod_bands.indexOf(band_selection) >= 0){
    var bandVis = rmseProVis;
  }else if (type_selection == 'rsqr'){
    var bandVis = rsqrVis;
  }else if (type_selection == 'rsqrA'){
    var bandVis = rsqrAdjVis;
  }else if (type_selection == 'Fconf'){
    var bandVis = confVis;
  }else if (type_selection == 'coeff' && cover_bands.indexOf(band_selection) >= 0){
    var bandVis = covcoeffVis;
  }else if (type_selection == 'coeff' && prod_bands.indexOf(band_selection) >= 0){
    var bandVis = procoeffVis;
  }else if (type_selection == 'Tconf'){
    var bandVis = confVis;
  }else if (type_selection.slice(0, 3) == 'cov'){
    var bandVis = covVis;
  }else if (type_selection.slice(0, 3) == 'pro'){
    var bandVis = proVis;
  }else if (type_selection.slice(0, 7) == 'covpred'){
    var bandVis = covVis;
  }else if (type_selection.slice(0, 7) == 'propred'){
    var bandVis = proVis;
  }
  Map.addLayer(im_to_show, bandVis);
}

var variable_dropdown = ui.Select({
  items:cover_bands.concat(prod_bands), 
  placeholder:'Select Variable', 
  onChange:renderVariable,
  style:widgetStyle
});

main_panel.add(variable_dropdown);

/////////////////////
//Change Metric Layer

function renderMetric(metric_selection){
  Map.layers().reset();
  type_selection = metric_selection;
  var im_to_show = main_fn(band_selection, ic_selection, type_selection);
  if (type_selection == 'rmse' && cover_bands.indexOf(band_selection) >= 0){
    var bandVis = rmseCovVis;
  }else if (type_selection == 'rmse' && prod_bands.indexOf(band_selection) >= 0){
    var bandVis = rmseProVis;
  }else if (type_selection == 'rsqr'){
    var bandVis = rsqrVis;
  }else if (type_selection == 'rsqrA'){
    var bandVis = rsqrAdjVis;
  }else if (type_selection == 'Fconf'){
    var bandVis = confVis;
  }else if (type_selection == 'coeff' && cover_bands.indexOf(band_selection) >= 0){
    var bandVis = covcoeffVis;
  }else if (type_selection == 'coeff' && prod_bands.indexOf(band_selection) >= 0){
    var bandVis = procoeffVis;
  }else if (type_selection == 'Tconf'){
    var bandVis = confVis;
  }
  Map.addLayer(im_to_show, bandVis);
}

var metric_dropdown = ui.Select({
  items:metric_strs, 
  placeholder:'Select Metric', 
  onChange:renderMetric,
  style:widgetStyle
});

main_panel.add(metric_dropdown);

/////////////////////
//Render RAP Layer

function renderYearA(year_selection){
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

var year_dropdownA = ui.Select({
  items:years_str_list, 
  placeholder:'Select RAP Year', 
  onChange:renderYearA,
  style:widgetStyle
});

main_panel.add(year_dropdownA);

////////////////////////
//Render Predicted Layer

function renderYearB(year_selection){
  Map.layers().reset();
  if (cover_bands.indexOf(band_selection) != -1){
    type_selection = 'predcov'.concat(year_selection);
  }else if (prod_bands.indexOf(band_selection) != -1){
    type_selection = 'predpro'.concat(year_selection);
  }
  var im_to_show = main_fn(band_selection, ic_selection, type_selection);
  if (type_selection.slice(0, 7) == 'predcov'){
    var bandVis = covVis;
  }else if (type_selection.slice(0, 7) == 'predpro'){
    var bandVis = proVis;
  }
  Map.addLayer(im_to_show, bandVis);
}

var year_dropdownB = ui.Select({
  items:years_str_list, 
  placeholder:'Select PRED. Year', 
  onChange:renderYearB,
  style:widgetStyle
});

main_panel.add(year_dropdownB);

/////////////////////////
//Render One-to-One Chart

function makePlotA(plot_ic, point_geo){
  var plot_ic_list = plot_ic.toList(999);
  var rap_ic = ee.ImageCollection(plot_ic_list.slice(0, n));
  var pred_ic = ee.ImageCollection(plot_ic_list.slice(n, plot_ic_list.size()));
  var rap_prop_list = rap_ic.getRegion(point_geo, scale).slice(1);
  var pred_prop_list = pred_ic.getRegion(point_geo, scale).slice(1);
  function get_sublist_fn(props){
    var value = ee.List(props).get(4);
  return value;
  }
  var rap_list = rap_prop_list.map(get_sublist_fn);
  var pred_list = pred_prop_list.map(get_sublist_fn);
  var min_value = ee.List([rap_list.reduce(ee.Reducer.min()), pred_list.reduce(ee.Reducer.min())]).reduce(ee.Reducer.min()).getInfo();
  var max_value = ee.List([rap_list.reduce(ee.Reducer.max()), pred_list.reduce(ee.Reducer.max())]).reduce(ee.Reducer.max()).getInfo();
  var oneone_chart = ui.Chart.array.values(pred_list, 0, rap_list)
  .setSeriesNames(['PredVsRap'])
  .setOptions({
    title:'One-to-One',
    titleTextStyle:{italic:false, bold:true, fontSize:26},
    legend:{position:'top-right'},
    hAxis:{viewWindow:{min:min_value*0.9, max:max_value*1.1}, title:'RAP', titleTextStyle:{italic:false, bold:true, fontSize:21}, gridlines:{count:6}},
    vAxis:{viewWindow:{min:min_value*0.9, max:max_value*1.1}, title:'PRED', titleTextStyle:{italic:false, bold:true, fontSize:21}, gridlines:{count:6}},
    colors:['#6a9f58'],
    pointSize:12,
    lineSize:0,
    chartArea:{height:500, width:500},
    trendlines:{0:{type:'linear', color:'green', lineWidth:3, opacity:0.3, visibleInLegend:true}}
  });
  return oneone_chart;
}

function clickCallbackA(clickInfo_obj){
  Map.remove(chart_panelA);
  var lat = clickInfo_obj.lat;
  var lon = clickInfo_obj.lon;
  var pt = ee.Geometry.Point([lon, lat]);
  var ic_to_show = main_fn(band_selection, ic_selection, 'OnetoOne');
  var chart_widget = makePlotA(ic_to_show, pt);
  chart_panelA = ui.Panel({
    widgets:chart_widget,
    layout:ui.Panel.Layout.Flow('horizontal')
  });
  Map.add(chart_panelA);
}

function renderOnetoOne(checkbox_bool){
  if (checkbox_bool === true){
    Map.onClick(clickCallbackA);
  }
  else{
    Map.unlisten();
    Map.remove(chart_panelA);
    chart_panelA = ui.Panel({style:chartPanelStyle});
  }
}

var OnetoOne_checkbox = ui.Checkbox({
  label:'One-to-One Plot on Click',
  onChange:renderOnetoOne,
  style:checkStyle
});

main_panel.add(OnetoOne_checkbox);

/////////////////////////
//Render Trend Chart

function makePlotB(rap_ic, point_geo){
  var rap_prop_list = rap_ic.getRegion(point_geo, scale).slice(1);
  function get_sublist_fn(props){
    var value = ee.List(props).get(4);
  return value;
  }
  var rap_list = rap_prop_list.map(get_sublist_fn);
  var min_value = rap_list.reduce(ee.Reducer.min()).getInfo();
  var max_value = rap_list.reduce(ee.Reducer.max()).getInfo();
  var trend_chart = ui.Chart.array.values(rap_list, 0, yr_idx_list)
  .setSeriesNames(['Rap'])
  .setOptions({
    title:'Trend',
    titleTextStyle:{italic:false, bold:true, fontSize:26},
    legend:{position:'top-right'},
    hAxis:{title:'Year', format:'####', titleTextStyle:{italic:false, bold:true, fontSize:21}, gridlines:{count:6}},
    vAxis:{title:'RAP', titleTextStyle:{italic:false, bold:true, fontSize:21}, gridlines:{count:6}},
    colors:['#6a9f58'],
    pointSize:12,
    lineSize:0,
    chartArea:{height:500, width:500},
    trendlines:{0:{type:'linear', color:'green', lineWidth:3, opacity:0.3, visibleInLegend:true}}
  });
  return trend_chart;
}

function clickCallbackB(clickInfo_obj){
  Map.remove(chart_panelB);
  var lat = clickInfo_obj.lat;
  var lon = clickInfo_obj.lon;
  var pt = ee.Geometry.Point([lon, lat]);
  var ic_to_show = main_fn(band_selection, ic_selection, 'Trend');
  var chart_widget = makePlotB(ic_to_show, pt);
  chart_panelB = ui.Panel({
    widgets:chart_widget,
    layout:ui.Panel.Layout.Flow('horizontal')
  });
  Map.add(chart_panelB);
}

function renderTrend(checkbox_bool){
  if (checkbox_bool === true){
    Map.onClick(clickCallbackB);
  }
  else{
    Map.unlisten();
    Map.remove(chart_panelB);
    chart_panelB = ui.Panel({style:chartPanelStyle});
  }
}

var trend_checkbox = ui.Checkbox({
  label:'RAP Trend on Click',
  onChange:renderTrend,
  style:checkStyle
});

main_panel.add(trend_checkbox);

ui.root.add(main_panel);





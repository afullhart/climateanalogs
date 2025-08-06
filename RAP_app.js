
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
var bio_bands = ['afgAGB', 'pfgAGB'];
var cover_ic = ee.ImageCollection('projects/rap-data-365417/assets/vegetation-cover-v3');
var prod_ic = ee.ImageCollection('projects/rap-data-365417/assets/npp-partitioned-v3');
var mat_ic = ee.ImageCollection('projects/rap-data-365417/assets/gridmet-MAT');
var metric_strs = ['avg', 'sdev', 'rmse', 'rsqr', 'rsqrA', 'Fconf', 'Tconf'];
var model_list = ['A1*pr + A2*tx + A3*tn + A4*td + A5*vx + A6*vn + K1', 'B1*pr + B2*tm + K2', 'C1*pr + K3'];
var coeff_list = ['Atrend', 'Ktrend', 'A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'B1', 'B2', 'C1', 'K1', 'K2', 'K3'];
var coeff_internal_list = ['c2', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'c2', 'c3', 'c2', 'c1', 'c1', 'c1'];

////////////////////////////
//Global main function args
var band_selection = 'PFG';
var ic_selection = cover_ic;
var type_selection = 'cov2000';
////////////////////////////

var model_selection = model_list[0];
var n = years_list.size();

//for F-statistic look-up table
//n - k = 39 - 6 = 33 //36 //37
//k - 1 = 6 - 1 = 5   //2   //1
//MODEL ZERO   //MODEL ONE   //MODEL TWO
//1% = 3.630   //1% = 5.248  //1% = 7.373
//5% = 2.503   //5% = 3.259  //5% = 4.105
//10% = 2.030  //10% = 2.456 //10% = 2.846
//for T-statistic look-up
//df = n - 1 = 39 - 1 = 38
//MODEL ZERO
//1% = 2.715
//5% = 2.026
//10% = 1.687

var first_im = cover_ic.first();
var proj = first_im.projection().getInfo();

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

var proj = first_im.projection();

var mat_ic_list = mat_ic.toList(999);
var prod_ic_list = prod_ic.toList(999);

function npp2biomass_fn(yr_idx){
  var mat_im = ee.Image(mat_ic_list.get(yr_idx));
  var rap_im = ee.Image(prod_ic_list.get(yr_idx));
  var rap_im = rap_im.select(['afgNPP', 'pfgNPP']);
  var year = ee.Date(rap_im.get('system:time_start')).format('YYYY');
  var fANPP_im = mat_im.multiply(0.0129).add(0.171); //Above-ground NPP fraction;
  var agb_im = rap_im.multiply(0.0001) // NPP scalar 
    .multiply(2.20462) // KgC to lbsC
    .multiply(4046.86) // m2 to acres
    .multiply(2.1276); // C to biomass
  var agb_im = ee.Image(agb_im).rename(['afgAGB', 'pfgAGB']);
  var agb_im = agb_im.multiply(fANPP_im);
  var im = ee.Image(agb_im.copyProperties(rap_im, rap_im.propertyNames()));
  var im = im.setDefaultProjection('EPSG:4326', transform_new);
  var im = im.reproject({crs:proj.crs(), crsTransform:transform_new});
  return im;
}

var agb_ic = ee.ImageCollection(yr_idx_list.map(npp2biomass_fn));

//PRISM STUFF

function yr_fn(yrobj){
  var year = ee.Number(yrobj);
  var start = ee.Date.fromYMD(year, 1, 1);
  var end = ee.Date.fromYMD(year.add(1), 1, 1);
  var precip_ic = prism_ic.filterDate(start, end).select('ppt');
  var precip_im = precip_ic.reduce(ee.Reducer.sum()).clip(area_shp);
  var tmean_ic = prism_ic.filterDate(start, end).select('tmean');
  var tmean_im = tmean_ic.reduce(ee.Reducer.mean()).clip(area_shp);
  var tmax_ic = prism_ic.filterDate(start, end).select('tmax');
  var tmax_im = tmax_ic.reduce(ee.Reducer.mean()).clip(area_shp);
  var tmin_ic = prism_ic.filterDate(start, end).select('tmin');
  var tmin_im = tmin_ic.reduce(ee.Reducer.mean()).clip(area_shp);
  var tdew_ic = prism_ic.filterDate(start, end).select('tdmean');
  var tdew_im = tdew_ic.reduce(ee.Reducer.mean()).clip(area_shp);
  var vmax_ic = prism_ic.filterDate(start, end).select('vpdmax');
  var vmax_im = vmax_ic.reduce(ee.Reducer.mean()).clip(area_shp);
  var vmin_ic = prism_ic.filterDate(start, end).select('vpdmin');
  var vmin_im = vmin_ic.reduce(ee.Reducer.mean()).clip(area_shp);
  var year_im = precip_im.addBands(tmean_im).addBands(tmax_im).addBands(tmin_im).addBands(tdew_im).addBands(vmax_im).addBands(vmin_im);
  var year_im = year_im.setDefaultProjection('EPSG:4326', transform_new);
  var year_im = year_im.reproject({crs:proj.crs(), crsTransform:transform_new});
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
  var rap_ic_list = rap_ic.toList(999);

  function nppscaled2npp_fn(yr_idx){
    var rap_im = ee.Image(rap_ic_list.get(yr_idx));
    var npp_im = rap_im.multiply(0.0001).multiply(2.20462).multiply(4046.86); // NPP scalar, KgC to lbsC, m2 to acre
    return npp_im;
  }

  if (band_selection.slice(3) == 'NPP'){
    var rap_ic = ee.ImageCollection(yr_idx_list.map(nppscaled2npp_fn));
  }

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
  
  //Climate model statistics
  if (model_selection == model_list[0]){
    var regr_ic = regr_ic.select(['constant', 'ppt_sum', 'tmax_mean', 'tmin_mean', 'tdmean_mean', 'vpdmax_mean', 'vpdmin_mean', band]);
    var result_im = regr_ic.reduce(ee.Reducer.linearRegression({numX:7, numY:1}));
    var rmsr_im = result_im.select('residuals').arrayProject([0]).arrayFlatten([['rmsr']]);
    var rss_im = rmsr_im.pow(2).multiply(n);
    var k = ee.Number(6);
    var dof = n.subtract(k).subtract(1);
    var sSquare_im = rss_im.divide(dof);
    var yVariance_im = regr_ic.select(band).reduce(ee.Reducer.sampleVariance());
    var rSquare_im = ee.Image(1).subtract(ee.Image(rss_im.divide(yVariance_im.multiply(n.subtract(1)))));
    var rSquareAdj_im = ee.Image(1).subtract(sSquare_im.divide(yVariance_im));
    var coeff_im = result_im.select('coefficients').arrayProject([0]).arrayFlatten([['c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7']]);
    var coeff_im = coeff_im.setDefaultProjection('EPSG:4326', transform_new);
    var coeff_im = coeff_im.reproject({crs:proj.crs(), crsTransform:transform_new});
    var top_im = rSquare_im.divide(k);
    var bot_im = ee.Image(ee.Image(1).subtract(rSquare_im)).divide(n.subtract(k).subtract(1));
    var f_im = top_im.divide(bot_im);
    var zero_im = rmsr_im.lt(0.0);
    var ninenine_im = zero_im.where(f_im.gte(3.630), 1);
    var ninefive_im = zero_im.where(f_im.gte(2.503), 1);
    var ninezero_im = zero_im.where(f_im.gte(2.030), 1);
    var conf_im = zero_im.add(ninenine_im).add(ninefive_im).add(ninezero_im);
    var conf_im = conf_im.setDefaultProjection('EPSG:4326', transform_new);
    var conf_im = conf_im.reproject({crs:proj.crs(), crsTransform:transform_new});
    function predictionA_fn(imobj){
      var regr_im = ee.Image(imobj);
      var c1 = coeff_im.select('c1');
      var c2 = regr_im.select('ppt_sum').multiply(coeff_im.select('c2'));
      var c3 = regr_im.select('tmax_mean').multiply(coeff_im.select('c3'));
      var c4 = regr_im.select('tmin_mean').multiply(coeff_im.select('c4'));
      var c5 = regr_im.select('tdmean_mean').multiply(coeff_im.select('c5'));
      var c6 = regr_im.select('vpdmax_mean').multiply(coeff_im.select('c6'));
      var c7 = regr_im.select('vpdmin_mean').multiply(coeff_im.select('c7'));
      var pred_im = c1.add(c2).add(c3).add(c4).add(c5).add(c6).add(c7);
      var pred_im = pred_im.setDefaultProjection('EPSG:4326', transform_new);
      var pred_im = pred_im.reproject({crs:proj.crs(), crsTransform:transform_new});
      return pred_im;
    }
    var predictionA_ic = ee.ImageCollection(regr_ic.map(predictionA_fn));
  }else if (model_selection == model_list[1]){
    var regr_ic = regr_ic.select(['constant', 'ppt_sum', 'tmean_mean', band]);
    var result_im = regr_ic.reduce(ee.Reducer.linearRegression({numX:3, numY:1}));
    var rmsr_im = result_im.select('residuals').arrayProject([0]).arrayFlatten([['rmsr']]);
    var rss_im = rmsr_im.pow(2).multiply(n);
    var k = ee.Number(2);
    var dof = n.subtract(k).subtract(1);
    var sSquare_im = rss_im.divide(dof);
    var yVariance_im = regr_ic.select(band).reduce(ee.Reducer.sampleVariance());
    var rSquare_im = ee.Image(1).subtract(ee.Image(rss_im.divide(yVariance_im.multiply(n.subtract(1)))));
    var rSquareAdj_im = ee.Image(1).subtract(sSquare_im.divide(yVariance_im));
    var coeff_im = result_im.select('coefficients').arrayProject([0]).arrayFlatten([['c1', 'c2', 'c3']]);
    var coeff_im = coeff_im.setDefaultProjection('EPSG:4326', transform_new);
    var coeff_im = coeff_im.reproject({crs:proj.crs(), crsTransform:transform_new});
    var top_im = rSquare_im.divide(k);
    var bot_im = ee.Image(ee.Image(1).subtract(rSquare_im)).divide(n.subtract(k).subtract(1));
    var f_im = top_im.divide(bot_im);
    var zero_im = rmsr_im.lt(0.0);
    var ninenine_im = zero_im.where(f_im.gte(5.248), 1);
    var ninefive_im = zero_im.where(f_im.gte(3.259), 1);
    var ninezero_im = zero_im.where(f_im.gte(2.456), 1);
    var conf_im = zero_im.add(ninenine_im).add(ninefive_im).add(ninezero_im);
    var conf_im = conf_im.setDefaultProjection('EPSG:4326', transform_new);
    var conf_im = conf_im.reproject({crs:proj.crs(), crsTransform:transform_new});
    function predictionB_fn(imobj){
      var regr_im = ee.Image(imobj);
      var c1 = coeff_im.select('c1');
      var c2 = regr_im.select('ppt_sum').multiply(coeff_im.select('c2'));
      var c3 = regr_im.select('tmean_mean').multiply(coeff_im.select('c3'));
      var pred_im = c1.add(c2).add(c3);
      var pred_im = pred_im.setDefaultProjection('EPSG:4326', transform_new);
      var pred_im = pred_im.reproject({crs:proj.crs(), crsTransform:transform_new});
      return pred_im;
    }
    var predictionB_ic = ee.ImageCollection(regr_ic.map(predictionB_fn));
  }else if (model_selection == model_list[2]){
    var regr_ic = regr_ic.select(['constant', 'ppt_sum', band]);
    var result_im = regr_ic.reduce(ee.Reducer.linearRegression({numX:2, numY:1}));
    var rmsr_im = result_im.select('residuals').arrayProject([0]).arrayFlatten([['rmsr']]);
    var rss_im = rmsr_im.pow(2).multiply(n);
    var k = ee.Number(1);
    var dof = n.subtract(k).subtract(1);
    var sSquare_im = rss_im.divide(dof);
    var yVariance_im = regr_ic.select(band).reduce(ee.Reducer.sampleVariance());
    var rSquare_im = ee.Image(1).subtract(ee.Image(rss_im.divide(yVariance_im.multiply(n.subtract(1)))));
    var rSquareAdj_im = ee.Image(1).subtract(sSquare_im.divide(yVariance_im));
    var coeff_im = result_im.select('coefficients').arrayProject([0]).arrayFlatten([['c1', 'c2']]);
    var coeff_im = coeff_im.setDefaultProjection('EPSG:4326', transform_new);
    var coeff_im = coeff_im.reproject({crs:proj.crs(), crsTransform:transform_new});
    var top_im = rSquare_im.divide(k);
    var bot_im = ee.Image(ee.Image(1).subtract(rSquare_im)).divide(n.subtract(k).subtract(1));
    var f_im = top_im.divide(bot_im);
    var zero_im = rmsr_im.lt(0.0);
    var ninenine_im = zero_im.where(f_im.gte(7.373), 1);
    var ninefive_im = zero_im.where(f_im.gte(4.105), 1);
    var ninezero_im = zero_im.where(f_im.gte(2.846), 1);
    var conf_im = zero_im.add(ninenine_im).add(ninefive_im).add(ninezero_im);
    var conf_im = conf_im.setDefaultProjection('EPSG:4326', transform_new);
    var conf_im = conf_im.reproject({crs:proj.crs(), crsTransform:transform_new});
    function predictionC_fn(imobj){
      var regr_im = ee.Image(imobj);
      var c1 = coeff_im.select('c1');
      var c2 = regr_im.select('ppt_sum').multiply(coeff_im.select('c2'));
      var pred_im = c1.add(c2);
      var pred_im = pred_im.setDefaultProjection('EPSG:4326', transform_new);
      var pred_im = pred_im.reproject({crs:proj.crs(), crsTransform:transform_new});
      return pred_im;
    }
    var predictionC_ic = ee.ImageCollection(regr_ic.map(predictionC_fn));
  }
  
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
  var res_im = reg_ic.reduce(ee.Reducer.linearRegression({numX:2, numY:1}));
  var rmse_im = res_im.select('residuals').arrayProject([0]).arrayFlatten([['rmsr']]);
  var top_im = ee.Image(rmse_im.pow(2).multiply(n).divide(n.subtract(2))).pow(0.5);
  var bot_im = reg_ic.select('constant_1').reduce(ee.Reducer.sampleVariance()).multiply(n.subtract(1)).pow(0.5);
  var sslp_im = top_im.divide(bot_im);
  var coef_im = res_im.select('coefficients').arrayProject([0]).arrayFlatten([['c1', 'c2']]);
  var coef_im = coef_im.setDefaultProjection('EPSG:4326', transform_new);
  var coef_im = coef_im.reproject({crs:proj.crs(), crsTransform:transform_new});
  var t_im = coef_im.select('c2').divide(sslp_im);
  var zer_im = rmse_im.lt(0.0);
  var ninenin_im = zer_im.where(t_im.gte(2.715), 1);
  var ninefiv_im = zer_im.where(t_im.gte(2.026), 1);
  var ninezer_im = zer_im.where(t_im.gte(1.687), 1);
  var con_im = zer_im.add(ninenin_im).add(ninefiv_im).add(ninezer_im);
  var con_im = con_im.setDefaultProjection('EPSG:4326', transform_new);
  var con_im = con_im.reproject({crs:proj.crs(), crsTransform:transform_new});

  if (out_im_type == 'avg'){
    var avg_im = rap_ic.reduce(ee.Reducer.mean());
    return avg_im;
  }else if (out_im_type == 'sdev'){
    var avg_im = rap_ic.reduce(ee.Reducer.mean());
    var sdev_im = rap_ic.reduce(ee.Reducer.sampleStdDev());
    var sdev_im = sdev_im.divide(avg_im).multiply(100.0);
    return sdev_im;
  }else if (out_im_type == 'rmse'){
    return rmsr_im;
  }else if (out_im_type == 'rsqr'){
    return rSquare_im;
  }else if (out_im_type == 'rsqrA'){
    return rSquareAdj_im;
  }else if (out_im_type == 'Fconf'){
    return conf_im;
  }else if (out_im_type == 'Tconf'){
    return con_im;    
  }else if (out_im_type.slice(0, 5) == 'coeff'){
    if (out_im_type.indexOf('trend') < 0){
      var coeff_str = out_im_type.slice(-2);
      var coeff_idx = coeff_list.indexOf(coeff_str);
      var coeff_str = coeff_internal_list[coeff_idx];
      var term_im = ee.Image(coeff_im.select(coeff_str));
      var term_im = term_im.setDefaultProjection('EPSG:4326', transform_new);
      var term_im = term_im.reproject({crs:proj.crs(), crsTransform:transform_new});      
      return term_im;
    }else{
      var coeff_str = out_im_type.slice(-6);
      var coeff_idx = coeff_list.indexOf(coeff_str);
      var coeff_str = coeff_internal_list[coeff_idx];
      var term_im = ee.Image(coef_im.select(coeff_str))
      var term_im = term_im.setDefaultProjection('EPSG:4326', transform_new);
      var term_im = term_im.reproject({crs:proj.crs(), crsTransform:transform_new});  
      return term_im;
    }
  }else if (out_im_type.slice(0, 7) == 'covpred' || out_im_type.slice(0, 7) == 'propred'|| out_im_type.slice(0, 7) == 'agbpred'){
    if (model_selection == model_list[0]){
      var prediction_ic = predictionA_ic;
    }else if (model_selection == model_list[1]){
      var prediction_ic = predictionB_ic;
    }else if (model_selection == model_list[2]){
      var prediction_ic = predictionC_ic;
    }
    var pred_ic_list = prediction_ic.toList(999);
    var year = ee.Number.parse(out_im_type.slice(7));
    var year_idx = years_list.indexOf(year);
    var pred_im = ee.Image(pred_ic_list.get(year_idx));
    return pred_im;
  }else if (out_im_type.slice(0, 3) == 'cov' || out_im_type.slice(0, 3) == 'pro' || out_im_type.slice(0, 3) == 'agb'){
    var year = ee.Number.parse(out_im_type.slice(3));
    var year_idx = years_list.indexOf(year);
    var rap_im = ee.Image(rap_ic_list.get(year_idx));
    return rap_im;
  }else if (out_im_type == 'OnetoOne'){
    if (model_selection == model_list[0]){
      var prediction_ic = predictionA_ic;
    }else if (model_selection == model_list[1]){
      var prediction_ic = predictionB_ic;
    }else if (model_selection == model_list[2]){
      var prediction_ic = predictionC_ic;
    }
    var oneone_ic = rap_ic.merge(prediction_ic);
    return oneone_ic;
  }else if (out_im_type == 'Trend'){
    return rap_ic;
  }else if (out_im_type == 'Debug'){
    print('DEBUG OUTPUT');
    return t_im;
  }
}

////////////////////////////////////////////////////
//START DeBugGing TEsTING deBUggiNg tESTiNG TEsTiNG
//var t_im = main_fn('PFG', cover_ic, 'Debug');
//Map.addLayer(t_im, {min:0, max:2})

// var regr_props = regr_ic.getRegion(geometry, scale);
// print(regr_props);
// var regr_props = regr_props.slice(1);

// function make_prop_feats_fnA(p_list){
//   var p_list = ee.List(p_list);
//   var a = p_list.get(4);
//   var b = p_list.get(5);
//   var c = p_list.get(6);
//   //var e = p_list.get(8);
//   //var f = p_list.get(9);
//   //var g = p_list.get(10);
//   //var ft = ee.Feature(null, {a:a, b:b, c:c, d:d, e:e, f:f, g:g});
//   var ft = ee.Feature(null, {a:a, b:b, c:c});
//   return ft;
// }

// var out_fc = ee.FeatureCollection(regr_props.map(make_prop_feats_fnA));

// Export.table.toDrive({
//   collection:out_fc,
//   description:'merge_ic'
// });

// var pred_ic = main_fn('PFG', cover_ic, 'Debug');
// var pred_props = pred_ic.getRegion(geometry, scale);
// var pred_props = pred_props.slice(1);

// function make_prop_feats_fnB(p_list){
//   var p_list = ee.List(p_list);
//   var a = p_list.get(4);
//   var ft = ee.Feature(null, {a:a});
//   return ft;
// }

// var pred_fc = ee.FeatureCollection(pred_props.map(make_prop_feats_fnB));

// Export.table.toDrive({
//   collection:pred_fc,
//   description:'vectorsPred'
// });
//END DeBugGing TEsTING dEbuGgINg tESTiNG TEsTiNG
/////////////////////////////////////////////////


///////////////////////
///////////////////////
//USER INTERFACE
///////////////////////
///////////////////////


//////////////////////
//Styles

//https://github.com/gee-community/ee-palettes
var palettes = require('users/gena/packages:palettes');
var rmseCovVis = {
  min:0,
  max:5,
  palette:palettes.kovesi.diverging_linear_bjr_30_55_c53[7]};

var palettes = require('users/gena/packages:palettes');
var rmseProVis = {
  min:0,
  max:200,
  palette:palettes.kovesi.diverging_linear_bjr_30_55_c53[7]};

var palettes = require('users/gena/packages:palettes');
var rmseBioVis = {
  min:0,
  max:250,
  palette:palettes.kovesi.diverging_linear_bjr_30_55_c53[7]};

var palettes = require('users/gena/packages:palettes');
var rsqrVis = {
  min:0,
  max:0.6,
  palette:palettes.kovesi.diverging_linear_bjr_30_55_c53[7].reverse()};

var palettes = require('users/gena/packages:palettes');
var rsqrAdjVis = {
  min:0,
  max:0.3,
  palette:palettes.kovesi.diverging_linear_bjr_30_55_c53[7].reverse()};

var palettes = require('users/gena/packages:palettes');
var covVis = {
  min:0,
  max:50,
  palette:palettes.niccoli.cubicl[7]};

var palettes = require('users/gena/packages:palettes');
var proVis = {
  min:0,
  max:500,
  palette:palettes.niccoli.cubicl[7]};

var palettes = require('users/gena/packages:palettes');
var bioVis = {
  min:0,
  max:500,
  palette:palettes.niccoli.cubicl[7]};

var confVis = {
  min:0,
  max:3,
  palette:['#FFFFFF', '#d7481d', '#59f720', '#800080']};

var widgetStyle = {
  position:'bottom-center',
  padding:'0px 0px',
  backgroundColor:'rgba(255, 255, 255, 0.7)'};

var checkStyle = {
  position:'bottom-center',
  padding:'0px 0px',
  backgroundColor:'rgba(255, 255, 255, 0.7)',
  fontSize:'12px'};

var mainPanelStyle = {
  position:'top-right',
  padding:'0px 0px', 
  backgroundColor:'rgba(255, 255, 255, 0.7)', 
  border:'1px solid black'};

var headerStyle = {
  padding:'0px 0px',
  backgroundColor:'rgba(255, 255, 255, 0.7)', 
  fontSize:'12px'}

var chartPanelStyle = {
  position:'bottom-center', 
  stretch:'vertical',
  height:'400px',
  width:'400px',
  margin:'10px 10px'};

var pixelLabelStyle = {
  height:'50px',
  width:'50px',
  padding:'0px',
  margin:'0px',
  position:'top-center',
  fontSize:'16px',
  fontWeight: 'bold'};

var pixelPanelStyle = {
  height:'100px',
  width:'100px',
  position:'bottom-center',
  margin:'10px 10px'};

var infoLabelStyle = {
  height:'1400px',
  width:'1000px',
  position:'bottom-center',
  whiteSpace:'preserve nowrap',
  padding:'1px',
  margin:'2px',
  textAlign:'left',
  fontSize:'12px'};

var textPanelStyle = {
  height:'500px',
  width:'700px',
  position:'bottom-center', 
  stretch:'vertical',
  margin:'10px 10px'};

var covcoeff_map = {};

covcoeff_map[coeff_list[0]] = {min:-0.3, max:0.3};
covcoeff_map[coeff_list[1]] = {min:-10, max:50};
covcoeff_map[coeff_list[2]] = {min:-0.02, max:0.02};
covcoeff_map[coeff_list[3]] = {min:-10, max:10};
covcoeff_map[coeff_list[4]] = {min:-20, max:20};
covcoeff_map[coeff_list[5]] = {min:-10, max:10};
covcoeff_map[coeff_list[6]] = {min:-5, max:5};
covcoeff_map[coeff_list[7]] = {min:-50, max:50};
covcoeff_map[coeff_list[8]] = {min:-0.1, max:0.1};
covcoeff_map[coeff_list[9]] = {min:-10, max:10};
covcoeff_map[coeff_list[10]] = {min:-0.05, max:0.05};
covcoeff_map[coeff_list[11]] = {min:-50, max:50};
covcoeff_map[coeff_list[12]] = {min:-50, max:50};
covcoeff_map[coeff_list[13]] = {min:-50, max:50};

var procoeff_map = {};

procoeff_map[coeff_list[0]] = {min:-5, max:5};
procoeff_map[coeff_list[1]] = {min:-100, max:500};
procoeff_map[coeff_list[2]] = {min:-0.5, max:0.5};
procoeff_map[coeff_list[3]] = {min:-100, max:100};
procoeff_map[coeff_list[4]] = {min:-200, max:200};
procoeff_map[coeff_list[5]] = {min:-100, max:100};
procoeff_map[coeff_list[6]] = {min:-50, max:50};
procoeff_map[coeff_list[7]] = {min:-500, max:500};
procoeff_map[coeff_list[8]] = {min:-0.5, max:0.5};
procoeff_map[coeff_list[9]] = {min:-100, max:100};
procoeff_map[coeff_list[10]] = {min:-0.5, max:0.5};
procoeff_map[coeff_list[11]] = {min:-1000, max:1000};
procoeff_map[coeff_list[12]] = {min:-1000, max:1000};
procoeff_map[coeff_list[13]] = {min:-1000, max:1000};

var biocoeff_map = {};

biocoeff_map[coeff_list[0]] = {min:-5, max:5};
biocoeff_map[coeff_list[1]] = {min:-100, max:500};
biocoeff_map[coeff_list[2]] = {min:-1, max:1};
biocoeff_map[coeff_list[3]] = {min:-100, max:100};
biocoeff_map[coeff_list[4]] = {min:-200, max:200};
biocoeff_map[coeff_list[5]] = {min:-100, max:100};
biocoeff_map[coeff_list[6]] = {min:-50, max:50};
biocoeff_map[coeff_list[7]] = {min:-500, max:500};
biocoeff_map[coeff_list[8]] = {min:-1, max:1};
biocoeff_map[coeff_list[9]] = {min:-100, max:100};
biocoeff_map[coeff_list[10]] = {min:-1, max:1};
biocoeff_map[coeff_list[11]] = {min:-1000, max:1000};
biocoeff_map[coeff_list[12]] = {min:-1000, max:1000};
biocoeff_map[coeff_list[13]] = {min:-1000, max:1000};

/////////////////////////////////////////
//Global Widget Vars and Initial Display

Map.setCenter(-109, 37, 5);

var palettes = require('users/gena/packages:palettes');
var im_to_show = main_fn(band_selection, ic_selection, type_selection);
var bandVis = covVis;
Map.addLayer(im_to_show, bandVis);

var main_panel = ui.Panel({
  layout:ui.Panel.Layout.flow('vertical'),
  style:mainPanelStyle
});

main_panel.add(ui.Label({value:'RAP Climate Sensitivity Viewer', style:{padding:'0px 0px', fontWeight:'bold'}}));

var legend_panel = ui.Panel({
  style:{position:'bottom-left', padding:'8px 15px'}
});

var chart_panelA = ui.Panel({style:chartPanelStyle});
var chart_panelB = ui.Panel({style:chartPanelStyle});
var pixel_panel = ui.Panel({style:pixelPanelStyle});

/////////////////
//Legend function

function makeLegend(){

  Map.remove(legend_panel);
  legend_panel.clear();
  
  if (type_selection == 'Fconf' || type_selection == 'Tconf'){
    var pos = 'middle-left';
  }else{
    var pos = 'bottom-left';
  }

  legend_panel = ui.Panel({
    style:{position:pos, padding:'8px 15px'}
  });

  if ((type_selection.slice(0, 3) == 'cov') || (type_selection.slice(0, 3) == 'avg' && cover_bands.indexOf(band_selection) >= 0) || (type_selection == 'rmse' && cover_bands.indexOf(band_selection) >= 0)){
    var lTitle = 'frac. %';
  }else if ((type_selection.slice(0, 3) == 'pro') || (type_selection.slice(0, 3) == 'avg' && prod_bands.indexOf(band_selection) >= 0) || (type_selection == 'rmse' && prod_bands.indexOf(band_selection) >= 0)){
    var lTitle = 'lbs/acre';
  }else if ((type_selection.slice(0, 3) == 'agb') || (type_selection.slice(0, 3) == 'avg' && bio_bands.indexOf(band_selection) >= 0) || (type_selection == 'rmse' && bio_bands.indexOf(band_selection) >= 0)){
    var lTitle = 'lbs/acre';
  }else if (type_selection.slice(0, 3) == 'sde'){
    var lTitle = '% of Avg.';
  }else if (type_selection == 'rsqr'){
    var lTitle = 'R²';
  }else if (type_selection == 'rsqrA'){
    var lTitle = 'R² Adj.';
  }else if (type_selection.indexOf('coeffK') >= 0){
    var lTitle = 'Intercept (y)';
  }else{
    var lTitle = 'Slope (y/x)';
  }

  var legendTitle = ui.Label({
    value:lTitle,
    style:{fontWeight:'bold', fontSize:'16px', margin:'0 0 4px 0', padding:'0'}
  });

  var panel = ui.Panel({
    widgets:[ui.Label('>= '.concat(String(bandVis['max'])))]
  });
  
  var lat = ee.Image.pixelLonLat().select('latitude');
  var gradient = lat.multiply((bandVis['max']-bandVis['min'])/100.0).add(bandVis['min']);
  var legendImage = gradient.visualize(bandVis);

  var thumbnail = ui.Thumbnail({
    image:legendImage,
    params:{bbox:'0,0,10,100', dimensions:'10x200'}, 
    style:{position:'bottom-center', stretch:'vertical', margin:'0px 8px', maxHeight:'200px'}
  });
  
  if (['cov', 'pro', 'agb', 'rms', 'avg', 'sde'].indexOf(type_selection.slice(0, 3)) < 0){
    var panel2 = ui.Panel({
      widgets:[ui.Label('<= '.concat(String(bandVis['min'])))]
    });
  }else{
    var panel2 = ui.Panel({
      widgets:[ui.Label(String(bandVis['min']))]
    });
  }

  if (type_selection == 'Fconf' || type_selection == 'Tconf'){

    var typeLabels = ['<90% conf.', '>=90% conf.', '>=95% conf.', '>=99% conf.'];
  
    for (var i = 0; i < [0, 1, 2, 3].length; i++){
      var colorBox = ui.Label({
        style:{
          backgroundColor:bandVis.palette[i],
          padding:'10px',
          margin:'0px',
          border:'1px solid black',
          fontSize:'10px'}
      });
    
      var label = ui.Label({
        value:typeLabels[i],
        style:{
          padding:'6px',
          margin:'0px',
          position:'middle-left',
          fontSize:'12px'}
      });
    
      var row = ui.Panel({
        widgets:[colorBox, label],
        layout:ui.Panel.Layout.Flow('horizontal')
      });
    
      legend_panel.add(row);
    }
  }
  
  if (type_selection != 'Fconf' && type_selection != 'Tconf'){
    legend_panel.add(legendTitle);
    legend_panel.add(panel);
    legend_panel.add(thumbnail);
    legend_panel.add(panel2);
  }

  Map.add(legend_panel);
}

//////////////////////
//Change Climate Model

function renderModel(model_string){
  Map.layers().reset();
  model_selection = model_string;
  im_to_show = main_fn(band_selection, ic_selection, type_selection);
  Map.addLayer(im_to_show, bandVis);
  makeLegend();
}

var model_dropdown = ui.Select({
  items:model_list, 
  placeholder:'Select Climate Regression Model', 
  onChange:renderModel,
  style:widgetStyle
});

main_panel.add(ui.Label({value:'Select Climate Regression Model', style:headerStyle}))
main_panel.add(model_dropdown);

/////////////////////
//Change RAP Variable

function renderVariable(var_selection){
  Map.layers().reset();
  band_selection = var_selection;
  if (cover_bands.indexOf(band_selection) >= 0){
    ic_selection = cover_ic;
    bandVis = covVis;
    if (type_selection.slice(0, 3) == 'pro'){
      var year = type_selection.slice(-4);
      type_selection = type_selection.replace('pro', 'cov').slice(0, -4).concat(year);
    }   
    if (type_selection.slice(0, 3) == 'agb'){
      var year = type_selection.slice(-4);
      type_selection = type_selection.replace('agb', 'cov').slice(0, -4).concat(year);
    }    
  }else if (prod_bands.indexOf(band_selection) >= 0){
    ic_selection = prod_ic;
    bandVis = proVis;
    if (type_selection.slice(0, 3) == 'cov'){
      var year = type_selection.slice(-4);
      type_selection = type_selection.replace('cov', 'pro').slice(0, -4).concat(year);
    }
    if (type_selection.slice(0, 3) == 'agb'){
      var year = type_selection.slice(-4);
      type_selection = type_selection.replace('agb', 'pro').slice(0, -4).concat(year);
    }
  }else if (bio_bands.indexOf(band_selection) >= 0){
    ic_selection = agb_ic;
    bandVis = bioVis;
    if (type_selection.slice(0, 3) == 'cov'){
      var year = type_selection.slice(-4);
      type_selection = type_selection.replace('cov', 'agb').slice(0, -4).concat(year);
    }
    if (type_selection.slice(0, 3) == 'pro'){
      var year = type_selection.slice(-4);
      type_selection = type_selection.replace('pro', 'agb').slice(0, -4).concat(year);
    }
  }
  im_to_show = main_fn(band_selection, ic_selection, type_selection);
  if (type_selection.slice(0, 3) == 'cov'){
    bandVis = covVis;
  }else if (type_selection.slice(0, 3) == 'pro'){
    bandVis = proVis;
  }else if (type_selection.slice(0, 3) == 'agb'){
    bandVis = bioVis;
  }else if (type_selection.slice(0, 7) == 'covpred'){
    bandVis = covVis;
  }else if (type_selection.slice(0, 7) == 'propred'){
    bandVis = proVis;
  }else if (type_selection.slice(0, 7) == 'agbpred'){
    bandVis = agbVis;
  }else if (type_selection == 'rmse' && cover_bands.indexOf(band_selection) >= 0){
    bandVis = rmseCovVis;
  }else if (type_selection == 'rmse' && prod_bands.indexOf(band_selection) >= 0){
    bandVis = rmseProVis;
  }else if (type_selection == 'rmse' && bio_bands.indexOf(band_selection) >= 0){
    bandVis = rmseBioVis;
  }else if (type_selection == 'rsqr'){
    bandVis = rsqrVis;
  }else if (type_selection == 'rsqrA'){
    bandVis = rsqrAdjVis;
  }else if (type_selection == 'Fconf'){
    bandVis = confVis;
  }else if (type_selection == 'Tconf'){
    bandVis = confVis;
  }else if (type_selection.indexOf('coeff') >= 0 && cover_bands.indexOf(band_selection) >= 0){
    var coeff_str = type_selection.slice(5);
    var im_max = covcoeff_map[coeff_str]['max'];
    var im_min = covcoeff_map[coeff_str]['min'];
    var palettes = require('users/gena/packages:palettes');
    bandVis = {min:im_min, max:im_max, palette:palettes.colorbrewer.BrBG[7]};
  }else if (type_selection.indexOf('coeff') >= 0 && prod_bands.indexOf(band_selection) >= 0){
    var coeff_str = type_selection.slice(5);
    var im_max = procoeff_map[coeff_str]['max'];
    var im_min = procoeff_map[coeff_str]['min'];
    var palettes = require('users/gena/packages:palettes');
    bandVis = {min:im_min, max:im_max, palette:palettes.colorbrewer.BrBG[7]};
  }else if (type_selection.indexOf('coeff') >= 0 && bio_bands.indexOf(band_selection) >= 0){
    var coeff_str = type_selection.slice(5);
    var im_max = biocoeff_map[coeff_str]['max'];
    var im_min = biocoeff_map[coeff_str]['min'];
    var palettes = require('users/gena/packages:palettes');
    bandVis = {min:im_min, max:im_max, palette:palettes.colorbrewer.BrBG[7]};
  }
  Map.addLayer(im_to_show, bandVis);
  makeLegend();
}

var variable_dropdown = ui.Select({
  items:cover_bands.concat(prod_bands).concat(bio_bands), 
  placeholder:'Select Variable', 
  onChange:renderVariable,
  style:widgetStyle
});

main_panel.add(ui.Label({value:'Selected RAP Variable', style:headerStyle}))
main_panel.add(variable_dropdown);

/////////////////////
//Change Metric Layer

function renderMetric(metric_selection){
  Map.layers().reset();
  type_selection = metric_selection;
  im_to_show = main_fn(band_selection, ic_selection, type_selection);

  if (type_selection == 'avg' && cover_bands.indexOf(band_selection) >= 0){
    bandVis = covVis;
  }else if (type_selection == 'avg' && prod_bands.indexOf(band_selection) >= 0){
    bandVis = proVis;
  }else if (type_selection == 'avg' && bio_bands.indexOf(band_selection) >= 0){
    bandVis = bioVis;
  }else if (type_selection == 'sdev' && cover_bands.indexOf(band_selection) >= 0){
    bandVis = covVis;
  }else if (type_selection == 'sdev' && prod_bands.indexOf(band_selection) >= 0){
    bandVis = proVis;
  }else if (type_selection == 'sdev' && bio_bands.indexOf(band_selection) >= 0){
    bandVis = bioVis;
  }else if (type_selection == 'rmse' && cover_bands.indexOf(band_selection) >= 0){
    bandVis = rmseCovVis;
  }else if (type_selection == 'rmse' && prod_bands.indexOf(band_selection) >= 0){
    bandVis = rmseProVis;
  }else if (type_selection == 'rmse' && bio_bands.indexOf(band_selection) >= 0){
    bandVis = rmseBioVis;
  }else if (type_selection == 'rsqr'){
    bandVis = rsqrVis;
  }else if (type_selection == 'rsqrA'){
    bandVis = rsqrAdjVis;
  }else if (type_selection == 'Fconf'){
    bandVis = confVis;
  }else if (type_selection == 'Tconf'){
    bandVis = confVis;
  }
  Map.addLayer(im_to_show, bandVis);
  makeLegend();
}

var metric_dropdown = ui.Select({
  items:metric_strs, 
  placeholder:'Select Metric', 
  onChange:renderMetric,
  style:widgetStyle
});

main_panel.add(ui.Label({value:'View Metric Map', style:headerStyle}))
main_panel.add(metric_dropdown);

///////////////////////////
//Render Coefficient Layer

function renderCoeff(coeff_str){
  Map.layers().reset();
  type_selection = 'coeff'.concat(coeff_str);
  if (cover_bands.indexOf(band_selection) >= 0){
    var coeff_map = covcoeff_map;
  }else if (prod_bands.indexOf(band_selection) >= 0){
    var coeff_map = procoeff_map;
  }else if (bio_bands.indexOf(band_selection) >= 0){
    var coeff_map = biocoeff_map;    
  }
  for (var i = 0; i < model_list.length; i++){
    if (model_list[i].indexOf(coeff_str) >= 0){
      model_selection = model_list[i]
    }
  }
  var palettes = require('users/gena/packages:palettes');
  var im_max = coeff_map[coeff_str]['max'];
  var im_min = coeff_map[coeff_str]['min'];
  bandVis = {min:im_min, max:im_max, palette:palettes.colorbrewer.BrBG[7]};
  im_to_show = main_fn(band_selection, ic_selection, type_selection);
  Map.addLayer(im_to_show, bandVis);
  makeLegend();
}

var coeff_dropdown = ui.Select({
  items:coeff_list, 
  placeholder:'Select Coefficient Term', 
  onChange:renderCoeff,
  style:widgetStyle
});

main_panel.add(ui.Label({value:'View Regression Coefficient Map', style:headerStyle}));
main_panel.add(coeff_dropdown);

///////////////////
//Render RAP Layer

function renderYearA(year_selection){
  Map.layers().reset();
  if (cover_bands.indexOf(band_selection) != -1){
    type_selection = 'cov'.concat(year_selection);
  }else if (prod_bands.indexOf(band_selection) != -1){
    type_selection = 'pro'.concat(year_selection);
  }else if (bio_bands.indexOf(band_selection) != -1){
    type_selection = 'agb'.concat(year_selection);
  }
  im_to_show = main_fn(band_selection, ic_selection, type_selection);
  if (type_selection.slice(0, 3) == 'cov'){
    bandVis = covVis;
  }else if (type_selection.slice(0, 3) == 'pro'){
    bandVis = proVis;
  }else if (type_selection.slice(0, 3) == 'agb'){
    bandVis = bioVis;
  }
  Map.addLayer(im_to_show, bandVis);
  makeLegend();
}

var year_dropdownA = ui.Select({
  items:years_str_list, 
  placeholder:'Select RAP Year', 
  onChange:renderYearA,
  style:widgetStyle
});

main_panel.add(ui.Label({value:'View RAP year', style:headerStyle}));
main_panel.add(year_dropdownA);

/////////////////////////
//Render Predicted Layer

function renderYearB(year_selection){
  Map.layers().reset();
  if (cover_bands.indexOf(band_selection) != -1){
    type_selection = 'covpred'.concat(year_selection);
  }else if (prod_bands.indexOf(band_selection) != -1){
    type_selection = 'propred'.concat(year_selection);
  }else if (bio_bands.indexOf(band_selection) != -1){
    type_selection = 'agbpred'.concat(year_selection);
  }
  im_to_show = main_fn(band_selection, ic_selection, type_selection);
  if (type_selection.slice(0, 7) == 'covpred'){
    bandVis = covVis;
  }else if (type_selection.slice(0, 7) == 'propred'){
    bandVis = proVis;
  }else if (type_selection.slice(0, 7) == 'agbpred'){
    bandVis = bioVis;
  }
  Map.addLayer(im_to_show, bandVis);
  makeLegend();
}

var year_dropdownB = ui.Select({
  items:years_str_list, 
  placeholder:'Select Predicted Year', 
  onChange:renderYearB,
  style:widgetStyle
});

main_panel.add(ui.Label({value:'View Predicted RAP year', style:headerStyle}));
main_panel.add(year_dropdownB);

//////////////////////////
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
    title:'Pred. vs. True',
    titleTextStyle:{italic:false, bold:true, fontSize:26},
    legend:{position:'top-right'},
    hAxis:{viewWindow:{min:min_value*0.9, max:max_value*1.1}, title:'RAP', titleTextStyle:{italic:false, bold:true, fontSize:21}, gridlines:{count:6}},
    vAxis:{viewWindow:{min:min_value*0.9, max:max_value*1.1}, title:'PRED', titleTextStyle:{italic:false, bold:true, fontSize:21}, gridlines:{count:6}},
    colors:['#6a9f58'],
    pointSize:12,
    lineSize:0,
    chartArea:{height:400, width:400},
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
  label:'Pred. vs. True Plot on Click',
  onChange:renderOnetoOne,
  style:checkStyle
});

main_panel.add(OnetoOne_checkbox);

/////////////////////
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
    hAxis:{title:'Years (1986-2024)', format:'####', titleTextStyle:{italic:false, bold:true, fontSize:21}, gridlines:{count:6}},
    vAxis:{title:'RAP', titleTextStyle:{italic:false, bold:true, fontSize:21}, gridlines:{count:6}},
    colors:['#6a9f58'],
    pointSize:12,
    lineSize:0,
    chartArea:{height:400, width:400},
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

/////////////////////
//Render Pixel Value

function makePlotC(point_geo){
  var point_fc = im_to_show.sample(point_geo, scale, proj);
  var prop_ft = point_fc.first();
  var prop_names = prop_ft.propertyNames();
  var sys_index_idx = prop_names.getInfo().indexOf('system:index');
  var prop_idx = 1 - sys_index_idx;
  var prop_val = prop_ft.get(prop_names.get(prop_idx));
  var pixel_label = ui.Label({
    value:'Pixel Value:'.concat('\n').concat(prop_val.getInfo()),
    style:pixelLabelStyle
  });
  return pixel_label;
}

function clickCallbackC(clickInfo_obj){
  Map.remove(pixel_panel);
  var lat = clickInfo_obj.lat;
  var lon = clickInfo_obj.lon;
  var pt = ee.Geometry.Point([lon, lat]);
  var pixel_label = makePlotC(pt);
  pixel_panel = ui.Panel({
    widgets:pixel_label,
    layout:ui.Panel.Layout.Flow('horizontal'),
    style:pixelPanelStyle
  });
  Map.add(pixel_panel);
}

function renderInspector(checkbox_bool){
  if (checkbox_bool === true){
    Map.onClick(clickCallbackC);
  }
  else{
    Map.unlisten();
    Map.remove(pixel_panel);
    pixel_panel = ui.Panel({style:pixelPanelStyle});
  }
}

var inspect_checkbox = ui.Checkbox({
  label:'Pixel Value on Click',
  onChange:renderInspector,
  style:checkStyle
});

main_panel.add(inspect_checkbox);

//////////////////
//Render info box

var info_str = 'OVERVIEW: \n' +
              'This app is built using the Google Earth Engine cloud platform and publically available datasets. \n' +
              'The Rangeland Assessment Platform (RAP) dataset of ground cover is visualized for a southwestern US  \n' +
              'coverage area and connections are identified between RAP and climate variables of the PRISM dataset. \n' +
              'Maps are generated at the ~800 m resolution of PRISM, such that RAP is spatially averaged from its \n' +
              'original resolution of 30 m. The sensitivity of RAP to climate is quantified based on the success of \n' +
              'multiple linear regressions that predict individual RAP variables using only climate predictors. \n' +
              'Each pixel is fitted with its own independent regression model. Trend analysis based on simple linear \n' +
              'regression of RAP annual time series can also be visualized for each pixel. This application is meant \n' +
              'to supplement features available on the offical RAP website and to additionally explore climate \n' +
              'sensitivity on rangeland in the southwestern US. \n' +
              '\n' +
              'DEFINITIONS: \n' +
              'RAP: Rangeland Assessment Platform is a dataset with rangeland-specific vegetation growth and cover. \n' + 
              'PRISM: A US gridded observational climate dataset. In this case, monthly ~800 m data is annually averaged. \n' +               
              'Ground Cover: The fraction of surface area covered by plants or other cover forms when viewed from above. \n' +
              'NPP: Net primary production, units of pounds carbon per acre (lbs / acre). \n' + 
              'AGB: Above ground biomass, units of pounds biomass per acre (lbs / acre). \n' + 
              'AFG: Annual forbs and grass ground cover (Frac. %). \n' + 
              'PFG: Perennial forbs and grass ground cover (Frac. %). \n' +
              'BRG: Bare ground cover (Frac. %).\n' +
              'SHR: Shurb ground cover (Frac. %). \n' +               
              'TRE: Tree ground cover (Frac. %). \n' +
              'afgNPP: NPP of annual forbs and grass (lbs/acre). \n' + 
              'pfgNPP: NPP of perennial forbs and grass (lbs/acre). \n' +
              'shrNPP: NPP of shrubs (lbs/acre). \n' +               
              'treNPP: NPP of trees (lbs/acre). \n' +
              'afgAGB: AGB of annual forbs and grass estimated from NPP (lbs/acre). \n' + 
              'pfgAGB: AGB of perennial forbs and grass estimated from NPP (lbs/acre). \n' +
              'LR: linear regression is a statistical model of the relationship between a dependend variable and one \n' + 
              'or more independent variables. \n' + 
              'avg: Average value of the selected RAP variable. \n' +
              'sdev: Standard deviation of the selected RAP variable. Expressed as a percentage fraction of the avg. \n' + 
              'RMSE: Root mean square error of the selected climate regression model. \n' + 
              'Rsqr: R-square of the selected climate regression model. \n' + 
              'RsqrA: Adjusted R-square of the selected climate regression model that penalizes the model for having \n' + 
              'added parameters. \n' + 
              'F-statistic: A metric used to calculate the p-value for overall significance of climate LR models. \n' +
              'Selecting the Fconf metric shows the p-value expressed as a confidence percentage. \n' +
              'T-statistic: A metric used to calculate the p-value for the slope coefficient of simple LR trends for \n' +
              'annual RAP time series. Selecting the Tconf metric shows the p-value expressed as a confidence percentage. \n' +
              '\n' + 
              'METHODOLOGY: \n' +
              'Models based on LR are fitted to predict individual RAP variables using the 39-year annual record \n' + 
              'given by RAP (1986-2024) and annually averaged PRISM variables for the same time period. \n' + 
              'Seven PRISM predictor variables were used: precipitation (pr), mean max/min temperature (tx, tn), \n' + 
              'mean temperature (tm), mean dewpoint temperature (td), and mean max/min vapor pressure deficit (vx, vn). \n' + 
              'Similarly, trend analysis is done using simple LR. The slope coefficient of the trend analysis is \n' +
              'a negative value when there is a downward trend, and positive if the trend is upwards. The magnitude \n' +
              'determines the steepness of the trend. The selectable metric maps allow the performance of the LR models \n' +
              'to be judged. The Fconf and Tconf metric maps indicate whether the fitted LR model parameters are \n' +
              'statistically meaningful and are expressed as a percent confidence. The Fconf maps are generated for \n' +
              'the climate regression models and are based on the F-statistic, which indicates whether the overall \n' +
              'LR model is statistically significant. The Tconf map is determined for the trend analysis LR model \n' +
              'and indicates whether the slope cofficient of the LR is statistically significant. Estimated AGB \n' +
              'for grass and forbs uses the same methodology as the RAP website. \n' +
              '\n' + 
              'USAGE: \n' +
              'Choosing any selection option will render a new map. For options that are not chosen, placeholder options \n' + 
              'are used until a selection is made. Rendered maps of metrics, coeffiecents, RAP years, and predicted \n' + 
              'RAP years depend on the selection of climate regression model and/or RAP variable. If you are unsure \n' + 
              'of what options were used to generate the current map, render a new map with the intended selections. \n' + 
              'Some maps will load slower because all processing is done on-the-fly. Layer transparency can be \n' + 
              'adjusted using the slider bar feature in the layer list (top-right). Rendered graphs can be expanded \n' + 
              'to a larger size by clicking the icon in the top-right corner of the graph. This opens an interactive \n' + 
              'version of the graph in a separate browser tab. Trends in RAP can be visualized by viewing coefficient  \n' + 
              'maps of Atrend and Ktrend, the slope and intercept, respectively, of the LR trend of annual RAP time series.  \n' + 
              '\n' + 
              'CITATIONS: \n' +
              'Daly, C., Halbleib, M., Smith, J. I., Gibson, W. P., Doggett, M. K., Taylor, G. H., ... & \n' + 
              'Pasteris, P. P. (2008). Physiographically sensitive mapping of climatological temperature and precipitation \n' +
              'across the conterminous United States. International Journal of Climatology: a Journal of the \n' + 
              'Royal Meteorological Society, 28(15), 2031-2064. \n' +
              '\n' + 
              'Jones, M. O., Robinson, N. P., Naugle, D. E., Maestas, J. D., Reeves, M. C., Lankston, R. W., &  \n' +
              'Allred, B. W. (2021). Annual and 16-day rangeland production estimates for the western United States. \n' +
              'Rangeland Ecology & Management, 77, 112-117. \n' +
              '\n' + 
              'Kleinhesselink, A. R., Kachergis, E. J., McCord, S. E., Shirley, J., Hupp, N. R., Walker, J., ... & \n' + 
              'Naugle, D. E. (2023). Long-term trends in vegetation on Bureau of Land Management rangelands in the \n' + 
              'western United States. Rangeland Ecology & Management, 87, 1-12. \n' + 
              '\n' + 
              'ADDITIONAL NOTES: \n' +
              'The official Rangeland Assessment Platform website is found at https://rangelands.app \n' + 
              'The PRISM dataset used is available from https://support.climateengine.org/article/80-prism \n' +
              'The code for this application can be found at www.github.com'
              'This app is coded by Andrew Fullhart and Gerardo Armendariz (U of Arizona SNRE, USDA-ARS-SWRC).';

var text_box = ui.Label({value:info_str, style:infoLabelStyle});
var text_panel = ui.Panel({widgets:null, layout:null, style:textPanelStyle});

function render_infobox(bool_obj){
  if (bool_obj === true){
    text_panel.add(text_box);
    Map.add(text_panel);
  }else{
    Map.remove(text_panel);
    text_panel.clear();
  }
}

var info_checkbox = ui.Checkbox({
  label:'Show ReadMe',
  onChange:render_infobox,
  style:checkStyle
});

main_panel.add(info_checkbox);

ui.root.add(main_panel);

//Transparent panels only works with Map.add?
//Map.add(main_panel);

var ic = ee.ImageCollection("NASA/NEX-DCP30");

var model = ee.String('CCSM4');

var ndays_months = ee.List([31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]);
var order_months = ee.List([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]);
var summr_months = ee.List([4, 5, 6, 7, 8, 9]);
var wintr_months = ee.List([1, 2, 3, 10, 11, 12]);

var start_year = 2000;
var end_year = 2029;

var modelfilter = ee.Filter.or(
  ee.Filter.eq('scenario', 'historical'),
  ee.Filter.eq('scenario', 'rcp45'));
var ic = ic.filter(modelfilter);
var ic = ic.filter(ee.Filter.eq('model', model));

var start = ee.Date.fromYMD(start_year, 1, 1);
var end = ee.Date.fromYMD(end_year+1, 1, 1);
var year_ic = ic.filterDate(start, end);

function make_p_ic_fn(month){
  var month = ee.Number(month);
  var ndays = ee.Number(ndays_months.get(month.subtract(1)));
  var mo_ic = year_ic.filter(ee.Filter.calendarRange(month, month,'month'));
  var p_im = mo_ic.select('pr').reduce(ee.Reducer.mean()).multiply(86400).multiply(ndays);
  return p_im; 
}

var p_ic = ee.ImageCollection(order_months.map(make_p_ic_fn));

function make_t_ic_fn(month){
  var month = ee.Number(month);
  var mo_ic = year_ic.filter(ee.Filter.calendarRange(month, month,'month'));
  var tmax_im = mo_ic.select('tasmax').reduce(ee.Reducer.mean()).subtract(273.15);
  var tmin_im = mo_ic.select('tasmin').reduce(ee.Reducer.mean()).subtract(273.15);
  var t_im = ee.Image(tmax_im.add(tmin_im)).divide(2.0);
  return t_im; 
}

var t_ic = ee.ImageCollection(order_months.map(make_t_ic_fn));

function weight_temps_fn(month){
  var month = ee.Number(month);
  var ndays = ee.Number(ndays_months.get(month.subtract(1)));
  var mo_im = ee.Image(t_ic.toList(12).get(month.subtract(1)));
  var wmo_im = mo_im.multiply(ndays).divide(365.25);
  return wmo_im;
}

var wt_ic = ee.ImageCollection(order_months.map(weight_temps_fn));

var tann_im = wt_ic.reduce(ee.Reducer.sum());
var pann_im = p_ic.reduce(ee.Reducer.sum());
var tw_im = t_ic.reduce(ee.Reducer.max());
var tc_im = t_ic.reduce(ee.Reducer.min());
var zero_im = pann_im.lt(0.0);

function make_p_seasn_fn(month){
  var month = ee.Number(month);
  var mo_im = ee.Image(p_ic.toList(12).get(month.subtract(1)));
  return mo_im; 
}

var pwintr_ic = ee.ImageCollection(wintr_months.map(make_p_seasn_fn));
var psummr_ic = ee.ImageCollection(summr_months.map(make_p_seasn_fn));
var pwintr_im = pwintr_ic.reduce(ee.Reducer.sum());
var psummr_im = psummr_ic.reduce(ee.Reducer.sum());
var pwintrw_im = pwintr_ic.reduce(ee.Reducer.max());
var pwintrc_im = pwintr_ic.reduce(ee.Reducer.min());
var psummrw_im = psummr_ic.reduce(ee.Reducer.max());
var psummrc_im = psummr_ic.reduce(ee.Reducer.min());

var test_im = ee.Image(pann_im.multiply(0.70));
var conA_im = pwintr_im.gte(test_im);
var conB_im = psummr_im.gte(test_im);
var conAB_im = conA_im.add(conB_im);
var conC_im = conAB_im.eq(0.0);

var pthrA_im = conA_im.where(conA_im, tann_im.multiply(2.0));
var pthrB_im = conB_im.where(conB_im, ee.Image(tann_im.multiply(2.0)).add(28.0));
var pthrC_im = conC_im.where(conC_im, ee.Image(tann_im.multiply(2.0)).add(14.0));
var pthr_im = pthrA_im.add(pthrB_im).add(pthrC_im);

//E
var e_im = tw_im.lt(10.0);

var con_et_im = tw_im.gte(0.0);
var mix_im = e_im.add(con_et_im);
var et_im = mix_im.eq(2.0);

var con_ef_im = tw_im.lt(0.0);
var mix_im = e_im.add(con_ef_im);
var ef_im = mix_im.eq(2.0);

//B
var sin_e_im = tw_im.gte(10.0);
var con_b_im = zero_im.where(pann_im.lt(pthr_im.multiply(10.0)), 1);
var mix_im = con_b_im.add(sin_e_im);
var b_im = mix_im.eq(2.0);

var con_bs_im = zero_im.where(pann_im.gt(pthr_im.multiply(5.0)), 1);
var mix_im = b_im.add(con_bs_im);
var bs_im = mix_im.eq(2.0);

var con_bw_im = zero_im.where(pann_im.lte(pthr_im.multiply(5.0)), 1);
var mix_im = b_im.add(con_bw_im);
var bw_im = mix_im.eq(2.0);

var con_bsh_im = zero_im.where(tann_im.gte(18.0), 1);
var mix_im = bs_im.add(con_bsh_im);
var bsh_im = mix_im.eq(2.0);

var con_bsk_im = zero_im.where(tann_im.lt(18.0), 1);
var mix_im = bs_im.add(con_bsk_im);
var bsk_im = mix_im.eq(2.0);

var con_bwh_im = zero_im.where(tann_im.gte(18.0), 1);
var mix_im = bw_im.add(con_bwh_im);
var bwh_im = mix_im.eq(2.0);

var con_bwk_im = zero_im.where(tann_im.lt(18.0), 1);
var mix_im = bw_im.add(con_bwk_im);
var bwk_im = mix_im.eq(2.0);

//D
var mix_im = e_im.add(b_im);
var sin_e_sin_b_im = mix_im.eq(0);
var con_d_im = zero_im.where(tc_im.lte(-3.0), 1);
var mix_im = sin_e_sin_b_im.add(con_d_im);
var d_im = mix_im.eq(2.0);



Map.addLayer(e_im, null, 'e_im');
Map.addLayer(b_im, null, 'b_im');
Map.addLayer(d_im, null, 'd_im');


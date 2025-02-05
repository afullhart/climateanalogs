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

//Binary test images, etc.
var pwintr_ic = ee.ImageCollection(wintr_months.map(make_p_seasn_fn));
var psummr_ic = ee.ImageCollection(summr_months.map(make_p_seasn_fn));
var pwintr_im = pwintr_ic.reduce(ee.Reducer.sum());
var psummr_im = psummr_ic.reduce(ee.Reducer.sum());
var pwintrw_im = pwintr_ic.reduce(ee.Reducer.max());
var pwintrd_im = pwintr_ic.reduce(ee.Reducer.min());
var psummrw_im = psummr_ic.reduce(ee.Reducer.max());
var psummrd_im = psummr_ic.reduce(ee.Reducer.min());

var test_im = ee.Image(pann_im.multiply(0.70));
var conA_im = pwintr_im.gte(test_im);
var conB_im = psummr_im.gte(test_im);
var conAB_im = conA_im.add(conB_im);
var conC_im = conAB_im.eq(0.0);

var pthrA_im = conA_im.where(conA_im, tann_im.multiply(2.0));
var pthrB_im = conB_im.where(conB_im, ee.Image(tann_im.multiply(2.0)).add(28.0));
var pthrC_im = conC_im.where(conC_im, ee.Image(tann_im.multiply(2.0)).add(14.0));
var pthr_im = pthrA_im.add(pthrB_im).add(pthrC_im);

var dry_summrA_im = zero_im.where(psummrd_im.lt(pwintrd_im), 1);
var dry_summrB_im = zero_im.where(pwintrw_im.gt(psummrd_im.multiply(3.0)), 1);
var dry_summrC_im = zero_im.where(psummrd_im.lt(40.0), 1);
var mix_im = dry_summrA_im.add(dry_summrB_im).add(dry_summrC_im);
var dry_summr_im = mix_im.eq(3.0);

var dry_wintrA_im = zero_im.where(pwintrd_im.lt(psummrd_im), 1);
var dry_wintrB_im = zero_im.where(psummrw_im.gt(pwintrd_im.multiply(10.0)), 1);
var mix_im = dry_wintrA_im.add(dry_wintrB_im);
var dry_wintr_im = mix_im.eq(2.0);

var hot_summr_im = zero_im.where(tw_im.gte(22.0), 1);
var sin_hot_summr_im = hot_summr_im.eq(0); 

function count_warm_months_fn(t_im){
  var warm_im = ee.Image(t_im.gte(10));
  return warm_im;
}

var warm_ic = ee.ImageCollection(t_ic.map(count_warm_months_fn));
var warm_mo_ct_im = warm_ic.reduce(ee.Reducer.sum());
var warm_mo_im = warm_mo_ct_im.gte(4);


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
var sin_e_b_im = mix_im.eq(0.0);

var con_d_im = zero_im.where(tc_im.lte(-3.0), 1);
var mix_im = sin_e_b_im.add(con_d_im);
var d_im = mix_im.eq(2.0);

var mix_im = d_im.add(dry_summr_im);
var ds_im = mix_im.eq(2.0);

var mix_im = d_im.add(dry_wintr_im);
var dw_im = mix_im.eq(2.0);

var mix_im = d_im.add(ds_im).add(dw_im);
var df_im = mix_im.eq(1.0);

//Dsa
var con_dsa = zero_im.where(tw_im.gte(22.0), 1);
var mix_im = ds_im.add(con_dsa);
var dsa_im = mix_im.eq(2.0);

//Dsb
var sin_dsa = dsa_im.eq(0.0);
var mix_im = sin_dsa.add(ds_im).add(warm_mo_im);
var dsb_im = mix_im.eq(3.0);

//Dsc
var sin_dsa = dsa_im.eq(0.0);
var sin_dsb = dsb_im.eq(0.0);
var mix_im = sin_dsa.add(sin_dsb).add(ds_im);
var sin_dsa_dsb_im = mix_im.eq(3.0);
var con_dsc_im = zero_im.where(tc_im.gt(-38.0), 1);
var mix_im = con_dsc_im.add(sin_dsa_dsb_im);
var dsc_im = mix_im.eq(2.0);

//Dsd
var sin_dsa = dsa_im.eq(0.0);
var sin_dsb = dsb_im.eq(0.0);
var mix_im = sin_dsa.add(sin_dsb).add(ds_im);
var sin_dsa_dsb_im = mix_im.eq(3.0);
var con_dsd_im = zero_im.where(tc_im.lte(-38.0), 1);
var mix_im = con_dsd_im.add(sin_dsa_dsb_im);
var dsd_im = mix_im.eq(2.0);


//Dwa
var con_dwa = zero_im.where(tw_im.gte(22.0), 1);
var mix_im = dw_im.add(con_dwa);
var dwa_im = mix_im.eq(2.0);

//Dwb
var sin_dwa = dwa_im.eq(0.0);
var mix_im = sin_dwa.add(dw_im).add(warm_mo_im);
var dwb_im = mix_im.eq(3.0);

//Dwc
var sin_dwa = dwa_im.eq(0.0);
var sin_dwb = dwb_im.eq(0.0);
var mix_im = sin_dwa.add(sin_dwb).add(dw_im);
var sin_dwa_dwb_im = mix_im.eq(3.0);
var con_dwc_im = zero_im.where(tc_im.gt(-38.0), 1);
var mix_im = con_dwc_im.add(sin_dwa_dwb_im);
var dwc_im = mix_im.eq(2.0);

//Dwd
var sin_dwa = dwa_im.eq(0.0);
var sin_dwb = dwb_im.eq(0.0);
var mix_im = sin_dwa.add(sin_dwb).add(dw_im);
var sin_dwa_dwb_im = mix_im.eq(3.0);
var con_dwd_im = zero_im.where(tc_im.lte(-38.0), 1);
var mix_im = con_dwd_im.add(sin_dwa_dwb_im);
var dwd_im = mix_im.eq(2.0);


//Dfa
var con_dfa = zero_im.where(tw_im.gte(22.0), 1);
var mix_im = df_im.add(con_dfa);
var dfa_im = mix_im.eq(2.0);

//Dfb
var sin_dfa = dfa_im.eq(0.0);
var mix_im = sin_dfa.add(df_im).add(warm_mo_im);
var dfb_im = mix_im.eq(3.0);

//Dfc
var sin_dfa = dfa_im.eq(0.0);
var sin_dfb = dfb_im.eq(0.0);
var mix_im = sin_dfa.add(sin_dfb).add(df_im);
var sin_dfa_dfb_im = mix_im.eq(3.0);
var con_dfc_im = zero_im.where(tc_im.gt(-38.0), 1);
var mix_im = con_dfc_im.add(sin_dfa_dfb_im);
var dfc_im = mix_im.eq(2.0);

//Dfd
var sin_dfa = dfa_im.eq(0.0);
var sin_dfb = dfb_im.eq(0.0);
var mix_im = sin_dfa.add(sin_dfb).add(df_im);
var sin_dfa_dfb_im = mix_im.eq(3.0);
var con_dfd_im = zero_im.where(tc_im.lte(-38.0), 1);
var mix_im = con_dfd_im.add(sin_dfa_dfb_im);
var dfd_im = mix_im.eq(2.0);



//C
var mix_im = e_im.add(b_im).add(d_im);
var sin_e_b_d_im = mix_im.eq(0.0);

var con_c_im = zero_im.where(tc_im.lt(18.0), 1);
var mix_im = sin_e_b_d_im.add(con_c_im);
var c_im = mix_im.eq(2.0);

var mix_im = c_im.add(dry_summr_im);
var cs_im = mix_im.eq(2.0);

var mix_im = c_im.add(dry_wintr_im);
var cw_im = mix_im.eq(2.0);

var mix_im = c_im.add(cs_im).add(cw_im);
var cf_im = mix_im.eq(1.0);


//Csa
var con_csa = zero_im.where(tw_im.gte(22.0), 1);
var mix_im = cs_im.add(con_csa);
var csa_im = mix_im.eq(2.0);

//Csb
var sin_csa = csa_im.eq(0.0);
var mix_im = sin_csa.add(cs_im).add(warm_mo_im);
var csb_im = mix_im.eq(3.0);

//Csc
var sin_csa = csa_im.eq(0.0);
var sin_csb = csb_im.eq(0.0);
var mix_im = sin_csa.add(sin_csb).add(cs_im);
var sin_csa_csb_im = mix_im.eq(3.0);
var con_csc_im = zero_im.where(tc_im.gt(-38.0), 1);
var mix_im = con_dsc_im.add(sin_csa_csb_im);
var csc_im = mix_im.eq(2.0);

//Csd
var sin_csa = csa_im.eq(0.0);
var sin_csb = csb_im.eq(0.0);
var mix_im = sin_csa.add(sin_csb).add(cs_im);
var sin_csa_csb_im = mix_im.eq(3.0);
var con_csd_im = zero_im.where(tc_im.lte(-38.0), 1);
var mix_im = con_csd_im.add(sin_csa_csb_im);
var csd_im = mix_im.eq(2.0);


//Cwa
var con_cwa = zero_im.where(tw_im.gte(22.0), 1);
var mix_im = cw_im.add(con_cwa);
var cwa_im = mix_im.eq(2.0);

//Cwb
var sin_cwa = cwa_im.eq(0.0);
var mix_im = sin_cwa.add(cw_im).add(warm_mo_im);
var cwb_im = mix_im.eq(3.0);

//Cwc
var sin_cwa = cwa_im.eq(0.0);
var sin_cwb = cwb_im.eq(0.0);
var mix_im = sin_cwa.add(sin_cwb).add(cw_im);
var sin_cwa_cwb_im = mix_im.eq(3.0);
var con_cwc_im = zero_im.where(tc_im.gt(-38.0), 1);
var mix_im = con_cwc_im.add(sin_cwa_cwb_im);
var cwc_im = mix_im.eq(2.0);

//Cwd
var sin_cwa = cwa_im.eq(0.0);
var sin_cwb = cwb_im.eq(0.0);
var mix_im = sin_cwa.add(sin_cwb).add(cw_im);
var sin_cwa_cwb_im = mix_im.eq(3.0);
var con_cwd_im = zero_im.where(tc_im.lte(-38.0), 1);
var mix_im = con_cwd_im.add(sin_cwa_cwb_im);
var cwd_im = mix_im.eq(2.0);


//Cfa
var con_dfa = zero_im.where(tw_im.gte(22.0), 1);
var mix_im = df_im.add(con_dfa);
var dfa_im = mix_im.eq(2.0);

//Cfb
var sin_dfa = dfa_im.eq(0.0);
var mix_im = sin_dfa.add(df_im).add(warm_mo_im);
var dfb_im = mix_im.eq(3.0);

//Cfc
var sin_dfa = dfa_im.eq(0.0);
var sin_dfb = dfb_im.eq(0.0);
var mix_im = sin_dfa.add(sin_dfb).add(df_im);
var sin_dfa_dfb_im = mix_im.eq(3.0);
var con_dfc_im = zero_im.where(tc_im.gt(-38.0), 1);
var mix_im = con_dfc_im.add(sin_dfa_dfb_im);
var dfc_im = mix_im.eq(2.0);

//Cfd
var sin_dfa = dfa_im.eq(0.0);
var sin_dfb = dfb_im.eq(0.0);
var mix_im = sin_dfa.add(sin_dfb).add(df_im);
var sin_dfa_dfb_im = mix_im.eq(3.0);
var con_dfd_im = zero_im.where(tc_im.lte(-38.0), 1);
var mix_im = con_dfd_im.add(sin_dfa_dfb_im);
var dfd_im = mix_im.eq(2.0);









Map.addLayer(ds_im, null, 'ds_im');
Map.addLayer(dsa_im, null, 'dsa_im');
Map.addLayer(dsb_im, null, 'dsb_im');
Map.addLayer(dsc_im, null, 'dsc_im');
Map.addLayer(dsd_im, null, 'dsd_im');

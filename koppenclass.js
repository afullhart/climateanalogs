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
  var p_im = mo_ic.select('pr').reduce(ee.Reducer.mean()).multiply(86400.0).multiply(ndays);
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
var pd_im = p_ic.reduce(ee.Reducer.min());
var zero_im = pann_im.lt(0.0);

//Binary test images, etc.
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
var pwintrd_im = pwintr_ic.reduce(ee.Reducer.min());
var psummrw_im = psummr_ic.reduce(ee.Reducer.max());
var psummrd_im = psummr_ic.reduce(ee.Reducer.min());
var pd_in_summr_im = psummrd_im.lt(pwintrd_im);
var pd_in_wintr_im = pwintrd_im.lt(psummrd_im);

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
  var warm_im = ee.Image(t_im.gte(10.0));
  return warm_im;
}

var warm_ic = ee.ImageCollection(t_ic.map(count_warm_months_fn));
var warm_mo_ct_im = warm_ic.reduce(ee.Reducer.sum());
var warm_mo_im = warm_mo_ct_im.gte(4);



//E
var e_im = tw_im.lt(10.0);
 
//Et
var con_et_im = tw_im.gte(0.0);
var mix_im = e_im.add(con_et_im);
var et_im = mix_im.eq(2.0);

//Ef
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

//Bsh
var con_bsh_im = zero_im.where(tann_im.gte(18.0), 1);
var mix_im = bs_im.add(con_bsh_im);
var bsh_im = mix_im.eq(2.0);

//Bsk
var con_bsk_im = zero_im.where(tann_im.lt(18.0), 1);
var mix_im = bs_im.add(con_bsk_im);
var bsk_im = mix_im.eq(2.0);

//Bwh
var con_bwh_im = zero_im.where(tann_im.gte(18.0), 1);
var mix_im = bw_im.add(con_bwh_im);
var bwh_im = mix_im.eq(2.0);

//Bwk
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
var con_cfa = zero_im.where(tw_im.gte(22.0), 1);
var mix_im = cf_im.add(con_cfa);
var cfa_im = mix_im.eq(2.0);

//Cfb
var sin_cfa = cfa_im.eq(0.0);
var mix_im = sin_cfa.add(cf_im).add(warm_mo_im);
var cfb_im = mix_im.eq(3.0);

//Cfc
var sin_cfa = cfa_im.eq(0.0);
var sin_cfb = cfb_im.eq(0.0);
var mix_im = sin_cfa.add(sin_cfb).add(cf_im);
var sin_cfa_cfb_im = mix_im.eq(3.0);
var con_cfc_im = zero_im.where(tc_im.gt(-38.0), 1);
var mix_im = con_cfc_im.add(sin_cfa_cfb_im);
var cfc_im = mix_im.eq(2.0);

//Cfd
var sin_cfa = cfa_im.eq(0.0);
var sin_cfb = cfb_im.eq(0.0);
var mix_im = sin_cfa.add(sin_cfb).add(cf_im);
var sin_cfa_cfb_im = mix_im.eq(3.0);
var con_cfd_im = zero_im.where(tc_im.lte(-38.0), 1);
var mix_im = con_cfd_im.add(sin_cfa_cfb_im);
var cfd_im = mix_im.eq(2.0);

//A
var sin_b_im = b_im.eq(0.0);
var con_a_im = zero_im.where(tc_im.gte(18.0), 1);
var mix_im = con_a_im.add(sin_b_im);
var a_im = mix_im.eq(2.0);

//Am
var con_am_im = zero_im.where(pann_im.gte(ee.Image(ee.Image(pd_im.multiply(-1.0)).add(100.0)).multiply(25.0)), 1);
var mix_im = con_am_im.add(a_im);
var am_im = mix_im.eq(2.0);

//Af
var sin_am_im = con_am_im.eq(0.0);
var con_af_im = zero_im.where(pd_im.gte(60.0), 1);
var mix_im = con_af_im.add(sin_am_im).add(a_im);
var af_im = mix_im.eq(3.0);

//As
var con_as_im = zero_im.where(pd_im.lt(60.0), 1);
var mix_im = con_as_im.add(sin_am_im).add(a_im).add(pd_in_summr_im);
var as_im = mix_im.eq(4.0);

//Aw
var con_aw_im = zero_im.where(pd_im.lt(60.0), 1);
var mix_im = con_aw_im.add(sin_am_im).add(a_im).add(pd_in_wintr_im);
var aw_im = mix_im.eq(4.0);



//Type value assignments
var af_im = af_im.where(af_im.eq(1.0), 1);
var am_im = am_im.where(am_im.eq(1.0), 2);
//As not present?
var aw_im = aw_im.where(aw_im.eq(1.0), 3);

var bwh_im = bwh_im.where(bwh_im.eq(1.0), 4);
var bwk_im = bwk_im.where(bwk_im.eq(1.0), 5);
var bsh_im = bsh_im.where(bsh_im.eq(1.0), 6);
var bsk_im = bsk_im.where(bsk_im.eq(1.0), 7);

var csa_im = csa_im.where(csa_im.eq(1.0), 8);
var csb_im = csb_im.where(csb_im.eq(1.0), 9);
var csc_im = csc_im.where(csc_im.eq(1.0), 10);
//csd not present?
var cwa_im = cwa_im.where(cwa_im.eq(1.0), 11);
var cwb_im = cwb_im.where(cwb_im.eq(1.0), 12);
var cwc_im = cwc_im.where(cwc_im.eq(1.0), 13);
//cwd not present?
var cfa_im = cfa_im.where(cfa_im.eq(1.0), 14);
var cfb_im = cfb_im.where(cfb_im.eq(1.0), 15);
var cfc_im = cfc_im.where(cfc_im.eq(1.0), 16);
//cfd not present?

var dsa_im = dsa_im.where(dsa_im.eq(1.0), 17);
var dsb_im = dsb_im.where(dsb_im.eq(1.0), 18);
var dsc_im = dsc_im.where(dsc_im.eq(1.0), 19);
var dsd_im = dsd_im.where(dsd_im.eq(1.0), 20);
var dwa_im = dwa_im.where(dwa_im.eq(1.0), 21);
var dwb_im = dwb_im.where(dwb_im.eq(1.0), 22);
var dwc_im = dwc_im.where(dwc_im.eq(1.0), 23);
var dwd_im = dwd_im.where(dwd_im.eq(1.0), 24);
var dfa_im = dfa_im.where(dfa_im.eq(1.0), 25);
var dfb_im = dfb_im.where(dfb_im.eq(1.0), 26);
var dfc_im = dfc_im.where(dfc_im.eq(1.0), 27);
var dfd_im = dfd_im.where(dfd_im.eq(1.0), 28);

var et_im = et_im.where(et_im.eq(1.0), 29);
var ef_im = ef_im.where(ef_im.eq(1.0), 30);

var type_ic = ee.ImageCollection([af_im, am_im, aw_im, bwh_im, bwk_im, bsh_im, bsk_im, csa_im, csb_im, csc_im, cwa_im, cwb_im, cwc_im, cfa_im, cfb_im, cfc_im, dsa_im, dsb_im, dsc_im, dsd_im, dwa_im, dwb_im, dwc_im, dwd_im, dfa_im, dfb_im, dfc_im, dfd_im, et_im, ef_im]);

function change_band_name_fn(im){
  var bLabel = im.bandNames().get(0);
  return im.select([bLabel],['B1']);
}

var type_ic = ee.ImageCollection(type_ic.map(change_band_name_fn));
var type_ic = ee.ImageCollection(type_ic.cast({B1:'int64'}, ['B1']));
var type_im = type_ic.reduce(ee.Reducer.sum());

var typePalette = [
  '#0000FF', // Af
  '#0078FF', //Am
  '#46FAAA', //Aw
  '#FF0000', //Bwh
  '#FF9696', //Bwk
  '#F5A500', //Bsh
  '#FFDC64', //Bsk
  '#FFFF00', //Csa
  '#C8C800', //Csb
  '#969600', //Csc
  '#96FF96', //Cwa
  '#64C864', //Cwb
  '#329632', //Cwc
  '#C8FF50', //Cfa
  '#64FF50', //Cfb
  '#32C800', //Cfc
  '#FF00FF', //Dsa
  '#C800C8', //Dsb
  '#963296', //Dsc
  '#966496', //Dsd
  '#AAAFFF', //Dwa
  '#5A78DC', //Dwb
  '#4B50B4', //Dwc
  '#320087', //Dwd
  '#00FFFF', //Dfa
  '#37C8FF', //Dfb
  '#007D7D', //Dfc
  '#00465F', //Dfd
  '#B2B2B2', //Et
  '#666666' //Ef
];

var singleBandVis = {
  min: 1,
  max: 30,
  palette: typePalette
};

Map.addLayer(type_im, singleBandVis);



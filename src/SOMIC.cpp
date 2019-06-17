// SOMic Model ver 1.00
// Author Dominic Woolf
// d.woolf@cornell.edu

#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

double calc_atsmd (double prec, double pet, double max_tsmd, int cover, double init_tsmd){
  // Calculate acquired topsoil moisture deficit, based on pet (potential evapotranspiration) and precip (precipitation)
  // Only used when SOMic is called with use_atsmd = TRUE
  // otherwise, SOMic uses h2O and sat vectors (soil water content and saturation capacity)
  double atsmd = 0.0;
  double excess = prec - pet;
  if (cover) atsmd = std::max(init_tsmd + excess, max_tsmd);
  else if (init_tsmd < max_tsmd/1.8) atsmd = init_tsmd + std::max(excess, 0.0);
  return(std::min(atsmd, 0.0));
}

double calc_clayfact (double clay, double mclay) {
  return(mclay*clay + 1 - 23*mclay);
}

double calc_cue (double temp, double cue_0, double mcue) {
  //http://onlinelibrary.wiley.com.proxy.library.cornell.edu/doi/10.1890/15-2110.1/full
  //http://onlinelibrary.wiley.com/doi/10.1111/ele.12113/full
  //http://link.springer.com/article/10.1007/s10533-013-9948-8
  return(cue_0 - (temp-15)*mcue);
}

double doc_leached (double doc, double thickness, double velocity, double kdoc) {
  double leached;
  // for now use constant velocity
  // at later time we may want to adjust leaching rate with water content
  leached = doc * (velocity / thickness) * std::exp(-kdoc); 
  leached = std::min(std::max(0.0, leached), doc); // we can't leach more doc than is remaining or less than none (upward capillary flow not yet included) 
  return (leached);
}

double wmean (NumericVector vals, NumericVector wts){
  // Weighted mean
  double sv=0.0;
  double sw=0.0;
  unsigned int n = vals.size();
  for (unsigned int i = 0; i<n; i++){
    sv = sv + vals[i]*wts[i];
    sw = sw + wts[i];
  }
  return(sv/sw);
}

// [[Rcpp::export]]
List somic(
    DataFrame soc_data,
    double mic_vmax,
    double mic_km,
    double kdissolution,
    double kdepoly,
    double kdeath_and_exudates,
    double kdesorb,
    double ksorb,
    double kmicrobial_uptake,
    double cue_0,
    double mcue,
    double mclay,
    double clay,
    double thickness = 23.0,
    double sat = 1.0,
    double max_tsmd = 0.0,
    bool use_atsmd = 0,
    double init_spm_14c_age = 0.0,
    double init_ipm_14c_age = 0.0,
    double init_doc_14c_age = 0.0,
    double init_mb_14c_age = 0.0,
    double init_mac_14c_age = 0.0,
    double init_soc_14c_age = NA_REAL) {
  
  // extract input dataframe columns into C++ vectors
  IntegerVector time = soc_data["time"];
  NumericVector spm = soc_data["spm"];
  NumericVector ipm = soc_data["ipm"];
  NumericVector doc = soc_data["doc"];
  NumericVector mb = soc_data["mb"];
  NumericVector mac = soc_data["mac"];
  NumericVector soc = soc_data["soc"];
  NumericVector add_spm = soc_data["added.spm"];
  NumericVector add_ipm = soc_data["added.ipm"];
  NumericVector add_doc = soc_data["added.doc"];
  NumericVector add_mb = soc_data["added.mb"];
  NumericVector add_mac = soc_data["added.mac"];
  IntegerVector cover = soc_data["cover"];
  NumericVector atsmd = soc_data["atsmd"];
  NumericVector temp = soc_data["temp"];
  NumericVector precip = soc_data["precip"];
  NumericVector pet = soc_data["pet"];
  NumericVector h2o = soc_data["h2o"];
  NumericVector a = soc_data["a"];
  NumericVector c = soc_data["c"];
  NumericVector add_d13c = soc_data["added.d13c"];
  NumericVector soc_d13c = soc_data["soc.d13c"];
  NumericVector add_14c_age = soc_data["add_14c_age"];
  NumericVector velocity = soc_data["velocity"];
  
  // define new vectors
  unsigned int n = time.size();
  NumericVector dec_spm(n), dec_ipm(n), dec_doc(n), dec_mb(n), dec_mac(n), dec_tot(n);
  NumericVector sorption(n), microbial_uptake(n), growth(n);
  NumericVector min_doc(n), min_cum(n);
  NumericVector b(n), mic(n), cue(n);
  NumericVector spm_d13c(n), ipm_d13c(n), doc_d13c(n), mb_d13c(n), mac_d13c(n), co2_d13c(n);
  NumericVector spm_14c(n), ipm_14c(n), doc_14c(n), mb_14c(n), mac_14c(n), soc_14c(n);
  NumericVector leached_doc(n), leached_14C(n);
  // New internal variables
  double clayfact = calc_clayfact(clay, mclay);
  double ksorb_altered, kmicrobial_uptake_altered, fsorb, fmic, kdoc;
  
  //====================================================================================================================================
  // Initialisations
  // ---------------
  // initialise soc pools d13C to bulk soc d13C
  spm_d13c[0] = ipm_d13c[0] = doc_d13c[0] = mb_d13c[0] = mac_d13c[0] = co2_d13c[0] = soc_d13c[0]; 
  // initialise 14C ages 
  if (NumericVector::is_na(init_soc_14c_age)) { // if bulk soil 14C not provided, then calculate from 14C of pools (weighted mean)
    init_soc_14c_age = (init_spm_14c_age * spm[0] + init_ipm_14c_age * ipm[0] + init_doc_14c_age * doc[0] + 
      init_mb_14c_age * mb[0] + init_mac_14c_age * mac[0]) / soc[0];
  }
  // initialise 14C ages of pools equal to bulk soil
  spm_14c[0] = ipm_14c[0] = doc_14c[0] = mb_14c[0] = mac_14c[0] = soc_14c[0] = leached_14C[0] = init_soc_14c_age; 
  
  //====================================================================================================================================
  // Calculate pool sizes and isoptope values
  // ----------------------------------------	
  for (unsigned int i = 0; i < n; i++) {
    
    if (i > 0) {  // first month has pool sizes and delta 13C values already initialized. Calculate only for subsequent months
      // pools
      spm[i] = spm[i-1] - dec_spm[i-1] + add_spm[i-1];
      ipm[i] = ipm[i-1] - dec_ipm[i-1] + add_ipm[i-1];
      doc[i] = doc[i-1] - dec_doc[i-1] + add_doc[i-1] + dec_spm[i-1] + dec_ipm[i-1] + dec_mb[i-1] + dec_mac[i-1];
      leached_doc[i] = doc_leached(doc[i], thickness, velocity[i], kdoc);
      doc[i] = doc[i] - leached_doc[i];
      mb[i] = mb[i-1] - dec_mb[i-1] + add_mb[i-1] + growth[i-1];
      mac[i] = mac[i-1] - dec_mac[i-1] + add_mac[i-1] + sorption[i-1];
      soc[i] = spm[i] + ipm[i] + doc[i] + mb[i] + mac[i];
      
      //d13C of soil carbon pools
      spm_d13c[i] = (spm[i] > 0) ? ((spm[i-1] - dec_spm[i-1]) * spm_d13c[i-1] + add_spm[i-1] *add_d13c[i-1]) / spm[i] : spm_d13c[i-1]; // if there's no spm remaining, just use previous value for spm 13C
      ipm_d13c[i] = (ipm[i] > 0) ? ((ipm[i-1] - dec_ipm[i-1]) * ipm_d13c[i-1] + add_ipm[i-1] *add_d13c[i-1]) / ipm[i] : ipm_d13c[i-1];
      doc_d13c[i] = ((doc[i-1] - dec_doc[i-1]) * doc_d13c[i-1] + add_doc[i-1] *add_d13c[i-1] + dec_spm[i-1] * spm_d13c[i-1] + dec_ipm[i-1] * ipm_d13c[i-1] + dec_mb[i-1] * mb_d13c[i-1] + dec_mac[i-1] * mac_d13c[i-1]) / doc[i];
      mb_d13c[i] = ((mb[i-1] - dec_mb[i-1]) * mb_d13c[i-1] + growth[i-1]  * doc_d13c[i-1]) / mb[i];
      mac_d13c[i] = ((mac[i-1] - dec_mac[i-1]) * mac_d13c[i-1] + add_mac[i-1] * add_d13c[i-1] + sorption[i-1] * doc_d13c[i-1]) / mac[i];
      
      //14C age of soil carbon pools
      spm_14c[i] = (spm[i] > 0) ? ((spm[i-1] - dec_spm[i-1]) * spm_14c[i-1] + add_spm[i-1] * add_14c_age[i-1]) / spm[i] + 1 : spm_14c[i-1]; // if there's no spm remaining, just use previous value for spm 14C
      ipm_14c[i] = (ipm[i] > 0) ? ((ipm[i-1] - dec_ipm[i-1]) * ipm_14c[i-1] + add_ipm[i-1] * add_14c_age[i-1]) / ipm[i] + 1 : ipm_14c[i-1];
      doc_14c[i] = wmean(NumericVector::create(doc_14c[i-1],            add_14c_age[i-1], spm_14c[i-1], ipm_14c[i-1], mb_14c[i-1], mac_14c[i-1]),
                         NumericVector::create(doc[i-1] - dec_doc[i-1], add_doc[i-1],     dec_spm[i-1], dec_ipm[i-1], dec_mb[i-1], dec_mac[i-1]));
      doc_14c[i] = doc_14c[i] + 1.0;
      mb_14c[i] = ((mb[i-1] - dec_mb[i-1]) * mb_14c[i-1] + growth[i-1]  * doc_14c[i-1]) / mb[i] + 1;
      mac_14c[i] = ((mac[i-1] - dec_mac[i-1]) * mac_14c[i-1] + add_mac[i-1] * add_14c_age[i-1] + sorption[i-1] * doc_14c[i-1]) / mac[i] + 1;
    }
    //====================================================================================================================================
    // Calculate rate modifiers
    // ------------------------	
    // rate modifiers for moisture and microbial biomass
    if (use_atsmd) {
      // note topsoil moisture deficit (TSMD) is measured in mm
      atsmd[i] = (i == 0) ? 0 : calc_atsmd(precip[i], pet[i], max_tsmd, cover[i], atsmd[i-1]);
      b[i] = (atsmd[i] > 0.444*max_tsmd) ? 1 : (0.2 + 0.8*(max_tsmd-atsmd[i]) / (0.556*max_tsmd));
    }
    else {
      // h2o is dimensionless, m/m water content
      b[i] =  (h2o[i] > 0.556*sat) ? 1 : 0.2 + 0.8 * h2o[i] / (0.556*sat);
    }
    mic[i] = mic_vmax*mb[i] / (mic_km + mb[i]); // Use Michaelis-Menten dynamics to relate MB to rate
    
    //====================================================================================================================================
    // Calculate partition fractions
    // -----------------------------	
    // partitioning of doc between competing processes of sorption and microbial uptake
    ksorb_altered = ksorb * clayfact * a[i] * b[i] * c[i]; // sorption rate multiplied by rate modifiers
    kmicrobial_uptake_altered = kmicrobial_uptake * mic[i] * a[i] * b[i] * c[i]; // microbial uptake multiplied by modifiers
    fsorb = ksorb_altered / (ksorb_altered + kmicrobial_uptake_altered); // fraction of doc turnover sorbed
    fmic = 1 - fsorb; // fraction of doc turnover taken up by microbes
    kdoc = (fsorb * ksorb_altered) + (fmic * kmicrobial_uptake_altered); // rate constant for doc loss is weighted mean of sorption and microbial uptake
    dec_doc[i] = doc[i] * (1 - std::exp(-kdoc));
    sorption[i] = dec_doc[i] * fsorb;
    microbial_uptake[i] = dec_doc[i] * fmic;
    // partitioning of microbial uptake between growth and respiration
    cue[i] = calc_cue(temp[i], cue_0, mcue);
    growth[i] = cue[i] * microbial_uptake[i];
    min_doc[i] = (1- cue[i]) * microbial_uptake[i];               // Mineralization in time period
    min_cum[i] = i==0 ? min_doc[i] : min_cum[i-1] + min_doc[i];   // cumulative mineralization
    
    //====================================================================================================================================
    // Calculate decomposition rates
    // -----------------------------	
    // Decomposition of non-doc pools
    dec_spm[i] = spm[i] * (1 - std::exp(-kdissolution * a[i] * b[i] * c[i] * mic[i]));
    dec_ipm[i] = ipm[i] * (1 - std::exp(-kdepoly * a[i] * b[i] * c[i] * mic[i]));
    dec_mb[i] = mb[i] * (1 - std::exp(-kdeath_and_exudates * a[i] * b[i] * c[i] * mic[i]));
    dec_mac[i] = mac[i] * (1 - std::exp(-kdesorb * a[i] * b[i] * c[i] * mic[i]));
    dec_tot[i] = dec_spm[i] + dec_ipm[i] + dec_doc[i] + dec_mb[i] + dec_mac[i];
    
    //====================================================================================================================================
    // Isotopes of respired CO2 and total soil organic carbon
    // ------------------------------------------------------
    co2_d13c[i] = doc_d13c[i];
    soc_d13c[i] = (spm[i]*spm_d13c[i] + ipm[i]*ipm_d13c[i] + doc[i]*doc_d13c[i] + mb[i]*mb_d13c[i] + mac[i]*mac_d13c[i]) / soc[i];
    soc_14c[i] = (spm[i]*spm_14c[i] + ipm[i]*ipm_14c[i] + doc[i]*doc_14c[i] + mb[i]*mb_14c[i] + mac[i]*mac_14c[i]) / soc[i];
    leached_14C[i] = doc_14c[i];
  }
  
  //export values back to R as a list
  List output;
  output["time"] = time;
  output["spm"] = spm;
  output["ipm"] = ipm;
  output["doc"] = doc;
  output["mb"] = mb;
  output["mac"] = mac;
  output["soc"] = soc;
  output["atsmd"] = atsmd;
  output["cue"] = cue;
  output["b"] = b;
  output["mic"] = mic;
  output["sorption"] = sorption;
  output["microbial_uptake"] = microbial_uptake;
  output["growth"] = growth;
  output["decomp.spm"] = dec_spm;
  output["decomp.ipm"] = dec_ipm;
  output["decomp.doc"] = dec_doc;
  output["decomp.mb"] = dec_mb;
  output["decomp.mac"] = dec_mac;
  output["decomp.tot"] = dec_tot;
  output["mineralise.doc"] = min_doc;
  output["mineralise.cum"] = min_cum;
  output["spm.d13c"] = spm_d13c;
  output["ipm.d13c"] = ipm_d13c;
  output["doc.d13c"] = doc_d13c;
  output["mb.d13c"] = mb_d13c;
  output["mac.d13c"] = mac_d13c;
  output["soc.d13c"] = soc_d13c;
  output["co2.d13c"] = co2_d13c;
  output["spm.14c"] = spm_14c;
  output["ipm.14c"] = ipm_14c;
  output["doc.14c"] = doc_14c;
  output["mb.14c"] = mb_14c;
  output["mac.14c"] = mac_14c;
  output["soc.14c"] = soc_14c;
  output["leached_doc"] = leached_doc;
  output["leached_14C"] = leached_14C;
  return output;
}

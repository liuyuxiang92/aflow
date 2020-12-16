// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow COREY OSES - Duke University 2013-2020                  *
// *                                                                         *
// ***************************************************************************
// Written by Corey Oses
// corey.oses@duke.edu

#ifndef _AFLOW_ML_CPP_
#define _AFLOW_ML_CPP_

#include "aflow.h"
#include "aflow_ml.h"
#include "aflow_cce.h"

#define _DEBUG_ML_ false

#define _ENERGIES_IONIZATION_MAX_ 5
#define _INJECT_ELEMENTAL_COMBINATIONS_ true
#define _TM_ONLY_ true

namespace aflowML {
  double getLIB0EnergyPAWPBE(const string& species_pp){
    if(species_pp=="Ac") return -0.330411;
    if(species_pp=="Ac_s") return -0.985025;
    if(species_pp=="Ag") return -0.337763;
    if(species_pp=="Al") return -0.313186;
    if(species_pp=="Al_h") return -0.645972;
    if(species_pp=="Ar") return -0.051018;
    if(species_pp=="As") return -1.6833;
    if(species_pp=="Au") return -0.286055;
    if(species_pp=="B_h") return -0.444479;
    if(species_pp=="B_s") return -0.485521;
    if(species_pp=="Ba_sv") return -0.034562;
    if(species_pp=="Be") return -0.025499;
    if(species_pp=="Be_sv") return -0.017215;
    if(species_pp=="Bi") return -1.32564;
    if(species_pp=="Bi_d") return -1.51376;
    if(species_pp=="Br") return -0.261686;
    if(species_pp=="C") return -1.37179;
    if(species_pp=="C_h") return -1.3699;
    if(species_pp=="C_s") return -1.35209;
    if(species_pp=="Ca_pv") return -0.061212;
    if(species_pp=="Ca_sv") return -0.093014;
    if(species_pp=="Cd") return -0.166094;
    if(species_pp=="Cl") return -0.373905;
    if(species_pp=="Cl_h") return -0.35806;
    if(species_pp=="Co") return -1.99048;
    if(species_pp=="Cr") return -5.48369;
    if(species_pp=="Cr_pv") return -5.59721;
    if(species_pp=="Cs_sv") return -0.137636;
    if(species_pp=="Cu") return -0.241525;
    if(species_pp=="Cu_pv") return -0.599875;
    if(species_pp=="Dy_3") return -0.387234;
    if(species_pp=="Er_3") return -0.390044;
    if(species_pp=="Eu") return -8.36334;
    if(species_pp=="Eu_2") return -0.045949;
    if(species_pp=="F") return -0.704599;
    if(species_pp=="F_h") return -0.820093;
    if(species_pp=="F_s") return -0.663234;
    if(species_pp=="Fe") return -3.45752;
    if(species_pp=="Fe_pv") return -3.46819;
    if(species_pp=="Ga") return -0.277948;
    if(species_pp=="Ga_d") return -0.397361;
    if(species_pp=="Ga_h") return -0.282091;
    if(species_pp=="Ge") return -0.771177;
    if(species_pp=="Ge_d") return -0.885122;
    if(species_pp=="Ge_h") return -0.78229;
    if(species_pp=="H") return -1.11416;
    if(species_pp=="H_h") return -1.11179;
    if(species_pp=="He") return -0.000538;
    if(species_pp=="Hf") return -3.47402;
    if(species_pp=="Hf_pv") return -3.53246;
    if(species_pp=="Hg") return -0.122356;
    if(species_pp=="Ho_3") return -0.381068;
    if(species_pp=="I") return -0.212174;
    if(species_pp=="In") return -0.22805;
    if(species_pp=="In_d") return -0.407416;
    if(species_pp=="Ir") return -1.62733;
    if(species_pp=="K_pv") return -0.158564;
    if(species_pp=="K_sv") return -0.227694;
    if(species_pp=="La") return -0.5631;
    if(species_pp=="La_s") return -0.516331;
    if(species_pp=="Li") return -0.293945;
    if(species_pp=="Li_sv") return -0.297272;
    if(species_pp=="Lu") return -0.479887;
    if(species_pp=="Lu_3") return -0.382234;
    if(species_pp=="Mg") return -0.038135;
    if(species_pp=="Mg_pv") return -0.094786;
    if(species_pp=="Mn") return -5.15867;
    if(species_pp=="Mn_pv") return -5.32776;
    if(species_pp=="Mo") return -4.59003;
    if(species_pp=="Mo_pv") return -4.60432;
    if(species_pp=="N") return -3.12444;
    if(species_pp=="N_h") return -3.11812;
    if(species_pp=="N_s") return -3.13766;
    if(species_pp=="Na") return -0.21998;
    if(species_pp=="Na_pv") return -0.228379;
    if(species_pp=="Na_sv") return -0.235575;
    if(species_pp=="Nb_sv") return -3.21181;
    if(species_pp=="Nd") return -4.15081;
    if(species_pp=="Ne") return -0.012381;
    if(species_pp=="O") return -1.90885;
    if(species_pp=="O_h") return -1.96399;
    if(species_pp=="O_s") return -1.89585;
    if(species_pp=="Os") return -2.90915;
    if(species_pp=="Os_pv") return -2.92055;
    if(species_pp=="P") return -1.88731;
    if(species_pp=="P_h") return -1.88281;
    if(species_pp=="Pb") return -0.58418;
    if(species_pp=="Pb_d") return -0.765497;
    if(species_pp=="Pd") return -1.47315;
    if(species_pp=="Pd_pv") return -1.74771;
    if(species_pp=="Pm_3") return -0.499376;
    if(species_pp=="Pr") return -2.48756;
    if(species_pp=="Pr_3") return -0.558962;
    if(species_pp=="Pt") return -0.604572;
    if(species_pp=="Rb_pv") return -0.156737;
    if(species_pp=="Rb_sv") return -0.188051;
    if(species_pp=="Re") return -4.60422;
    if(species_pp=="Re_pv") return -4.69957;
    if(species_pp=="Rh") return -1.64118;
    if(species_pp=="Rh_pv") return -1.7322;
    if(species_pp=="Ru") return -2.52896;
    if(species_pp=="Ru_pv") return -2.57847;
    if(species_pp=="S") return -1.07971;
    if(species_pp=="S_h") return -1.06433;
    if(species_pp=="Sb") return -1.41051;
    if(species_pp=="Sc_sv") return -2.18958;
    if(species_pp=="Se") return -0.878847;
    if(species_pp=="Si") return -0.871529;
    if(species_pp=="Si_h") return -0.878813;
    if(species_pp=="Sm_3") return -0.478431;
    if(species_pp=="Sn") return -0.643434;
    if(species_pp=="Sn_d") return -0.832202;
    if(species_pp=="Sr_sv") return -0.068251;
    if(species_pp=="Ta") return -3.66606;
    if(species_pp=="Ta_pv") return -3.74451;
    if(species_pp=="Tc") return -3.44648;
    if(species_pp=="Tc_pv") return -3.53766;
    if(species_pp=="Te") return -0.730289;
    if(species_pp=="Th_s") return -0.948771;
    if(species_pp=="Ti") return -2.23219;
    if(species_pp=="Ti_pv") return -2.61834;
    if(species_pp=="Ti_sv") return -2.64228;
    if(species_pp=="Tl") return -0.196234;
    if(species_pp=="Tl_d") return -0.360034;
    if(species_pp=="Tm_3") return -0.287414;
    if(species_pp=="W") return -4.53997;
    if(species_pp=="W_pv") return -4.65963;
    if(species_pp=="Xe") return -0.008491;
    if(species_pp=="Y_sv") return -2.2892;
    if(species_pp=="Yb") return -0.026493;
    if(species_pp=="Yb_2") return -0.062556;
    if(species_pp=="Zn") return -0.163774;
    if(species_pp=="Zr") return -2.11254;
    if(species_pp=="Zr_sv") return -2.30965;
    return NNN;
  }

  double getLIB1EnergyPAWPBE(const string& species_pp){
    if(species_pp=="Ac") return -4.09393;
    if(species_pp=="Ac_s") return -4.03196;
    if(species_pp=="Ag") return -2.82746;
    if(species_pp=="Al") return -3.74356;
    if(species_pp=="Al_h") return -3.80021;
    if(species_pp=="Ar") return -0.065575;
    if(species_pp=="As") return -4.65243;
    if(species_pp=="Au") return -3.27212;
    if(species_pp=="Ba_sv") return -1.924;
    if(species_pp=="Be") return -3.75366;
    if(species_pp=="Be_sv") return -3.74102;
    if(species_pp=="Bi") return -3.87274;
    if(species_pp=="Bi_d") return -4.03716;
    if(species_pp=="C") return -9.22034;
    if(species_pp=="C_h") return -9.19568;
    if(species_pp=="C_s") return -9.19508;
    if(species_pp=="Ca_pv") return -1.97638;
    if(species_pp=="Ca_sv") return -2.00107;
    if(species_pp=="Cd") return -0.906233;
    if(species_pp=="Cl") return -1.78721;
    if(species_pp=="Cl_h") return -1.77694;
    if(species_pp=="Co") return -7.10876;
    if(species_pp=="Cr") return -9.51285;
    if(species_pp=="Cr_pv") return -9.6294;
    if(species_pp=="Cs_sv") return -0.852451;
    if(species_pp=="Cu") return -3.71935;
    if(species_pp=="Cu_pv") return -4.09701;
    if(species_pp=="Dy_3") return -4.58813;
    if(species_pp=="Eu") return -10.2377;
    if(species_pp=="F") return -1.85953;
    if(species_pp=="F_h") return -1.90899;
    if(species_pp=="F_s") return -1.79087;
    if(species_pp=="Fe") return -8.31138;
    if(species_pp=="Fe_pv") return -8.45502;
    if(species_pp=="Ga_h") return -2.90153;
    if(species_pp=="Ge") return -4.49261;
    if(species_pp=="Ge_d") return -4.62213;
    if(species_pp=="Ge_h") return -4.50372;
    if(species_pp=="H") return -3.38635;
    if(species_pp=="H_h") return -3.36923;
    if(species_pp=="He") return 0.236864;
    if(species_pp=="Hf") return -9.95711;
    if(species_pp=="Hf_pv") return -9.95294;
    if(species_pp=="Ho_3") return -4.56866;
    if(species_pp=="I") return -1.51735;
    if(species_pp=="In_d") return -2.72115;
    if(species_pp=="Ir") return -8.85711;
    if(species_pp=="K_pv") return -1.02692;
    if(species_pp=="K_sv") return -1.09647;
    if(species_pp=="La") return -4.91782;
    if(species_pp=="La_s") return -4.86223;
    if(species_pp=="Li") return -1.89739;
    if(species_pp=="Mg") return -1.54148;
    if(species_pp=="Mg_pv") return -1.59348;
    if(species_pp=="Mn") return -9.028;
    if(species_pp=="Mo") return -10.9465;
    if(species_pp=="Mo_pv") return -10.8439;
    if(species_pp=="N") return -8.3187;
    if(species_pp=="N_h") return -8.3242;
    if(species_pp=="Na") return -1.3064;
    if(species_pp=="Na_sv") return -1.3139;
    if(species_pp=="Nb_pv") return -10.083;
    if(species_pp=="Nb_sv") return -10.2253;
    if(species_pp=="Ne") return -0.032434;
    if(species_pp=="Ni") return -5.57108;
    if(species_pp=="Ni_pv") return -5.77783;
    if(species_pp=="O") return -4.93114;
    if(species_pp=="O_h") return -5.01797;
    if(species_pp=="O_s") return -4.70103;
    if(species_pp=="Os") return -11.244;
    if(species_pp=="Os_pv") return -11.2193;
    if(species_pp=="P") return -5.32404;
    if(species_pp=="Pb") return -3.57056;
    if(species_pp=="Pb_d") return -3.7047;
    if(species_pp=="Pd") return -5.17847;
    if(species_pp=="Pd_pv") return -5.38197;
    if(species_pp=="Pt") return -6.05454;
    if(species_pp=="Re") return -12.4113;
    if(species_pp=="Re_pv") return -12.4325;
    if(species_pp=="Rh") return -7.27044;
    if(species_pp=="Rh_pv") return -7.34058;
    if(species_pp=="Ru") return -9.20414;
    if(species_pp=="Ru_pv") return -9.27135;
    if(species_pp=="S") return -4.12636;
    if(species_pp=="Sc_sv") return -6.33212;
    if(species_pp=="Se") return -3.48266;
    if(species_pp=="Si") return -5.42373;
    if(species_pp=="Si_h") return -5.44184;
    if(species_pp=="Sm_3") return -4.71227;
    if(species_pp=="Sn") return -3.79514;
    if(species_pp=="Sn_d") return -3.96372;
    if(species_pp=="Sr_sv") return -1.68354;
    if(species_pp=="Ta") return -11.86;
    if(species_pp=="Ta_pv") return -11.8489;
    if(species_pp=="Tc") return -10.3047;
    if(species_pp=="Tc_pv") return -10.3597;
    if(species_pp=="Te") return -3.14141;
    if(species_pp=="Ti") return -7.76396;
    if(species_pp=="Ti_pv") return -7.89056;
    if(species_pp=="Ti_sv") return -7.93863;
    if(species_pp=="Tl") return -2.2435;
    if(species_pp=="Tl_d") return -2.36274;
    if(species_pp=="V") return -8.94402;
    if(species_pp=="V_pv") return -9.07822;
    if(species_pp=="V_sv") return -9.11496;
    if(species_pp=="W") return -13.0124;
    if(species_pp=="W_pv") return -12.9546;
    if(species_pp=="Y_sv") return -6.46317;
    if(species_pp=="Yb") return -1.67034;
    if(species_pp=="Yb_2") return -1.51978;
    if(species_pp=="Zn") return -1.26581;
    if(species_pp=="Zr") return -8.47756;
    if(species_pp=="Zr_sv") return -8.54365;
    return NNN;
  }

  void insertElementalProperties(const vector<string>& vproperties,const xelement::xelement& xel,vector<string>& vitems) {
    uint i=0,j=0;
    int index=0;
    for(i=0;i<vproperties.size();i++){
      if(vproperties[i]=="lattice_constants"||vproperties[i]=="lattice_angles"){
        const xvector<double>& xvec=xel.getPropertyXVectorDouble(vproperties[i]);
        for(index=xvec.lrows;index<=xvec.urows;index++){
          vitems.push_back(aurostd::utype2string(xvec[index],_DOUBLE_WRITE_PRECISION_));
        }
      }
      else if(vproperties[i]=="energies_ionization"){
        //vitems.push_back(xel.getPropertyString(vproperties[i],",",2)); //only go to second ionization
        const vector<double> vec=xel.energies_ionization;
        for(j=0;j<_ENERGIES_IONIZATION_MAX_;j++){
          if(j>vec.size()-1){vitems.push_back(aurostd::utype2string(NNN,_DOUBLE_WRITE_PRECISION_));continue;}
          vitems.push_back(aurostd::utype2string(vec[j],_DOUBLE_WRITE_PRECISION_));
        }
      }
      else{vitems.push_back(xel.getPropertyString(vproperties[i],","));}
    }
  }
  void insertElementalPropertiesCCE(const vector<string>& vproperties,
      const xelement::xelement& xel,
      double M_X_bonds,
      double natoms_per_fu,
      vector<string>& vitems) {
    uint i=0,j=0;
    string units="";
    double d=0.0;
    for(i=0;i<vproperties.size();i++){
      units=xel.getUnits(vproperties[i]);
      if(units=="J/mol" && vproperties[i]=="energies_ionization"){
        const vector<double> vec=xel.energies_ionization;
        for(j=0;j<_ENERGIES_IONIZATION_MAX_;j++){
          if(j>vec.size()-1){vitems.push_back(aurostd::utype2string(NNN,_DOUBLE_WRITE_PRECISION_));continue;}
          d=vec[j];
          if(aurostd::isNaN(d)){vitems.push_back(aurostd::utype2string(NNN,_DOUBLE_WRITE_PRECISION_));continue;}
          d*=natoms_per_fu/M_X_bonds;
          vitems.push_back(aurostd::utype2string( d ,_DOUBLE_WRITE_PRECISION_));
        }
      }
      else if(units=="J/mol"||units=="J"){
        d=xel.getPropertyDouble(vproperties[i]);
        if(aurostd::isNaN(d)){vitems.push_back(aurostd::utype2string(NNN,_DOUBLE_WRITE_PRECISION_));continue;}
        d*=natoms_per_fu/M_X_bonds;
        vitems.push_back(aurostd::utype2string( d ,_DOUBLE_WRITE_PRECISION_));
      }
      else if(units=="K"){
        d=xel.getPropertyDouble(vproperties[i]);
        if(aurostd::isNaN(d)){vitems.push_back(aurostd::utype2string(NNN,_DOUBLE_WRITE_PRECISION_));continue;}
        d*=KBOLTZEV*natoms_per_fu/M_X_bonds;
        vitems.push_back(aurostd::utype2string( d ,_DOUBLE_WRITE_PRECISION_));
      }
    }
  }
  void insertCrystalProperties(const string& structure_path,const string& anion,const vector<string>& vheaders,vector<string>& vitems,aflowlib::_aflowlib_entry& entry,const string& e_props) {
    bool LDEBUG=(TRUE || _DEBUG_ML_ || XHOST.DEBUG);
    string soliloquy=XPID+"aflowML::insertCrystalProperties():";
    uint i=0,index_cation=0,index_anion=0;
    string entry_path="/common/LIB2/RAW/"+structure_path+"/aflowlib.out";
    if(!aurostd::FileExist(entry_path)){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,entry_path+" not found",_FILE_NOT_FOUND_);}
    entry.file2aflowlib(entry_path);
    //
    //stoich features
    vector<string> vheaders_stoich;
    vector<double> vfeatures_stoich;
    entry.getStoichFeatures(vheaders_stoich,vfeatures_stoich,false,e_props);  //_AFLOW_XELEMENT_PROPERTIES_ALL_
    for(i=0;i<vfeatures_stoich.size();i++){vitems.push_back(aurostd::utype2string(vfeatures_stoich[i],_DOUBLE_WRITE_PRECISION_));}
    //
    for(i=0;i<entry.vgeometry.size();i++){vitems.push_back(aurostd::utype2string(entry.vgeometry[i],_DOUBLE_WRITE_PRECISION_));}
    vitems.push_back(aurostd::utype2string(entry.vgeometry[0]/entry.vgeometry[1],_DOUBLE_WRITE_PRECISION_));
    vitems.push_back(aurostd::utype2string(entry.vgeometry[1]/entry.vgeometry[2],_DOUBLE_WRITE_PRECISION_));
    vitems.push_back(aurostd::utype2string(entry.vgeometry[0]/entry.vgeometry[2],_DOUBLE_WRITE_PRECISION_));
    vitems.push_back(aurostd::utype2string(entry.natoms));
    index_cation=index_anion=AUROSTD_MAX_UINT;
    if(entry.vspecies.size()!=2){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"non-binary found",_FILE_CORRUPT_);}
    if(!(entry.vspecies[0]==anion||entry.vspecies[1]==anion)){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"no "+anion+" found",_FILE_CORRUPT_);}
    if(entry.vspecies[0]==anion){index_cation=1;index_anion=0;}else{index_cation=0;index_anion=1;}
    vitems.push_back(aurostd::utype2string(entry.vcomposition[index_cation]));
    vitems.push_back(aurostd::utype2string(entry.vcomposition[index_anion]));
    vitems.push_back(aurostd::utype2string(entry.vstoichiometry[index_cation]));
    vitems.push_back(aurostd::utype2string(entry.vstoichiometry[index_anion]));
    vitems.push_back(aurostd::utype2string(entry.density,_DOUBLE_WRITE_PRECISION_));
    vitems.push_back(aurostd::utype2string(entry.vspecies_pp_ZVAL[index_cation]));
    vitems.push_back(aurostd::utype2string(entry.vspecies_pp_ZVAL[index_anion]));
    vitems.push_back(aurostd::utype2string(entry.valence_cell_iupac));
    vitems.push_back(aurostd::utype2string(entry.valence_cell_std));
    vitems.push_back(aurostd::utype2string(entry.volume_cell,_DOUBLE_WRITE_PRECISION_));
    vitems.push_back(aurostd::utype2string(entry.volume_atom,_DOUBLE_WRITE_PRECISION_));
    vitems.push_back(aurostd::utype2string(entry.energy_cell,_DOUBLE_WRITE_PRECISION_));
    vitems.push_back(aurostd::utype2string(entry.energy_atom,_DOUBLE_WRITE_PRECISION_));
    vitems.push_back(aurostd::utype2string(entry.energy_cutoff,_DOUBLE_WRITE_PRECISION_));
    vitems.push_back(aurostd::utype2string(entry.enthalpy_formation_cell,_DOUBLE_WRITE_PRECISION_));
    vitems.push_back(aurostd::utype2string(entry.enthalpy_formation_atom,_DOUBLE_WRITE_PRECISION_));
    if(entry.vnbondxx.size()!=3){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"non-binary found (nbondxx)",_FILE_CORRUPT_);}
    vitems.push_back(aurostd::utype2string(entry.vnbondxx[ (index_cation==0?0:2) ],_DOUBLE_WRITE_PRECISION_));
    vitems.push_back(aurostd::utype2string(entry.vnbondxx[ 1 ],_DOUBLE_WRITE_PRECISION_));
    vitems.push_back(aurostd::utype2string(entry.vnbondxx[ (index_cation==0?2:0) ],_DOUBLE_WRITE_PRECISION_));
    //for(i=0;i<entry.vnbondxx.size();i++){vitems.push_back(aurostd::utype2string(entry.vnbondxx[i],_DOUBLE_WRITE_PRECISION_));}
    vitems.push_back(aurostd::utype2string(entry.spacegroup_relax));
    vitems.push_back(aurostd::utype2string(entry.point_group_order));
    vitems.push_back(entry.Bravais_lattice_relax);
    //
    //patch oxidation_anion
    bool found=false;
    double oxidation_cation=0;
    for(i=0;i<vheaders.size()&&!found;i++){
      if(vheaders[i]=="oxidation_cation"){
        oxidation_cation=aurostd::string2utype<double>(vitems[i]);
        found=true;
      }
    }
    if(!found){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"oxidation_cation not found",_RUNTIME_ERROR_);}
    double oxidation_anion=-(entry.vcomposition[index_cation]*oxidation_cation)/entry.vcomposition[index_anion];
    if(LDEBUG){
      cerr << soliloquy << " oxidation_cation=" << oxidation_cation << endl;
      cerr << soliloquy << " oxidation_anion=" << oxidation_anion << endl;
    }
    found=false;
    for(i=0;i<vheaders.size()&&!found;i++){
      if(vheaders[i]=="oxidation_anion"){
        vitems[i]=aurostd::utype2string(oxidation_anion);
        found=true;
      }
    }
    if(!found){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"oxidation_anion not found",_RUNTIME_ERROR_);}
    //
  }
  double getStatistic(const xvector<double>& xvec,const string& stat){
    if(stat=="min"){return aurostd::min(xvec);}
    if(stat=="max"){return aurostd::max(xvec);}
    if(stat=="range"){return aurostd::max(xvec)-aurostd::min(xvec);}
    if(stat=="sum"){return aurostd::sum(xvec);}
    if(stat=="mean"){return aurostd::mean(xvec);}
    if(stat=="std"){return aurostd::stddev(xvec);}
    if(stat=="mode"){return aurostd::mode(xvec);}
    throw aurostd::xerror(_AFLOW_FILE_NAME_,"aflowML::getStatistic():","Unknown statistic: "+stat,_FILE_CORRUPT_);
    return 0.0;
  }
  void insertElementalCombinations(const vector<string>& vproperties,vector<string>& vheaders){
    xelement::xelement xel1;
    xelement::xelement xel2;
    aflowlib::_aflowlib_entry entry;
    vector<double> vfeatures;
    return insertElementalCombinations(vproperties,xel1,xel2,entry,1.0,1.0,1.0,vheaders,vfeatures,true);
  }
  void insertElementalCombinations(
      const vector<string>& vproperties,
      const xelement::xelement& xel_cation,
      const xelement::xelement& xel_anion,
      const aflowlib::_aflowlib_entry& entry,
      double M_X_bonds,
      double natoms_per_fu_cation,
      double natoms_per_fu_anion,
      vector<string>& vheaders,
      vector<double>& vfeatures,
      bool vheaders_only,
      uint count_vcols){
    bool LDEBUG=(FALSE || _DEBUG_ML_ || XHOST.DEBUG);
    string soliloquy=XPID+"aflowML::insertElementalCombinations():";
    vheaders.clear();vfeatures.clear();
    if(!_INJECT_ELEMENTAL_COMBINATIONS_){return;}
    if(count_vcols!=AUROSTD_MAX_UINT){vheaders.reserve(count_vcols);vfeatures.reserve(count_vcols);}
    else{
      if(vheaders_only==false){
        //get vheaders.size() and resize vfeatures
        insertElementalCombinations(vproperties,xel_cation,xel_anion,entry,M_X_bonds,natoms_per_fu_cation,natoms_per_fu_anion,vheaders,vfeatures,true);
        uint vheaders_size=vheaders.size();
        vheaders.clear();vfeatures.clear();
        vheaders.reserve(vheaders_size);vfeatures.reserve(vheaders_size);
      }
    }
    
    vector<string> vlattice_constants_variants;aurostd::string2tokens("a,b,c",vlattice_constants_variants,",");
    vector<string> vlattice_angles_variants;aurostd::string2tokens("alpha,beta,gamma",vlattice_angles_variants,",");
    vector<string> vions;aurostd::string2tokens("cation,anion",vions,",");
    vector<string> venvs;aurostd::string2tokens("crystal,atomenv",venvs,",");
    vector<string> vstats;aurostd::string2tokens("min,max,range,sum,mean,std,mode",vstats,",");
    
    uint index_cation=0,index_anion=0;
    if(!vheaders_only){
      if(entry.vspecies.size()!=2){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"non-binary found",_FILE_CORRUPT_);}
      if(!(entry.vspecies[0]==xel_anion.symbol||entry.vspecies[1]==xel_anion.symbol)){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"no "+xel_anion.symbol+" found",_FILE_CORRUPT_);}
      if(entry.vspecies[0]==xel_anion.symbol){index_cation=1;index_anion=0;}else{index_cation=0;index_anion=1;}
    }

    vector<int> vsizes;
    aurostd::xcombos xc;
    uint count=0;
    bool divide_by_adjacency=false;
    bool multiply_by_natoms=false;
    string property="";
    int index_property=0;
    int index_vec=0;
    int index_env=0;
    int index_ion=0;
    int index_adjacency=0;
    int index_stat=0;
    xvector<double> xvec_env;
    uint ncation=0,nanion=0,i=0,j=0;
    bool has_NaN=false;
    double d=0.0;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //SINGLE
    //last is fastest iterator
    vsizes.clear();
    vsizes.push_back((int)vproperties.size());
    vsizes.push_back((int)venvs.size());
    vsizes.push_back(_ENERGIES_IONIZATION_MAX_);  //index for (x)vector properties (lattice_constants, lattice_angles, energies_ionization), max is 5
    vsizes.push_back((int)vions.size());
    vsizes.push_back(3);  //normal, divide by M-X_bonds, multiply by natoms and dvided by M-X_bonds
    vsizes.push_back((int)vstats.size());
    
    xc.reset(vsizes,'E');
    count=0;
    divide_by_adjacency=false;
    multiply_by_natoms=false;
    property="";
    index_property=0;
    index_vec=0;
    index_env=0;
    index_ion=0;
    index_adjacency=0;
    index_stat=0;
    xvec_env.clear();
    ncation=0;nanion=0;i=0;j=0;
    has_NaN=false;
    d=0.0;
    while(xc.increment()){
      //
      const vector<int>& indices=xc.getCombo();
      //
      index_property=indices[0];
      index_env=indices[1];
      index_vec=indices[2]; //order this way so we resize xvector as little as possible
      index_ion=indices[3];
      index_adjacency=indices[4];
      index_stat=indices[5];
      //
      const string& property_element=vproperties[index_property];
      const string& env=venvs[index_env];
      const string& ion=vions[index_ion];
      divide_by_adjacency=(index_adjacency!=0);
      multiply_by_natoms=(index_adjacency==2);
      const string& stat=vstats[index_stat];
      //
      if(env=="crystal"&&index_ion>0){continue;}  //crystal does not need cation/anion  //ion!=vions[0]
      if(!(
            property_element=="lattice_constants"||
            property_element=="lattice_angles"||
            property_element=="energies_ionization"||
            false)
          &&index_vec>0){
        continue;
      }
      if((property_element=="lattice_constants"||property_element=="lattice_angles")&&index_vec>2){continue;}
      //
      property="";
      if(divide_by_adjacency){property+="(";}
      property+=property_element;
      if(property_element=="lattice_constants"){property+="_"+vlattice_constants_variants[index_vec];}
      else if(property_element=="lattice_angles"){property+="_"+vlattice_angles_variants[index_vec];}
      else if(property_element=="energies_ionization"){property+="_"+aurostd::utype2string(index_vec+1);}
      if(multiply_by_natoms){property+="_*_natoms_per_f.u.";}
      if(divide_by_adjacency){property+=")_/_(M-X_bonds)";}
      property+="_"+env;
      if(env=="atomenv"){property+="_"+ion;}
      property+="_"+stat;
      //
      if(LDEBUG){cerr << soliloquy << " count_single=" << count++ << " " << property << endl;}
      vheaders.push_back(property);
      if(vheaders_only==false){
        if(index_stat==0){  //do not waste cycles on recreating xvec over and over again
          if(env=="crystal"){xvec_env.resize(entry.natoms);ncation=entry.vcomposition[index_cation];nanion=entry.vcomposition[index_anion];}
          else if(env=="atomenv"){
            xvec_env.resize((int)std::ceil(M_X_bonds)+1);
            if(ion=="cation"){ncation=1;nanion=(uint)std::ceil(M_X_bonds);}
            else if(ion=="anion"){ncation=(uint)std::ceil(M_X_bonds);nanion=1;}
            else{throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"Unknown ion: "+ion,_RUNTIME_ERROR_);}
          }
          else{throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"Unknown env: "+env,_RUNTIME_ERROR_);}
          i=0;
          has_NaN=false;
          //cation
          d=xel_cation.getPropertyDouble(property_element,index_vec);
          has_NaN=(has_NaN || aurostd::isNaN(d));
          if(has_NaN==false && multiply_by_natoms){d*=natoms_per_fu_cation;}
          for(j=0;j<ncation;j++){xvec_env[xvec_env.lrows+(i++)]=d;}
          //anion
          d=xel_anion.getPropertyDouble(property_element,index_vec);
          has_NaN=(has_NaN || aurostd::isNaN(d));
          if(has_NaN==false && multiply_by_natoms){d*=natoms_per_fu_anion;}
          for(j=0;j<nanion;j++){xvec_env[xvec_env.lrows+(i++)]=d;}
          //
          if(has_NaN==false && divide_by_adjacency){xvec_env/=std::ceil(M_X_bonds);}
          if(LDEBUG){
            cerr << soliloquy << " xvec[\""+property+"\"]=" << xvec_env << endl;
            cerr << soliloquy << " has_NaN=" << has_NaN << endl;
          }
        }
        if(has_NaN){d=NNN;}
        else{d=getStatistic(xvec_env,stat);}
        if(LDEBUG){cerr << soliloquy << " " << stat << "=" << d << endl;}
        vfeatures.push_back(d);
      }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    int index_property2=0;
    int index_vec2=0;
    int index_operation=0;
    double d2=0.0;
    double d3=0.0;

    //convert to SI units for combinations
    if(false){
      xelement::xelement xel_cation_SI=xel_cation;xel_cation_SI.convertUnits();
      xelement::xelement xel_anion_SI=xel_anion;xel_anion_SI.convertUnits();
    }
    const xelement::xelement xel_cation_SI=xel_cation;
    const xelement::xelement xel_anion_SI=xel_anion;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //DOUBLE
    //last is fastest iterator
    vsizes.clear();
    vsizes.push_back((int)vproperties.size());
    vsizes.push_back((int)vproperties.size());
    vsizes.push_back((int)venvs.size());
    vsizes.push_back(_ENERGIES_IONIZATION_MAX_);  //index for (x)vector properties (lattice_constants, lattice_angles, energies_ionization), max is 5
    vsizes.push_back(_ENERGIES_IONIZATION_MAX_);  //index for (x)vector properties (lattice_constants, lattice_angles, energies_ionization), max is 5
    vsizes.push_back((int)vions.size());
    vsizes.push_back(2);  //multiplication/division
    vsizes.push_back(3);  //normal, divide by M-X_bonds, multiply by natoms and dvided by M-X_bonds
    vsizes.push_back((int)vstats.size());
    
    xc.reset(vsizes,'E');
    count=0;
    divide_by_adjacency=false;
    multiply_by_natoms=false;
    property="";
    index_property=0;
    index_property2=0;
    index_vec=0;
    index_vec2=0;
    index_env=0;
    index_ion=0;
    index_operation=0;
    index_adjacency=0;
    index_stat=0;
    xvec_env.clear();
    ncation=0;nanion=0;i=0;j=0;
    has_NaN=false;
    d=0.0;
    d2=0.0;
    d3=0.0;
    while(xc.increment()){
      //
      const vector<int>& indices=xc.getCombo();
      //
      index_property=indices[0];
      index_property2=indices[1];
      if(index_property2<=index_property){continue;}
      index_env=indices[2];
      index_vec=indices[3]; //order this way so we resize xvector as little as possible
      index_vec2=indices[4]; //order this way so we resize xvector as little as possible
      index_ion=indices[5];
      index_operation=indices[6];
      index_adjacency=indices[7];
      index_stat=indices[8];
      //
      const string& property_element=vproperties[index_property];
      const string& property_element2=vproperties[index_property2];
      const string& env=venvs[index_env];
      const string& ion=vions[index_ion];
      divide_by_adjacency=(index_adjacency!=0);
      multiply_by_natoms=(index_adjacency==2);
      const string& stat=vstats[index_stat];
      //
      if(env=="crystal"&&index_ion>0){continue;}  //crystal does not need cation/anion  //ion!=vions[0]
      if(!(
            property_element=="lattice_constants"||
            property_element=="lattice_angles"||
            property_element=="energies_ionization"||
            false)
          &&index_vec>0){
        continue;
      }
      if(!(
            property_element2=="lattice_constants"||
            property_element2=="lattice_angles"||
            property_element2=="energies_ionization"||
            false)
          &&index_vec2>0){
        continue;
      }
      if((property_element=="lattice_constants"||property_element=="lattice_angles")&&index_vec>2){continue;}
      if((property_element2=="lattice_constants"||property_element2=="lattice_angles")&&index_vec2>2){continue;}
      //
      property="";
      property+="(";
      //
      property+=property_element;
      if(property_element=="lattice_constants"){property+="_"+vlattice_constants_variants[index_vec];}
      else if(property_element=="lattice_angles"){property+="_"+vlattice_angles_variants[index_vec];}
      else if(property_element=="energies_ionization"){property+="_"+aurostd::utype2string(index_vec+1);}
      //
      if(index_operation==0){property+="_*_";}
      else{property+="_/_";}
      //
      property+=property_element2;
      if(property_element2=="lattice_constants"){property+="_"+vlattice_constants_variants[index_vec2];}
      else if(property_element2=="lattice_angles"){property+="_"+vlattice_angles_variants[index_vec2];}
      else if(property_element2=="energies_ionization"){property+="_"+aurostd::utype2string(index_vec2+1);}
      //
      if(multiply_by_natoms){property+="_*_natoms_per_f.u.";}
      //
      property+=")";
      if(divide_by_adjacency){property+="_/_(M-X_bonds)";}
      property+="_"+env;
      if(env=="atomenv"){property+="_"+ion;}
      property+="_"+stat;
      //
      if(LDEBUG){cerr << soliloquy << " count_double=" << count++ << " " << property << endl;}
      vheaders.push_back(property);
      if(vheaders_only==false){
        if(index_stat==0){  //do not waste cycles on recreating xvec over and over again
          if(env=="crystal"){xvec_env.resize(entry.natoms);ncation=entry.vcomposition[index_cation];nanion=entry.vcomposition[index_anion];}
          else if(env=="atomenv"){
            xvec_env.resize((int)std::ceil(M_X_bonds)+1);
            if(ion=="cation"){ncation=1;nanion=(uint)std::ceil(M_X_bonds);}
            else if(ion=="anion"){ncation=(uint)std::ceil(M_X_bonds);nanion=1;}
            else{throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"Unknown ion: "+ion,_RUNTIME_ERROR_);}
          }
          else{throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"Unknown env: "+env,_RUNTIME_ERROR_);}
          i=0;
          has_NaN=false;
          //cation
          d=xel_cation_SI.getPropertyDouble(property_element,index_vec);
          d2=xel_cation_SI.getPropertyDouble(property_element2,index_vec2);
          has_NaN=(has_NaN || aurostd::isNaN(d));
          has_NaN=(has_NaN || aurostd::isNaN(d2));
          if(index_operation==0){d3=d*d2;} //multiplication
          else{d3=d/d2;}  //division
          if(has_NaN==false && multiply_by_natoms){d3*=natoms_per_fu_cation;}
          //check issues with multiplication/division
          if(!std::isfinite(d3)){
            if(LDEBUG){
              if(index_operation==0){cerr << soliloquy << " multiplication yields inf: " << d << " * " << d2 << endl;}
              else{cerr << soliloquy << " division yields inf: " << d << " / " << d2 << endl;}
            }
            has_NaN=true;
          }
          //
          for(j=0;j<ncation;j++){xvec_env[xvec_env.lrows+(i++)]=d3;}
          //anion
          d=xel_anion_SI.getPropertyDouble(property_element,index_vec);
          d2=xel_anion_SI.getPropertyDouble(property_element2,index_vec2);
          has_NaN=(has_NaN || aurostd::isNaN(d));
          has_NaN=(has_NaN || aurostd::isNaN(d2));
          if(index_operation==0){d3=d*d2;} //multiplication
          else{d3=d/d2;}  //division
          if(has_NaN==false && multiply_by_natoms){d3*=natoms_per_fu_anion;}
          //check issues with multiplication/division
          if(!std::isfinite(d3)){
            if(LDEBUG){
              if(index_operation==0){cerr << soliloquy << " multiplication yields inf: " << d << " * " << d2 << endl;}
              else{cerr << soliloquy << " division yields inf: " << d << " / " << d2 << endl;}
            }
            has_NaN=true;
          }
          //
          for(j=0;j<nanion;j++){xvec_env[xvec_env.lrows+(i++)]=d3;}
          //
          if(has_NaN==false && divide_by_adjacency){xvec_env/=std::ceil(M_X_bonds);}
          if(LDEBUG){
            cerr << soliloquy << " xvec[\""+property+"\"]=" << xvec_env << endl;
            cerr << soliloquy << " has_NaN=" << has_NaN << endl;
          }
        }
        if(has_NaN){d=NNN;}
        else{d=getStatistic(xvec_env,stat);}
        if(LDEBUG){cerr << soliloquy << " " << stat << "=" << d << endl;}
        vfeatures.push_back(d);
      }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  }
  void getColumn(const vector<vector<string> >& table,uint icol,vector<string>& column,bool& isfloat,bool& isinteger,bool include_header) {
    isfloat=false;
    isinteger=false;
    column.clear();
    for(uint j=(include_header?0:1);j<table.size();j++){
      isfloat=( isfloat || aurostd::isfloat(table[j][icol]) );
      isinteger=( isinteger || (isfloat && aurostd::isinteger(aurostd::string2utype<double>(table[j][icol]))) );
      column.push_back(table[j][icol]);
    }
  }
  void delColumn(vector<vector<string> >& table,uint icol){
    for(uint j=0;j<table.size();j++){
      table[j].erase(table[j].begin()+icol);
    }
  }
  void oneHotFeatures(vector<vector<string> >& table,const string& features_categories) {
    bool LDEBUG=(FALSE || _DEBUG_ML_ || XHOST.DEBUG);
    string soliloquy=XPID+"aflowML::oneHotFeatures():";
    stringstream message;

    if(LDEBUG){cerr << soliloquy << " BEGIN" << endl;}
    
    message << "creating one-hot features for " << features_categories;
    pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,cout,_LOGGER_MESSAGE_);
    
    uint i=0,j=0,k=0;

    if(LDEBUG){
      cerr << soliloquy << " table_orig=" << endl;
      for(i=0;i<table.size();i++){cerr << aurostd::joinWDelimiter(table[i],",") << endl;}
    }

    if(table.size()<2){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"table.size()<2",_RUNTIME_ERROR_);}

    vector<string> vfeatures_categories;
    aurostd::string2tokens(features_categories,vfeatures_categories,",");
    if(LDEBUG){cerr << soliloquy << " vfeatures_categories.size()=" << vfeatures_categories.size() << endl;}

    //table[0] are headers
    vector<uint> vicol;
    bool found=false;
    for(i=0;i<vfeatures_categories.size();i++){
      found=false;
      for(j=0;j<table[0].size()&&!found;j++){
        if(vfeatures_categories[i]==table[0][j]){
          vicol.push_back(j);
          found=true;
        }
      }
    }
    if(vfeatures_categories.size()!=vicol.size()){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"vfeatures_categories.size()!=vicol.size()",_RUNTIME_ERROR_);}
    
    //sort and reverse, so we can erase/insert
    std::sort(vicol.rbegin(),vicol.rend());

    bool isfloat=false,isinteger=false;
    vector<string> column_orig,column;
    vector<int> column_int;
    string header="",header_new="";
    for(i=0;i<vicol.size();i++){
      column.clear();
      header=table[0][vicol[i]];
      getColumn(table,vicol[i],column,isfloat,isinteger,false);
      if(LDEBUG){
        cerr << soliloquy << " column[" << header << "]=" << aurostd::joinWDelimiter(column,",") << endl;
        cerr << soliloquy << " isinteger=" << isinteger << endl;
      }
      column_orig.clear();for(j=0;j<column.size();j++){column_orig.push_back(column[j]);} //column_orig=column;
      if(isinteger){
        column_int.clear();
        for(j=0;j<column.size();j++){column_int.push_back(aurostd::string2utype<int>(column[j]));}
        std::sort(column_int.begin(),column_int.end());column_int.erase( std::unique( column_int.begin(), column_int.end() ), column_int.end() );  //get unique values
        column.clear();
        for(j=0;j<column_int.size();j++){column.push_back(aurostd::utype2string(column_int[j]));}
      }else{
        std::sort(column.begin(),column.end());column.erase( std::unique( column.begin(), column.end() ), column.end() );  //get unique values
      }
      if(LDEBUG){cerr << soliloquy << " unique=" << aurostd::joinWDelimiter(column,",") << endl;}
      delColumn(table,vicol[i]);
      for(k=0;k<column.size();k++){
        header_new=header+"_"+column[k];
        if(LDEBUG){cerr << soliloquy << " header_new=" << header_new << endl;}
        table[0].insert(table[0].begin()+vicol[i]+k,header_new);
        for(j=1;j<table.size();j++){
          table[j].insert(table[j].begin()+vicol[i]+k, column_orig[j-1]==column[k]?"1":"0" );
        }
      }
    }

    if(LDEBUG){
      cerr << soliloquy << " table_new=" << endl;
      for(i=0;i<table.size();i++){cerr << aurostd::joinWDelimiter(table[i],",") << endl;}
    }

    //check we didn't mess up
    uint ncols=table[0].size();
    for(i=0;i<table.size();i++){
      if(table[i].size()!=ncols){
        throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"table[i="+aurostd::utype2string(i)+"].size()!=ncols",_RUNTIME_ERROR_);
      }
    }
  }
  void replaceNaN(xvector<double>& xvec,double val){for(int i=xvec.lrows;i<=xvec.urows;i++){if(aurostd::isNaN(xvec[i])){xvec[i]=val;}}}
  void removeNaN(const xvector<double>& xvec,xvector<double>& xvec_new){
    int count=0;
    int i=0,j=0;
    for(i=xvec.lrows;i<=xvec.urows;i++){if(!aurostd::isNaN(xvec[i])){count++;}}
    if(count==0){xvec_new=aurostd::null_xv<double>();return;}
    xvec_new.resize(xvec.lrows,xvec.lrows+count-1);
    j=xvec_new.lrows;
    for(i=xvec.lrows;i<=xvec.urows;i++){if(!aurostd::isNaN(xvec[i])){xvec_new[j++]=xvec[i];}}
  }
  void MinMaxScale(xvector<double>& xvec){
    double min=aurostd::min(xvec);
    double max=aurostd::max(xvec);
    double denom=max-min;
    if(denom==0.0){xvec.reset();} //set it all to 0
    else{for(int i=xvec.lrows;i<=xvec.urows;i++){xvec[i]=(xvec[i]-min)/denom;}}
  }
  void reduceFeatures(vector<vector<string> >& table,const string& yheader,double var_threshold,double ycorr_threshold,double selfcorr_threshold) {
    vector<uint> vicol2skip;
    return reduceFeatures(table,yheader,vicol2skip,var_threshold,ycorr_threshold,selfcorr_threshold);
  }
  void reduceFeatures(vector<vector<string> >& table,const string& yheader,const string& header2skip,double var_threshold,double ycorr_threshold,double selfcorr_threshold) {
    vector<string> vheaders2skip;aurostd::string2tokens(header2skip,vheaders2skip,",");
    return reduceFeatures(table,yheader,vheaders2skip,var_threshold,ycorr_threshold,selfcorr_threshold);
  }
  void reduceFeatures(vector<vector<string> >& table,const string& yheader,const vector<string>& vheaders2skip,double var_threshold,double ycorr_threshold,double selfcorr_threshold) {
    bool LDEBUG=(FALSE || _DEBUG_ML_ || XHOST.DEBUG);
    string soliloquy=XPID+"aflowML::reduceFeatures():";
    
    if(LDEBUG){cerr << soliloquy << " BEGIN" << endl;}
    
    uint i=0,j=0;

    //table[0] are headers
    vector<uint> vicol2skip;
    bool found=false;
    for(i=0;i<vheaders2skip.size();i++){
      found=false;
      for(j=0;j<table[0].size()&&!found;j++){
        if(vheaders2skip[i]==table[0][j]){
          vicol2skip.push_back(j);
          found=true;
        }
      }
    }
    if(vheaders2skip.size()!=vicol2skip.size()){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"vheaders2skip.size()!=vicol2skip.size()",_RUNTIME_ERROR_);}

    return reduceFeatures(table,yheader,vicol2skip,var_threshold,ycorr_threshold,selfcorr_threshold);
  }
  void reduceFeatures(vector<vector<string> >& table,const string& yheader,uint icol2skip,double var_threshold,double ycorr_threshold,double selfcorr_threshold) {
    vector<uint> vicol2skip;vicol2skip.push_back(icol2skip);
    return reduceFeatures(table,yheader,vicol2skip,var_threshold,ycorr_threshold,selfcorr_threshold);
  }
  void reduceFeatures(vector<vector<string> >& table,const string& yheader,const vector<uint>& vicol2skip,double var_threshold,double ycorr_threshold,double selfcorr_threshold) {
    bool LDEBUG=(FALSE || _DEBUG_ML_ || XHOST.DEBUG);
    string soliloquy=XPID+"aflowML::reduceFeatures():";
    stringstream message;

    if(LDEBUG){cerr << soliloquy << " BEGIN" << endl;}
    
    uint i=0,j=0,ncols=0,ncols_orig=0;

    if(LDEBUG){
      cerr << soliloquy << " table_orig=" << endl;
      for(i=0;i<table.size();i++){cerr << aurostd::joinWDelimiter(table[i],",") << endl;}
    }

    if(table.size()<2){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"table.size()<2",_RUNTIME_ERROR_);}
    
    ncols=ncols_orig=table[0].size();
    message << "ncols_orig=" << ncols;
    pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,cout,_LOGGER_MESSAGE_);
    
    message << "converting string table to xvector table";
    pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,cout,_LOGGER_MESSAGE_);

    vector<xvector<double> > xvtable;
    bool isfloat=false,isinteger=false;
    vector<string> column;
    xvector<double> xv,xv_clean;
    xvector<double> nullxv=aurostd::null_xv<double>();
    vector<uint> xvindices;
    vector<double> vmeans,vstddevs;
    uint yiheader=AUROSTD_MAX_UINT;
    for(i=0;i<ncols;i++){
      const string& header=table[0][i];
      if(LDEBUG){cerr << soliloquy << " looking at column " << header << endl;}
      if(header==yheader){yiheader=i;}
      if(aurostd::WithinList(vicol2skip,i)){
        if(LDEBUG){cerr << soliloquy << " skipping " << header << " as requested" << endl;}
        xvtable.push_back(nullxv);vmeans.push_back(0.0);vstddevs.push_back(0.0);continue;
      }
      getColumn(table,i,column,isfloat,isinteger,false);
      if(!isfloat){
        if(LDEBUG){cerr << soliloquy << " skipping " << header << ": not floats" << endl;}
        xvtable.push_back(nullxv);vmeans.push_back(0.0);vstddevs.push_back(0.0);continue;
      }
      xv=aurostd::vector2xvector<double>(column);
      if(LDEBUG){cerr << soliloquy << " xv=" << xv << endl;}
      xvtable.push_back(xv);vmeans.push_back(aurostd::mean(xv));vstddevs.push_back(aurostd::stddev(xv));
      xvindices.push_back(i);
    }
    if(xvtable.size()!=ncols){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"xvtable.size()!=ncols",_RUNTIME_ERROR_);}
    if(vmeans.size()!=ncols){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"vmeans.size()!=ncols",_RUNTIME_ERROR_);}
    if(vstddevs.size()!=ncols){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"vstddevs.size()!=ncols",_RUNTIME_ERROR_);}
    if(yiheader==AUROSTD_MAX_UINT){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"yiheader==AUROSTD_MAX_UINT",_INPUT_MISSING_);}
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    //variance
    
    message << "identifying and removing null / low-variance features";
    pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,cout,_LOGGER_MESSAGE_);

    //table[0] are headers
    //vector<uint> vicol2remove;
    xvector<int> xvicol2remove(ncols);  //if xvicol2remove[i]==1, then remove that column
    double var=0.0;
    uint index=0;
    for(i=0;i<xvindices.size();i++){
      index=xvindices[i];
      const string& header=table[0][index];
      if(LDEBUG){cerr << soliloquy << " looking at column " << header << endl;}
      const xvector<double>& xvec=xvtable[index];
      if(LDEBUG){cerr << soliloquy << " xvec(orig  )=" << xvec << endl;}
      removeNaN(xvec,xv_clean);
      if(LDEBUG){cerr << soliloquy << " xvec(clean )=" << xv_clean << endl;}
      if(xv_clean.rows==0||xv_clean.rows==1){
        if(LDEBUG){cerr << soliloquy << " removing " << header << ": " << (xv_clean.rows==0?"null-":"1-") << "vector" << endl;}
        //vicol2remove.push_back(index);
        xvicol2remove[xvicol2remove.lrows+index]=1;
        continue;
      }
      if(xv_clean.rows<(xvec.rows/2)){
        if(LDEBUG){cerr << soliloquy << " removing " << header << ": " << " vector more than half filled with NaNs (" << xv_clean.rows << " out of " << xvec.rows << ")" << endl;}
        //vicol2remove.push_back(index);
        xvicol2remove[xvicol2remove.lrows+index]=1;
        continue;
      }
      MinMaxScale(xv_clean);
      if(LDEBUG){cerr << soliloquy << " xvec(scaled)=" << xv_clean << endl;}
      var=aurostd::var(xv_clean,1); //sample variance
      if(LDEBUG){cerr << soliloquy << " var(xvec)=" << var << endl;}
      if(var<var_threshold){
        if(LDEBUG){cerr << soliloquy << " no variance in column " << header << endl;}
        //vicol2remove.push_back(index);
        xvicol2remove[xvicol2remove.lrows+index]=1;
        continue;
      }
    }

    //std::sort(vicol2remove.begin(),vicol2remove.end());vicol2remove.erase( std::unique( vicol2remove.begin(), vicol2remove.end() ), vicol2remove.end() );  //get unique values
    //message << "found " << vicol2remove.size() << " null / low-variance columns";
    message << "found " << aurostd::sum(xvicol2remove) << " null / low-variance columns";
    pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,cout,_LOGGER_MESSAGE_);

    //[ERASE IS TOO SLOW]if(0){  //erase is TOO slow
    //[ERASE IS TOO SLOW]  //sort and reverse, so we can erase/insert
    //[ERASE IS TOO SLOW]  std::sort(vicol2remove.rbegin(),vicol2remove.rend());
    //[ERASE IS TOO SLOW]
    //[ERASE IS TOO SLOW]  for(i=0;i<vicol2remove.size();i++){
    //[ERASE IS TOO SLOW]    if(LDEBUG){cerr << soliloquy << " removing column " << table[0][vicol2remove[i]] << endl;}
    //[ERASE IS TOO SLOW]    delColumn(table,vicol2remove[i]);
    //[ERASE IS TOO SLOW]  }
    //[ERASE IS TOO SLOW]}
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    //correlated with y
    
    message << "identifying and removing low-correlated features";
    pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,cout,_LOGGER_MESSAGE_);

    //get xvec of y
    getColumn(table,yiheader,column,isfloat,isinteger,false);
    if(!isfloat){
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"yvec is not float",_INPUT_ILLEGAL_);
    }
    xvector<double> yvec=aurostd::vector2xvector<double>(column);
    if(LDEBUG){cerr << soliloquy << " yvec[" << yheader << "]=" << yvec << endl;}
    double ymean=aurostd::mean(yvec);
    double ystddev=aurostd::stddev(yvec);

    double corr=0.0;
    for(i=0;i<xvindices.size();i++){
      index=xvindices[i];
      if(LDEBUG){cerr << soliloquy << " index=" << index << endl;}
      const string& header=table[0][index];
      //if(aurostd::WithinList(vicol2remove,index))
      if(xvicol2remove[xvicol2remove.lrows+index]==1){
        if(LDEBUG){cerr << soliloquy << " skipping " << header << ": to be removed" << endl;}
        continue;
      }
      const xvector<double>& xvec=xvtable[index]; //no need to remove NNN, correlation would be the same with 0 or NNN
      if(LDEBUG){cerr << soliloquy << " xvec[" << header << "]=" << xvec << endl;}
      
      corr=std::pow(aurostd::correlation_Pearson_fast(xvec,vmeans[index],vstddevs[index],yvec,ymean,ystddev),2.0);
      if(LDEBUG){cerr << soliloquy << " corr(xvec,yvec)^2=" << corr << endl;}
      if(corr<ycorr_threshold){
        if(LDEBUG){cerr << soliloquy << " low correlation between columns " << header << " and " <<  yheader << endl;}
        //vicol2remove.push_back(index);
        xvicol2remove[xvicol2remove.lrows+index]=1;
      }
    }

    //std::sort(vicol2remove.begin(),vicol2remove.end());vicol2remove.erase( std::unique( vicol2remove.begin(), vicol2remove.end() ), vicol2remove.end() );  //get unique values
    //message << "found " << vicol2remove.size() << " null / low-variance / low-correlated columns";
    message << "found " << aurostd::sum(xvicol2remove) << " null / low-variance / low-correlated columns";
    pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,cout,_LOGGER_MESSAGE_);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    //self-correlated
    
    message << "identifying and removing self-correlated features";
    pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,cout,_LOGGER_MESSAGE_);

    //uint ncpus=1;
    uint index2=0;

    //if(ncpus==1){
    for(i=0;i<xvindices.size()-1;i++){
      index=xvindices[i];
      if(1||LDEBUG){cerr << soliloquy << " index=" << index << endl;}
      const string& header=table[0][index];
      //if(aurostd::WithinList(vicol2remove,index))
      if(xvicol2remove[xvicol2remove.lrows+index]==1){
        if(LDEBUG){cerr << soliloquy << " skipping " << header << ": to be removed" << endl;}
        continue;
      }
      const xvector<double>& xvec=xvtable[index]; //no need to remove NNN, correlation would be the same with 0 or NNN
      if(LDEBUG){cerr << soliloquy << " xvec[" << header << "]=" << xvec << endl;}
      
      for(j=i+1;j<xvindices.size();j++){
        index2=xvindices[j];
        if(LDEBUG){cerr << soliloquy << " index2=" << index2 << endl;}
        const string& header2=table[0][index2];
        //if(aurostd::WithinList(vicol2remove,index2))
        if(xvicol2remove[xvicol2remove.lrows+index2]==1){
          if(LDEBUG){cerr << soliloquy << " skipping " << header2 << ": to be removed" << endl;}
          continue;
        }
        const xvector<double>& xvec2=xvtable[index2]; //no need to remove NNN, correlation would be the same with 0 or NNN
        if(LDEBUG){cerr << soliloquy << " xvec2[" << header2 << "]=" << xvec2 << endl;}

        corr=std::pow(aurostd::correlation_Pearson_fast(xvec,vmeans[index],vstddevs[index],xvec2,vmeans[index2],vstddevs[index2]),2.0);
        if(LDEBUG){cerr << soliloquy << " corr(xvec,xvec2)^2=" << corr << endl;}
        if(corr>selfcorr_threshold){
          if(LDEBUG){cerr << soliloquy << " high correlation between columns " << header << " and " <<  header2 << endl;}
          //vicol2remove.push_back(index2);
          xvicol2remove[xvicol2remove.lrows+index2]=1;
        }
      }
    }
    //}else{
    //  vector<vector<uint> > vinputs;
    //  for(i=0;i<xvindices.size()-1;i++){
    //    index=xvindices[i];
    //    for(j=i+1;j<xvindices.size();j++){
    //      index2=xvindices[j];
    //      vinputs.push_back(vector<uint>(0));
    //      vinputs.back().push_back(index);
    //      vinputs.back().push_back(index2);
    //    }
    //  }
    //  cerr << vinputs.size() << endl;
    //  exit(0);
    //}
    
    //std::sort(vicol2remove.begin(),vicol2remove.end());vicol2remove.erase( std::unique( vicol2remove.begin(), vicol2remove.end() ), vicol2remove.end() );  //get unique values
    //message << "removing " << vicol2remove.size() << " null / low-variance / low-correlated / self-correlated columns";
    message << "removing " << aurostd::sum(xvicol2remove) << " null / low-variance / low-correlated / self-correlated columns";
    pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,cout,_LOGGER_NOTICE_);

    //[ERASE IS TOO SLOW]if(0){  //erase is TOO slow
    //[ERASE IS TOO SLOW]  //sort and reverse, so we can erase/insert
    //[ERASE IS TOO SLOW]  std::sort(vicol2remove.rbegin(),vicol2remove.rend());
    //[ERASE IS TOO SLOW]
    //[ERASE IS TOO SLOW]  for(i=0;i<vicol2remove.size();i++){
    //[ERASE IS TOO SLOW]    if(LDEBUG){cerr << soliloquy << " removing column " << table[0][vicol2remove[i]] << endl;}
    //[ERASE IS TOO SLOW]    delColumn(table,vicol2remove[i]);
    //[ERASE IS TOO SLOW]  }
    //[ERASE IS TOO SLOW]}
    
    //two options
    //flip row and column
    //get list of columns to keep
    
    //vector<uint> vicol2keep;
    //for(i=0;i<ncols;i++){
    //  if(aurostd::WithinList(vicol2remove,i,true)){continue;}
    //  vicol2keep.push_back(i);
    //}

    message << "adding back Z_cation and Mendeleev_number_cation";
    pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,cout,_LOGGER_MESSAGE_);
    for(i=0;i<ncols;i++){
      if( table[0][i]=="Z_cation" || table[0][i]=="Mendeleev_number_cation" ){  //keep these
        xvicol2remove[xvicol2remove.lrows+i]=0;
      }
    }

    vector<vector<string> > table_new;
    for(i=0;i<table.size();i++){table_new.push_back(vector<string>(ncols-aurostd::sum(xvicol2remove)));}
    index=0;
    for(i=0;i<ncols;i++){
      //if(!aurostd::WithinList(vicol2keep,i)){continue;}
      if(xvicol2remove[xvicol2remove.lrows+i]==1){continue;}
      for(j=0;j<table.size();j++){
        table_new[j][index]=table[j][i];
      }
      index++;
    }
    if(LDEBUG){cerr << soliloquy << " done creating table_new" << endl;}

    //for(i=0;i<table.size();i++){
    //  table_new.push_back(vector<string>(0));
    //  for(j=0;j<table[i].size();j++){
    //    if(!aurostd::WithinList(vicol2keep,j)){continue;}
    //    table_new.back().push_back(table[i][j]);
    //  }
    //}
    table.clear();table=table_new;  //overwrite, should be faster than erase

    //check we didn't mess up
    ncols=table[0].size();
    //if((ncols_orig-ncols)!=vicol2remove.size()){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"(ncols_orig-ncols)!=vicol2remove.size()",_RUNTIME_ERROR_);}
    if((int)(ncols_orig-ncols)!=aurostd::sum(xvicol2remove)){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"(ncols_orig-ncols)!=aurostd::sum(xvicol2remove)",_RUNTIME_ERROR_);}
    for(i=0;i<table.size();i++){
      if(table[i].size()!=ncols){
        throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"table[i="+aurostd::utype2string(i)+"].size()!=ncols",_RUNTIME_ERROR_);
      }
    }
    
    message << "ncols_new=" << ncols;
    pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,cout,_LOGGER_MESSAGE_);
    
    if(1||LDEBUG){
      cerr << soliloquy << " table_new=" << endl;
      for(i=0;i<table.size();i++){cerr << aurostd::joinWDelimiter(table[i],",") << endl;}
    }
    
  }

  string reduceEProperties(double var_threshold,double selfcorr_threshold) {
    bool LDEBUG=(TRUE || _DEBUG_ML_ || XHOST.DEBUG);
    string soliloquy=XPID+"aflowML::reduceEProperties():";
    stringstream message;

    uint i=0,j=0,ielement=0,index=0,index2=0;
    xelement::xelement xel("N");  //dummy to start

    vector<string> vproperties_elements_full;aurostd::string2tokens(_AFLOW_XELEMENT_PROPERTIES_ALL_,vproperties_elements_full,",");
    vector<string> vproperties_elements_numbers;
    //number properties next
    for(j=0;j<vproperties_elements_full.size();j++){
      if(xel.getType(vproperties_elements_full[j])=="number"){  //ignore type==="numbers", too much work to rewrite other methods to eliminate single component
        if(vproperties_elements_full[j]!="oxidation_states" && vproperties_elements_full[j]!="oxidation_states_preferred"){ //exclude these entirely
          vproperties_elements_numbers.push_back(vproperties_elements_full[j]);
        }
      }
    }

    if(LDEBUG){cerr << soliloquy << " vproperties_elements_numbers=" << aurostd::joinWDelimiter(vproperties_elements_numbers,",") << endl;}
    
    if(!_TM_ONLY_){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"going beyond transition metals requires rewriting this function",_FILE_CORRUPT_);}

    vector<vector<string> > vvtable;
    for(ielement=0;ielement<100;ielement++){
      xel.populate(ielement); //,ioxidation
      xel.convertUnits();
      if(_TM_ONLY_){
        if(!(xel.group>=3 && xel.group<=12 && xel.period>=4 && xel.period<=6)){continue;}
      }
      vvtable.push_back(vector<string>(0));
      insertElementalProperties(vproperties_elements_numbers,xel,vvtable.back());
      if(vvtable.back().size()!=vproperties_elements_numbers.size()){
        throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"vvtable.back().size()["+aurostd::utype2string(vvtable.back().size())+"]!=vproperties_elements_numbers.size()["+aurostd::utype2string(vproperties_elements_numbers.size())+"]",_FILE_CORRUPT_);
      }
    }

    vector<xvector<double> > vcols;
    for(i=0;i<vproperties_elements_numbers.size();i++){
      vcols.push_back(xvector<double>(vvtable.size()));
    }

    for(i=0;i<vvtable.size();i++){
      for(j=0;j<vvtable[i].size();j++){
        vcols[j][vcols[j].lrows+i]=aurostd::string2utype<double>(vvtable[i][j]);
      }
    }

    if(LDEBUG){
      for(i=0;i<vcols.size();i++){
        cerr << soliloquy << " vcols[" << vproperties_elements_numbers[i] << "]=" << vcols[i] << endl;
      }
    }

    vector<uint> vicol2remove;
    xvector<double> xv_clean;
    vector<double> vmeans,vstddevs;
    double var=0.0;
    for(i=0;i<vcols.size();i++){
      index=i;
      const string& header=vproperties_elements_numbers[index];
      if(LDEBUG){cerr << soliloquy << " looking at column " << header << endl;}
      const xvector<double>& xvec=vcols[index];
      vmeans.push_back(aurostd::mean(xvec));vstddevs.push_back(aurostd::stddev(xvec));
      if(LDEBUG){cerr << soliloquy << " xvec(orig  )=" << xvec << endl;}
      removeNaN(xvec,xv_clean);
      if(LDEBUG){cerr << soliloquy << " xvec(clean )=" << xv_clean << endl;}
      if(xv_clean.rows==0||xv_clean.rows==1){
        if(LDEBUG){cerr << soliloquy << " removing " << header << ": " << (xv_clean.rows==0?"null-":"1-") << "vector" << endl;}
        vicol2remove.push_back(index);
        continue;
      }
      if(xv_clean.rows<(xvec.rows/2)){
        if(LDEBUG){cerr << soliloquy << " removing " << header << ": " << " vector more than half filled with NaNs (" << xv_clean.rows << " out of " << xvec.rows << ")" << endl;}
        vicol2remove.push_back(index);
        continue;
      }
      MinMaxScale(xv_clean);
      if(LDEBUG){cerr << soliloquy << " xvec(scaled)=" << xv_clean << endl;}
      var=aurostd::var(xv_clean,1); //sample variance
      if(LDEBUG){cerr << soliloquy << " var(xvec)=" << var << endl;}
      if(var<var_threshold){
        if(LDEBUG){cerr << soliloquy << " no variance in column " << header << endl;}
        vicol2remove.push_back(index);
        continue;
      }
    }
    
    std::sort(vicol2remove.begin(),vicol2remove.end());vicol2remove.erase( std::unique( vicol2remove.begin(), vicol2remove.end() ), vicol2remove.end() );  //get unique values
    message << "found " << vicol2remove.size() << " null / low-variance columns";
    pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,cout,_LOGGER_MESSAGE_);

    double corr=0.0;
    for(i=0;i<vcols.size();i++){
      index=i;
      if(1||LDEBUG){cerr << soliloquy << " index=" << index << endl;}
      const string& header=vproperties_elements_numbers[index];
      if(aurostd::WithinList(vicol2remove,index)){
        if(LDEBUG){cerr << soliloquy << " skipping " << header << ": to be removed" << endl;}
        continue;
      }
      const xvector<double>& xvec=vcols[index]; //no need to remove NNN, correlation would be the same with 0 or NNN
      if(LDEBUG){cerr << soliloquy << " xvec[" << header << "]=" << xvec << endl;}
      
      for(j=i+1;j<vcols.size();j++){
        index2=j;
        if(LDEBUG){cerr << soliloquy << " index2=" << index2 << endl;}
        const string& header2=vproperties_elements_numbers[index2];
        if(aurostd::WithinList(vicol2remove,index2)){
          if(LDEBUG){cerr << soliloquy << " skipping " << header2 << ": to be removed" << endl;}
          continue;
        }
        const xvector<double>& xvec2=vcols[index2]; //no need to remove NNN, correlation would be the same with 0 or NNN
        if(LDEBUG){cerr << soliloquy << " xvec2[" << header2 << "]=" << xvec2 << endl;}

        corr=std::pow(aurostd::correlation_Pearson_fast(xvec,vmeans[index],vstddevs[index],xvec2,vmeans[index2],vstddevs[index2]),2.0);
        if(LDEBUG){cerr << soliloquy << " corr(xvec,xvec2)^2=" << corr << endl;}
        if(corr>selfcorr_threshold){
          if(LDEBUG){cerr << soliloquy << " high correlation between columns " << header << " and " <<  header2 << endl;}
          vicol2remove.push_back(index2);
        }
      }
    }
    
    std::sort(vicol2remove.begin(),vicol2remove.end());vicol2remove.erase( std::unique( vicol2remove.begin(), vicol2remove.end() ), vicol2remove.end() );  //get unique values
    message << "found " << vicol2remove.size() << " null / low-variance / self-correlated columns";
    pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,cout,_LOGGER_MESSAGE_);

    vector<string> vcol2remove;
    for(i=0;i<vicol2remove.size();i++){vcol2remove.push_back(vproperties_elements_numbers[vicol2remove[i]]);}

    vector<string> vproperties_elements_full_new;
    for(i=0;i<vproperties_elements_full.size();i++){
      if(aurostd::WithinList(vcol2remove,vproperties_elements_full[i])){continue;}
      vproperties_elements_full_new.push_back(vproperties_elements_full[i]);
    }
    
    if(LDEBUG){
      cerr << soliloquy << " reduced " << vproperties_elements_full.size() << " elemental properties down to " << vproperties_elements_full_new.size() << " (reduced by " << (vproperties_elements_full.size()-vproperties_elements_full_new.size()) << ")" << endl;
      cerr << soliloquy << " vproperties_elements_full_new=" << aurostd::joinWDelimiter(vproperties_elements_full_new,",") << endl;
    }

    return aurostd::joinWDelimiter(vproperties_elements_full_new,",");
  }
  void writeCCECSV() {
    bool LDEBUG=(TRUE || _DEBUG_ML_ || XHOST.DEBUG);
    string soliloquy=XPID+"aflowML::writeCCECSV():";

    if(LDEBUG){cerr << soliloquy << " BEGIN" << endl;}

    uint ielement=0,ioxidation=0,i=0,j=0,k=0,ipp=0;
    string correction_line="",bader_line="",input_pre="",input="";
    xelement::xelement xel_cation,xel_N("N"),xel_O("O");
    xel_N.convertUnits();xel_O.convertUnits();
    //cce::get_corrections_line_O("Al_+3_N");
    vector<vector<vector<string> > > vvlines;
    vector<string> vitems_pre,vitems,vtokens;
    vvlines.resize(2); //N,O

    string eproperties_full=reduceEProperties();
    
    vector<string> vheaders,vheaders_additional;
    vector<double> vfeatures;
    vector<string> vions;aurostd::string2tokens("cation,anion",vions,",");
    vector<string> vproperties_elements;aurostd::string2tokens("symbol",vproperties_elements,","); //string properties first
    vector<string> vproperties_elements_full;aurostd::string2tokens(eproperties_full,vproperties_elements_full,","); //_AFLOW_XELEMENT_PROPERTIES_ALL_
    vector<string> vproperties_elements_numbers;
    string type="",units="";

    //number properties next
    for(j=0;j<vproperties_elements_full.size();j++){
      type=xel_N.getType(vproperties_elements_full[j]);
      if(type=="number"||type=="numbers"){
        if(vproperties_elements_full[j]!="oxidation_states" && vproperties_elements_full[j]!="oxidation_states_preferred"){ //exclude these entirely
          vproperties_elements_numbers.push_back(vproperties_elements_full[j]);
        }
      }
    }
    vproperties_elements.insert(vproperties_elements.end(),vproperties_elements_numbers.begin(),vproperties_elements_numbers.end());  //insert numbers properties
    //
    vector<string> vlattice_constants_variants;aurostd::string2tokens("a,b,c",vlattice_constants_variants,",");
    vector<string> vlattice_angles_variants;aurostd::string2tokens("alpha,beta,gamma",vlattice_angles_variants,",");
    for(i=0;i<vions.size();i++){
      for(j=0;j<vproperties_elements.size();j++){
        if(vproperties_elements[j]=="lattice_constants"){
          for(k=0;k<vlattice_constants_variants.size();k++){
            vheaders.push_back(vproperties_elements[j]+"_"+vlattice_constants_variants[k]+"_"+vions[i]);
          }
        }
        else if(vproperties_elements[j]=="lattice_angles"){
          for(k=0;k<vlattice_angles_variants.size();k++){
            vheaders.push_back(vproperties_elements[j]+"_"+vlattice_angles_variants[k]+"_"+vions[i]);
          }
        }
        else if(vproperties_elements[j]=="energies_ionization"){
          for(k=0;k<_ENERGIES_IONIZATION_MAX_;k++){vheaders.push_back(vproperties_elements[j]+"_"+aurostd::utype2string(k+1)+"_"+vions[i]);}
        }
        else{vheaders.push_back(vproperties_elements[j]+"_"+vions[i]);}
      }
      vheaders.push_back("oxidation_"+vions[i]);
      vheaders.push_back("energy_groundstate_PBE_"+vions[i]);
      vheaders.push_back("EATOM_PBE_"+vions[i]);
    }
    //
    vheaders.push_back("charge_bader_cation_PBE");
    //
    vheaders.push_back("PBE_298.15K");
    vheaders.push_back("PBE_0K");
    vheaders.push_back("LDA_298.15K");
    vheaders.push_back("LDA_0K");
    vheaders.push_back("SCAN_298.15K");
    vheaders.push_back("SCAN_0K");
    vheaders.push_back("PBE+U_298.15K");
    vheaders.push_back("PBE+U_0K");
    vheaders.push_back("exp_298.15K");
    vheaders.push_back("M-X_bonds");
    vheaders.push_back("natoms_per_f.u._cation");
    vheaders.push_back("natoms_per_f.u._anion");
    vheaders.push_back("enthalpy_formation_atom_exp");
    //
    //cce-like
    string tmp_str="";
    for(i=0;i<vions.size();i++){
      for(j=0;j<vproperties_elements.size();j++){
        units=xel_N.getUnits(vproperties_elements[j]);
        if(units=="J/mol" && vproperties_elements[j]=="energies_ionization"){
          for(k=0;k<_ENERGIES_IONIZATION_MAX_;k++){vheaders.push_back(vproperties_elements[j]+"_"+aurostd::utype2string(k+1)+"_"+vions[i]+"_per_bond");}
        }
        else if(units=="J/mol"||units=="J"){
          vheaders.push_back(vproperties_elements[j]+"_"+vions[i]+"_per_bond");
        }
        else if(units=="K"){
          tmp_str=vproperties_elements[j];
          aurostd::StringSubst(tmp_str,"temperature","energy");
          vheaders.push_back(tmp_str+"_"+vions[i]+"_per_bond");
        }
      }
    }
    //
    //
    aflowlib::_aflowlib_entry entry;
    entry.getStoichFeatures(vheaders_additional,eproperties_full); //just get headers //_AFLOW_XELEMENT_PROPERTIES_ALL_
    vheaders.insert(vheaders.end(),vheaders_additional.begin(),vheaders_additional.end());
    //
    vheaders.push_back("geometry_a_crystal");
    vheaders.push_back("geometry_b_crystal");
    vheaders.push_back("geometry_c_crystal");
    vheaders.push_back("geometry_alpha_crystal");
    vheaders.push_back("geometry_beta_crystal");
    vheaders.push_back("geometry_gamma_crystal");
    vheaders.push_back("geometry_a_/_b_crystal");
    vheaders.push_back("geometry_b_/_c_crystal");
    vheaders.push_back("geometry_a_/_c_crystal");
    vheaders.push_back("natoms_crystal");
    vheaders.push_back("natoms_cation_crystal");
    vheaders.push_back("natoms_anion_crystal");
    vheaders.push_back("stoich_cation_crystal");
    vheaders.push_back("stoich_anion_crystal");
    vheaders.push_back("density_crystal");
    vheaders.push_back("ZVAL_cation_crystal");  //not Z
    vheaders.push_back("ZVAL_anion_crystal"); //not Z
    vheaders.push_back("ZVAL_iupac_crystal");  //IUPAC
    vheaders.push_back("ZVAL_std_crystal");
    vheaders.push_back("volume_cell_crystal");
    vheaders.push_back("volume_atom_crystal");
    vheaders.push_back("energy_cell_crystal");
    vheaders.push_back("energy_atom_crystal");
    vheaders.push_back("energy_cutoff_crystal");
    vheaders.push_back("enthalpy_formation_cell_crystal");
    vheaders.push_back("enthalpy_formation_atom_crystal");
    vheaders.push_back("distance_cation_cation_crystal");
    vheaders.push_back("distance_cation_anion_crystal");
    vheaders.push_back("distance_anion_anion_crystal");
    vheaders.push_back("spacegroup_crystal");
    vheaders.push_back("point_group_order_crystal");
    vheaders.push_back("Bravais_lattice_crystal");

    uint vheaders_size=vheaders.size();

    //create environmental features

    //single features
    insertElementalCombinations(vproperties_elements_numbers,vheaders_additional);
    vheaders.insert(vheaders.end(),vheaders_additional.begin(),vheaders_additional.end());
    uint count_vcols_ecombocce=vheaders_additional.size();

    if(LDEBUG){
      for(i=0;i<vheaders.size();i++){cerr << soliloquy << " vheaders[i=" << i << "]=" << vheaders[i] << endl;}
    }

    for(i=0;i<vvlines.size();i++){vvlines[i].push_back(vheaders);}
    
    string species_pp="";
    bool found_pp=false;
    string structure_path="";
    double M_X_bonds=0.0;
    double natoms_per_fu_cation=0.0,natoms_per_fu_anion=0.0;
    for(ielement=0;ielement<100;ielement++){
      for(ioxidation=0;ioxidation<10;ioxidation++){
        xel_cation.populate(ielement,ioxidation);
        xel_cation.convertUnits();
        //transition metals only
        //df[(df["group_cation"]>=3) & (df["group_cation"]<=12) & (df["period_cation"]>=4) & (df["period_cation"]<=6)]
        if(_TM_ONLY_){
          if(!(xel_cation.group>=3 && xel_cation.group<=12 && xel_cation.period>=4 && xel_cation.period<=6)){continue;}
        }
        input_pre=xel_cation.symbol;
        input_pre+="_+"+aurostd::utype2string(ioxidation);
        if(LDEBUG){cerr << soliloquy << " input=" << input_pre << endl;}
        vitems_pre.clear();
        //cation
        insertElementalProperties(vproperties_elements,xel_cation,vitems_pre);
        //
        vitems_pre.push_back(aurostd::utype2string(ioxidation,_DOUBLE_WRITE_PRECISION_)); //ioxidation
        try{species_pp=AVASP_Get_PseudoPotential_PAW_PBE(xel_cation.symbol);}
        catch(aurostd::xerror& excpt){continue;}
        found_pp=false;
        for(ipp=0;ipp<vxpseudopotential.size()&&found_pp==false;ipp++) {
          if(vxpseudopotential[ipp].species_pp_type[0]=="PAW_PBE" && vxpseudopotential[ipp].species_pp[0]==species_pp){
            found_pp=true;
            //vitems_pre.push_back(aurostd::utype2string(vxpseudopotential[ipp].species_Z[0],_DOUBLE_WRITE_PRECISION_)); //Z
            vitems_pre.push_back(aurostd::utype2string(vxpseudopotential[ipp].species_pp_groundstate_energy[0],_DOUBLE_WRITE_PRECISION_)); //groundstate_energy
            vitems_pre.push_back(aurostd::utype2string(vxpseudopotential[ipp].vEATOM[0],_DOUBLE_WRITE_PRECISION_)); //EATOM
          }
        }

        //N
        input=input_pre+"_N";
        correction_line=cce::get_corrections_line_N(input);
        if(!correction_line.empty()){
          vitems.clear();vitems.reserve(vheaders_size);for(i=0;i<vitems_pre.size();i++){vitems.push_back(vitems_pre[i]);}
          //anion
          insertElementalProperties(vproperties_elements,xel_N,vitems);
          //
          vitems.push_back(aurostd::utype2string(-3,_DOUBLE_WRITE_PRECISION_)); //ioxidation - fix later if needed
          try{species_pp=AVASP_Get_PseudoPotential_PAW_PBE(xel_N.symbol);}
          catch(aurostd::xerror& excpt){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"cannot run AVASP_Get_PseudoPotential_PAW_PBE() for nitrogen",_FILE_CORRUPT_);}
          found_pp=false;
          for(ipp=0;ipp<vxpseudopotential.size()&&found_pp==false;ipp++) {
            if(vxpseudopotential[ipp].species_pp_type[0]=="PAW_PBE" && vxpseudopotential[ipp].species_pp[0]==species_pp){
              found_pp=true;
              if(vxpseudopotential[ipp].species_pp_groundstate_structure[0]!="diatom"){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"pp_groundstate_structure for "+xel_N.symbol+" is NOT diatom",_INPUT_ILLEGAL_);}
              //vitems.push_back(aurostd::utype2string(vxpseudopotential[ipp].species_Z[0],_DOUBLE_WRITE_PRECISION_)); //Z
              vitems.push_back(aurostd::utype2string(vxpseudopotential[ipp].species_pp_groundstate_energy[0],_DOUBLE_WRITE_PRECISION_)); //groundstate_energy
              vitems.push_back(aurostd::utype2string(vxpseudopotential[ipp].vEATOM[0],_DOUBLE_WRITE_PRECISION_)); //EATOM
            }
          }
          //
          vitems.push_back("0.0");  //no bader yet

          //
          if(LDEBUG){cerr << soliloquy << " correction_line=\"" << correction_line << "\"" << endl;}
          aurostd::string2tokens(correction_line,vtokens," ");
          if(vtokens.size()<16){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"vtokens.size()<16",_FILE_CORRUPT_);}
          for(i=2;i<=11;i++){vitems.push_back(vtokens[i]);}
          M_X_bonds=aurostd::string2utype<double>(vtokens[11]);
          if(LDEBUG){cerr << soliloquy << " M-X_bonds=" << M_X_bonds << endl;}
          //
          structure_path=vtokens[12];
          if(LDEBUG){cerr << soliloquy << " structure_path=" << structure_path << endl;}
          //
          vitems.push_back(vtokens[13]);  //ncations_per_f.u.
          vitems.push_back(vtokens[14]);  //nanions_per_f.u.
          vitems.push_back(vtokens[15]);  //H_f^exp_298.15K
          //
          natoms_per_fu_cation=aurostd::string2utype<double>(vtokens[13]);
          natoms_per_fu_anion=aurostd::string2utype<double>(vtokens[14]);
          //
          insertElementalPropertiesCCE(vproperties_elements,xel_cation,M_X_bonds,natoms_per_fu_cation,vitems);
          insertElementalPropertiesCCE(vproperties_elements,xel_N,M_X_bonds,natoms_per_fu_anion,vitems);
          //
          insertCrystalProperties(structure_path,"N",vheaders,vitems,entry,eproperties_full);
          //
          insertElementalCombinations(vproperties_elements_numbers,xel_cation,xel_N,entry,M_X_bonds,natoms_per_fu_cation,natoms_per_fu_anion,vheaders_additional,vfeatures,false,count_vcols_ecombocce);
          for(i=0;i<vfeatures.size();i++){vitems.push_back(aurostd::utype2string(vfeatures[i],_DOUBLE_WRITE_PRECISION_));}
          //
          if(vitems.size()!=vheaders.size()){
            for(uint ii=0;ii<vheaders.size()&&ii<vitems.size();ii++){
              cerr << vheaders[ii] << "=" << vitems[ii] << endl;
            }
            throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"vitems.size()["+aurostd::utype2string(vitems.size())+"]!=vheaders.size()["+aurostd::utype2string(vheaders.size())+"]",_FILE_CORRUPT_);
          }
          vvlines[0].push_back(vitems);
        }

        //O
        input=input_pre+"_O";
        correction_line=cce::get_corrections_line_O(input);
        if(!correction_line.empty()){
          vitems.clear();vitems.reserve(vheaders_size);for(i=0;i<vitems_pre.size();i++){vitems.push_back(vitems_pre[i]);}
          //anion
          insertElementalProperties(vproperties_elements,xel_O,vitems);
          //
          vitems.push_back(aurostd::utype2string(-2,_DOUBLE_WRITE_PRECISION_)); //ioxidation - fix later if needed
          try{species_pp=AVASP_Get_PseudoPotential_PAW_PBE(xel_O.symbol);}
          catch(aurostd::xerror& excpt){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"cannot run AVASP_Get_PseudoPotential_PAW_PBE() for oxygen",_FILE_CORRUPT_);}
          found_pp=false;
          for(ipp=0;ipp<vxpseudopotential.size()&&found_pp==false;ipp++) {
            if(vxpseudopotential[ipp].species_pp_type[0]=="PAW_PBE" && vxpseudopotential[ipp].species_pp[0]==species_pp){
              found_pp=true;
              if(vxpseudopotential[ipp].species_pp_groundstate_structure[0]!="diatom"){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"pp_groundstate_structure for "+xel_O.symbol+" is NOT diatom",_INPUT_ILLEGAL_);}
              //vitems.push_back(aurostd::utype2string(vxpseudopotential[ipp].species_Z[0],_DOUBLE_WRITE_PRECISION_)); //Z
              vitems.push_back(aurostd::utype2string(vxpseudopotential[ipp].species_pp_groundstate_energy[0],_DOUBLE_WRITE_PRECISION_)); //groundstate_energy
              vitems.push_back(aurostd::utype2string(vxpseudopotential[ipp].vEATOM[0],_DOUBLE_WRITE_PRECISION_)); //EATOM
            }
          }
          //
          bader_line=cce::get_Bader_templates(xel_cation.symbol);
          if(!bader_line.empty()){
            if(LDEBUG){cerr << soliloquy << " bader_line=\"" << bader_line << "\"" << endl;}
            aurostd::string2tokens(bader_line,vtokens," ");
            if(vtokens.size()<6){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"vtokens.size()<13",_FILE_CORRUPT_);}
            vitems.push_back(vtokens[2]);
          }

          //
          if(LDEBUG){cerr << soliloquy << " correction_line=\"" << correction_line << "\"" << endl;}
          aurostd::string2tokens(correction_line,vtokens," ");
          if(vtokens.size()<16){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"vtokens.size()<16",_FILE_CORRUPT_);}
          for(i=2;i<=11;i++){vitems.push_back(vtokens[i]);}
          M_X_bonds=aurostd::string2utype<double>(vtokens[11]);
          if(LDEBUG){cerr << soliloquy << " M-X_bonds=" << M_X_bonds << endl;}
          //
          structure_path=vtokens[12];
          if(LDEBUG){cerr << soliloquy << " structure_path=" << structure_path << endl;}
          //
          vitems.push_back(vtokens[13]);  //ncations_per_f.u.
          vitems.push_back(vtokens[14]);  //nanions_per_f.u.
          vitems.push_back(vtokens[15]);  //H_f^exp_298.15K
          //
          natoms_per_fu_cation=aurostd::string2utype<double>(vtokens[13]);
          natoms_per_fu_anion=aurostd::string2utype<double>(vtokens[14]);
          //
          insertElementalPropertiesCCE(vproperties_elements,xel_cation,M_X_bonds,natoms_per_fu_cation,vitems);
          insertElementalPropertiesCCE(vproperties_elements,xel_O,M_X_bonds,natoms_per_fu_anion,vitems);
          //
          insertCrystalProperties(structure_path,"O",vheaders,vitems,entry,eproperties_full);
          //
          insertElementalCombinations(vproperties_elements_numbers,xel_cation,xel_O,entry,M_X_bonds,natoms_per_fu_cation,natoms_per_fu_anion,vheaders_additional,vfeatures,false,count_vcols_ecombocce);
          for(i=0;i<vfeatures.size();i++){vitems.push_back(aurostd::utype2string(vfeatures[i],_DOUBLE_WRITE_PRECISION_));}
          //
          if(vitems.size()!=vheaders.size()){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"vitems.size()!=vheaders.size()",_FILE_CORRUPT_);}
          vvlines[1].push_back(vitems);
        }
      }
    }
    
    for(i=0;i<vvlines.size();i++){
      oneHotFeatures(vvlines[i],"Bravais_lattice_crystal,spacegroup_crystal,spacegroup_number_cation,spacegroup_number_anion");
      reduceFeatures(vvlines[i],"PBE_0K","PBE_298.15K,PBE_0K,LDA_298.15K,LDA_0K,SCAN_298.15K,SCAN_0K,PBE+U_298.15K,PBE+U_0K,exp_298.15K"); //skip y
    }

    stringstream file;
    //N
    aurostd::StringstreamClean(file);
    for(i=0;i<vvlines[0].size();i++){file << aurostd::joinWDelimiter(vvlines[0][i],",") << endl;}
    aurostd::stringstream2file(file,"cce_data_N.csv");
    //O
    aurostd::StringstreamClean(file);
    for(i=0;i<vvlines[1].size();i++){file << aurostd::joinWDelimiter(vvlines[1][i],",") << endl;}
    aurostd::stringstream2file(file,"cce_data_O.csv");
    //total
    if(0){  //we've removed anion properties in reduceFeatures()
      aurostd::StringstreamClean(file);
      for(i=0;i<vvlines[0].size();i++){file << aurostd::joinWDelimiter(vvlines[0][i],",") << endl;}
      for(i=1;i<vvlines[1].size();i++){file << aurostd::joinWDelimiter(vvlines[1][i],",") << endl;} //skip header
      aurostd::stringstream2file(file,"cce_data_NO.csv");
    }

  }
} // namespace aflowML

namespace aflowlib {
  void insertStoichStats(const vector<string> vstats,const xvector<double>& nspecies_xv,const xvector<double>& stoich_xv,vector<double>& vfeatures){
    bool LDEBUG=(FALSE || _DEBUG_ML_ || XHOST.DEBUG);
    string soliloquy=XPID+"aflowlib::insertStoichStats():";

    uint k=0,l=0;
    int index=0;
    double d_tmp;
    vector<uint> vi_tmp;

    //check for NNN or AUROSTD_NAN
    bool has_NaN=false;
    for(index=nspecies_xv.lrows;index<=nspecies_xv.urows&&!has_NaN;index++){
      if(aurostd::isNaN(nspecies_xv[index])){has_NaN=true;}
    }

    if(LDEBUG){cerr << soliloquy << " nspecies_xv=" << nspecies_xv << endl;}
    for(k=0;k<vstats.size();k++){
      if(vstats[k]=="min"){if(has_NaN){vfeatures.push_back(NNN);continue;} vfeatures.push_back( aurostd::min(nspecies_xv) );}
      else if(vstats[k]=="max"){if(has_NaN){vfeatures.push_back(NNN);continue;} vfeatures.push_back( aurostd::max(nspecies_xv) );}
      else if(vstats[k]=="range"){if(has_NaN){vfeatures.push_back(NNN);continue;} vfeatures.push_back( aurostd::max(nspecies_xv) - aurostd::min(nspecies_xv) );}
      else if(vstats[k]=="mean"){if(has_NaN){vfeatures.push_back(NNN);continue;} vfeatures.push_back( aurostd::scalar_product(stoich_xv,nspecies_xv) );}
      else if(vstats[k]=="dev"){
        if(has_NaN){vfeatures.push_back(NNN);continue;}
        d_tmp=aurostd::scalar_product(stoich_xv,nspecies_xv);  //mean
        vfeatures.push_back( aurostd::scalar_product(stoich_xv,aurostd::abs(nspecies_xv-d_tmp)) );
      }
      else if(vstats[k]=="mode"){ //property of most promiment species
        if(has_NaN){vfeatures.push_back(NNN);continue;}
        d_tmp=aurostd::max(stoich_xv);  //stoich_max
        vi_tmp.clear();
        for(index=stoich_xv.lrows;index<=stoich_xv.urows;index++){
          if(aurostd::isequal(stoich_xv[index],d_tmp)){vi_tmp.push_back(index);}
        }
        if(vi_tmp.size()==1){vfeatures.push_back( nspecies_xv[vi_tmp[0]] );}  //easy case
        else{
          //take average
          d_tmp=0;
          for(l=0;l<vi_tmp.size();l++){d_tmp+=nspecies_xv[vi_tmp[l]];}
          d_tmp/=(double)vi_tmp.size();
          vfeatures.push_back( d_tmp );
        }
      }
      else{throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"Unknown statistic type: "+vstats[k],_RUNTIME_ERROR_);}
    }
  }

  void _aflowlib_entry::getStoichFeatures(vector<string>& vheaders,const string& e_props){
    vector<double> vfeatures; //dummy
    return getStoichFeatures(vheaders,vfeatures,true,e_props);
  }
  void _aflowlib_entry::getStoichFeatures(vector<string>& vheaders,vector<double>& vfeatures,bool vheaders_only,const string& e_props){
    //follows supplementary of 10.1038/npjcompumats.2016.28
    bool LDEBUG=(FALSE || _DEBUG_ML_ || XHOST.DEBUG);
    string soliloquy=XPID+"_aflowlib_entry::getStoichFeatures():";
    vheaders.clear();vfeatures.clear();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //headers
    stringstream tmp_ss;
    xelement::xelement xel;
    vector<xelement::xelement> vxel;
    uint i=0,j=0,k=0;

    //L^p norms
    vector<uint> vp;
    aurostd::string2tokens("0,2,3,4,5,6,7,8,9,10",vp,",");
    for(i=0;i<vp.size();i++){
      aurostd::StringstreamClean(tmp_ss);
      tmp_ss << "stoich_norm_p_" << std::setfill('0') << std::setw(2) << vp[i];
      vheaders.push_back(tmp_ss.str());
    }

    //element-property-based
    //get which properties to average
    xel.populate(1);  //dummy to get properties
    vector<string> vproperties_full,vproperties;
    vector<string> vstats;
    aurostd::string2tokens(e_props,vproperties_full,",");
    aurostd::string2tokens("min,max,range,mean,dev,mode",vstats,",");
    //load up vheaders
    for(i=0;i<vproperties_full.size();i++){
      if(vproperties_full[i]=="oxidation_states"){continue;} //skip this
      if(xel.getType(vproperties_full[i])=="number"||xel.getType(vproperties_full[i])=="numbers"){
        vproperties.push_back(vproperties_full[i]);
        
        if(xel.getType(vproperties.back())=="number"){
          if(LDEBUG){cerr << soliloquy << " " << vproperties.back() << " is a number" << endl;}
          for(j=0;j<vstats.size();j++){vheaders.push_back(vproperties.back()+"_stoich_"+vstats[j]);}
        }
        else if(xel.getType(vproperties.back())=="numbers"){
          if(LDEBUG){cerr << soliloquy << " " << vproperties.back() << " are numbers" << endl;}
          if(vproperties.back()=="lattice_constants"){
            for(j=0;j<vstats.size();j++){vheaders.push_back(vproperties.back()+"_a_stoich_"+vstats[j]);}
            for(j=0;j<vstats.size();j++){vheaders.push_back(vproperties.back()+"_b_stoich_"+vstats[j]);}
            for(j=0;j<vstats.size();j++){vheaders.push_back(vproperties.back()+"_c_stoich_"+vstats[j]);}
          }
          else if(vproperties.back()=="lattice_angles"){
            for(j=0;j<vstats.size();j++){vheaders.push_back(vproperties.back()+"_alpha_stoich_"+vstats[j]);}
            for(j=0;j<vstats.size();j++){vheaders.push_back(vproperties.back()+"_beta_stoich_"+vstats[j]);}
            for(j=0;j<vstats.size();j++){vheaders.push_back(vproperties.back()+"_gamma_stoich_"+vstats[j]);}
          }
          else if(vproperties.back()=="oxidation_states_preferred"){
            for(j=0;j<vstats.size();j++){
              vheaders.push_back(vproperties.back()+"_stoich_"+vstats[j]);  //only use 0th oxidation_state_preferred
            }
          }
          else if(vproperties.back()=="energies_ionization"){
            for(k=0;k<_ENERGIES_IONIZATION_MAX_;k++){
              for(j=0;j<vstats.size();j++){vheaders.push_back(vproperties.back()+"_"+aurostd::utype2string(k+1)+"_stoich_"+vstats[j]);}
            }
            //for(j=0;j<vstats.size();j++){vheaders.push_back(vproperties.back()+"_2_stoich_"+vstats[j]);}
            //for(j=0;j<vstats.size();j++){vheaders.push_back(vproperties.back()+"_3_stoich_"+vstats[j]);}
            //for(j=0;j<vstats.size();j++){vheaders.push_back(vproperties.back()+"_4_stoich_"+vstats[j]);}
            //for(j=0;j<vstats.size();j++){vheaders.push_back(vproperties.back()+"_5_stoich_"+vstats[j]);}
          }
          else{throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"Unknown numbers type: "+vproperties.back(),_RUNTIME_ERROR_);}
        }
      }
    }

    //valence (un)occupation
    vector<string> vorbitals;
    aurostd::string2tokens("s,p,d,f",vorbitals,",");
    for(i=0;i<vorbitals.size();i++){
      vheaders.push_back("valence_fraction_occupied_"+vorbitals[i]);
      vheaders.push_back("valence_fraction_unoccupied_"+vorbitals[i]);
    }

    //ionic character
    vheaders.push_back("formability_ionic");
    vector<string> vEN;
    for(i=0;i<vproperties.size();i++){
      if(vproperties[i].find("electronegativity")!=string::npos && xel.getUnits(vproperties[i]).empty()){ //must have NO units (goes in exp)
        vEN.push_back(vproperties[i]);
      }
    }
    for(i=0;i<vEN.size();i++){
      vheaders.push_back("character_ionic_"+vEN[i]+"_max");
      vheaders.push_back("character_ionic_"+vEN[i]+"_mean");
    }

    if(LDEBUG){
      for(i=0;i<vheaders.size();i++){cerr << soliloquy << " vheaders[i=" << i << "]=\"" << vheaders[i] << "\"" << endl;}
    }

    if(vheaders_only) return;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //features
    
    uint nspecies=vspecies.size();
    if(nspecies==0){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"nspecies==0",_RUNTIME_ERROR_);}
    xvector<double> nspecies_xv(nspecies);
    
    //L^p norms
    xvector<double> stoich_xv(nspecies);
    for(j=0;j<nspecies;j++){stoich_xv[stoich_xv.lrows+j]=vcomposition[j]/natoms;}
    if(LDEBUG){cerr << soliloquy << " stoich_xv=" << stoich_xv << endl;}
    for(i=0;i<vp.size();i++){
      for(j=0;j<nspecies;j++){
        nspecies_xv[nspecies_xv.lrows+j]=std::pow(stoich_xv[stoich_xv.lrows+j],(double)vp[i]);
      }
      if(LDEBUG){cerr << soliloquy << " nspecies_xv[\"stoich_norm_p_"+aurostd::utype2string(vp[i])+"\"]=" << nspecies_xv << endl;}
      vfeatures.push_back( std::pow(sum(nspecies_xv),(vp[i]==0?1.0:1.0/vp[i])) );
    }
    
    //element-property-based
    //load up vxel
    int index=0,index_min=0,index_max=0;
    for(j=0;j<nspecies;j++){
      vxel.push_back(xelement::xelement(vspecies[j]));
    }
    for(i=0;i<vproperties.size();i++){
      if(xel.getType(vproperties[i])=="number"){
        for(j=0;j<nspecies;j++){
          nspecies_xv[nspecies_xv.lrows+j]=vxel[j].getPropertyDouble(vproperties[i]);
        }
        if(LDEBUG){cerr << soliloquy << " nspecies_xv[\""+vproperties[i]+"\"]=" << nspecies_xv << endl;}
        insertStoichStats(vstats,nspecies_xv,stoich_xv,vfeatures);
      }
      if(xel.getType(vproperties[i])=="numbers"){
        if(vproperties[i]=="lattice_constants"||vproperties[i]=="lattice_angles"){
          index_min=1;index_max=3;
          for(index=index_min;index<=index_max;index++){
            for(j=0;j<nspecies;j++){
              const xvector<double>& xvec=vxel[j].getPropertyXVectorDouble(vproperties[i]);
              nspecies_xv[nspecies_xv.lrows+j]=xvec[index];
            }
            if(LDEBUG){cerr << soliloquy << " nspecies_xv[\""+vproperties[i]+"_index_"+aurostd::utype2string(index)+"\"]=" << nspecies_xv << endl;}
            insertStoichStats(vstats,nspecies_xv,stoich_xv,vfeatures);
          }
        }
        else if(vproperties[i]=="oxidation_states_preferred"||vproperties[i]=="energies_ionization"){
          if(vproperties[i]=="oxidation_states_preferred"){index_min=0;index_max=0;}
          else if(vproperties[i]=="energies_ionization"){index_min=0;index_max=_ENERGIES_IONIZATION_MAX_-1;}
          else{throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"Unknown numbers property (vector): "+vproperties[i],_RUNTIME_ERROR_);}
          for(index=index_min;index<=index_max;index++){
            for(j=0;j<nspecies;j++){
              const vector<double>& vec=vxel[j].getPropertyVectorDouble(vproperties[i]);
              nspecies_xv[nspecies_xv.lrows+j]=(index<(int)vec.size()?vec[index]:NNN);
            }
            if(LDEBUG){cerr << soliloquy << " nspecies_xv[\""+vproperties[i]+"_index="+aurostd::utype2string(index)+"\"]=" << nspecies_xv << endl;}
            insertStoichStats(vstats,nspecies_xv,stoich_xv,vfeatures);
          }
        }
        else{throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"Unknown numbers property: "+vproperties[i],_RUNTIME_ERROR_);}
      }
    }
    
    //valence (un)occupation
    for(j=0;j<nspecies;j++){nspecies_xv[nspecies_xv.lrows+j]=vxel[j].getPropertyDouble("valence_std");}
    double denom=aurostd::scalar_product(stoich_xv,nspecies_xv); //same for all quantities
    vector<double> vval_total;
    aurostd::string2tokens("2,6,10,14",vval_total,",");
    for(i=0;i<vorbitals.size();i++){
      for(j=0;j<nspecies;j++){nspecies_xv[nspecies_xv.lrows+j]=vxel[j].getPropertyDouble("valence_"+vorbitals[i]);}  //populate with orbital occupation
      vfeatures.push_back( aurostd::scalar_product(stoich_xv,nspecies_xv)/denom );  //occupied
      for(j=0;j<nspecies;j++){nspecies_xv[nspecies_xv.lrows+j]=vval_total[i]-nspecies_xv[nspecies_xv.lrows+j];} //has occupied inside already
      vfeatures.push_back( aurostd::scalar_product(stoich_xv,nspecies_xv)/denom );  //unoccupied
    }
    
    //ionic character
    //ionic formability
    bool formability_ionic=false;
    bool has_NaN=false;
    vector<int> nspecies_v;
    aurostd::xcombos xc;
    for(j=0;j<nspecies&&!has_NaN;j++){
      const vector<double>& oxidation_states=vxel[j].getPropertyVectorDouble("oxidation_states");
      if(oxidation_states.size()==0){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"No oxidation states found for element: "+vxel[j].symbol,_RUNTIME_ERROR_);}  //should have NNN
      if(oxidation_states.size()==1 && aurostd::isNaN(oxidation_states[0])){has_NaN=true;}
      nspecies_v.push_back((int)oxidation_states.size());
    }
    if(has_NaN){cerr << soliloquy << " has NaN" << endl;}
    if(!has_NaN){
      xvector<double> natoms_xv(natoms);
      if(LDEBUG){cerr << soliloquy << " oxidation_states_count=" << aurostd::joinWDelimiter(nspecies_v,",") << endl;}
      xc.reset(nspecies_v,'E');
      while(xc.increment()&&!formability_ionic){
        const vector<int>& indices=xc.getCombo();
        if(LDEBUG){cerr << soliloquy << " indices=" << aurostd::joinWDelimiter(indices,",") << endl;}
        i=0;
        for(j=0;j<nspecies&&!has_NaN;j++){
          const vector<double>& oxidation_states=vxel[j].getPropertyVectorDouble("oxidation_states");
          for(k=0;k<vcomposition[j];k++){
            if(aurostd::isNaN(oxidation_states[indices[j]])){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"Found NaN among populated oxidation_states",_RUNTIME_ERROR_);}
            natoms_xv[natoms_xv.lrows+(i++)]=oxidation_states[indices[j]];
          }
        }
        if(LDEBUG){cerr << soliloquy << " natoms_xv[\"oxidation_states\"]=" << natoms_xv << endl;}
        if(aurostd::isequal(aurostd::sum(natoms_xv),0.0)){formability_ionic=true;}
      }
    }
    vfeatures.push_back( (!has_NaN&&formability_ionic?1:0) );
    //character_ionic_max and _mean
    xvector<double> pairs_xv(aurostd::nCk((int)nspecies,2));  //electronegativities
    xvector<double> pairs2_xv(aurostd::nCk((int)nspecies,2)); //stoich
    for(i=0;i<vEN.size();i++){
      if(LDEBUG){cerr << soliloquy << " EN=" << vEN[i] << endl;}
      xc.reset(nspecies,2);
      k=0;
      while(xc.increment()){
        const vector<int>& indices=xc.getIndices();
        if(LDEBUG){cerr << soliloquy << " indices=" << aurostd::joinWDelimiter(indices,",") << endl;}
        pairs_xv[pairs_xv.lrows+k]=1.0-std::exp(-0.25*(std::pow(vxel[indices[0]].getPropertyDouble(vEN[i])-vxel[indices[1]].getPropertyDouble(vEN[i]),2.0)));
        pairs2_xv[pairs2_xv.lrows+k]=stoich_xv[stoich_xv.lrows+indices[0]]*stoich_xv[stoich_xv.lrows+indices[1]];
        k++;
      }
      if(LDEBUG){cerr << soliloquy << " pairs_xv[\""+vEN[i]+"\"]=" << pairs_xv << endl;}
      vfeatures.push_back( aurostd::max(pairs_xv) );
      vfeatures.push_back( aurostd::scalar_product(pairs_xv,pairs2_xv) );
    }

    if(vheaders.size()!=vfeatures.size()){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"vheaders.size()!=vfeatures.size()",_RUNTIME_ERROR_);}

    if(LDEBUG){
      for(i=0;i<vheaders.size();i++){cerr << soliloquy << " " << vheaders[i] << "=" << vfeatures[i] << endl;}
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  }
} // namespace aflowlib

#endif  // _AFLOW_ML_CPP_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow COREY OSES - Duke University 2013-2020                  *
// *                                                                         *
// ***************************************************************************

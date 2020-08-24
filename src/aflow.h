// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *                                                                         *
// ***************************************************************************

#ifndef _AFLOW_H_
#define _AFLOW_H_

// --------------------------------------------------------------------------
// Standard stuff

#include "AUROSTD/aurostd.h"

// #define  _AFLOW_TEMP_PRESERVE_  // to preseve /tmp files for debug

#define NNN   -123456
//[CO20200502 - moved to aurostd.h]#define GCC_VERSION (__GNUC__ * 10000  + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#define _ANRL_NOWEB_ //DX

//ZERO PRECISION DEFINITIONS - TIGHT (DEFAULT) AND LOOSE
#define _ZERO_PRECISION_ 10
#define _ZERO_TOL_ std::pow(10,-_ZERO_PRECISION_) //DX
#define _ZERO_PRECISION_LOOSE_ 3
#define _ZERO_TOL_LOOSE_ std::pow(10,-_ZERO_PRECISION_LOOSE_) //DX
#define _DOUBLE_PRECISION_ 8
#define _DOUBLE_TOL_ std::pow(10,-_DOUBLE_PRECISION_) 
#define _FLOAT_PRECISION_ 6
#define _FLOAT_TOL_ std::pow(10,-_FLOAT_PRECISION_)  //ME20200519 - tolerance for float precision
//PRECISION and TOLERANCE definitions
#define _DOUBLE_WRITE_PRECISION_MAX_ 14  //CO20180509 - used for xstrctures
#define _DOUBLE_WRITE_PRECISION_ 12  //CO20180509 - used in writing doubles in qmvasp
#define _AFLOWLIB_STOICH_PRECISION_ _DOUBLE_PRECISION_ //[CO20200731 - too many different precisions... just use default]9  //CO20200731
#define _AFLOWLIB_DATA_DOUBLE_PREC_ _DOUBLE_PRECISION_ //[CO20200731 - too many different precisions... just use default]6 //CO20200731
#define _AFLOWLIB_DATA_GEOMETRY_PREC_ _DOUBLE_PRECISION_ //[CO20200731 - too many different precisions... just use default]7 //CO20200731
#define _AFLOW_POCC_PRECISION_ _DOUBLE_PRECISION_ //8 //must be less than _DOUBLE_WRITE_PRECISION_MAX_, which is currently set to 14
#define _AFLOW_POCC_ZERO_TOL_ std::pow(10,-_AFLOW_POCC_PRECISION_) 
#define _XPROTO_TOO_CLOSE_ERROR_ 0.60 // was 0.75
#define _XPROTO_ZERO_VOL_ _FLOAT_TOL_  //CO20190218

//ME20190628 BEGIN - moved from CHULL for broader access
// Output formats
enum filetype {   //CO20190629
  //GENERAL FILE TYPES
  txt_ft,         //general plain text
  json_ft,
  csv_ft,
  latex_ft,
  gnuplot_ft,
  jupyter2_ft,    //python 2 jupyter
  jupyter3_ft,    //python 3 jupyter
  //CHULL SPECIFIC
  chull_apool_ft,
  chull_web_ft,
};
//[CO20190629 - obsolete with enum filetype]#define _apool_         'a'  // apool
//[CO20190629 - obsolete with enum filetype]#define _json_          'j'  // standard json
//[CO20190629 - obsolete with enum filetype]#define _pdf_           'p'  // pdf
//[CO20190629 - obsolete with enum filetype]#define _txt_           't'  // plain text
//[CO20190629 - obsolete with enum filetype]#define _web_           'w'  // web json
//[CO20190629 - obsolete with enum filetype]#define _latex_         'l'  // latex
//[CO20190629 - obsolete with enum filetype]#define _gnuplot_       'g'  // gnuplot
//[CO20190629 - obsolete with enum filetype]#define _jupyterthree_  'y'  // jupyter python 3
//[CO20190629 - obsolete with enum filetype]#define _jupytertwo_    'z'  // jupyter python 2

// Vector reduction types
enum vector_reduction_type {   //CO20190629
  frac_vrt,   //reduce to fractions (normalized to 1)
  gcd_vrt,    //reduce by gcd
  no_vrt,     //no reduction
};
//[CO20190629 - obsolete with enum vector_reduction_type]#define _frac_          'f'  //fractional
//[CO20190629 - obsolete with enum vector_reduction_type]#define _gcd_           'g'  //gcd
//[CO20190629 - obsolete with enum vector_reduction_type]#define _none_          'n'  //none
//ME20190628 END

//compound specification is how a compound is specified
//composition (Mn2Pt3) is ORTHOGONAL to pseudopotential string (Mn_pvPt)
//for instance, H1.25 can be a pseudopotential and NOT a composition
enum compound_designation {
  composition_string,
  pp_string,
};

//CO20190712 - see VASP_PseudoPotential_CleanName_InPlace() in aflow_ivasp.cpp
const string CAPITAL_LETTERS_PP_LIST="_GW2"    //CO20190712 - potpaw_LDA/potpaw_LDA.20100505/Li_AE_GW2
",_GW"    //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/As_GW
",_ZORA"  //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Pt_ZORA
",_LDApU" //CO20190712 - potpaw_LDA/potpaw_LDA.20100505/Zn_sv_LDApU
",_AE"    //CO20190712 - potpaw_LDA/potpaw_LDA.20100505/Li_AE_GW2
",_NC2"   //CO20190712 - potpaw_LDA/potpaw_LDA.20100505/As_NC2
",_200eV"
"";

//MESSAGE defaults - CO20200502
#define _AFLOW_MESSAGE_DEFAULTS_ "user,host,pid,time" //tid //CO20200624 - only depends on XHOST (not aflags)

//CO20200731 START
static const string POCC_TITLE_TAG=":POCC_";
static const string POCC_TITLE_TOL_TAG=":TOL_";
static const string ARUN_TITLE_TAG=":ARUN.";
static const string POCC_ARUN_TITLE_TAG=ARUN_TITLE_TAG+"POCC_";
//CO20200731 END

//XRD
#define XRAY_RADIATION_COPPER_Kalpha 1.5418   //Angstroms     //CO20190622

//moved from avasp.cpp for broader access (chull.cpp)
#define SPECIE_TRANSITION_METALS string("Ag,Au,Cd,Co,Cr_pv,Cu_pv,Fe_pv,Hf_pv,Hg,Ir,La,Mn_pv,Mo_pv,Nb_sv,Ni_pv,Os_pv,Pd_pv,Pt,Re_pv,Rh_pv,Ru_pv,Sc_sv,Ta_pv,Tc_pv,Ti_sv,V_sv,W_pv,Y_sv,Zn,Zr_sv")
#define SPECIE_RAW_LIB2U SPECIE_TRANSITION_METALS
#define SPECIE_RAW_LIB2 string("Ag,Al,As,Au,B_h,Ba_sv,Be_sv,Bi_d,Br,Ca_sv,Cd,Cl,Co,Cr_pv,Cu_pv,Fe_pv,Ga_h,Ge_h,Hf_pv,Hg,In_d,Ir,K_sv,La,Li_sv,Mg_pv,Mn_pv,Mo_pv,Na_pv,Nb_sv,Ni_pv,Os_pv,P,Pb_d,Pd_pv,Pt,Re_pv,Rh_pv,Ru_pv,Sb,Sc_sv,Se,Si,Sn,Sr_sv,Ta_pv,Tc_pv,Te,Ti_sv,Tl_d,V_sv,W_pv,Y_sv,Zn,Zr_sv")

#define SPECIE_RAW_LIB3 string("Ag,Al,As,Au,B_h,Ba_sv,Be_sv,Bi_d,Br,Ca_sv,Cd,Cl,Co,Cr_pv,Cu_pv,Fe_pv,Ga_h,Ge_h,Hf_pv,Hg,In_d,Ir,K_sv,La,Li_sv,Mg_pv,Mn_pv,Mo_pv,Na_sv,Nb_sv,Ni_pv,Os_pv,P,Pb_d,Pd_pv,Pt,Re_pv,Rh_pv,Ru_pv,Sb,Sc_sv,Se,Si,Sn,Sr_sv,Ta_pv,Tc_pv,Te,Ti_sv,Tl_d,V_sv,W_pv,Y_sv,Zn,Zr_sv")
//#define SPECIE_RAW_LIB3 string("Ag,Au,Cd,Co,Cr_pv,Cu_pv,Fe_pv,Hf_pv,Hg,Ir,La,Mn_pv,Mo_pv,Nb_sv,Ni_pv,Os_pv,Pd_pv,Pt,Re_pv,Rh_pv,Ru_pv,Sc_sv,Ta_pv,Tc_pv,Ti_sv,V_sv,W_pv,Y_sv,Zn,Zr_sv")
//#define SPECIE_RAW_LIB3 string("Ag,Al,As,Au,B_h,Bi_d,Cd,Co,Cr_pv,Cu_pv,Fe_pv,Ga_h,Ge_h,Hf_pv,Hg,In_d,Ir,La,Mg_pv,Mn_pv,Mo_pv,Nb_sv,Ni_pv,Os_pv,P,Pb_d,Pd_pv,Pt,Re_pv,Rh_pv,Ru_pv,Sb,Sc_sv,Se,Si,Sn,Ta_pv,Te,Tc_pv,Ti_sv,V_sv,W_pv,Y_sv,Zn,Zr_sv")

//Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,La,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg - not in order!
#define SPECIE_RAW_LIB4 SPECIE_TRANSITION_METALS

//CO20180729 - OBSOLETE - use xerror
//[OBSOLETE]//CO20180419 - global exception handling - START
//[OBSOLETE]class AFLOWRuntimeError : public std::runtime_error {
//[OBSOLETE]  public:
//[OBSOLETE]    AFLOWRuntimeError(const std::string& function,const std::string& message);
//[OBSOLETE]    AFLOWRuntimeError(const std::string& function,std::stringstream& message);
//[OBSOLETE]    string where();
//[OBSOLETE]    ~AFLOWRuntimeError() throw() {};
//[OBSOLETE]  private:
//[OBSOLETE]    string f_name;  //cannot be const &, as it goes out of scope //const string& f_name;
//[OBSOLETE]};
//[OBSOLETE]class AFLOWLogicError : public std::logic_error {
//[OBSOLETE]  public:
//[OBSOLETE]    AFLOWLogicError(const std::string& function,const std::string& message);
//[OBSOLETE]    AFLOWLogicError(const std::string& function,std::stringstream& message);
//[OBSOLETE]    string where();
//[OBSOLETE]    ~AFLOWLogicError() throw() {};
//[OBSOLETE]  private:
//[OBSOLETE]    string f_name;  //cannot be const &, as it goes out of scope //const string& f_name;
//[OBSOLETE]};
//[OBSOLETE]//CO20180419 - global exception handling - STOP

// --------------------------------------------------------------------------
// definitions for MULTHREADS
//#define MAX_ALLOCATABLE_PTHREADS     1024
#define MAX_ALLOCATABLE_PTHREADS     256
#define PTHREADS_DEFAULT 8
namespace AFLOW_PTHREADS {
  extern bool FLAG;         // run pthread YES/NO
  extern int MAX_PTHREADS;  // how many MAX threads I can use  default or --np
  extern int RUNNING;       // how many threads are actually running
  extern pthread_t vpthread[MAX_ALLOCATABLE_PTHREADS];  // the actual thread
  extern int viret[MAX_ALLOCATABLE_PTHREADS];          // the thread runnings
  extern bool vpthread_busy[MAX_ALLOCATABLE_PTHREADS];  // is the thread busy
}

extern string _AFLOWIN_; 
extern string _AFLOWLOCK_; 

// --------------------------------------------------------------------------
// definitions for aflow
// aflow2 default definitions
#define AFLOW_MATERIALS_SERVER_DEFAULT        string("materials.duke.edu")
#define AFLOW_WEB_SERVER_DEFAULT              string("nietzsche.mems.duke.edu")
#define AFLOWLIB_SERVER_DEFAULT               string("aflowlib.duke.edu")
#define AFLOWLIB_MATERIALS_SERVER             string("aflow.org")
#define AFLOWLIB_CONSORTIUM_STRING            string("aflow.org")
#define _XENTRY_ string("index.php")

// [OBSOLETE] #define DEFAULT_KZIP_BIN              string("bzip2")           // moved to aflow_aflowrc.cpp in V3.1.194
// [OBSOLETE] #define DEFAULT_KZIP_EXT              string(".bz2")            // moved to aflow_aflowrc.cpp in V3.1.194
#define DEFAULT_KBIN_ALIEN_BIN        string("ls -las")
#define DEFAULT_KBIN_MATLAB_BIN       string("/usr/local/bin/matlab -nodesktop -nosplash -nodisplay ")
// [OBSOLETE] #define DEFAULT_KBIN_CONVERT_BIN      string("convert")
// [OBSOLETE] #define DEFAULT_KBIN_EPSTOPDF_BIN     string("epstopdf")

#define QSUB_COMMAND_DEFAULT          "qsub"
#define QSUB_PARAMS_DEFAULT           " "

#define KBIN_SYMMETRY_SGROUP_RADIUS_DEFAULT 3.0
#define KBIN_SYMMETRY_SGROUP_MAX_NUMBER 1000000

#define KBIN_SUBDIRECTORIES           string("ARUN.")

#define KBIN_NEIGHBOURS_MAX_NUMBER      30000
#define KBIN_NEIGHBOURS_RADIUS_DEFAULT 3.0
#define KBIN_NEIGHBOURS_DRADIUS_DEFAULT 0.1

#define ALIEN_INPUT_FILE_NAME_DEFAULT  "./input"
#define ALIEN_EXTERNAL_INPUT_DEFAULT   "../input_external"
#define ALIEN_OUTPUT_FILE_NAME_DEFAULT  "./output"

// aflow1 definitions (soon to be obsolete)
#define _MPI_NP_STRINGS_ "MPI_NP","mpi_np","-MPI_NP","-mpi_np"
#define _MPI_NCPUS_DEF_ 4
#define VASP_OPTIONS_MPI_DEFAULT         ""
#define VASPLS_BIN_POSTFIX_DEFAULT       "LS"
#define GRND_BIN_DEFAULT                 "./grnd_intel"

#define _VASP_POSCAR_MODE_EXPLICIT_START_ "[VASP_POSCAR_MODE_EXPLICIT]START."  //CO20200624
#define _VASP_POSCAR_MODE_EXPLICIT_STOP_ "[VASP_POSCAR_MODE_EXPLICIT]STOP."  //CO20200624

// --------------------------------------------------------------------------
// definitions for projects
// [OBSOLETE] #define SERVER_PROJECT_GNDSTATE       string("/common/GNDSTATE")
#define LIBRARY_NOTHING 256
extern uint LIBRARY_AUID,LIBRARY_ICSD,LIBRARY_LIB0,LIBRARY_LIB1,LIBRARY_LIB2,LIBRARY_LIB3,LIBRARY_LIB4,LIBRARY_LIB5,LIBRARY_LIB6,LIBRARY_LIB7,LIBRARY_LIB8,LIBRARY_LIB9;  // not in order.. will be nailed by init.cpp

#define LIBRARY_ALL         100

// [OBSOLETE] #define DEFAULT_FILE_AFLOWLIB_ENTRY_OUT        string("aflowlib.out")      // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_AFLOWLIB_ENTRY_JSON       string("aflowlib.json")     // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_EDATA_ORIG_OUT            string("edata.orig.out")    // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_EDATA_RELAX_OUT           string("edata.relax.out")   // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_EDATA_BANDS_OUT           string("edata.bands.out")   // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_DATA_ORIG_OUT             string("data.orig.out")     // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_DATA_RELAX_OUT            string("data.relax.out")    // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_DATA_BANDS_OUT            string("data.bands.out")    // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_EDATA_ORIG_JSON           string("edata.orig.json")   // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_EDATA_RELAX_JSON          string("edata.relax.json")  // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_EDATA_BANDS_JSON          string("edata.bands.json")  // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_DATA_ORIG_JSON            string("data.orig.json")    // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_DATA_RELAX_JSON           string("data.relax.json")   // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_DATA_BANDS_JSON           string("data.bands.json")   // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_TIME_OUT                  string("time")              // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_SPACEGROUP1_OUT           string("SpaceGroup")        // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_SPACEGROUP2_OUT           string("SpaceGroup2")       // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_VOLDISTPARAMS_OUT         string("VOLDISTParams")     // moved to aflow_aflowrc.cpp in V3.1.178
// [OBSOLETE] #define DEFAULT_FILE_VOLDISTEVOLUTION_OUT      string("VOLDISTEvolution")  // moved to aflow_aflowrc.cpp in V3.1.178 
#define _AFLOWLIB_ENTRY_SEPARATOR_       string(" | ")
#define _APENNSY_STYLE_OLD_              FALSE

// --------------------------------------------------------------------------
// definition for frozsl files
#define _FROZSL_VASPSETUP_FILE_ "./aflow.frozsl_vaspsetup_file"
// --------------------------------------------------------------------------
// definitions for WEB PHP
#define AFLOW_PHP_APOOL_REFERENCES       string("19,20,49,50,51,53,54,55,56,57,59,61,62,63,65,66,67,70,71,74,75,76,81,87,99")

// --------------------------------------------------------------------------
// Definitions
//#define DEFAULT_AFLOW_FIND_PARAMETERS "-follow"
#define DEFAULT_AFLOW_FIND_PARAMETERS_NORMAL     string("-follow")
#define DEFAULT_AFLOW_FIND_PARAMETERS_NOLEAF     string("-noleaf -follow")
#define BUFFER_MAXLEN 1024
// [OBSOLETE] #define DEFAULT_AFLOW_PRESCRIPT_OUT            string("aflow.prescript.out")  // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_PRESCRIPT_COMMAND        string("aflow.prescript.command")  // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_POSTSCRIPT_OUT           string("aflow.postscript.out")  // moved to aflow_aflowrc.cpp in V3.1.189
// [OBSOLETE] #define DEFAULT_AFLOW_POSTSCRIPT_COMMAND       string("aflow.postscript.command")  // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_PGROUP_OUT               string("aflow.pgroup.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_PGROUP_XTAL_OUT          string("aflow.pgroup_xtal.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_PGROUPK_OUT              string("aflow.pgroupk.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_PGROUPK_XTAL_OUT         string("aflow.pgroupk_xtal.out") // moved to aflow_aflowrc.cpp in V3.1.189  //DX20171205 - Added pgroupk_xtal
// [OBSOLETE] #define DEFAULT_AFLOW_FGROUP_OUT               string("aflow.fgroup.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_SGROUP_OUT               string("aflow.sgroup.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_AGROUP_OUT               string("aflow.agroup.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_IATOMS_OUT               string("aflow.iatoms.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_PGROUP_JSON              string("aflow.pgroup.json") // moved to aflow_aflowrc.cpp in V3.1.189       //DX20170802 - Add JSON
// [OBSOLETE] #define DEFAULT_AFLOW_PGROUP_XTAL_JSON         string("aflow.pgroup_xtal.json") // moved to aflow_aflowrc.cpp in V3.1.189  //DX20170802 - Add JSON
// [OBSOLETE] #define DEFAULT_AFLOW_PGROUPK_JSON             string("aflow.pgroupk.json")  // moved to aflow_aflowrc.cpp in V3.1.189     //DX20170802 - Add JSON
// [OBSOLETE] #define DEFAULT_AFLOW_PGROUPK_XTAL_JSON        string("aflow.pgroupk_xtal.json") // moved to aflow_aflowrc.cpp in V3.1.189 //DX20170802 - Add JSON //DX20171205 - Added pgroupk_xtal
// [OBSOLETE] #define DEFAULT_AFLOW_FGROUP_JSON              string("aflow.fgroup.json") // moved to aflow_aflowrc.cpp in V3.1.189       //DX20170802 - Add JSON
// [OBSOLETE] #define DEFAULT_AFLOW_SGROUP_JSON              string("aflow.sgroup.json")  // moved to aflow_aflowrc.cpp in V3.1.189      //DX20170802 - Add JSON
// [OBSOLETE] #define DEFAULT_AFLOW_AGROUP_JSON              string("aflow.agroup.json")  // moved to aflow_aflowrc.cpp in V3.1.189      //DX20170802 - Add JSON
// [OBSOLETE] #define DEFAULT_AFLOW_IATOMS_JSON              string("aflow.iatoms.json")  // moved to aflow_aflowrc.cpp in V3.1.189      //DX20170802 - Add JSON
// [OBSOLETE] #define DEFAULT_AFLOW_PHONON_FILE  string("aflow.phonons.out") // abandoned
// [OBSOLETE] #define DEFAULT_AFLOW_ICAGES_OUT               string("aflow.icages.out")  // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_SURFACE_OUT              string("aflow.surface.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_QMVASP_OUT               string("aflow.qmvasp.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_ERVASP_OUT               string("aflow.error.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_IMMISCIBILITY_OUT        string("aflow.immiscibility.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_MEMORY_OUT               string("aflow.memory.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_FROZSL_INPUT_OUT         string("aflow.frozsl_input.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_FROZSL_POSCAR_OUT        string("aflow.frozsl_poscar.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_FROZSL_MODES_OUT         string("aflow.frozsl_energies.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_FROZSL_EIGEN_OUT         string("aflow.frozsl_eigen.out") // moved to aflow_aflowrc.cpp in V3.1.189 
// [OBSOLETE] #define DEFAULT_AFLOW_END_OUT                  string("aflow.end.out") // moved to aflow_aflowrc.cpp in V3.1.189 

// --------------------------------------------------------------------------
// include all prototypes for aflow
#define SWAP(a,b)      {temp=(a);(a)=(b);(b)=temp;}
#define RCYCLIC(a,b,c) {temp=(c);(b)=(a);(c)=(b);a=temp;}
#define LCYCLIC(a,b,c) {temp=(a);(a)=(b);(b)=(c);c=temp;}

#define NANOPARTICLE_RADIUS_DEFAULT   10.0
#define NANOPARTICLE_DISTANCE_DEFAULT 10.0
using aurostd::min;
using aurostd::max;
using aurostd::mod;
using aurostd::_isodd;
using aurostd::_iseven;
using aurostd::_isfloat;
using aurostd::_iscomplex;
using aurostd::sign;
using aurostd::nint;
using aurostd::xcomplex;
using aurostd::xmatrix;
using aurostd::clear;
using aurostd::xvector;
//[OBSOLETE ME20180705]using aurostd::xtensor3;
//[OBSOLETE ME20180705]using aurostd::xtensor4;
//[OBSOLETE ME20180705]using aurostd::xtensor5;
//[OBSOLETE ME20180705]using aurostd::xtensor6;
//[OBSOLETE ME20180705]using aurostd::xtensor7;
//[OBSOLETE ME20180705]using aurostd::xtensor8;
using aurostd::xoption;
using aurostd::xcombos;

//BANDGAP  //CO20191110
#define _METALGAP_ -1.0*AUROSTD_NAN
#define _METALEDGE_ -1.0

// --------------------------------------------------------------------------
// this is a container of general global choices
class _XHOST {
  public:
    // constructor destructor                         // constructor/destructor
    _XHOST();                                         // default, just allocate
    ~_XHOST();                                        // kill everything
    // _XHOST(const _XHOST& b);                          // constructor copy
    const _XHOST& operator=(const _XHOST &b);         // copy
    // BOOT
    int PID,TID;                // aflow_init.cpp  PID/TID number  //CO20200508
    ostringstream ostrPID,ostrTID; // aflow_init.cpp  PID/TID in ostringstream... //CO20200508
    string sPID,sTID;           // aflow_init.cpp  [PID=12345678]  [TID=12345678]
    bool showPID,showTID;       // aflow_init.cpp  check if --showPID
    // machinery
    bool QUIET,TEST,DEBUG,MPI;
    bool GENERATE_AFLOWIN_ONLY; //CT20180719
    bool POSTPROCESS; //CO20200624 - generic postprocessing, including --lib2raw and --lib2lib
    bool ARUN_POSTPROCESS; //CT20181212 - this is for the --postprocess flag needed for AEL/AGL, can be extended to other modules too
    bool AVOID_RUNNING_VASP; //CO20200624
    bool PSEUDOPOTENTIAL_GENERATOR; //SC20200327
    // HARDWARE/SOFTWARE
    string hostname,machine_type,tmpfs,user,group,home,shell,progname;
    string Find_Parameters;
    bool sensors_allowed;
    // ARGUMENTS
    vector<string> argv;          // argv of line command
    // SERVERS
    string AFLOW_MATERIALS_SERVER,AFLOW_WEB_SERVER;
    long double RAM,RAM_MB,RAM_GB;
    int CPU_Cores;
    string CPU_Model;
    string CPU_MHz;
    vector<double> vTemperatureCore;
    long double Time_starting,Time_now;
    long int Date;
    string Day,Month,Year;
    string Copyright_Years; // =string("2003-YEAR_FROM_DATE");
    // MULTHREADS
    // bool PTHREADS_FLAG;        // run pthread YES/NO 
    // int  PTHREADS_MAX;         // how many MAX threads I can use  default or --np 
    // int PTHREADS_RUNNING;      // how many threads are actually running 
    // vector<pthread_t> thread;  // the actual thread 
    // vector<int> iret;          // the thread runnings 
    // vector<bool> thread_busy;  // is the thread busy 
    // COMMANDS
    vector<string> vcmd;  
    // RAM CHECK
    double maxmem;
    // FUNCTIONS
    string command(const string& command);
    bool is_command(const string& command);
    // AFLOW STUFF
    // vflag_aflow.flag("LOOP");
    // vflag_aflow.flag("CLEAN");
    // vflag_aflow.isflag*"XCLEAN");
    bool AFLOW_RUNDIRflag;
    bool AFLOW_MULTIflag;
    bool AFLOW_RUNXflag;
    uint AFLOW_RUNXnumber;
    // QUEQUE STUFF
    bool is_PBS;int PBS_NUM_PPN,PBS_NNODES;
    bool is_SLURM;int SLURM_CPUS_ON_NODE,SLURM_NNODES,SLURM_NTASKS;
    bool is_MACHINE_FULTON_MARYLOU; // some flags
    // APENNST stuff
    bool APENNSY_USE_SERVER;
    bool APENNSY_USE_LIBRARY;
    bool APENNSY_SERVER_AFLOWLIB_ORG;
    // Library_CALCULATED*
    vector<uint>   vGlobal_uint;      // parameters uint
    vector<string> vGlobal_string;    // parameters as strings
    vector<vector<string> > vvGlobal_string; // parameters as vector strings
    vector<vector<string> > vvLIBS; // parameters as vector strings
    // vector<string> vLibrary_ICSD;     // ordered by #species (needs to be allocated)
    // vector<string> vLibrary_ICSD_ALL; // line by line
    // string Library_ICSD_ALL;          // the complete library
    // vector<string> vVASP_POTCAR_DIRECTORIES;
    // vector<string> vAFLOW_LIBRARY_DIRECTORIES;
    // vector<string> vAFLOW_PROJECTS_DIRECTORIES;
    // AFLOW flags/options
    aurostd::xoption vflag_aflow;  // argv/argc options following the xoption structure
    aurostd::xoption vflag_pflow;  // argv/argc options following the xoption structure
    aurostd::xoption vflag_apennsy;  // argv/argc options following the xoption structure
    aurostd::xoption vflag_outreach;  // argv/argc options following the xoption structure
    aurostd::xoption vflag_control;  // argv/argc options following the xoption structure
    aurostd::xoption vschema;        // keywords, names, units etc etc
    // USUAL COMMANDS
    vector<string> vcat; //     cat, bzcat, xzcat, gzcat
    vector<string> vext; //      "",  .bz2,   .xz,   .gz
    vector<string> vzip; //      "", bzip2,    xz,  gzip
    // AFLOWRC
    string aflowrc_filename;
    string aflowrc_content;
    vector<string> vaflowrc;  
    xoption adefault;            // default  xoption
    // AFLOWSYM
    bool SKEW_TEST; //DX20171019
    double SKEW_TOL; //DX20171019
    // WEB MODE
    //[CO20200404 - overload with --www]bool WEB_MODE;  //CO20190401
  private:                                                //
    void free();                                           // free space
    void copy(const _XHOST& b);                            //
    void clear();                                          // free space
};

#define XPID XHOST.sPID
#define XTID XHOST.sTID

#define XHOST_vGlobal_MAX                              256
#define XHOST_Library_HTQC                             XHOST.vGlobal_string.at(0)
#define XHOST_aflowlib_icsd                            XHOST.vGlobal_string.at(23)
#define XHOST_aflowlib_lib0                            XHOST.vGlobal_string.at(24)
#define XHOST_aflowlib_lib1                            XHOST.vGlobal_string.at(25)
#define XHOST_aflowlib_lib2                            XHOST.vGlobal_string.at(26)
#define XHOST_aflowlib_lib3                            XHOST.vGlobal_string.at(27)
#define XHOST_aflowlib_lib4                            XHOST.vGlobal_string.at(28)
#define XHOST_aflowlib_lib5                            XHOST.vGlobal_string.at(29)
#define XHOST_aflowlib_lib6                            XHOST.vGlobal_string.at(30)
#define XHOST_aflowlib_lib7                            XHOST.vGlobal_string.at(31)
#define XHOST_aflowlib_lib8                            XHOST.vGlobal_string.at(32)
#define XHOST_aflowlib_lib9                            XHOST.vGlobal_string.at(33)
//#define XHOST_AUID                                     XHOST.vGlobal_string.at(34)
//#define XHOST_AURL                                     XHOST.vGlobal_string.at(35)
//#define XHOST_LOOP                                     XHOST.vGlobal_string.at(36)
#define XHOST_Library_ICSD_ALL                         XHOST.vGlobal_string.at(38)
//  string Library_ICSD_ALL;          // the complete library

#define XHOST_vLIBS XHOST.vvLIBS
#define XHOST_vAURL XHOST.vvLIBS.at(0)
#define XHOST_vAUID XHOST.vvLIBS.at(1)
#define XHOST_vLOOP XHOST.vvLIBS.at(2)
#define XHOST_LIBRARY_JSONL                            XHOST.vGlobal_string.at(3)

#define vVASP_POTCAR_DIRECTORIES                       XHOST.vvGlobal_string.at(4)
#define vAFLOW_LIBRARY_DIRECTORIES                     XHOST.vvGlobal_string.at(5)
#define vAFLOW_PROJECTS_DIRECTORIES                    XHOST.vvGlobal_string.at(6)
#define XHOST_vLibrary_ICSD                            XHOST.vvGlobal_string.at(7)
#define XHOST_vLibrary_ICSD_ALL                        XHOST.vvGlobal_string.at(8)
#define XHOST_Library_CALCULATED_ICSD_LIB              XHOST.vvGlobal_string.at(9)
#define XHOST_Library_CALCULATED_ICSD_RAW              XHOST.vvGlobal_string.at(10)
#define XHOST_Library_CALCULATED_LIB0_LIB              XHOST.vvGlobal_string.at(11)
#define XHOST_Library_CALCULATED_LIB0_RAW              XHOST.vvGlobal_string.at(12)
#define XHOST_Library_CALCULATED_LIB1_LIB              XHOST.vvGlobal_string.at(13)
#define XHOST_Library_CALCULATED_LIB1_RAW              XHOST.vvGlobal_string.at(14)
#define XHOST_Library_CALCULATED_LIB2_LIB              XHOST.vvGlobal_string.at(15)
#define XHOST_Library_CALCULATED_LIB2_RAW              XHOST.vvGlobal_string.at(16)
#define XHOST_Library_CALCULATED_LIB3_LIB              XHOST.vvGlobal_string.at(17)
#define XHOST_Library_CALCULATED_LIB3_RAW              XHOST.vvGlobal_string.at(18)
#define XHOST_Library_CALCULATED_LIB4_LIB              XHOST.vvGlobal_string.at(19)
#define XHOST_Library_CALCULATED_LIB4_RAW              XHOST.vvGlobal_string.at(20)
#define XHOST_Library_CALCULATED_LIB5_LIB              XHOST.vvGlobal_string.at(21)
#define XHOST_Library_CALCULATED_LIB5_RAW              XHOST.vvGlobal_string.at(22)
#define XHOST_Library_CALCULATED_LIB6_LIB              XHOST.vvGlobal_string.at(23)
#define XHOST_Library_CALCULATED_LIB6_RAW              XHOST.vvGlobal_string.at(24)
#define XHOST_Library_CALCULATED_LIB7_LIB              XHOST.vvGlobal_string.at(25)
#define XHOST_Library_CALCULATED_LIB7_RAW              XHOST.vvGlobal_string.at(26)
#define XHOST_Library_CALCULATED_LIB8_LIB              XHOST.vvGlobal_string.at(27)
#define XHOST_Library_CALCULATED_LIB8_RAW              XHOST.vvGlobal_string.at(28)
#define XHOST_Library_CALCULATED_LIB9_LIB              XHOST.vvGlobal_string.at(29)
#define XHOST_Library_CALCULATED_LIB9_RAW              XHOST.vvGlobal_string.at(30)

//  vector<string> vLibrary_ICSD;     // ordered by #species
//  vector<string> vLibrary_ICSD_ALL; // line by line

#define XHOST_README_AFLOW_LICENSE_GPL3_TXT            XHOST.vGlobal_string.at(40)
#define XHOST_README_AFLOW_TXT                         XHOST.vGlobal_string.at(41)
#define XHOST_README_AFLOW_VERSIONS_HISTORY_TXT        XHOST.vGlobal_string.at(42)
#define XHOST_README_AFLOW_PFLOW_TXT                   XHOST.vGlobal_string.at(43)
#define XHOST_README_AFLOW_APENNSY_TXT                 XHOST.vGlobal_string.at(44)
#define XHOST_README_AFLOW_SCRIPTING_TXT               XHOST.vGlobal_string.at(45)
#define XHOST_README_AFLOW_FROZSL_TXT                  XHOST.vGlobal_string.at(46)
#define XHOST_README_AFLOW_POCC_TXT                    XHOST.vGlobal_string.at(47)
#define XHOST_README_AFLOW_APL_TXT                     XHOST.vGlobal_string.at(48)
#define XHOST_README_AFLOW_QHA_SCQHA_QHA3P_TXT         XHOST.vGlobal_string.at(49)
#define XHOST_README_AFLOW_AGL_TXT                     XHOST.vGlobal_string.at(50)
#define XHOST_README_AFLOW_AEL_TXT                     XHOST.vGlobal_string.at(51)
#define XHOST_README_AFLOW_ANRL_TXT                    XHOST.vGlobal_string.at(52)
#define XHOST_README_AFLOW_COMPARE_TXT                 XHOST.vGlobal_string.at(53)
#define XHOST_README_AFLOW_GFA_TXT                     XHOST.vGlobal_string.at(53)  //CO20190401
#define XHOST_README_AFLOW_SYM_TXT                     XHOST.vGlobal_string.at(54)
#define XHOST_README_AFLOW_CCE_TXT                     XHOST.vGlobal_string.at(55)  //CO20190620
#define XHOST_README_AFLOW_CHULL_TXT                   XHOST.vGlobal_string.at(56)  //CO20190620
#define XHOST_README_AFLOW_EXCEPTIONS_TXT              XHOST.vGlobal_string.at(57)  //ME20180705
#define XHOST_README_PROTO_TXT                         XHOST.vGlobal_string.at(58)
#define XHOST_README_AFLOW_XAFLOW_TXT                  XHOST.vGlobal_string.at(59)
#define XHOST_README_AFLOW_AFLOWRC_TXT                 XHOST.vGlobal_string.at(60)

#define XHOST_FINDSYM_data_space_txt                   XHOST.vGlobal_string.at(70)
#define XHOST_FINDSYM_data_wyckoff_txt                 XHOST.vGlobal_string.at(71)
#define XHOST_FROZSL_data_space_txt                    XHOST.vGlobal_string.at(72)
#define XHOST_FROZSL_data_wyckoff_txt                  XHOST.vGlobal_string.at(73)
#define XHOST_FROZSL_data_images_txt                   XHOST.vGlobal_string.at(74)
#define XHOST_FROZSL_data_irreps_txt                   XHOST.vGlobal_string.at(75)
#define XHOST_FROZSL_data_isotropy_txt                 XHOST.vGlobal_string.at(76)
#define XHOST_FROZSL_data_little_txt                   XHOST.vGlobal_string.at(77)
#define XHOST_FROZSL_symmetry2_dat                     XHOST.vGlobal_string.at(78)
#define XHOST_FROZSL_const_dat                         XHOST.vGlobal_string.at(79)
#define XHOST_FROZSL_phvaspsetup_AFLOW                 XHOST.vGlobal_string.at(80)
#define XHOST_FROZSL_phvaspsetup_POSCAR                XHOST.vGlobal_string.at(81)
#define XHOST_ElectronStoppingPower_txt                XHOST.vGlobal_string.at(82)
#define XHOST_PhotonCrossSection_txt                   XHOST.vGlobal_string.at(83)
#define XHOST_PhotonStoppingPower_txt                  XHOST.vGlobal_string.at(84)
#define XHOST_ICSD_List_txt                            XHOST.vGlobal_string.at(85)
#define XHOST_AFLOW_PSEUDOPOTENTIALS                   XHOST.vGlobal_string.at(86)
#define XHOST_AFLOW_PSEUDOPOTENTIALS_TXT               XHOST.vGlobal_string.at(87)
#define XHOST_AFLOW_PSEUDOPOTENTIALS_LIST_TXT          XHOST.vGlobal_string.at(88)
#define XHOST_f144468a7ccc2d3a72ba44000715efdb         XHOST.vGlobal_string.at(90)
#define XHOST_d0f1b0e47f178ae627a388d3bf65d2d2         XHOST.vGlobal_string.at(91)
#define XHOST_decf00ca3ad2fe494eea8e543e929068         XHOST.vGlobal_string.at(92)
// [OBSOLETE] #define XHOST_AFLOW_BinaryRead           XHOST.vGlobal_string.at(93)
// [OBSOLETE] #define XHOST_AFLOW_Binary_Angle_Read    XHOST.vGlobal_string.at(94)


// LOADENTRIES DEFAULTS
#define _AFLOW_LIB_MAX_ 10                             //LIB11 does not exist yet, modify accordingly

#define XHOST_LIBRARY_LIB0                             XHOST.vGlobal_uint.at(0)
#define XHOST_LIBRARY_LIB1                             XHOST.vGlobal_uint.at(1)
#define XHOST_LIBRARY_LIB2                             XHOST.vGlobal_uint.at(2)
#define XHOST_LIBRARY_LIB3                             XHOST.vGlobal_uint.at(3)
#define XHOST_LIBRARY_LIB4                             XHOST.vGlobal_uint.at(4)
#define XHOST_LIBRARY_LIB5                             XHOST.vGlobal_uint.at(5)
#define XHOST_LIBRARY_LIB6                             XHOST.vGlobal_uint.at(6)
#define XHOST_LIBRARY_LIB7                             XHOST.vGlobal_uint.at(7)
#define XHOST_LIBRARY_LIB8                             XHOST.vGlobal_uint.at(8)
#define XHOST_LIBRARY_LIB9                             XHOST.vGlobal_uint.at(9)
#define XHOST_LIBRARY_ICSD                             XHOST.vGlobal_uint.at(10)
#define XHOST_LIBRARY_AUID                             XHOST.vGlobal_uint.at(11)


// max is 128
extern _XHOST XHOST; // this will be global

//DX20180131 - add symmetry definitions - START
// symmetry 
#define SG_SETTING_1    1
#define SG_SETTING_2    2
#define SG_SETTING_ANRL 3
//DX20180131 - add symmetry definitions - END

//DX20191122 START
// atom environment modes
#define ATOM_ENVIRONMENT_MODE_1    1 // minimum coordination shell
#define ATOM_ENVIRONMENT_MODE_2    2 // [FUTURE] out to a given radius
#define ATOM_ENVIRONMENT_MODE_3    3 // [FUTURE] largest gap in radial distribution function (GFA)
//DX20191122 END

// ----------------------------------------------------------------------------
// aflow_aflowrc.cpp
#define _AFLOW_AFLOWRC_H_
#define _AFLOW_AFLOWRC_CPP_
#include "aflow_aflowrc.cpp"
#undef _AFLOW_AFLOWRC_CPP_
#undef _AFLOW_AFLOWRC_H_

// --------------------------------------------------------------------------

// Structures for flags and properties to share FAST !
// STRUCTURES
#define AFLOWIN_SEPARATION_LINE  string("[AFLOW] ************************************************************************************************************************** ")

// --------------------------------------------------------------------------
// general flags to run aflow
class _aflags {
  public:
    // trivial constructurs/destuctors/operators
    _aflags();                                          // default, just allocate
    ~_aflags();                                         // kill everything
    _aflags(const _aflags& b);                          // constructor copy
    const _aflags& operator=(const _aflags &b);         // copy
    void clear(void);                                   // clear
    // CONTENT
    bool QUIET;
    int  AFLOW_PTHREADS_NUMBER;                         // cant be GLOBAL as this is a local run stuff
    // particular
    string LocalDirectory;                              // where is aflow now
    string Directory;                                   // where aflow must run
    bool AFLOW_FORCE_RUN;                               // Force run also in database
    bool AFLOW_PERFORM_DIRECTORY;                       // Directory is specified (sometimes it is useful).
    bool AFLOW_PERFORM_FILE;                            // File is specified (sometimes it is useful).
    bool AFLOW_PERFORM_ORDER_SORT;                      // Sorts the _AFLOWIN_ in the list
    bool AFLOW_PERFORM_ORDER_REVERSE;                   // Reverse the _AFLOWIN_ in the list
    bool AFLOW_PERFORM_ORDER_RANDOM;                    // Randomize the _AFLOWIN_ in the list
    bool AFLOW_MODE_GENERATE;                           // TODO OVVERRIDE all _AFLOWIN_
    bool AFLOW_MODE_QSUB_MODE1;                         // TODO OVVERRIDE all _AFLOWIN_
    bool AFLOW_MODE_QSUB_MODE2;                         // TODO OVVERRIDE all _AFLOWIN_
    bool AFLOW_MODE_QSUB_MODE3;                         // TODO OVVERRIDE all _AFLOWIN_
    // general flags to operate in the directory
    bool KBIN_RUN_AFLOWIN;
    bool KBIN_GEN_GENERAL; //CO20180409
    bool KBIN_GEN_VASP_FROM_AFLOWIN;
    bool KBIN_GEN_AIMS_FROM_AFLOWIN;
    bool KBIN_GEN_AFLOWIN_FROM_VASP;
    //DX START
    bool KBIN_GEN_SYMMETRY_OF_AFLOWIN;
    //DX END 
    bool KBIN_DELETE_AFLOWIN;
    bool AFLOW_FORCE_MPI;     // not yet implemented
    bool AFLOW_FORCE_SERIAL;  // not yet implemented
    int  AFLOW_GLOBAL_NCPUS;         // Forced CPUS
    // Perform TASKS
    bool AFLOW_PERFORM_CLEAN;       // to clean a directory
    // host related things
    xoption AFLOW_MACHINE_GLOBAL;
    xoption AFLOW_MACHINE_LOCAL;       // flag for duke_beta_mpich
    // APENNSY
    aurostd::xoption vflag;
    string APENNSY_LATTICE_flag;                        // APENNSY flags
    string APENNSY_GNUPLOT_FONT_str;                    // APENNSY content
    string APENNSY_GNUPLOT_FONT_BOLD_str;               // APENNSY content
    string APENNSY_GNUPLOT_FONT_ITALICS_str;            // APENNSY content
    // [OBSOLETE] "APENNSY::HELP;                                   // APENNSY flags
    // [OBSOLETE] "APENNSY::VERBOSE_flag;                           // APENNSY flags
    // [OBSOLETE] "APENNSY_LIST_flag;                               // APENNSY flags
    // [OBSOLETE] "APENNSY::SERVER_AFLOWLIB_ORG_flag;               // APENNSY flags
    // [OBSOLETE] "APENNSY::LATEX_SNAPSHOT;                         // APENNSY flags
    // [OBSOLETE] "APENNSY::LATEX_OUTPUT;                           // APENNSY flags
    // [OBSOLETE] "APENNSY::LATEX_CITE;                             // APENNSY flags
    // [OBSOLETE] "APENNSY::ENTHALPY_TOT;                           // APENNSY flags
    // [OBSOLETE] "APENNSY::ENTHALPY_ATOM;                          // APENNSY flags
    // [OBSOLETE] "APENNSY::ENTHALPY_FORMATION_ATOM;                // APENNSY flags
    // [OBSOLETE] "APENNSY::LOAD_LIB2;                              // APENNSY flags
    // [OBSOLETE] "APENNSY::LOAD_LIB2U;                             // APENNSY flags
    // [OBSOLETE] "APENNSY::LOAD_LIB2PGM;                           // APENNSY flags
    // [OBSOLETE] "APENNSY::LOAD_LIB2X;                             // APENNSY flags
    // [OBSOLETE] "APENNSY::LOAD_ALLOY;                             // APENNSY flags
    // [OBSOLETE] "APENNSY::APOOL_PUBLIC;                           // APENNSY flags
    // [OBSOLETE] "APENNSY::APOOL_PRIVATE;                          // APENNSY flags
    // [OBSOLETE] "APENNSY::APOOL_TEST;                             // APENNSY flags
    // [OBSOLETE] "APENNSY::DATA;                                   // APENNSY flags
    // [OBSOLETE] "APENNSY::UNCLE;                                  // APENNSY flags
    // [OBSOLETE] "APENNSY::WEB;                                    // APENNSY flags
    // [OBSOLETE] "APENNSY::ALL;                                    // APENNSY flags
    // [OBSOLETE] "APENNSY::FCC;                                    // APENNSY flags
    // [OBSOLETE] "APENNSY::BCC;                                    // APENNSY flags
    // [OBSOLETE] "APENNSY::HCP;                                    // APENNSY flags
    // [OBSOLETE] "APENNSY::COUT;                                   // APENNSY flags
    // [OBSOLETE] "APENNSY::CERR;                                   // APENNSY flags
    // [OBSOLETE] "APENNSY::ENTHALPY_LIST"                          // APENNSY flags
    // [OBSOLETE] "APENNSY::PS_ENERGY_LIST"                         // APENNSY flags
    // [OBSOLETE] "APENNSY::CONVEX_HULL"                            // APENNSY flags
    // [OBSOLETE] "APENNSY::MATLAB"                                 // APENNSY flags
    // [OBSOLETE] "APENNSY::GNUPLOT"                                // APENNSY flags
    // [OBSOLETE] "APENNSY::SMALL_CONVEX_HULL_MATLAB"               // APENNSY flags
    // [OBSOLETE] "APENNSY::HISTOGRAM_LIST"                         // APENNSY flags
    // [OBSOLETE] "APENNSY::MATLAB_LIB"                             // APENNSY flags
    // [OBSOLETE] "APENNSY::RULES"                                  // APENNSY flags
    // [OBSOLETE] "APENNSY::STRUCTURES"                             // APENNSY flags
    // [OBSOLETE] "APENNSY::VASPIN;                                 // APENNSY flags
    // [OBSOLETE] "APENNSY::ORDER;                                  // APENNSY flags
    // [OBSOLETE] "APENNSY::INFO;                                   // APENNSY flags
    // [OBSOLETE] "APENNSY::MISCIBILITY;                            // APENNSY flags
    // [OBSOLETE] "APENNSY::MISCIBILITY_EXPERIMENTS;                // APENNSY flags
    // [OBSOLETE] "APENNSY::MISCIBILITY_MIEDEMA;                    // APENNSY flags
    // [OBSOLETE] "APENNSY::MISCIBILITY_HUMEROTHERY;                // APENNSY flags
    // [OBSOLETE] "APENNSY::MISCIBILITY_TABLE;                      // APENNSY flags
    // [OBSOLETE] "APENNSY::MISCIBILITY_STATISTICS;                 // APENNSY flags
    // [OBSOLETE] "APENNSY::STRUCTURE_VOLUMES;                      // APENNSY flags
    // [OBSOLETE] "APENNSY::REFERENCE;                              // APENNSY flags
    // [OBSOLETE] "APENNSY_PROTOCHECK;                              // APENNSY flags
    // [OBSOLETE] "APENNSY::NEGLECT_STRUCTURES;                     // APENNSY flags
    // [OBSOLETE] "APENNSY::UPDATE;                                 // APENNSY flags
    // [OBSOLETE] "APENNSY::CSWAP;                                  // APENNSY flags
    vector<string> APENNSY_NEGLECT_STRUCTURES_vstrs;    // APENNSY content
  private:                                              //
    void free();                                        // free space
    void copy(const _aflags& b);                        //
};

//ME20181026 - Container for APL options
struct _moduleOptions {
  // APL
  vector<aurostd::xoption> aplflags;

  // AAPL
  vector<aurostd::xoption> aaplflags;

  // QHA
  vector<aurostd::xoption> qhaflags;

  // AEL
  vector<aurostd::xoption> aelflags;

  // AGL
  vector<aurostd::xoption> aglflags;
};

// --------------------------------------------------------------------------
// general flags for kbinary (all)
class _kflags {
  public:
    // trivial constructurs/destuctors/operators
    _kflags();                                          // default, just allocate
    ~_kflags();                                         // kill everything
    _kflags(const _kflags& b);                          // constructor copy
    const _kflags& operator=(const _kflags &b);         // copy
    void clear(void);                                   // clear
    // CONTENT
    // in this struct we put all the flags which are used on LOCAL DIRECTORIES in KBIN MODE
    //
    bool AFLOW_MODE_ALIEN;
    //
    bool AFLOW_MODE_MATLAB;
    bool AFLOW_MATLAB_MODE_EXPLICIT;
    bool AFLOW_MATLAB_MODE_EXPLICIT_START_STOP;
    bool AFLOW_MATLAB_MODE_IMPLICIT;
    bool AFLOW_MATLAB_MODE_EXTERNAL;
    bool AFLOW_MATLAB_FILE;
    bool AFLOW_MATLAB_FILE_FILE;
    bool AFLOW_MATLAB_FILE_COMMAND;
    //
    bool AFLOW_MODE_VASP;
    bool AFLOW_MODE_AIMS;
    //
    bool AFLOW_MODE_PRESCRIPT_EXPLICIT;
    bool AFLOW_MODE_PRESCRIPT_EXPLICIT_START_STOP;
    stringstream AFLOW_MODE_PRESCRIPT;
    bool AFLOW_MODE_POSTSCRIPT_EXPLICIT;
    bool AFLOW_MODE_POSTSCRIPT_EXPLICIT_START_STOP;
    stringstream AFLOW_MODE_POSTSCRIPT;
    //
    bool AFLOW_MODE_EMAIL;
    // normal binary
    string KBIN_BIN;
    string KBIN_SERIAL_BIN; //ME20190107
    string KZIP_BIN;
    bool   KZIP_COMPRESS;
    // MPI binaries and flags
    bool   KBIN_MPI;
    int    KBIN_MPI_NCPUS;
    string KBIN_MPI_NCPUS_STRING; //ME20181216
    int    KBIN_MPI_NCPUS_BUFFER;
    string KBIN_MPI_START;
    string KBIN_MPI_STOP;
    string KBIN_MPI_COMMAND;
    bool   KBIN_MPI_AUTOTUNE;
    string KBIN_MPI_BIN;
    string KBIN_MPI_OPTIONS;
    // QSUB
    bool   KBIN_QSUB;
    bool   KBIN_QSUB_MODE1;
    bool   KBIN_QSUB_MODE2;
    bool   KBIN_QSUB_MODE3;
    string KBIN_QSUB_COMMAND;
    string KBIN_QSUB_PARAMS;
    bool   KBIN_QSUB_MODE_EXPLICIT;
    bool   KBIN_QSUB_MODE_EXPLICIT_START_STOP;
    bool   KBIN_QSUB_MODE_IMPLICIT;
    bool   KBIN_QSUB_FILE;
    // symmetry operation lists
    bool  KBIN_SYMMETRY_CALCULATION;
    //DX START
    bool  KBIN_SYMMETRY_NO_SCAN;
    double KBIN_SYMMETRY_EPS;
    bool  KBIN_SYMMETRY_CALCULATE_PGROUP;       //DX20170814 - Specify what to calculate/verify
    bool  KBIN_SYMMETRY_CALCULATE_PGROUPK;      //DX20170814 - Specify what to calculate/verify
    bool  KBIN_SYMMETRY_CALCULATE_FGROUP;       //DX20170814 - Specify what to calculate/verify
    bool  KBIN_SYMMETRY_CALCULATE_PGROUP_XTAL;  //DX20170814 - Specify what to calculate/verify
    bool  KBIN_SYMMETRY_CALCULATE_PGROUPK_XTAL; //DX20171205 - Specify what to calculate/verify; Added pgroupk_xtal
    bool  KBIN_SYMMETRY_CALCULATE_PGROUPK_PATTERSON;  //DX20200129 - Specify what to calculate/verify
    bool  KBIN_SYMMETRY_CALCULATE_IATOMS;       //DX20170814 - Specify what to calculate/verify
    bool  KBIN_SYMMETRY_CALCULATE_AGROUP;       //DX20170814 - Specify what to calculate/verify
    bool  KBIN_SYMMETRY_CALCULATE_SGROUP;       //DX20170814 - Specify what to calculate/verify
    //DX END
    bool  KBIN_SYMMETRY_PGROUP_WRITE;      // taken TRUE by default
    bool  KBIN_SYMMETRY_PGROUPK_WRITE;     // taken TRUE by default
    bool  KBIN_SYMMETRY_PGROUP_XTAL_WRITE; // taken TRUE by default
    bool  KBIN_SYMMETRY_PGROUPK_XTAL_WRITE; //DX20171205 - Added pgroupk_xtal
    bool  KBIN_SYMMETRY_PGROUPK_PATTERSON_WRITE;  //DX20200129 - taken TRUE by default
    bool  KBIN_SYMMETRY_FGROUP_WRITE;      // taken TRUE by default
    bool  KBIN_SYMMETRY_SGROUP_WRITE;
    bool  KBIN_SYMMETRY_AGROUP_WRITE;      // taken TRUE by default
    bool  KBIN_SYMMETRY_IATOMS_WRITE;      // taken TRUE by default
    double KBIN_SYMMETRY_SGROUP_RADIUS;
    // neighbours operation lists
    bool  KBIN_NEIGHBOURS_CALCULATION;
    bool  KBIN_NEIGHBOURS_WRITE;
    double KBIN_NEIGHBOURS_RADIUS;
    double KBIN_NEIGHBOURS_DRADIUS;
    // pocc operation lists
    bool   KBIN_POCC;
    bool   KBIN_POCC_CALCULATION;
    string KBIN_POCC_TEMPERATURE_STRING;  //CO20191114
    string KBIN_POCC_ARUNS2SKIP_STRING;   //CO20200627
    // frozsl operation lists
    bool   KBIN_FROZSL;
    bool   KBIN_FROZSL_DOWNLOAD;
    bool   KBIN_FROZSL_FILE;
    string KBIN_FROZSL_FILE_NAME;
    // [OBSOLETE]  bool   KBIN_FROZSL_PRESCRIPT_MODE_EXPLICIT;
    // [OBSOLETE]  bool   KBIN_FROZSL_PRESCRIPT_MODE_EXPLICIT_START_STOP;
    // [OBSOLETE]  string KBIN_FROZSL_PRESCRIPT_STRING;
    // [OBSOLETE]  bool   KBIN_FROZSL_POSTSCRIPT_MODE_EXPLICIT;
    // [OBSOLETE]  bool   KBIN_FROZSL_POSTSCRIPT_MODE_EXPLICIT_START_STOP;
    // [OBSOLETE]  string KBIN_FROZSL_POSTSCRIPT_STRING;
    bool   KBIN_FROZSL_STRUCTURE_MODE_FILE;
    bool   KBIN_FROZSL_STRUCTURE_MODE_EXPLICIT_START_STOP;
    string KBIN_FROZSL_STRUCTURE_STRING;
    bool   KBIN_FROZSL_DIELECTRIC_MODE_FILE;
    bool   KBIN_FROZSL_DIELECTRIC_MODE_EXPLICIT_START_STOP;
    bool   KBIN_FROZSL_DIELECTRIC_ZEFF;
    string KBIN_FROZSL_DIELECTRIC_STRING;
    // phonons operation lists
    bool   KBIN_PHONONS_CALCULATION_APL;
    bool   KBIN_PHONONS_CALCULATION_QHA;  //CO20170601
    bool   KBIN_PHONONS_CALCULATION_AAPL; //CO20170601
    bool   KBIN_PHONONS_CALCULATION_AGL;
    bool   KBIN_PHONONS_CALCULATION_AEL;
    bool   KBIN_PHONONS_CALCULATION_FROZSL;
    string KBIN_PHONONS_CALCULATION_FROZSL_output;
    string KBIN_PHONONS_CALCULATION_FROZSL_poscars;
    _moduleOptions KBIN_MODULE_OPTIONS;  //ME20181027
  private:                                             //
    void free();                                        // free space
    void copy(const _kflags& b);                        //
};

// --------------------------------------------------------------------------
// general flags for vasp mode
class xstructure; // prototype of structure, just to compile
class _vflags {
  public:
    // trivial constructurs/destuctors/operators
    _vflags();                                            // default, just allocate
    ~_vflags();                                           // kill everything
    _vflags(const _vflags& b);                            // constructor copy
    const _vflags& operator=(const _vflags &b);           // copy
    void clear(void);                                     // clear
    // CONTENT
    // in this struct we put all the flags which are used on LOCAL DIRECTORIES in VASP MODE
    xoption AFLOW_SYSTEM;                         //ME20181121
    int KBIN_VASP_RUN_NRELAX;
    xoption KBIN_VASP_RUN;                        // GENERATE, STATIC, KPOINTS, RELAX, RELAX_STATIC, RELAX_STATIC_BANDS, STATIC_BANDS, DIELECTRIC_STATIC, DIELECTRIC_DYNAMIC, DSCF
    xoption KBIN_VASP_REPEAT;                     // REPEAT_BANDS REPEAT_STATIC_BANDS REPEAT_DELSOL
    xoption KBIN_VASP_FORCE_OPTION_NOTUNE;        // NOTUNE
    xoption KBIN_VASP_FORCE_OPTION_SYSTEM_AUTO;   // SYSTEM_AUTO
    xoption KBIN_VASP_FORCE_OPTION_RELAX_MODE;    // RELAX_MODE  forces/energy
    xoption KBIN_VASP_FORCE_OPTION_RELAX_TYPE;    // RELAX_TYPE  STATIC, ALL, IONS, CELL_SHAPE, CELL_VOLUME, IONS_CELL_VOLUME
    xoption KBIN_VASP_FORCE_OPTION_PREC;          // PREC 
    xoption KBIN_VASP_FORCE_OPTION_ALGO;          // ALGO 
    xoption KBIN_VASP_FORCE_OPTION_METAGGA;       // METAGGA 
    xoption KBIN_VASP_FORCE_OPTION_IVDW;          // IVDW
    xoption KBIN_VASP_FORCE_OPTION_ABMIX;         // ABMIX
    xoption KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS; // AUTO_PSEUDOPOTENTIALS
    // ENMAX_MULTIPLY
    xoption KBIN_VASP_FORCE_OPTION_ENMAX_MULTIPLY_EQUAL;  // isentry and content_int
    // NBANDS
    xoption KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL;  // isentry and content_int
    bool KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry;
    // POTIM
    xoption KBIN_VASP_FORCE_OPTION_POTIM_EQUAL;   // isentry and content_double
    // PSTRESS
    xoption KBIN_VASP_FORCE_OPTION_PSTRESS_EQUAL; // isentry and content_double
    // EDIFFG
    xoption KBIN_VASP_FORCE_OPTION_EDIFFG_EQUAL; // isentry and content_double
    // ISMEAR
    xoption KBIN_VASP_FORCE_OPTION_ISMEAR_EQUAL;  // isentry and content_double //CO20181129
    // SIGMA
    xoption KBIN_VASP_FORCE_OPTION_SIGMA_EQUAL;  // isentry and content_double //CO20181129
    // RWIGS
    bool KBIN_VASP_FORCE_OPTION_RWIGS_STATIC;  
    xoption KBIN_VASP_FORCE_OPTION_SKIP_NOMIX;    // SKIP_NOMIX
    xoption KBIN_VASP_FORCE_OPTION_SPIN;          // SPIN 
    bool KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1;
    bool KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2;
    // xoption KBIN_VASP_FORCE_OPTION_TRISTATE;      //  SYM 
    xoption KBIN_VASP_FORCE_OPTION_BADER;         // BADER=ON | OFF | NONE
    xoption KBIN_VASP_FORCE_OPTION_ELF;           // ELF=ON | OFF | NONE
    xoption KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM;   // AUTO_MAGMOM
    xoption KBIN_VASP_FORCE_OPTION_SYM;           // SYM
    xoption KBIN_VASP_FORCE_OPTION_WAVECAR;       // WAVECAR
    xoption KBIN_VASP_FORCE_OPTION_CHGCAR;        // CHGCAR
    xoption KBIN_VASP_FORCE_OPTION_CHGCAR_FILE;   //ME20191028
    xoption KBIN_VASP_FORCE_OPTION_LSCOUPLING;    // LSCOUPLING
    xoption KBIN_VASP_FORCE_OPTION_LDAU0;         // LDAU0
    xoption KBIN_VASP_FORCE_OPTION_LDAU1;         // LDAU1
    xoption KBIN_VASP_FORCE_OPTION_LDAU2;         // LDAU2
    xoption KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC;// LDAU_ADIABATIC
    xoption KBIN_VASP_FORCE_OPTION_LDAU_CUTOFF;   // LDAU_CUTOFF
    string KBIN_VASP_LDAU_SPECIES;
    string KBIN_VASP_LDAU_PARAMETERS;
    bool KBIN_VASP_LDAU_AFLOW_AUTO_flag;
    // FORCE_OPTION
    xoption KBIN_VASP_FORCE_OPTION_TYPE;          // TYPE 
    bool KBIN_VASP_FORCE_OPTION_NSW_EQUAL;
    int  KBIN_VASP_FORCE_OPTION_NSW_EQUAL_VALUE;

    xoption KBIN_VASP_FORCE_OPTION_IGNORE_AFIX;   // AFIX
    // xoption kopts;
    xoption KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL;   // CONVERT_UNIT_CELL
    xoption KBIN_VASP_FORCE_OPTION_VOLUME;        // EQUAL_EQUAL, MULTIPLY_EQUAL,PLUS_EQUAL
    xoption KBIN_VASP_FORCE_OPTION_KPOINTS;       // KPOINTS 
    xoption KBIN_VASP_INCAR_MODE;                 // EXPLICIT, EXPLICIT_START_STOP, IMPLICIT, EXTERNAL;
    // RELAX
    xoption KBIN_VASP_KPOINTS_MODE;               // EXPLICIT, EXPLICIT_START_STOP, IMPLICIT, EXTERNAL;
    xoption KBIN_VASP_KPOINTS_KMODE;              // isentry and content_int
    xoption KBIN_VASP_KPOINTS_KPPRA;              // isentry and content_int
    xoption KBIN_VASP_KPOINTS_KSCHEME;            // isentry and content_string
    xoption KBIN_VASP_KPOINTS_KSHIFT;             // isentry and content_string
    // STATIC
    xoption KBIN_VASP_KPOINTS_STATIC_KMODE;       // isentry and content_int
    xoption KBIN_VASP_KPOINTS_STATIC_KPPRA;       // isentry and content_int
    xoption KBIN_VASP_KPOINTS_STATIC_KSCHEME;     // isentry and content_string
    xoption KBIN_VASP_KPOINTS_STATIC_KSHIFT;      // isentry and content_string
    // PHONONS
    xoption KBIN_VASP_KPOINTS_PHONONS_KPPRA;      // isentry and content_int
    xoption KBIN_VASP_KPOINTS_PHONONS_KSCHEME;    // isentry and content_string
    xoption KBIN_VASP_FORCE_OPTION_KPOINTS_PHONONS_PARITY;  // EVEN ODD
    xoption KBIN_VASP_KPOINTS_PHONONS_GRID;       //ME20200427
    xoption KBIN_VASP_KPOINTS_PHONONS_SHIFT;      //ME20200427
    // BANDS
    xoption KBIN_VASP_KPOINTS_BANDS_LATTICE;
    //  bool KBIN_VASP_KPOINTS_BANDS_LATTICE_FLAG;
    //  string KBIN_VASP_KPOINTS_BANDS_LATTICE_VALUE;
    bool KBIN_VASP_KPOINTS_BANDS_LATTICE_AUTO_FLAG;
    bool KBIN_VASP_KPOINTS_BANDS_GRID_FLAG;
    uint KBIN_VASP_KPOINTS_BANDS_GRID_VALUE;
    bool KBIN_VASP_WRITE_KPOINTS;

    xoption KBIN_VASP_POSCAR_MODE; // EXPLICIT, EXPLICIT_START_STOP, EXPLICIT_START_STOP_POINT, IMPLICIT, EXTERNAL;
    std::vector<string> KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRING;
    std::vector<xstructure> KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE;
    xoption KBIN_VASP_POTCAR_MODE;  // EXPLICIT, IMPLICIT, EXTERNAL;
    bool KBIN_VASP_INCAR_VERBOSE;   // VERBOSITY
    xoption KBIN_VASP_INCAR_FILE;   // KEYWORD, SYSTEM_AUTO, FILE, COMMAND
    stringstream KBIN_VASP_INCAR_EXPLICIT;  //ME20181127
    stringstream KBIN_VASP_INCAR_EXPLICIT_START_STOP;  //ME20181127
    xoption KBIN_VASP_KPOINTS_FILE; // KEYWORD, FILE, COMMAND
    stringstream KBIN_VASP_KPOINTS_EXPLICIT;  //ME20181127
    stringstream KBIN_VASP_KPOINTS_EXPLICIT_START_STOP;  //ME20181127
    xoption KBIN_VASP_POSCAR_FILE;  // KEYWORD, PROTOTYPE, FILE, COMMAND
    xoption KBIN_VASP_POSCAR_FILE_VOLUME; // EQUAL_EQUAL, MULTIPLY_EQUAL PLUS_EQUAL
    xoption KBIN_VASP_POTCAR_FILE; // KEYWORD, SYSTEM_AUTO, PREFIX, SUFFIX, FILE, COMMAND, WRITE
    stringstream KBIN_VASP_POTCAR_EXPLICIT;  //CO20181226
  private:                                                   //
    void free();                                              // free space
    void copy(const _vflags& b);                              //
};

// --------------------------------------------------------------------------
// general flags for aims mode
class _aimsflags {
  public:
    _aimsflags();
    ~_aimsflags();
    _aimsflags(const _aimsflags& b);
    const _aimsflags& operator=(const _aimsflags& b);
    void clear();
    // CONTENT
    // in this struct we put all the flags which are used on LOCAL DIRECTORIES in AIMS MODE
    xoption KBIN_AIMS_FORCE_OPTION_NOTUNE;
    xoption KBIN_AIMS_RUN;
    xoption KBIN_AIMS_GEOM_MODE;
    xoption KBIN_AIMS_GEOM_FILE;
    std::vector<string> KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRING;
    std::vector<xstructure> KBIN_AIMS_GEOM_MODE_EXPLICIT_VSTRUCTURE;
    xoption KBIN_AIMS_GEOM_FILE_VOLUME;
    xoption KBIN_AIMS_FORCE_OPTION_VOLUME;
    xoption KBIN_AIMS_FORCE_OPTION_CONVERT_UNIT_CELL;
    xoption KBIN_AIMS_CONTROL_MODE;
    xoption KBIN_AIMS_CONTROL_FILE;
    bool KBIN_AIMS_CONTROL_VERBOSE;   // VERBOSITY
  private:
    void free();                                              // free space
    void copy(const _aimsflags& b);                           //
};

// --------------------------------------------------------------------------
// general flags for alien mode
class _alienflags {
  public:
    // trivial constructurs/destuctors/operators
    _alienflags();                                            // default, just allocate
    ~_alienflags();                                           // kill everything
    _alienflags(const _alienflags& b);                        // constructor copy
    const _alienflags& operator=(const _alienflags &b);       // copy
    void clear(void);                                         // clear
    // CONTENT
    bool KBIN_ALIEN_COMMAND_BINARY_FLAG;
    string KBIN_ALIEN_COMMAND_BINARY_VALUE;
    bool KBIN_ALIEN_COMMAND_BINARY_START_STOP_FLAG;
    // in this struct we put all the flags which are used on LOCAL DIRECTORIES in ALIEN MODE
    bool KBIN_ALIEN_FORCE_OPTION_NOTUNE;
    bool KBIN_ALIEN_FORCE_OPTION_SOMETHING;                   // SOMETHING

    bool KBIN_ALIEN_INPUT_MODE_EXPLICIT;
    bool KBIN_ALIEN_INPUT_MODE_EXPLICIT_START_STOP;
    bool KBIN_ALIEN_INPUT_MODE_IMPLICIT;
    bool KBIN_ALIEN_INPUT_MODE_EXTERNAL;
    bool KBIN_ALIEN_INPUT_FILE;
    bool KBIN_ALIEN_INPUT_FILE_FILE_FLAG;
    string KBIN_ALIEN_INPUT_FILE_FILE_VALUE;
    bool KBIN_ALIEN_INPUT_FILE_COMMAND_FLAG;
    string KBIN_ALIEN_INPUT_FILE_COMMAND_VALUE;
    bool KBIN_ALIEN_INPUT_MODE_INPUT_FLAG;
    string KBIN_ALIEN_INPUT_MODE_INPUT_VALUE;
    bool KBIN_ALIEN_OUTPUT_MODE_OUTPUT_FLAG;
    string KBIN_ALIEN_OUTPUT_MODE_OUTPUT_VALUE;
  private:                                              //
    void free();                                             // free space
    void copy(const _alienflags& b);                         //
};

// --------------------------------------------------------------------------
// general container for any set of flags
class _xflags {
  public:
    _xflags();
    _xflags(_vflags& vflags);
    _xflags(_aimsflags& aimsflags);
    _xflags(_alienflags& alienflags);
    ~_xflags();
    _xflags(const _xflags& b);
    const _xflags& operator=(const _xflags& b);
    void clear();
    bool AFLOW_MODE_VASP;
    _vflags vflags;
    bool AFLOW_MODE_AIMS;
    _aimsflags aimsflags;
    bool AFLOW_MODE_ALIEN;
    _alienflags alienflags;
    void setVFlags(_vflags& vflags);
    void setAIMSFlags(_aimsflags& aimsflags);
    //add qe and others here eventually
    void setALIENFlags(_alienflags& alienflags);
  private:
    void free();                                              // free space
    void copy(const _xflags& b);                           //
};

// --------------------------------------------------------------------------
// aflow_init.cpp
namespace init {
  int GetCPUCores();
  int InitMachine(bool INIT_VERBOSE,vector<string>& argv,vector<string>& cmds,std::ostream& outf);  //ME20200724 - changed to int
  string InitLoadString(string string2load,bool=FALSE);
  string InitGlobalObject(string string2load,string="",bool=FALSE);
  string InitLibraryObject(string string2load,bool=FALSE);
  string AFLOW_Projects_Directories(string string2load);
  long GetRAM(void);
  uint GetTEMP(void);
  double WaitTEMP(double TRESHOLD=AFLOWRC_AFLOW_CORE_TEMPERATURE_HALT,ostream& oss=cout,bool LVERBOSE=FALSE,vector<string> vmessage=vector<string>(0));
  uint InitSchema(bool INIT_VERBOSE);
} // namespace init

uint AFLOW_getTEMP(vector<string> argv);
uint AFLOW_monitor(vector<string> argv);
double AFLOW_checkMEMORY(string="",double=102.0);
bool CheckMaterialServer(const string& message);  //CO20200624
bool CheckMaterialServer(void);
string aflow_get_time_string(void);
string aflow_convert_time_ctime2aurostd(const string& time_LOCK); //CO20200624
string aflow_get_time_string_short(void);
// [OBSOLETE] string strPID(void);

string Message(const string& list2print="");  //CO20200713
string Message(const string& list2print,const string& filename);  //CO20200713
string Message(const _aflags& aflags,const string& list2print="",const string& filename="");  //CO20200713
bool AFLOW_BlackList(const string& h);  //CO20200713
namespace init {
  void MessageOption(const string& options, const string& routine,vector<string> vusage);  //CO20200624 - should go to cerr for web //DX20200724 - bool to void
  void MessageOption(const string& options, const string& routine,string vusage);  //CO20200624 - should go to cerr for web //DX20200724 - bool to void
  void ErrorOption(const string& options, const string& routine,vector<string> vusage);  //CO20200624 - should go to cerr for web //DX20200724 - bool to void
  void ErrorOption(const string& options, const string& routine,string vusage);  //CO20200624 - should go to cerr for web //DX20200724 - bool to void
}

// --------------------------------------------------------------------------
// aflow_aflowrc.cpp
namespace aflowrc {
  bool is_available(std::ostream& oss,bool AFLOWRC_VERBOSE);
  bool read(std::ostream& oss,bool AFLOWRC_VERBOSE);
  bool write_default(std::ostream& oss,bool AFLOWRC_VERBOSE);
  bool print_aflowrc(std::ostream& oss,bool AFLOWRC_VERBOSE);
} // namespace aflowrc

// --------------------------------------------------------------------------
// aflow_arguments  
uint PflowARGs(vector<string> &argv,vector<string> &cmds,aurostd::xoption &vpflow); // called inside Init::InitMachine coded in aflow_pflow_main.cpp
uint ApennsyARGs(vector<string> &argv,vector<string> &cmds,aurostd::xoption &vflag); // called inside Init::InitMachine coded in aflow_apennsy_main.cpp

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
// aflow_xatom.cpp

#define _COORDS_FRACTIONAL_ 0
#define _COORDS_CARTESIAN_  1
#define _UPDATE_LATTICE_VECTORS_TO_ABCANGLES_   2
#define _UPDATE_LATTICE_ABCANGLES_TO_VECTORS_   3

#define _PGROUP_ 0             // for point group lattice
#define _PGROUPK_ 5            // for point group klattice
#define _PGROUP_XTAL_ 6        // for point group crystal
#define _PGROUPK_XTAL_ 7       // for point group kcrystal
#define _PGROUPK_PATTERSON_ 8   // for point group Patterson //DX20200129
#define _FGROUP_ 1             // for factor group
#define _SGROUP_ 2             // for space group
#define _AGROUP_ 3             // for site positions point group
#define _IATOMS_ 4             // for equivalent atoms

// --------------------------------------------------------------------------
//DX+CO START
//DX20190214 [OBSOLETE] extern thread_local double _SYM_TOL_; // tolerance control for isequal_RHT in atom class (RHT)
//DX+CO END

class _atom { // simple class.. nothing fancy
  public:
    // constructor destructor                              // constructor/destructor
    _atom();                                               // default, just allocate
    ~_atom();                                              // kill everything
    _atom(const _atom& b);                                 // constructor copy
    const _atom& operator=(const _atom &b);                // copy
    void clear();
    // content                                             // content
    xvector<double> fpos;                                  // positions are with respect to ijk lattice cell
    xvector<double> cpos;                                  // so if fpos/cpos is outside a cell, you can shift
    xvector<double> corigin;                               // origin for convasp purposes
    xvector<double> coord;                                 // general coordinate for symmetry routines (RHT)
    vector<string> fpos_equation;                          //DX20180607 - lattice equation for atomic position 
    vector<string> cpos_equation;                          //DX20180607 - Cartesian equation for atomic position 
    double spin;                                           // spin along z in VASP MODE
    bool spin_is_given;                                    // TRUE if spin has been set //DX20170921
    xvector<double> noncoll_spin;                          // non-collinear spin                //DX20171205
    bool noncoll_spin_is_given;                            // TRUE if noncoll_spin has been set //DX20171205
    double mass;                                           // mass 
    int    type;                                           // with bringincell, which adjust cpos/fpos and ijk as well
    string name;                                           // the name read from the INPUT
    bool   name_is_given;                                  // TRUE is atom name has been given
    string cleanname;                                      // a chemical clean version of the name
    int info;                                              // container for misc. information  //RHT
    int    atomic_number;                                  // 0 by defauls
    int    number;                                         // atom number reference for convasp, from zero to the sky
    string sd;                                             // ?
    xvector<int> ijk;                                      // xvector identifier of the lattice (but you must give str)
    bool   isincell;                                       // is in cell ? (i==j==k==0 ?)
    int    basis;                                          // identifier of position in the basis, from zero to the sky
    double reference;                                      // reference/measure for ordering
    int    ireference;                                     // sort of number in the list
    // for symmetry
    int    equivalent;                                     // points to the equivalent atom in the cell (-1 if irreducible)
    bool   is_inequivalent;                                // if atom is irreducible
    uint   num_equivalents;                                // say how many they are (only for is_inequivalent)
    uint   index_iatoms;                                   // if calculated on the xstructure, the index within iatoms for the identical atoms
    // for order parameter                                 // order parameter
    int    order_parameter_value;                          // order parameter
    bool   order_parameter_atom;                           // order parameter
    // for partial occupation                              // partial occupation
    double partial_occupation_value;                       // partial occupation
    bool   partial_occupation_flag;                        // partial occupation
    int shell;                                             // neighbour shell number
    // printing
    bool   verbose;                                        // verbose in printing
    bool   print_RHT;                                      // a printer for coord and name (general position)   //RHT
    bool   print_cartesian;                                // print frac or cartesian
    // operators/functions                                 // operator/functions
    friend ostream& operator<<(ostream &,const _atom&);    // print
    void CleanName(void);                                  // function to clean up the name
    void CleanSpin(void);                                  // function to clean up the spin from EZ vasp script
    void ClearSymmetry(void);                              // clears symmetry //CO20190219
  private:                                                //
    void free();                                           // free space
    void copy(const _atom& b); //
};

class _atom_reference_cmp {                              // sorting through reference
  public:
    bool operator()(const _atom& atom1,const _atom& atom2) const {
      return (bool) (atom1.reference<atom2.reference);}
};
class _atom_type_cmp {                                   // sorting through type
  public:
    bool operator()(const _atom& atom1,const _atom& atom2) const {
      return (bool) (atom1.type<atom2.type);}
};

#define NUM_ELEMENTS (103+1)  // up to Uranium
extern std::vector<string> vatom_symbol;             // store starting from ONE
extern std::vector<string> vatom_name;               // store starting from ONE
extern std::vector<double> vatom_mass;               // store starting from ONE
extern std::vector<double> vatom_volume;             // store starting from ONE
extern std::vector<int> vatom_valence_iupac;         // store starting from ONE http://en.wikipedia.org/wiki/Valence_(chemistry)
extern std::vector<int> vatom_valence_std;           // store starting from ONE http://en.wikipedia.org/wiki/Valence_(chemistry)
extern std::vector<double> vatom_miedema_phi_star;       // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
extern std::vector<double> vatom_miedema_nws;            // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
extern std::vector<double> vatom_miedema_Vm;             // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
extern std::vector<double> vatom_miedema_gamma_s;        // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
extern std::vector<double> vatom_miedema_BVm;            // store starting from ONE Miedema Rule Table 1a Physica 100B (1980) 1-28
extern std::vector<double> vatom_radius;             // store starting from ONE - Saxena
extern std::vector<double> vatom_radius_covalent;    // store starting from ONE - Codero, Covalent radii revisited, DOI: 10.1039/b801115j //DX+CO20170904
extern std::vector<double> vatom_electronegativity;  // store starting from ONE - Saxena
extern std::vector<string> vatom_crystal;            // store starting from ONE - Ashcroft Mermin
extern std::vector<double> vatom_xray_scatt;              // store starting from ONE
extern std::vector<double> vatom_pettifor_scale;              // store starting from ONE - Chemical Scale Pettifor Solid State Communications 51 31-34 1984
extern std::vector<double> vatom_pearson_coefficient;         //ME20181020

void atoms_initialize(void);
uint GetAtomNumber(const string& symbol);
std::string GetAtomName(const string& symbol);
std::string GetAtomName(const uint& atnum);
std::string GetAtomSymbol(const string& symbol);
std::string GetAtomSymbol(const uint& atnum);
double GetAtomMass(const string& symbol,bool clean=true);  // in Kg //CO20181129
double GetAtomMass(const uint& atnum); // in Kg
double GetAtomComptonCrossSection(const string& symbol); // barn (1 barn = 1e-28 m^2)
double GetAtomComptonCrossSection(const uint& atnum); // barn (1 barn = 1e-28 m^2)
double GetAtomPhotoelectricCrossSection(const string& symbol);  // barn (1 barn = 1e-28 m^2)
double GetAtomPhotoelectricCrossSection(const uint& atnum);  // barn (1 barn = 1e-28 m^2)
double GetAtomVolume(const string& symbol,bool clean=true);  //CO20181129
double GetAtomVolume(const uint& atnum);
int GetAtomValenceIupac(const string& symbol);
int GetAtomValenceIupac(const uint& atnum);
int GetAtomValenceStd(const string& symbol);
int GetAtomValenceStd(const uint& atnum);
double GetAtomRadius(const string& symbol);
double GetAtomRadius(const uint& atnum);
double GetAtomRadiusCovalent(const string& symbol); //DX+CO20170904
double GetAtomRadiusCovalent(const uint& atnum); //DX+CO20170904
double GetAtomElectronegativity(const string& symbol);
double GetAtomElectronegativity(const uint& atnum);
string GetAtomCrystal(const string& symbol);
string GetAtomCrystal(const uint& atnum);
double GetAtomPettiforScale(const string& symbol);
double GetAtomPettiforScale(const uint& atnum);
bool GetAtomPettiforScale(const vector<string>& vsymbol,vector<double>& vvalue);
bool GetAtomPettiforScale(const vector<uint>& vatnum,vector<double>& vvalue);
bool GetAtomPettiforScale(const vector<string>& vsymbol,xvector<double>& vvalue);
bool GetAtomPettiforScale(const vector<uint>& vatnum,xvector<double>& vvalue);
bool SortAtomsPettiforScale(vector<string> &vsymbols,xvector<int> &vorders,xvector<double> &vvalues);
bool SortAtomsPettiforScale(vector<string> &vsymbols,vector<int> &vorders,vector<double> &vvalues);
bool SortAtomsPettiforScale(vector<string> &vsymbol,vector<int> &vorder);
bool SortAtomsPettiforScale(vector<string> &vsymbol,vector<double> &vvalue);
bool SortAtomsPettiforScale(vector<string> &vsymbol);
double GetPearsonCoefficient(const string&);
double GetPearsonCoefficient(const int&);
double GetAtomXrayScatt(const string& symbol);
double GetAtomXrayScatt(const uint& atnum);
vector<string> GetGroupOfAtoms(string& group_name); //DX20181220 
double GetCompoundAttenuationLength(const vector<string>& species,const vector<double>& composition,const double& density);  // density in g/cm^3, return in cm
double GetCompoundAttenuationLength(const deque<string>& _species,const deque<int>& _composition,const double& density);  // density in g/cm^3, return in cm
//DX+CO START
//DX20190214 [OBSOLETE]bool isequalRHT(const _atom& a, const _atom& b,double=_SYM_TOL_);       // bool equality only checks 'coord' and 'name' (RHT)  //RHT
//DX+CO END
// routines of general use
string XATOM_AlphabetizationSpecies(string speciesA,string speciesB);
string XATOM_AlphabetizationSpecies(vector<string> vspecies);
string XATOM_AlphabetizationSpecies(vector<string> vspecies,vector<double> vnumbers);
void XATOM_AlphabetizationSpecies(string& system, vector<string>& vspecies,vector<double>& vnumbers);
void XATOM_AlphabetizationCompound(string& system, vector<string>& vspecies,vector<double>& vnumbers);
void XATOM_AlphabetizationSpecies(string& system, vector<string>& vspecies);
void XATOM_AlphabetizationSpecies(string& system);
void XATOM_AlphabetizationCompound(string& system);
uint XATOM_SplitAlloySpecies(string alloy_in, vector<string> &speciesX);
uint XATOM_SplitAlloySpecies(string alloy_in, vector<string> &speciesX, vector<double> &natomsX);
uint XATOM_SplitAlloyPseudoPotentials(string alloy_in, vector<string> &species_ppX);
uint XATOM_SplitAlloyPseudoPotentials(string alloy_in, vector<string> &species_ppX, vector<double> &natomsX);
vector<uint> composition2stoichiometry(string& composition); //DX20181009
// neighbour things
void GetUnitCellRep(const xvector<double>& ppos,xvector<double>& p_cell0,xvector<int>& ijk,const xmatrix<double>& lattice,const bool coord_flag);

string xstructure2json(xstructure& xstr); //DX20170831 - xstructure2json
string atom2json(_atom& atom, int coord_flag, int poccupation); //DX20170831 - atom2json

// --------------------------------------------------------------------------
class _sym_op {
  public:
    // constructor destructor
    _sym_op();                                                    // default, just allocate
    _sym_op(const _sym_op& b);                                    // constructor copy
    ~_sym_op();                                                   // kill everything
    // content
    // for _PGROUP_
    xmatrix<double>  Uc;            // 3x3                        // uniques (not irreducible) operations on positions (Uc cartesian)
    xmatrix<double>  Uf;            // 3x3                        // uniques (not irreducible) operations on indices   (Uf fractional)
    xmatrix<double>  generator;     // 3x3                        // generator A, U=exp(A*theta)
    xvector<double>  generator_coefficients;                      // generator coefficients on Lx, Ly, Lz basis //DX20171206
    xmatrix<xcomplex<double> > SU2_matrix; // 2x2                 // SU(2) 2x2 complex matrix //DX20180115
    xvector<xcomplex<double> > su2_coefficients;                  // su(2) coefficients on sigma_1, sigma_2, sigma_3 basis (Pauli matrices) //DX20180115
    double           angle;                                       // angle axis
    xvector<double>  axis;          // 3                          // (1,2,3)=axis
    xvector<double>  quaternion_vector;				//GG
    xmatrix<double>  quaternion_matrix;				//GG
    string           str_type;                                    // generic type of the operation
    string           str_Hermann_Mauguin;                         // Hermann_Mauguin notation
    string           str_Schoenflies;                             // Schoenflies notation
    bool             flag_inversion;                              // flag if inversion
    bool             is_pgroup;                                   // bool is_pgroup
    // for _PGROUP_XTAL_
    bool             is_pgroup_xtal;                              // bool is_pgroup_xtal
    // for _PGROUPK_PATTERSON_
    bool             is_pgroupk_Patterson;                        // bool is_pgroup_Patterson //DX20200129
    // for _PGROUPK_
    bool             is_pgroupk;                                  // bool is_pgroupk
    // for _PGROUPK_XTAL_                  
    bool             is_pgroupk_xtal;                             // bool is_pgroupk_xtal //DX20171205 - Added pgroupk_xtal
    // for _FGROUP_
    xvector<double>  ctau;          // 3                          // translation in CARTESIAN       // FACTOR GROUP only, [0,1[
    xvector<double>  ftau;          // 3                          // translation in FRACTIONAL      // FACTOR GROUP only, [0,1[
    vector<int>      basis_atoms_map;                             // this is the vector that tell where the basis atom gets mapped by the operation
    vector<int>      basis_types_map;                             // this is the vector that tell where the basis species gets mapped by the operation
    bool             basis_map_calculated;                        //have we've calculated it?
    bool             is_fgroup;                                   // bool is_fgroup
    // for _SGROUP_
    xvector<double>  ctrasl;        // 3                          // translation in CARTESIAN       // SPACE GROUP only, [integers]
    xvector<double>  ftrasl;        // 3                          // translation in FRACTIONAL      // SPACE GROUP only, [ingegers]
    bool             is_sgroup;                                   // bool is_sgroup
    // operators
    // for _AGROUP_
    uint             site;                                        // uint site          // site index //DX20170803
    bool             is_agroup;                                   // bool is_agroup     // for site operation point group

    const _sym_op& operator=(const _sym_op& b);
    friend ostream& operator<<(ostream &,const _sym_op&);

    void setUc(const xmatrix<double>& Uc,const xmatrix<double>& lattice);  //CO20190321
    void setUf(const xmatrix<double>& Uf,const xmatrix<double>& lattice);  //CO20190321
    void setctau(const xmatrix<double>& Uc,const xmatrix<double>& lattice);  //CO20190321
    void setftau(const xmatrix<double>& Uf,const xmatrix<double>& lattice);  //CO20190321

  private:
    void free();
    void copy(const _sym_op& b);
};

//DX201801107 - add _kpoint class - START
// --------------------------------------------------------------------------
class _kpoint {
  public:
    // constructor destructor
    _kpoint();                                           // default, just allocate
    ~_kpoint();                                          // default, just allocate
    // content
    char iomode;                                         // store format (not used yet)
    xmatrix<double> klattice;                            // reciprocal lattice
    xvector<double> fpos;                                // fractional position of kpoint
    xvector<double> cpos;                                // Cartesian position of kpoint (not used yet)
    string label;                                        // kpoint label (i.e., high-symmetry point labels)
    bool is_transformed;                                 // indicates if kpoint is transformed from AFLOW standard
    const _kpoint& operator=(const _kpoint& b);          // assignment operator
    // operators/functions                               // operator/functions
    string str() const;                                  // prints "fpos ! label" (e.g., 0.0000 0.0000 0.0000 ! \\Gamma)
    void TransformKpoint(const xmatrix<double>& P);      // transforms kpoint via P matrix (k'=k*P) and klattice via Q matrix (L_recip'=Q*L_recip) (see ITC-A pg. 79)
    friend ostream& operator<<(ostream&,const _kpoint&); // ostream operator
  private:
    void free();
};
//DX201801107 - add _kpoint class - END

// --------------------------------------------------------------------------
class wyckoffsite_ITC { //Also for wyckoff sites
  public:
    wyckoffsite_ITC(void);
    wyckoffsite_ITC(const wyckoffsite_ITC& b);
    ~wyckoffsite_ITC(void);
    // OPERATORS                                                  // --------------------------------------
    const wyckoffsite_ITC& operator=(const wyckoffsite_ITC& b);             // some operators
    bool operator<(const wyckoffsite_ITC& b) const;                          // < operator //DX20190130 - so we can sort by Wyckoff letter, then by species
    friend ostream& operator<<(ostream&,const wyckoffsite_ITC&);       // ostream
    // CONTENT
    xvector<double> coord;
    string type; //chemical label etc
    string wyckoffSymbol;
    string letter;                                                //DX20190128 - add Wyckoff letter
    string site_symmetry;                                         //DX20190128 - add Wyckoff site symmetry
    uint multiplicity;                                            //DX20190128 - add Wyckoff multiplicity
    double site_occupation;                                       //DX20190128 - add Wyckoff site occupation
    vector<vector<string> > equations;                            //DX20190128 - add Wyckoff equations
  private:                                                       // ---------------------------------------
    void free();                                                  // to free everything
};

//DX20181010 - grouped Wyckoff class - START
// --------------------------------------------------------------------------
// ===== GroupedWyckoffPosition Class ===== //
class GroupedWyckoffPosition{
  public:
    GroupedWyckoffPosition();
    ~GroupedWyckoffPosition();
    friend ostream& operator<<(ostream& oss, const GroupedWyckoffPosition& GroupedWyckoffPosition);
    const GroupedWyckoffPosition& operator=(const GroupedWyckoffPosition& b);
    bool operator<(const GroupedWyckoffPosition& b) const;
    GroupedWyckoffPosition(const GroupedWyckoffPosition& b);
    uint type;
    string element;
    vector<string> site_symmetries;
    vector<uint> multiplicities;
    vector<string> letters;
  private:
    void free();
    void copy(const GroupedWyckoffPosition& b);
};
// --------------------------------------------------------------------------
//DX20181010 - grouped Wyckoff class - END

// --------------------------------------------------------------------------
// AtomEnvironment Class //DX20191120 
class AtomEnvironment{
  public:
    AtomEnvironment();                                                                      // constructor operator
    ~AtomEnvironment();                                                                     // destructor operator
    friend ostream& operator<<(ostream& oss, const AtomEnvironment& AtomEnvironment);       // stringstream operator (printing)
    const AtomEnvironment& operator=(const AtomEnvironment& b);                             // assignment operator
    AtomEnvironment(const AtomEnvironment& b);                                              // copy constructor
    string element_center;                                                                  // species/element at center of environment                                                                   
    uint type_center;                                                                       // type (uint) at center of environment
    vector<string> elements_neighbor;                                                       // species/element of atoms neighboring center atom
    vector<uint> types_neighbor;                                                            // types (uint) of atoms neighboring center atom
    vector<double> distances_neighbor;                                                      // distances to atoms neighboring atoms (typically put in a bin with small tolerance threshold)                                             
    vector<uint> coordinations_neighbor;                                                    // coordination of neighboring distance                                              
    vector<vector<xvector<double> > > coordinates_neighbor;                                 // coordinates of atoms neighboring atoms (center is assumed to be zero,i.e. coord=neighbor-origin)
    //functions
    void getAtomEnvironment(const xstructure& xstr, uint center_index, uint mode=ATOM_ENVIRONMENT_MODE_1);                                          // get environment around atom index                                               
    void getAtomEnvironment(const xstructure& xstr, uint center_index, const vector<string>& neighbor_elements, uint mode=ATOM_ENVIRONMENT_MODE_1); // get restricted environment (via specified elements) around atom index
  private:
    void free();                                                                            // free operator
    void copy(const AtomEnvironment& b);                                                    // copy constructor
};
// --------------------------------------------------------------------------

#define MAX_TITLE_SIZE 512

#define IOAFLOW_AUTO   0
#define IOVASP_AUTO    1
#define IOVASP_POSCAR  2
#define IOVASP_ABCCAR  3
#define IOVASP_WYCKCAR 4
#define IOQE_AUTO      5
#define IOQE_GEOM      6
#define IOABINIT_AUTO  7
#define IOABINIT_GEOM  8
#define IOAIMS_AUTO    9
#define IOAIMS_GEOM   10
#define IOCIF         11 //DX20180723
#define IOELK_AUTO    12 //DX20200310
#define IOELK_GEOM    13 //DX20200310

#define NOSG string("NNN #0")

#define _EQUIV_FPOS_EPS_    2.0e-5    // NOV 2009 Israel  used to be 1.0e-6 too small for ICSD
#define _pocc_no_sublattice_ -1
//[CO20181226 OBSOLETE]#define DEFAULT_PARTIAL_OCCUPATION_TOLERANCE 0.02

bool sortAtomsTypes(const _atom& a1,const _atom& a2);		// sort atoms by types
bool sortAtomsNames(const _atom& a1,const _atom& a2);		// sort atoms by names
bool sortAtomsDist(const _atom& a1,const _atom& a2);		// sort atoms by dist  //CO20180420
bool sortAtomsEquiv(const _atom& a1,const _atom& a2); // cluster by equivalent atoms //CO20190116

class xstructure {
  public:
    // constructors/destructors                                   // --------------------------------------
    xstructure(string="");                                        // constructor default
    xstructure(const xstructure& b);                              // constructor copy
    xstructure(istream& input,int=IOVASP_POSCAR);                 // constructor from istream
    xstructure(ifstream& input,int=IOVASP_POSCAR);                // constructor from ifstream
    xstructure(stringstream& input,int=IOVASP_POSCAR);            // constructor from stringstream
    xstructure(const string& input,int);                          // constructor from file
    xstructure(const string& url,const string& file,int=IOVASP_POSCAR); // constructor from URL
    ~xstructure();                                                // destructor
    // I/O, mutators                                              // --------------------------------------
    bool GetStoich(void);                                         // get stoich_each_type - CO20170724
    bool sortAtomsEquivalent(void);                               // sort by equivalent atoms - CO20190116
    bool FixLattices(void);                                       // Reciprocal/f2c/c2f
    void SetCoordinates(const int& mode);                         // change coordinates
    void MakeBasis(void);                                         // make basis for atoms (basis and number)
    void MakeTypes(void);                                         // refresh types based on num_each_type  //CO20180420
    void AddAtom(const _atom& atom);                              // adding an atom
    void AddAtom_POCC(const _atom& atom);                         // adding an atom FOR POCC ONLY
    void RemoveAtom(const uint& iat);                             // deleting an atom (index)
    void RemoveAtom(vector<uint>& v_atoms_to_remove);             // deleting many atoms (indices)
    void ReplaceAtoms(const deque<_atom>& new_atoms);             // replace all atoms SAFELY/CLEANLY
    void RemoveCopies(double=1.0e-3);                             // deleting atoms too close F/C
    void RemoveFractionalCopies(double=1.0e-3);                   // deleting atoms too close F
    void RemoveCartesianCopies(double=1.0e-3);                    // deleting atoms too close C
    void AddCorners(void);                                        // for picturing purpose
    void clear(void);                                             // clear everything //DX20191220 - uppercase to lowercase clear
    void clean(void);                                             // performs stringstream clean //DX20191220 - uppercase to lowercase clean
    void ClearSpecies(void);                                      // Clear all the symmetry
    void ShiftOriginToAtom(const int& iat);                       // Shift the origin to atom(iat)
    void IdenticalAtoms(void);                                    // Make identical atoms
    void SwapCoordinates(const uint& i,const uint& j);            // Permute Coordinates i with j
    string SpeciesLabel(const uint& A);                           // Returns the Label of the specie A (if available)
    void SpeciesSwap(const uint& A,const uint& B);                // Permute Species A with B (safe for species C).
    bool SpeciesGetAlphabetic(void);                              // Check is species are in alphabetic order
    bool SpeciesPutAlphabetic(void);                              // Put Species in alphabetic
    string SpeciesString(void);                                   // Gives a string with the list of all the species
    uint SetSpecies(const std::deque<string>& vspecies);          // Set the species
    void GetLatticeType(xstructure& sp,xstructure& sc);           // Get all lattices
    void GetLatticeType(void);                                    // Get all lattices
    void Standard_Primitive_UnitCellForm(void);                   // Reduce the Unit Cell to Standard Primitive Form
    void GetStandardPrimitive(void);                              // stub for void Standard_Primitive_UnitCellForm(void);
    void Standard_Conventional_UnitCellForm(void);                // Reduce the Unit Cell to Standard Conventional Form
    void GetStandardConventional(void);                           // stub for void Standard_Conventional_UnitCellForm(void);
    void NiggliUnitCellForm(void);                                // Reduce the Unit Cell to Niggli Form
    void MinkowskiBasisReduction(void);                           // Reduce the Basis to the max orthogonality (Minkowski)
    void LatticeReduction(void);                                  // Lattice Reduction to Max Orthogonality (MINK) and then Niggly Form
    //DX20190905 [OBSOLETE] void BringInCell(void);                                       // Bring all the atoms in the origin
    //DX20190905 [OBSOLETE] void BringInCell(double);                                     // Bring all the atoms in the origin
    void BringInCell(double tolerance=_ZERO_TOL_, double upper_bound=1.0, double lower_bound=0.0); //DX20190904
    void BringInCompact(void);                                    // Bring all the atoms near the origin
    void BringInWignerSeitz(void);                                // Bring all the atoms in the Wigner Seitz Cell
    void GetPrimitive(void);                                      // Make it primitive, if possible
    void GetPrimitive(double tol);                                // Make it primitive, if possible
    void GetPrimitive1(void);                                     // Make it primitive, if possible
    void GetPrimitive2(void);                                     // Make it primitive, if possible
    void GetPrimitive3(void);                                     // Make it primitive, if possible
    uint GetPrimitiveCell(void);                                  // Make it primitive, if possible. Returns 1 if routine fails (RHT)   //RHT
    double MinDist(void);                                         // get minimum interatomic distance -- CO20171024
    void ReScale(const double &in_scale);                         // Change scale but keep volume fixed
    void SetScale(const double &in_scale);                        // Change scale
    void UpdateCartesianCoordinates();                            //AS20200514
    void SetVolume(const double &in_volume);                      // Change volume
    void SetAutoVolume(bool use_AFLOW_defaults_in=false);         // Change volume to sum of atoms  //CO20191010
    void InflateLattice(const double &coefficient);               // Inflate lattice
    void InflateVolume(const double &coefficient);                // Inflate volume
    string platon2print(bool,bool,double,double,double,double);   // Create Platon input file >=51108
    void FakeNames(void);                                         // Fix names as fakes - useful for platon
    string platon2sg(bool P_EQUAL=DEFAULT_PLATON_P_EQUAL,
        bool P_EXACT=DEFAULT_PLATON_P_EXACT,
        double P_ang=DEFAULT_PLATON_P_ANG,
        double P_d1=DEFAULT_PLATON_P_D1,
        double P_d2=DEFAULT_PLATON_P_D2,
        double P_d3=DEFAULT_PLATON_P_D3);
    string findsym2sg(double tolerance=DEFAULT_FINDSYM_TOL);
    string findsym2execute(double tolerance=DEFAULT_FINDSYM_TOL);
    string findsym2print(double tolerance=DEFAULT_FINDSYM_TOL);
    //  string platon2sg(void);
    double GetVolume(void) const;                                 // Return volume  //CO20200201
    double Volume(void) const;                                    // Return volume  //CO20200201
    double GetZVAL(const vector<double>& vZVAL);                  // Given the ZVAL of each species, it returns total ZVAL of cell
    double GetPOMASS(const vector<double>& vPOMASS);              // Given the POMASS of each species, it returns total POMASS of cell
    void ClearSymmetry(void);                                     // Clear all the symmetry
    bool CalculateSymmetry(bool,double);                          // Calculate the symmetry
    bool CalculateSymmetry(void);                                 // Calculate the symmetry
    void CalculateSymmetryPointGroup(bool);                       // Calculate the symmetry
    void CalculateSymmetryPointGroup(void);                       // Calculate the symmetry
    void CalculateSymmetryFactorGroup(bool);                      // Calculate the symmetry
    void CalculateSymmetryFactorGroup(void);                      // Calculate the symmetry
    void CalculateSymmetryPointGroupCrystal(bool);                // Calculate the symmetry
    void CalculateSymmetryPointGroupCrystal(void);                // Calculate the symmetry
    void CalculateSymmetryPointGroupKLattice(bool);               // Calculate the symmetry
    void CalculateSymmetryPointGroupKLattice(void);               // Calculate the symmetry
    void CalculateSymmetryPointGroupKCrystal(bool);               // Calculate the symmetry  //ME20200114
    void CalculateSymmetryPointGroupKCrystal(void);               // Calculate the symmetry  //ME20200114
    void CalculateSymmetryPointGroupKPatterson(bool);             // Calculate the symmetry  //ME20200129
    void CalculateSymmetryPointGroupKPatterson(void);             // Calculate the symmetry  //ME20200129
    int  GenerateGridAtoms(int,int,int,int,int,int);              // generate grid of atoms
    int  GenerateGridAtoms(int,int,int);                          // generate grid of atoms
    int  GenerateGridAtoms(int);                                  // generate grid of atoms
    int  GenerateLIJK(double);                                    // generate lijk look up table
    // QUANTUM ESPRESSO AND ABINIT AND AIMS                       // --------------------------------------
    void buildGenericTitle(bool vasp_input=false,bool force_fix=false); // build a nice title with atoms
    void xstructure2qe(void);                                     // some wrap up IOs to convert format to QE
    void xstructure2vasp(void);                                   // some wrap up IOs to convert format to VASP
    void xstructure2abinit(void);                                 // some wrap up IOs to convert format to ABINIT
    void xstructure2aims(void);                                   // some wrap up IOs to convert format to AIMS
    void xstructure2cif(void);                                    // some wrap up IOs to convert format to CIF //DX20190123
    void xstructure2abccar(void);                                 // some wrap up IOs to convert format to ABCCAR //DX20190123
    void xstructure2elk(void);                                    // some wrap up IOs to convert format to ELK //DX20200313
    //[CO20180420 - moved outside of xstructure]bool sortAtomsTypes(const _atom& a1,const _atom& a2);		// sort atoms by types
    //[CO20180420 - moved outside of xstructure]bool sortAtomsNames(const _atom& a1,const _atom& a2);		// sort atoms by names
    // OPERATORS                                                  // --------------------------------------
    const xstructure& operator=(const xstructure& b);             // some operators
    friend istream& operator>>(istream&,xstructure&);             // istream
    friend ostream& operator<<(ostream&,const xstructure&);       // ostream
    // CONTENT                                                    // --------------------------------------
    string title;                                                 // Title of the structure
    string directory;                                             // Directory where xstructure came from //DX
    string prototype;                                             // Prototype of the structure
    string info;                                                  // Info of the structure
    int iomode;                                                   // IOVASP_POSCAR/IOXXXX
    // int num_types=num_each_type.size();                        // old useless stuff
    // int num_atoms=atoms.size();                                // old useless stuff
    bool neg_scale;                                               // flag for negative scale (for printing)
    double scale;                                                 // scale (always linear A)
    bool neg_scale_second;                                        // POCC (hnf vs. tol) //CO20180409
    double scale_second;                                          // POCC hnf/stoich tol/site tol //CO20180409
    aurostd::xoption scale_third;                                 // if there is a third scale number provided, we use isentry and content_double //CO20170803 - site tol
    char coord_type[2];                                           // type of coordinates
    bool coord_flag;                                              // _COORDS_FRACTIONAL_ (0) fractional, _COORDS_CARTESIAN_ (1) cartesian.
    bool isd;                                                     // TRUE=Selective dynamics; FALSE=no selective dynamics.
    xmatrix<double> lattice;                                      // LATTICE in REAL SPACE (meters)            // vector per RAW (must trasp per algebra)
    double a,b,c,alpha,beta,gamma;                                // LATTICE in a,b,c,alpha,beta,gamma
    xmatrix<double> klattice;                                     // LATTICE in MOMENTUM SPACE (1/meters)      // vevror per RAW (must trasp per algebra)
    xvector<double> origin;                                       // origin
    xmatrix<double> f2c;                                          // transformation matrix for F vector per COLUM f2c=trasp(lattice)
    xmatrix<double> c2f;                                          // transformation matrix for C vector per ROW   c2f=inv(trasp(lattice))
    double equiv_fpos_epsilon;                                    // when they are the same DEFAULT _EQUIV_FPOS_EPS_
    std::deque<int> num_each_type;                                // WARNING: we use starting from 0
    std::deque<double> comp_each_type;                            // WARNING: we use starting from 0
    std::deque<double> stoich_each_type;                          // WARNING: we use starting from 0 - 20170724
    std::deque<_atom> atoms;                                      // WARNING: we use starting from 0
    std::deque<string> species,species_pp,species_pp_type,species_pp_version; // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
    std::deque<double> species_pp_ZVAL; // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
    std::deque<std::deque<double> > species_pp_vLDAU;             // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
    std::deque<double> species_volume;                            // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
    std::deque<double> species_mass;                              // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
    //  ----------------------------------------------------------------------------------------
    // SYMBOLIC MATH stuff
    bool symbolic_math_representation_only;                       // print symbolic math representation only //DX20180618 
    bool constrained_symmetry_calculation;                        // append symbolic math representation for constrained symmetry calculation //DX20180618 
    vector<vector<string> > symbolic_math_lattice;                // symbolic math representation of lattice //DX20180618 
    uint num_parameters;                                          // number of parameters ANRL 20180618
    uint num_lattice_parameters;                                  // number of lattice parameters ANRL 20180618
    vector<string> prototype_parameter_list;                      // prototype parameter list ANRL 20180618
    vector<double> prototype_parameter_values;                    // prototype parameters values ANRL 20180618
    //  ----------------------------------------------------------------------------------------
    bool is_vasp4_poscar_format;                                  // flags for VASP4*
    bool is_vasp5_poscar_format;                                  // flags for VASP5*
    bool Niggli_calculated;                                       // flags for calculation
    bool Niggli_avoid;                                            // flags for avoiding the calculation
    bool Minkowski_calculated;                                    // flags for calculation
    bool Minkowski_avoid;                                         // flags for avoiding the calculation
    bool LatticeReduction_calculated;                             // flags for calculation
    bool LatticeReduction_avoid;                                  // flags for avoiding the calculation
    //  ----------------------------------------------------------------------------------------
    // PRINTING stuff
    string PrintSymbolicMathRepresentation(void);                 // Print symbolic math representation of structure //DX20180618
    string PrintUNCLE(void);                                      // Print in UNCLE format
    //  ----------------------------------------------------------------------------------------
    // LATTICE stuff
    bool Standard_Lattice_calculated;                             // flags for calculation
    bool Standard_Lattice_avoid;                                  // flags for avoiding the calculation
    bool Standard_Lattice_primitive;                              // flags for calculation
    bool Standard_Lattice_conventional;                           // flags for calculation
    bool Standard_Lattice_has_failed;                             // flags for Lattice has failed ?
    string bravais_lattice_type;                                  // lattice type as a string  (14)
    string bravais_lattice_variation_type;                        // lattice type as a string WSETYAWAN mod  (with the mods of WSETYAWAN)
    string bravais_lattice_system;                                // lattice system http://en.wikipedia.org/wiki/Bravais_lattice (7)
    string bravais_lattice_lattice_type;                          // lattice_lattice type as a string  (14)
    string bravais_lattice_lattice_variation_type;                // lattice_lattice type as a string WSETYAWAN mod  (with the mods of WSETYAWAN)
    string bravais_lattice_lattice_system;                        // lattice_lattice system http://en.wikipedia.org/wiki/Bravais_lattice (7)
    string pearson_symbol;                                        // pearson symbol as a string
    string reciprocal_lattice_type;                               // reciprocal lattice type as a string
    string reciprocal_lattice_variation_type;                     // reciprocal lattice type as a string WSETYAWAN mod
    //string reciprocal_conventional_lattice_type;                // reciprocal lattice type as a string
    string bravais_superlattice_type;                             // super lattice type as a string (identical atoms)
    string bravais_superlattice_variation_type;                   // super lattice type as a string (identical atoms) WSETYAWAN mod
    string bravais_superlattice_system;                           // lattice system http://en.wikipedia.org/wiki/Bravais_lattice (7)
    string pearson_symbol_superlattice;                           // pearson symbol of the superlattice (identical atoms)
    bool volume_changed_original2new;                             // flag for volume has changed between original and new (i.e., transformation won't work) //DX20181105
    xmatrix<double> transform_coordinates_original2new;           // transform coordinate system from original to new; (Q in ITC notation) //DX20181105
    xmatrix<double> transform_coordinates_new2original;           // transform coordinate system from new to original; (Q^-1 in ITC notation) //DX20181105
    xmatrix<double> rotate_lattice_original2new;                  // rotate from original to new lattice; (P in ITC notation) //DX20181105
    xmatrix<double> rotate_lattice_new2original;                  // rotate from new to original lattice; (P^-1 in ITC notation) //DX20181105
    //  ----------------------------------------------------------------------------------------
    // GENERAL PURPOSE LABELS                                     // general purpose label
    uint label_uint;                                              // general purpose label_uint
    int label_int;                                                // general purpose label_int
    double label_double;                                          // general purpose label_double
    // ----------------------------------------------------------------------------------------
    // ORDER PARAMETER                                            // order parameter for xstructure
    bool order_parameter_structure;                               // order parameter for xstructure
    std::vector<uint> order_parameter_atoms;                      // indices of atoms to be shuffled
    uint order_parameter_orbit;                                   // number of equivalent configurations with the factor group
    int order_parameter_sum;                                      // sum of all the order parameters
    // ----------------------------------------------------------------------------------------
    // PARTIAL OCCUPATION                                         // partial occupation for xstructure
    bool partial_occupation_flag;                                 // flags for partial occupation TRUE/FALSE
    double partial_occupation_site_tol;                           // tolerance for partial occupation site >=0.0 <=1.0   //CO20180409
    double partial_occupation_stoich_tol;                         // tolerance for partial occupation stoich >=0.0 <=1.0 //CO20180409
    int partial_occupation_HNF;                                   // volume HNF size
    vector<int> partial_occupation_sublattice;                    // contains the information about the specie# of the sublattice in the partial occupation otherwise _pocc_no_sublattice_
    // ----------------------------------------------------------------------------------------
    // GEOMETRY ENERGETICS after the QM calculations              // --------------------------------------
    void qm_clear(void);                                          // QM create/clean all the vectors
    void qm_recycle(void);                                        // QM shift data from QM to GEOM
    void qm_load(string directory,string suffix="",int=IOVASP_POSCAR);                    // QM results load from an ab-initio calculation
    bool qm_calculated;                                           // QM calculation
    double qm_scale;                                              // QM scale (always linear A)
    xmatrix<double> qm_lattice;                                   // QM LATTICE in REAL SPACE (meters)
    xmatrix<double> qm_klattice;                                  // QM LATTICE in MOMENTUM SPACE (1/meters)     // SAVED TRASP
    xvector<double> qm_origin;                                    // QM origin
    xmatrix<double> qm_f2c;                                       // QM transformation matrix for F vector per COLUM f2c=trasp(lattice)
    xmatrix<double> qm_c2f;                                       // QM transformation matrix for C vector per ROW   c2f=inv(trasp(lattice))
    std::deque<_atom> qm_atoms;                                  // QM WARNING: we use starting from 0
    std::vector<xvector<double> > qm_forces;                      // QM FORCES calculation
    bool qm_forces_write;                                         // QM FORCES calculation
    std::vector<xvector<double> > qm_positions;                   // QM POSITIONS calculation
    bool qm_positions_write;                                      // QM POSITIONS calculation
    double qm_E_cell,qm_dE_cell,qm_H_cell,qm_PV_cell,qm_P,qm_mag_cell;                     // QM energetics PER CELL
    double qm_E_atom,qm_dE_atom,qm_H_atom,qm_PV_atom,qm_mag_atom;                     // QM energetics ATOMIC
    // ----------------------------------------------------------------------------------------
    // KPOINTS                                                    // --------------------------------------
    int kpoints_mode;                                             // mode of kpoints
    int    kpoints_k1,kpoints_k2,kpoints_k3;                      // parameters that are plug during
    double kpoints_s1,kpoints_s2,kpoints_s3;                      // parameters that are plug during
    int kpoints_kmax,kpoints_kppra;                               // load/unload and calculations
    string kpoints_kscheme;                                       // of ab-initio
    // ---------------------- SYMMETRY --------------------------------------------------------
    // A=U*B but in A and B we plug vectors as columns watch lattice is per row Uc=A*inv(B)
    // A is the lattice (vectors per colum), B is the test lattice (epr column)
    // Uc is the point_group operation which operates AFTER the vector (row)
    // as: new_vector_row=old_vector_row*Uc  and point_group is the list of all the Uc !!!
    //DX+CO START
    double dist_nn_min; 
    // SYMMETRY TOLERANCE ----------------------------
    bool sym_eps_calculated;                                      // was it calculated automatically per symmetry operations (aflowSYM)?
    double sym_eps;                                               // universal tolerance for symmetry (dictates resolution and mapping tolerances)                     
    uint sym_eps_change_count;                                  // universal tolerance count for symmetry //DX20180223 - added count to xstructure
    //DX+CO END
    // POINT GROUP                                                // POINT GROUP LATTICE
    std::vector<_sym_op> pgroup;                                  // rotations/inversions operations
    bool pgroup_calculated;                                       // WARNING: we use starting from 0
    // POINT GROUP CRYSTAL                                        // POINT GROUP CRYSTAL
    std::vector<_sym_op> pgroup_xtal;                             // rotations/inversions operations
    bool pgroup_xtal_calculated;                                  // WARNING: we use starting from 0
    string crystal_family;                                        // crystal and point group properties
    string crystal_system;                                        // crystal and point group properties
    string point_group_crystal_class;                             // crystal and point group properties
    string point_group_Shoenflies;                                // crystal and point group properties
    string point_group_Hermann_Mauguin;                           // crystal and point group properties
    string point_group_orbifold;                                  // crystal and point group properties
    string point_group_type;                                      // crystal and point group properties
    string point_group_order;                                     // crystal and point group properties
    string point_group_structure;                                 // crystal and point group properties
    // POINT GROUP PATTERSON                                      // POINT GROUP PATTERSON //DX20200129
    std::vector<_sym_op> pgroupk_Patterson;                       // rotations/inversions operations
    bool pgroupk_Patterson_calculated;                            // WARNING: we use starting from 0
    // POINT GROUP KLATTICE                                       // POINT GROUP
    std::vector<_sym_op> pgroupk;                                 // rotations/inversions operations
    bool pgroupk_calculated;                                      // WARNING: we use starting from 0
    // POINT GROUP KCRYSTAL                                       // POINT GROUP
    std::vector<_sym_op> pgroupk_xtal;                            // rotations/inversions operations
    bool pgroupk_xtal_calculated;                                 // WARNING: we use starting from 0
    // FACTOR GROUP                                               // FACTOR GROUP
    std::vector<_sym_op> fgroup;                                  // rotations/inversions + incell_translations operations
    bool fgroup_calculated;                                       // WARNING: we use starting from 0
    // SPACE GROUP                                                // SPACE GROUP
    std::vector<_sym_op> sgroup;                                  // rotations/inversions + incell        //outcell_translations operations
    bool sgroup_calculated;                                       // WARNING: we use starting from 0
    double sgroup_radius;                                         // radius of application (all ops connecting objects inside the sphere)
    xvector<int> sgroup_radius_dims;                              // dimension of the radius (in +- integers)
    // SITE POINT GROUP                                           // SITE POINT GROUP
    bool agroup_calculated;                                       // WARNING: we use starting from 0
    std::vector<std::vector<_sym_op> > agroup;                    // rotations/inversions operations on sites, has one for each atom (all)
    // INEQUIVALENTE ATOMS                                        // --------------------------------------
    bool iatoms_calculated;                                       // given the symmetry, the atoms are mapped in inequivalent
    std::vector<std::vector<int> > iatoms;                        // equivalent/inequivalent atoms lookup table
    // SPACE GROUP CALCULATION WITH PLATON/FINDSYM                // with platon >= 51108
    string spacegroup;                                            // space group symbol
    string spacegrouplabel;                                       // the number with #
    string spacegroupoption;                                      // origin, axes and so on
    int    spacegroupnumber;                                      // the number
    int    spacegroupnumberoption;                                // the option as number
    bool is_spacegroup_platon,is_spacegroup_findsym,is_spacegroup_aflow; // got spacegroup
    // SPACE GROUP ITC
    //DX+CO START
    string crystal_system_ITC;                                    // aflow_symmetry_spacegroup.cpp (RHT)
    string point_group_ITC;                                       // aflow_symmetry_spacegroup.cpp (RHT)
    char bravais_label_ITC;                                       // aflow_symmetry_spacegroup.cpp (RHT)
    char lattice_label_ITC;                                       // aflow_symmetry_spacegroup.cpp (RHT)
    uint space_group_ITC;                                         // aflow_symmetry_spacegroup.cpp (RHT)
    string wyckoff_library_entry_ITC;                             // aflow_symmetry_spacegroup.cpp (RHT)
    int setting_ITC;                                              // aflow_symmetry_spacegroup.cpp (RHT) //DX20170830 - SGDATA
    xvector<double> origin_ITC;                                   // aflow_symmetry_spacegroup.cpp (RHT) //DX20170830 - SGDATA
    vector<string> general_position_ITC;                          // aflow_symmetry_spacegroup.cpp (RHT) //DX20170830 - SGDATA
    //DX+CO END
    // double volume;  USE double GetVolume(const xstructure& a);  and SetVolume //
    //  int number_of_atoms; USE (int) atoms.size();  looks you do not need.
    //DX+CO START
    vector<string> wyccar_ITC;                                    // aflow_symmetry_spacegroup.cpp (RHT)
    xmatrix<double> standard_lattice_ITC;                         // aflow_symmetry_spacegroup.cpp (RHT)
    deque<_atom> standard_basis_ITC;                              // aflow_symmetry_spacegroup.cpp (RHT)    
    vector<wyckoffsite_ITC> wyckoff_sites_ITC;                    // aflow_symmetry_spacegroup.cpp (RHT) //(x,y,z) XX
    vector<string> wyckoff_symbols_ITC;                           // aflow_symmetry_spacegroup.cpp (RHT)
    uint SpaceGroup_ITC(void);                                    // aflow_symmetry_spacegroup.cpp (RHT)
    uint SpaceGroup_ITC(double& use_tol);                         // aflow_symmetry_spacegroup.cpp (RHT)
    uint SpaceGroup_ITC(double& use_tol, bool& no_scan);          // aflow_symmetry_spacegroup.cpp (RHT)
    //uint SpaceGroup_ITC(double& use_tol, const int& manual_it);// aflow_symmetry_spacegroup.cpp (RHT)
    uint SpaceGroup_ITC(double& use_tol, const int& setting);// aflow_symmetry_spacegroup.cpp (RHT) //DX20180806
    //DX20170905 - [OBSOLETE] uint SpaceGroup_ITC(double& use_tol, double& orig_tolerance, const int& manual_it,int& change_sym_count,bool& no_scan);// aflow_symmetry_spacegroup.cpp (RHT)
    //uint SpaceGroup_ITC(double& use_tol,const int& manual_it,bool& no_scan);// aflow_symmetry_spacegroup.cpp (RHT)
    uint SpaceGroup_ITC(double& use_tol,const int& setting,bool& no_scan);// aflow_symmetry_spacegroup.cpp (RHT) //DX20180806
    uint SpaceGroup_ITC(double& use_tol,const int& manual_it,const int& setting,bool& no_scan);// aflow_symmetry_spacegroup.cpp (RHT) //DX20180806
    string aflow2sg(void);                                        // aflow_symmetry_spacegroup.cpp (DX)
    string aflow2sg(double& use_tol);                             // aflow_symmetry_spacegroup.cpp (DX)
    string aflow2sg(double& use_tol, const int& manual_it);       // aflow_symmetry_spacegroup.cpp (DX)
    //DX+CO END
    // ---------------------- FROZSL ---------------------- 
    // [OBSOLETE]  string prototype_label();
    ostringstream FROZSL_output(vector<string> Kvectors);
    // ---------------------- PHONONS ---------------------------------------------------------
    // based on the Maradudin harmonic and deformation analysis
    // LIJK OBEJCTS                                               // LIJK OBEJCTS   WORKING
    bool lijk_calculated;                                         // is calculated ?
    vector<xvector<int> >    lijk_table;                          // bravais lattice look up table l <-> i,j,k
    vector<xvector<double> > lijk_fpos;                           // bravais lattice look up table fpos same as lijk_table
    vector<xvector<double> > lijk_cpos;                           // bravais lattice look up table cpos
    xvector<int> lijk_dims;                                       // dimension
    // ----------------------------------------------------------------------------------------
    // GRID ATOMS                                                 // --------------------------------------
    bool grid_atoms_calculated;                                   // GRID ATOMS from dimsL to dimsH has been calculated ?
    xvector<int> grid_atoms_dimsL;                                // dims low of i,j,k (minus is NOT included)
    xvector<int> grid_atoms_dimsH;                                // dims high of i,j,k (plus is NOT included)
    std::deque<_atom> grid_atoms;                                // WARNING: we use starting from 0
    int grid_atoms_number;                                        // how many...
    vector<int> grid_atoms_sc2pcMap;                              //CO20170804 - mapping between grid_atoms (sc) and atoms (pc)
    vector<int> grid_atoms_pc2scMap;                              //CO20170804 - mapping between grid_atoms (sc) and atoms (pc)
    // ----------------------------------------------------------------------------------------
    // NEIGHBOURS OBEJCTS EXPERIMENTAL/UNFINISHED                 // NEIGHBOURS OBEJCTS   WORKING EXPERIMENTAL
    void neighbours_clear(void);                                  // NN create/clean all the vectors
    void neighbours_calculate(void);                              // NN shift data from QM to GEOM
    bool neighbours_calculated;                                   // NN calculation
    double neighbours_radius;                                     // radius of application
    double neighbours_dradius;                                    // delta radius of shell (bins)
    vector<vector<double> > neighbours_atoms_func_r_vs_nn;        // contains function distance vs neighbours (each atom)
    vector<vector<int> > neighbours_atoms_func_num_vs_nn;         // contains function distance vs neighbours (each atom)
    vector<double> neighbours_func_r_vs_nn;                       // contains function distance vs neighbours (all atoms)
    vector<int> neighbours_func_num_vs_nn;                        // contains function number vs neighbours (all atoms)
    // xvector<int> ndims;                                        // dimension of the radius (in +- integers)
    // std::deque<_atom> ashell;                                 // all the atoms in the shell
    // std::deque<deque<_atom> > natoms;                        // vector of vectors
    // std::vector<vector<double> > rshell;                       // vector of shells
    // std::vector<vector<int> > nshell;                          // vector of density in shells
    // int nbins;                                                 // number of bins
    // NEIGHBOURS OBEJCTS OLD-ACONVASP BUT WORKS                  // NEIGHBOURS OBEJCTS 
    // GetNeighData collects all the neighbor data between rmin and rmax and stores it for each atom in a vector of atom objects in order of increasing distance.  
    void GetNeighData(const deque<_atom>& in_atom_vec,const double& rmin, const double& rmax,deque<deque<_atom> >& neigh_mat);
    // GetStrNeighData collects all the neighbor data out to some cutoff and stores it for each atom in the structure.
    void GetStrNeighData(const double cutoff,deque<deque<_atom> >& neigh_mat);
    // ----------------------------------------------------------------------------------------
    // OUTPUT/ERROR FLAGS                                         // --------------------------------------
    bool Niggli_has_failed;                                       // Niggli has failed ?
    bool Minkowski_has_failed;                                    // Minkowski has failed ?
    bool LatticeReduction_has_failed;                             // LatticeReduction has failed ?
    bool write_lattice_flag;                                      // flag for OUTPUT printing
    bool write_klattice_flag;                                     // flag for OUTPUT printing
    bool write_inequivalent_flag;                                 // flag for OUTPUT printing
    bool write_DEBUG_flag;                                        // flag for OUTPUT printing
    bool error_flag;                                              // flag TRUE is error
    string error_string;                                          // contains type of error
    // END OF CONTENT                                             //
  private:                                                       // ---------------------------------------
    void free();                                                  // to free everything
    void copy(const xstructure& b);                               // the flag is necessary because sometimes you need to allocate the space.
};

void LightCopy(const xstructure&, xstructure&);  //ME20200220

//CO20180420
//for stream management with objects
class xStream {
  public:
    //NECESSARY PUBLIC CLASS METHODS - START
    //constructors - START
    xStream(ostream& oss=cout);
    xStream(ofstream& ofs,ostream& oss=cout);
    //constructors - STOP
    ~xStream();
    //NECESSARY PUBLIC CLASS METHODS - END

    //getters
    ostream* getOSS() const; //CO20191110
    ofstream* getOFStream() const; //CO20191110
  protected:
    //NECESSARY private CLASS METHODS - START
    void free();
    void copy(const xStream& b);
    void initialize(ostream& oss=cout);  //ME20200427
    void initialize(ofstream& ofs,ostream& oss=cout);  //ME20200427
    //NECESSARY END CLASS METHODS - END
    //logger variables
    ostream* p_oss;
    ofstream* p_FileMESSAGE;
    bool f_new_ofstream;  //for deletion later
    //general setters
    void setOFStream(ofstream& FileMESSAGE);
    void setOSS(ostream& oss);
};

// for queue
class _xqsub {
  public:
    // trivial constructurs/destuctors/operators
    _xqsub();                                                     // default, just allocate
    ~_xqsub();                                                    // kill everything
    _xqsub(const _xqsub& b);                                      // constructor copy
    const _xqsub& operator=(const _xqsub &b);                     // copy
    void clear(void);                                             // clear
    stringstream QSUB;                                            //
    stringstream QSUB_orig;                                       //
    bool         QSUB_generated;                                  //
    bool         QSUB_changed;                                    //
  private:                                                       //
    void free();                                                  // free space
    void copy(const _xqsub& b);                                   //
};

// for a vasp run
class _xvasp {
  public:
    // trivial constructurs/destuctors/operators
    _xvasp();                                                     // default, just allocate
    ~_xvasp();                                                    // kill everything
    _xvasp(const _xvasp& b);                                      // constructor copy
    const _xvasp& operator=(const _xvasp &b);                     // copy
    void clear(void);                                             // clear
    // aflow_xatom.cpp contains the code
    // CONTENT
    xstructure   str;
    string       Directory;
    string       AnalyzeLabel;
    _xqsub       xqsub;
    xoption      aopts;
    // --------------------------------
    // VASP INPUT CONTENT
    // [OBSOLETE] bool         AFLOWIN_FLAG::VASP;
    stringstream POSCAR;
    stringstream POSCAR_orig;
    // [OBSOLETE] bool         POSCAR_generated;
    // [OBSOLETE] bool         POSCAR_changed;
    uint         POSCAR_index;
    stringstream INCAR;
    stringstream INCAR_orig;
    // [OBSOLETE] bool         INCAR_generated;
    // [OBSOLETE] bool         INCAR_changed;
    stringstream KPOINTS;
    stringstream KPOINTS_orig;
    // [OBSOLETE] bool         KPOINTS_generated;
    // [OBSOLETE] bool         KPOINTS_changed;
    stringstream POTCAR;
    stringstream POTCAR_orig;
    // [OBSOLETE] bool         POTCAR_generated;
    // [OBSOLETE] bool         POTCAR_changed;
    string       POTCAR_TYPE;
    bool         POTCAR_TYPE_DATE_PRINT_flag;  // to add Zr:POT_PAW:01Apr2001 to directory...
    bool         POTCAR_TYPE_PRINT_flag;  // to add Zr:POT_PAW to directory... (no date) //CO20181226
    double       POTCAR_ENMAX;
    double       POTCAR_ENMIN;
    bool         POTCAR_PAW;
    stringstream POTCAR_POTENTIALS;
    deque<string> POTCAR_AUID; // md5sum one line by one line of the POTCAR
    // --------------------------------
    // VASP INPUT CONTENT
    stringstream OUTCAR;  // OUTPUT
    stringstream CONTCAR; // OUTPUT
    stringstream OSZICAR; // OUTPUT
    // --------------------------------
    // QE INPUT CONTENT
    // [OBSOLETE] bool         AFLOWIN_FLAG::QE;                // FUTURE
    stringstream QE_GEOM;                        // FUTURE
    stringstream QE_GEOM_orig;                   // FUTURE
    bool         QE_GEOM_generated;              // FUTURE
    bool         QE_GEOM_changed;                // FUTURE
    uint         QE_GEOM_index;                  // FUTURE
    // --------------------------------
    // ABINIT INPUT CONTENT
    // [OBSOLETE] bool         AFLOWIN_FLAG::ABINIT;            // FUTURE
    // --------------------------------
    // AIMS INPUT CONTENT
    // [OBSOLETE] bool         AFLOWIN_FLAG::AIMS;            // FUTURE
    // --------------------------------
    int          NCPUS;
    int          NRELAX; // job to do         // -1 (static) 0(error) 1,2,3,4.... (relaxes)  -2(run kpoints)
    int          NRELAXING; // job doing (to monitor odd/even)
    // --------------------------------
    // content for AVASP generation
    bool   AVASP_aflowin_only_if_missing;                          //
    bool AVASP_arun;  //ME20181019
    string AVASP_arun_mode;  //ME20181019
    string AVASP_arun_runname;  //ME20181019
    xoption aplopts;  //ME20181025
    xoption aaplopts;  //ME20181025
    xoption qhaopts; //AS20200302
    string AVASP_dirbase;
    string AVASP_libbase;
    string AVASP_label;
    string AVASP_parameters;
    string AVASP_pocc_parameters; //CO20181226
    string AVASP_pocc_tol; //CO20181226
    double AVASP_volume_in;
    string AVASP_potential;
    bool   AVASP_alpha_fix;
    uint   AVASP_prototype_mode;
    bool   AVASP_prototype_from_library_;
    bool   AVASP_directory_from_library_;
    // [OBSOLETE] xoption AVASP_flag_AFLOW_WRITE;
    // [OBSOLETE] bool   AVASP_flag_AUTO_PSEUDOPOTENTIALS;
    // [OBSOLETE] bool   AVASP_flag_NBANDS;
    // [OBSOLETE] bool   AVASP_flag_POTIM;
    // [OBSOLETE] double AVASP_value_POTIM;
    // [OBSOLETE] bool   AVASP_flag_PSTRESS;
    // [OBSOLETE] double AVASP_value_PSTRESS;
    // [OBSOLETE] bool   AVASP_flag_EDIFFG;
    // [OBSOLETE] double AVASP_value_EDIFFG;  
    // [OBSOLETE] bool   AVASP_flag_SKIP_NOMIX;
    // [OBSOLETE] bool   AVASP_flag_WAVECAR;
    // [OBSOLETE] bool   AVASP_flag_CHGCAR;
    // [OBSOLETE] bool   AVASP_flag_SPIN;
    // [OBSOLETE] bool   AVASP_flag_SPIN_REMOVE_RELAX_1;
    // [OBSOLETE] bool   AVASP_flag_SPIN_REMOVE_RELAX_2;
    // [OBSOLETE] bool   AVASP_flag_BADER;
    // [OBSOLETE] bool   AVASP_flag_ELF;
    // [OBSOLETE] bool   AVASP_flag_LSCOUPLING;
    // [OBSOLETE] bool   AVASP_flag_AUTO_MAGMOM;
    // [OBSOLETE] bool   AVASP_flag_RELAX_FORCES;
    int    AVASP_value_NSW;
    // [OBSOLETE] bool   AVASP_flag_KPPRA;
    int    AVASP_value_KPPRA;
    string AVASP_KSCHEME;
    int    AVASP_value_KPPRA_STATIC;
    string AVASP_STATIC_KSCHEME;
    string AVASP_KPOINTS; //ME20181226
    // [OBSOLETE] bool   AVASP_flag_PRECISION_flag;
    string AVASP_flag_PRECISION_scheme;
    // [OBSOLETE] bool   AVASP_flag_PRECISION_preserved;
    // [OBSOLETE] bool   AVASP_flag_ALGO_flag;
    string AVASP_flag_ALGO_scheme;
    // [OBSOLETE] bool   AVASP_flag_ALGO_preserved;
    // [OBSOLETE] string AVASP_flag_METAGGA_scheme;
    // [OBSOLETE] tring AVASP_flag_IVDW_scheme;
    // [OBSOLETE] bool   AVASP_flag_ABMIX_flag;
    string AVASP_flag_ABMIX_scheme;
    xoption AVASP_flag_TYPE;   // TYPE 
    // [OBSOLETE] bool   AVASP_flag_forceLDAU;
    // [OBSOLETE] bool   AVASP_flag_forceNOLDAU;  
    // [OBSOLETE] bool   AVASP_flag_LDAU1;
    // [OBSOLETE] bool   AVASP_flag_LDAU2;
    // [OBSOLETE] bool   AVASP_flag_LDAU_ADIABATIC;
    // [OBSOLETE] bool   AVASP_flag_LDAU_CUTOFF;
    string AVASP_LDAU_PARAMETERS_STRING;
    double AVASP_LDAU_PARAMETERS_UJSUM;
    // [OBSOLETE] xoption AVASP_flag_CONVERT_UNIT_CELL;   // CONVERT_UNIT_CELL
    // [OBSOLETE] bool   AVASP_flag_PRESERVE_VOLUME;
    // [OBSOLETE] bool   AVASP_flag_EXTRA_INCAR;
    stringstream AVASP_EXTRA_INCAR;
    vector<string> AVASP_INCAR_KEYWORD; //ME20181226
    stringstream AVASP_INCAR_EXPLICIT_START_STOP; //ME20181226
    vector<string> AVASP_KPOINTS_KEYWORD; //ME20181226
    stringstream AVASP_KPOINTS_EXPLICIT_START_STOP; //ME20181226
    vector<string> AVASP_POTCAR_KEYWORD; //ME20181226
    bool   AVASP_flag_MPI;
    bool   AVASP_flag_RUN_RELAX;
    bool   AVASP_flag_RUN_RELAX_STATIC;
    bool   AVASP_flag_RUN_RELAX_STATIC_BANDS;
    bool   AVASP_flag_RUN_STATIC_BANDS;
    string AVASP_path_BANDS;
    uint   AVASP_value_BANDS_GRID;
    bool   AVASP_flag_RUN_STATIC;
    bool   AVASP_flag_GENERATE;
    // [OBSOLETE] xoption AVASP_flag_preserve;   // PRESERVE
    // [OBSOLETE] //  bool   AVASP_flag_preserve_POSCAR;
    // [OBSOLETE] //  bool   AVASP_flag_preserve_KPOINTS;
    // [OBSOLETE] //  bool   AVASP_flag_preserve_CHGCAR;
    // [OBSOLETE] //  bool   AVASP_flag_preserve_WAVECAR;
    // [OBSOLETE] //  bool   AVASP_flag_preserve_WAVEDER;
    // -------------------------------- FUNCTIONS
    double GetZVAL(void);
    double GetCellAtomZVAL(string mode);  // CELL ATOM
    double GetPOMASS(void);
    double GetCellAtomPOMASS(string mode); // CELL ATOM
  private:                                                       //
    void free();                                                  // free space
    void copy(const _xvasp& b);                                   //
};

//for an aims run
class _xaims {
  public:
    _xaims();
    ~_xaims();
    _xaims(const _xaims& b);
    const _xaims& operator=(const _xaims& b);
    void clear();
    uint          GEOM_index;
    xstructure    str;
    string        Directory;
    _xqsub        xqsub;
    xoption       aopts;
    int           NCPUS;
    // --------------------------------
    // AIMS INPUT CONTENT
    stringstream  CONTROL;
    stringstream  CONTROL_orig;
    bool          CONTROL_generated;
    bool          CONTROL_changed;
    string        CONTROL_FILE_NAME;
    stringstream  GEOM;
    stringstream  GEOM_orig;
    bool          GEOM_generated;
    bool          GEOM_changed;
    string        GEOM_FILE_NAME;
    string        OUTPUT_FILE_NAME;
  private:                                                       //
    void free();                                                  // free space
    void copy(const _xaims& b);                                   //
};

// for a alien run
class _xalien {
  public:
    // trivial constructurs/destuctors/operators
    _xalien();                                                    // default, just allocate
    ~_xalien();                                                   // kill everything
    _xalien(const _xalien& b);                                    // constructor copy
    const _xalien& operator=(const _xalien &b);                   // copy
    void clear(void);                                             // clear
    // aflow_xatom.cpp contains the code
    // CONTENT
    string       Directory;
    _xqsub       xqsub;
    stringstream INPUT;
    stringstream INPUT_orig;
    bool         INPUT_generated;
    bool         INPUT_changed;
    string       INPUT_FILE_NAME;
    string       OUTPUT_FILE_NAME;
    // ----------------
    int          NCPUS;
    int          NRELAX;          // -1 (static) 0(error) 1,2,3,4.... (relaxes)  -2(run kpoints)
  private:                                                       //
    void free();                                                  // free space
    void copy(const _xalien& b);                                  //
};

// for a generic run
class _xinput {
  public:
    _xinput();
    _xinput(_xvasp& xvasp);
    _xinput(_xaims& xaims);
    _xinput(_xalien& xalien);
    ~_xinput();
    _xinput(const _xinput& b);
    const _xinput& operator=(const _xinput &b);
    void clear();
    bool AFLOW_MODE_VASP;
    _xvasp xvasp;
    bool AFLOW_MODE_AIMS;
    _xaims xaims;
    bool AFLOW_MODE_ALIEN;
    _xalien xalien;
    void setXVASP(_xvasp& xvasp);
    void setXAIMS(_xaims& xaims);
    void setXALIEN(_xalien& xalien);
    xstructure& getXStr();
    string& getDirectory();
    void setXStr(const xstructure& str,bool set_all=false);
    void setDirectory(const string Directory,bool set_all=false);
  private:
    void free();
    void copy(const _xinput& b);
};

// typedef struct { //
//   //  _xvasp   *pxvasp; //
//   _aflags  *paflags; //
//   // _kflags  *pkflags; //
//   // _vflags  *pvflags; //
//   // string   stringA; //
//   // string   stringB; //
//   // int      mode; //
//   // ofstream *pFileMESSAGE; //
//   // bool     *pQUIET; //
//   char  ***pargv; //
//   bool    *pbusy; //
// } _threaded_KBIN_params; //

//void xstructure::free();
//void xstructure::copy(const xstructure& b);
//xstructure::xstructure(string structure_title);
//xstructure::xstructure(const xstructure& b);
//xstructure::~xstructure();
//const xstructure& xstructure::operator=(const xstructure& b);
xstructure GetStructure(const int& iomode,ifstream& input);     // plug from cin
xstructure GetStructure(const int& iomode,const string& Directory); // plug from a directory
//void xstructure::SetCoordinates(const int& mode);
xstructure SetSDNumbers(const xstructure& a,const vector<string>& in_sd);
xstructure SetSDTypes(const xstructure& a,const vector<string>& in_sd);
vector<int> GetTypes(const xstructure& a);
vector<string> GetNames(const xstructure& a);
vector<string> GetCleanNames(const xstructure& a);
vector<double> GetSpins(const xstructure& a);
string GetElementName(string stringin);
string GetSpaceGroupName(int spacegroupnumber, string directory=""); //DX20180526 - add directory
int GetSpaceGroupNumber(const string& spacegroupsymbol, string directory=""); //DX20190708
string GetSpaceGroupLabel(int spacegroupnumber);
string GetSpaceGroupSchoenflies(int spacegroupnumber, string directory=""); //DX20170901 //DX20180526 - add directory
string GetSpaceGroupHall(int spacegroupnumber, int setting=1, string directory=""); //DX20170901 //DX20180526 - add directory //DX20180806 - added setting
string GetLaueLabel(string& point_group); //DX20170901 //DX20180526 - add directory

#define RADIANS 0
#define DEGREES  1
#define _calculate_symmetry_default_sgroup_radius_   2.0

xmatrix<double> MetricTensor(const xstructure& a); //CO20180409
xmatrix<double> MetricTensor(const xmatrix<double>& lattice,double scale=1.0); //CO20180409
xmatrix<double> ReciprocalLattice(const xstructure& a); //CO20180409
xmatrix<double> ReciprocalLattice(const xmatrix<double>& rlattice,double scale=1.0); //CO20180409
string KPPRA(int& k1,int& k2,int& k3,const xmatrix<double>& rlattice,const int& NK);
string KPPRA(xstructure& str,const int& _NK);
string KPPRA_DELTA(int& k1,int& k2,int& k3,const xmatrix<double>& rlattice,const double& DK);
string KPPRA_DELTA(xstructure& str,const double& DK);
int GetNBANDS(int electrons,int nions,int spineach,bool ispin);
double GetZVAL(const stringstream& sss,vector<double>& vZVAL);
double GetZVAL(const _xvasp& xvasp,vector<double>& vZVAL);
double GetZVAL(const string& directory,vector<double>& vZVAL);
double GetCellAtomZVAL(const stringstream& sss,vector<double>& vZVAL,const stringstream& sstr,vector<double>& sZVAL,string mode);  // sss sstr returns ZVAL cell, VAL and sZVAL
double GetCellAtomZVAL(const string& directory,vector<double>& vZVAL,vector<double>& sZVAL,string mode);  // from directory POT/POS returns total ZVAL cell, vZVAL and sZVAL
double GetPOMASS(const stringstream& sss,vector<double>& vPOMASS);
double GetPOMASS(const _xvasp& xvasp,vector<double>& vPOMASS);
double GetPOMASS(const string& directory,vector<double>& vPOMASS);
double GetCellAtomPOMASS(const stringstream& sss,vector<double>& vPOMASS,const stringstream& sstr,vector<double>& sPOMASS,string mode);  // sss sstr returns POMASS cell, VAL and sPOMASS
double GetCellAtomPOMASS(const string& directory,vector<double>& vPOMASS,vector<double>& sPOMASS,string mode);  // from directory POT/POS returns total POMASS cell, vPOMASS and sPOMASS
double GetVol(const xmatrix<double>& lat);
double det(const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3);
double GetVol(const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3);
double det(const double&,const double&,const double&,const double&,const double&,const double&,const double&,const double&,const double&);
//double getcos(const xvector<double>& a,const xvector<double>& b);  // removed and put in aurostd_xvector.h as cos(xvector,xvector) and sin(xvector,xvector)
xvector<double> Getabc_angles(const xmatrix<double>& lat,const int& mode);
xvector<long double> Getabc_angles(const xmatrix<long double>& lat,const int& mode);
xvector<double> Getabc_angles(const xmatrix<double>& lat,const xvector<int>& permut,const int& mode);
xvector<double> Getabc_angles(const xvector<double>& r1,const xvector<double>& r2,const xvector<double>& r3,const int& mode);
xvector<double> Getabc_angles(const xvector<double>& r1,const xvector<double>& r2,const xvector<double>& r3,const xvector<int>& permut,const int& mode);
#define _Getabc_angles Getabc_angles
//#define _Getabc_angles __NO_USE_Sortabc_angles
xvector<double> Sortabc_angles(const xmatrix<double>& lat,const int& mode);
xmatrix<double> GetClat(const xvector<double>& abc_angles);
xmatrix<double> GetClat(const double &a,const double &b,const double &c,const double &alpha,const double &beta,const double &gamma);
xstructure GetIntpolStr(xstructure strA,xstructure strB,const double& f,const string& path_flag);
double RadiusSphereLattice(const xmatrix<double>& lattice,double scale=1.0); //CO20180409
xvector<int> LatticeDimensionSphere(const xmatrix<double>& lattice,double radius,double scale=1.0); //CO20180409
xvector<int> LatticeDimensionSphere(const xstructure& str,double radius);
void resetLatticeDimensions(const xmatrix<double>& lattice, double radius, xvector<int>& dims,
    vector<xvector<double> >& l1, vector<xvector<double> >& l2, 
    vector<xvector<double> >& l3, vector<int>& a_index, 
    vector<int>& b_index, vector<int>& c_index); //DX20191122
xvector<double> F2C(const double& scale,const xmatrix<double>& lattice,const xvector<double>& fpos);    // fpos are F components per COLUMS !
xvector<double> F2C(const xmatrix<double>& lattice,const xvector<double>& fpos);                        // fpos are F components per COLUMS !
xvector<double> C2F(const double& scale,const xmatrix<double>& lattice,const xvector<double>& cpos);    // cpos are C components per COLUMS !
xvector<double> C2F(const xmatrix<double>& lattice,const xvector<double>& cpos);                        // cpos are C components per COLUMS !
xmatrix<double> F2C(const double& scale,const xmatrix<double>& lattice,const xmatrix<double>& fpos);    // fpos are F components per COLUMS !
xmatrix<double> F2C(const xmatrix<double>& lattice,const xmatrix<double>& fpos);                        // fpos are F components per COLUMS !
xmatrix<double> C2F(const double& scale,const xmatrix<double>& lattice,const xmatrix<double>& cpos);    // cpos are C components per COLUMS !
xmatrix<double> C2F(const xmatrix<double>& lattice,const xmatrix<double>& cpos);                        // cpos are C components per COLUMS !
_atom F2C(const double& scale,const xmatrix<double>& lattice,const _atom& iatom);                       // atom.fpos are F components per COLUMS !
_atom F2C(const xstructure& str,const _atom& iatom);                                                    // atom.fpos are F components per COLUMS !
_atom C2F(const double& scale,const xmatrix<double>& lattice,const _atom& iatom);                       // atom.cpos are C components per COLUMS !
_atom C2F(const xmatrix<double>& lattice,const _atom& iatom);                                           // atom.cpos are C components per COLUMS !
_atom F2C(const double& scale,const xstructure& str,const _atom& iatom);                                // atom.fpos are F components per COLUMS !
_atom F2C(const xstructure& str,const _atom& iatom);                                                    // atom.fpos are F components per COLUMS !
_atom C2F(const double& scale,const xstructure& str,const _atom& iatom);                                // atom.fpos are F components per COLUMS !
_atom C2F(const xstructure& str,const _atom& iatom);                                                    // atom.cpos are C components per COLUMS !
xmatrix<double> FF2CC(const double& scale,const xmatrix<double>& lattice,const xmatrix<double>& fmat);  // fmat is an operation in F coordinates
xmatrix<double> FF2CC(const xmatrix<double>& lattice,const xmatrix<double>& fmat);                      // fmat is an operation in F coordinates
xmatrix<double> CC2FF(const double& scale,const xmatrix<double>& lattice,const xmatrix<double>& cmat);  // cmat is an operation in C coordinates
xmatrix<double> CC2FF(const xmatrix<double>& lattice,const xmatrix<double>& cmat);                      // cmat is an operation in C coordinates
//DX20190905 START
//BringInCellInPlace() overloads
void BringInCellInPlace(double&, double=_ZERO_TOL_, double=1.0, double=0.0);  //ME+DX20190409
void BringInCellInPlace(xvector<double>&, double=_ZERO_TOL_, double=1.0, double=0.0);  //ME+DX20190409
void BringInCellInPlace(_atom& atom_in, const xmatrix<double>& lattice, double tolerance=_ZERO_TOL_, double upper_bound=1.0, double lower_bound=0.0); //DX20190904
void BringInCellInPlace(xstructure& xstr, double tolerance=_ZERO_TOL_, double upper_bound=1.0, double lower_bound=0.0); //DX20190904

//BringInCell() overloads
double BringInCell(double, double=_ZERO_TOL_, double=1.0, double=0.0);  //ME+DX20190409
xvector<double> BringInCell(const xvector<double>& fpos_in, double tolerance=_ZERO_TOL_, double upper_bound=1.0, double lower_bound=0.0); //DX20190904
_atom BringInCell(const _atom& atom_in, const xmatrix<double>& lattice, double tolerance=_ZERO_TOL_, double upper_bound=1.0, double lower_bound=0.0); //DX20190904
xstructure BringInCell(const xstructure& xstr_in, double tolerance=_ZERO_TOL_, double upper_bound=1.0, double lower_bound=0.0); //DX20190904

//BringInCellFPOS overloads
void BringInCellInPlaceFPOS(_atom& atom_in, double tolerance=_ZERO_TOL_, double upper_bound=1.0, double lower_bound=0.0); //DX20190904
_atom BringInCellFPOS(const _atom& atom_in, double tolerance=_ZERO_TOL_, double upper_bound=1.0, double lower_bound=0.0); //DX20190904
//DX20190905 END
//DX+CO START
//DX20190905 [OBSOLETE] double BringInCell(const double& x);
//DX20190905 [OBSOLETE] double BringInCell_20161115(const double& x);
//DX20190905 [OBSOLETE] double BringInCell_20160101(const double& x);
//DX20190905 [OBSOLETE] double BringInCell(const double& x);
//DX20190905 [OBSOLETE] double BringInCell_20160101(const double& x);
//DX20190905 [OBSOLETE] xvector<double> BringInCell(const xvector<double>& v_in,double epsilon);
//DX20190905 [OBSOLETE] xvector<double> BringInCell_20161115(const xvector<double>& v_in,double epsilon);
//DX20190905 [OBSOLETE] xvector<double> BringInCell_20160101(const xvector<double>& v_in,double epsilon);
//DX20190905 [OBSOLETE] xvector<double> BringInCell(const xvector<double>& v_in);
//DX20190905 [OBSOLETE] xvector<double> BringInCell2(const xvector<double>& v_in);
//DX20190905 [OBSOLETE] xvector<double> BringInCell2_20161115(const xvector<double>& v_in);
//DX20190905 [OBSOLETE] xvector<double> BringInCell2_20160101(const xvector<double>& v_in, double tolerance);
//DX+CO END
xstructure IdenticalAtoms(const xstructure& a);                                // Make identical atoms
//xstructure SwapSpecies(const xstructure& a,const uint& A,const uint& B);       // Permute Species A with B (safe for species C).
//xstructure SwapCoordinates(const xstructure& str,const uint& i,const uint& j); // Permute Coordinates i with j
//string SpeciesLabel(const xstructure& a,const uint& A);                        // Returns the Label of the specie A (if available)
//string SpeciesString(const xstructure& a);                                           // Gives a string with the list of all the species
bool GetNiggliCell(const xmatrix<double>& in_lat,xmatrix<double>& niggli_lat,xmatrix<double>& P,xmatrix<double>& Q);
bool GetNiggliCell_20180213(const xmatrix<double>& in_lat,xmatrix<double>& niggli_lat,xmatrix<double>& P,xmatrix<double>& Q); //DX20180213 - new dated function
bool GetNiggliCell_20180101(const xmatrix<double>& in_lat,xmatrix<double>& niggli_lat,xmatrix<double>& P,xmatrix<double>& Q); //DX20180213 - old dated function
// standard lattice reduction and type
string GetLatticeType(xmatrix<double> lattice);
string GetLatticeType(xvector<double> data);
xstructure Standard_Primitive_UnitCellForm(const xstructure& a);
xstructure GetStandardPrimitive(const xstructure& a);
xmatrix<double> GetStandardPrimitive(xmatrix<double> lattice);
xvector<double> GetStandardPrimitive(xvector<double> data);
xstructure Standard_Conventional_UnitCellForm(const xstructure& a);
xstructure GetStandardConventional(const xstructure& a);
xmatrix<double> GetStandardConventional(xmatrix<double> lattice);
xvector<double> GetStandardConventional(xvector<double> data);
// niggli
xstructure GetNiggliStr(const xstructure& in_str);
xmatrix<double> GetNiggliStr(const xmatrix<double>& lattice);
xstructure NiggliUnitCellForm(const xstructure& a);
xmatrix<double> NiggliUnitCellForm(const xmatrix<double>& lattice);
// minkowsky
xstructure MinkowskiBasisReduction(const xstructure& a);
xmatrix<double> MinkowskiBasisReduction(const xmatrix<double>& lattice);
// optimal lattice reduction
xstructure LatticeReduction(const xstructure& a);
xmatrix<double> LatticeReduction(const xmatrix<double>& lattice);
//CO20170807 START
//DX20190214 [OBSOLETE] deque<_atom> foldAtomsInCell(deque<_atom>& atoms, xmatrix<double>& c2f_new, xmatrix<double>& f2c_new, bool skew); //CO20190520 - removed pointers for bools and doubles, added const where possible
deque<_atom> foldAtomsInCell(const xstructure& a, const xmatrix<double>& lattice_new, bool skew, double tol, bool check_min_dists=true); //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - added check_min_dists bool
deque<_atom> foldAtomsInCell(const deque<_atom>& atoms, const xmatrix<double>& lattice_orig, const xmatrix<double>& lattice_new, bool skew, double tol, bool check_min_dists=true); //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 = added check_min_dists bool
xstructure GetPrimitiveVASP(const xstructure& a);
xstructure GetPrimitiveVASP(const xstructure& a,double tol);
//CO20170807 STOP
// bring cell in,compact, wigner seitz
//DX20190905 [OBSOLETE] _atom BringInCell(const _atom& atom_in,const xmatrix<double>& lattice,double epsilon);
//DX20190905 [OBSOLETE] //DX+CO START
//DX20190905 [OBSOLETE] _atom BringInCell_20161115(const _atom& atom_in,const xmatrix<double>& lattice,double epsilon); //DX
//DX20190905 [OBSOLETE] _atom BringInCell_20160101(const _atom& atom_in,const xmatrix<double>& lattice,double epsilon); //DX
//DX20190905 [OBSOLETE] _atom BringInCell(const _atom& atom_in,const xmatrix<double>& lattice);
//DX20190905 [OBSOLETE] _atom BringInCell_20161115(const _atom& atom_in,const xmatrix<double>& lattice); //DX
//DX20190905 [OBSOLETE] _atom BringInCell_20160101(const _atom& atom_in,const xmatrix<double>& lattice); //DX
//DX20190905 [OBSOLETE] xstructure BringInCell(const xstructure& a,double epsilon);
//DX20190905 [OBSOLETE] xstructure BringInCell_20161115(const xstructure& a,double epsilon); //DX
//DX20190905 [OBSOLETE] xstructure BringInCell_20160101(const xstructure& a,double epsilon); //DX
//DX20190905 [OBSOLETE] xstructure BringInCell(const xstructure& a);
//DX20190905 [OBSOLETE] xstructure BringInCell_20161115(const xstructure& a); //DX
//DX20190905 [OBSOLETE] xstructure BringInCell_20160101(const xstructure& a); //DX
//DX20190905 [OBSOLETE] //DX+CO END
xstructure BringInCompact(const xstructure& a);
xstructure BringInWignerSeitz(const xstructure& a);
// primitive stuff
xstructure GetPrimitive(const xstructure& a);
xstructure GetPrimitive(const xstructure& a,double tol);
xstructure GetPrimitive1(const xstructure& a);
xstructure GetPrimitive2(const xstructure& a);
xstructure GetPrimitive3(const xstructure& a);
bool IsTranslationFVector(const xstructure& a,const xvector<double>& ftvec);
bool IsTranslationCVector(const xstructure& a,const xvector<double>& ctvec);
// other eggs
xstructure ReScale(const xstructure& a,const double& in_scale);
xstructure SetScale(const xstructure& a,const double& in_scale);
xstructure SetVolume(const xstructure& a,const double& in_volume);
xstructure InflateLattice(const xstructure& a,const double& coefficient);
xstructure InflateVolume(const xstructure& a,const double& coefficient);
double GetVolume(const xstructure& a);
double Volume(const xstructure& a);
//DX20180726 START
_atom BringCloseToOrigin(_atom& atom, xmatrix<double>& f2c);
bool uniqueAtomInCell(_atom& atom, deque<_atom>& atoms);
bool alreadyInCell(_atom& atom, deque<_atom> atoms);
//DX20180726 END
//DX+CO START
bool atomInCell(const _atom& atom, double tolerance=_ZERO_TOL_); //DX20191125
bool inCell(const xvector<double>& pos_vec, double tolerance=_ZERO_TOL_); //DX20191125 - added tolerance
//DX+CO END
xstructure GetSuperCell(const xstructure& a,const xmatrix<double>& sc);
xstructure GetSuperCell(const xstructure& a,const xvector<double>& sc);
xstructure GetSuperCell(const xstructure& a,const xvector<int>& sc);
xstructure GetSuperCell(const xstructure& a, const int& sc11,const int& sc12,const int& sc13, const int& sc21,const int& sc22,const int& sc23, const int& sc31,const int& sc32,const int& sc33);
xstructure GetSuperCell(const xstructure& a,const int& sc1,const int& sc2,const int& sc3);
//CO START
xstructure GetSuperCell(const xstructure& a,const xmatrix<double>& sc,vector<int>& sc2pcMap,vector<int>& pc2scMap,bool get_symmetry, bool get_full_basis, bool force_supercell_matrix=false,bool force_strict_pc2scMap=false); //DX20190319 - added force_supercell_matrix //CO20190409 - added force_strict_pc2scMap
xstructure GetSuperCell(const xstructure& a,const xvector<double>& sc,vector<int>& sc2pcMap,vector<int>& pc2scMap,bool get_symmetry, bool get_full_basis, bool force_supercell_matrix=false,bool force_strict_pc2scMap=false); //DX20190319 - added force_supercell_matrix //CO20190409 - added force_strict_pc2scMap
xstructure GetSuperCell(const xstructure& a,const xvector<int>& sc,vector<int>& sc2pcMap,vector<int>& pc2scMap,bool get_symmetry, bool get_full_basis, bool force_supercell_matrix=false,bool force_strict_pc2scMap=false); //DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
xstructure GetSuperCell(const xstructure& a,const int& sc11,const int& sc12,const int& sc13, const int& sc21,const int& sc22,const int& sc23, const int& sc31,const int& sc32,const int& sc33,vector<int>& sc2pcMap,vector<int>& pc2scMap,bool get_symmetry, bool get_full_basis, bool force_supercell_matrix=false,bool force_strict_pc2scMap=false); //DX20190319 - added force_supercell_matrix //CO20190409 - added force_strict_pc2scMap
xstructure GetSuperCell(const xstructure& a,const int& sc1,const int& sc2,const int& sc3,vector<int>& sc2pcMap,vector<int>& pc2scMap,bool get_symmetry, bool get_full_basis, bool force_supercell_matrix=false,bool force_strict_pc2scMap=false); //DX20190319 - added force_supercell_matrix  //CO20190409 - added force_strict_pc2scMap
//CO END
bool CalculateSymmetry(xstructure& str,bool ossverbose,ostream& oss,bool fffverbose,double radius);
bool CalculateSymmetry(xstructure& str,bool ossverbose,ostream& oss,bool fffverbose);
bool CalculateSymmetry(xstructure& str,bool ossverbose,ostream& oss,double radius);
bool CalculateSymmetry(xstructure& str,bool ossverbose,double radius);
bool CalculateSymmetry(xstructure& str,double radius);
bool CalculateSymmetry(xstructure& str,bool ossverbose);
bool CalculateSymmetry(xstructure& str);
void CalculateSymmetryPointGroup(xstructure& str,bool ossverbose,ostream& oss,bool fffverbose);
void CalculateSymmetryPointGroup(xstructure& str,bool ossverbose,ostream& oss);
void CalculateSymmetryPointGroup(xstructure& str,bool ossverbose);
void CalculateSymmetryPointGroup(xstructure& str);
void CalculateSymmetryPointGroupCrystal(xstructure& str,bool ossverbose,ostream& oss,bool fffverbose);
void CalculateSymmetryPointGroupCrystal(xstructure& str,bool ossverbose,ostream& oss);
void CalculateSymmetryPointGroupCrystal(xstructure& str,bool ossverbose);
void CalculateSymmetryPointGroupCrystal(xstructure& str);
void CalculateSymmetryFactorGroup(xstructure& str,bool ossverbose,ostream& oss,bool fffverbose);
void CalculateSymmetryFactorGroup(xstructure& str,bool ossverbose,ostream& oss);
void CalculateSymmetryFactorGroup(xstructure& str,bool ossverbose);
void CalculateSymmetryFactorGroup(xstructure& str);
void CalculateSymmetryPointGroupKLattice(xstructure& str,bool ossverbose,ostream& oss,bool fffverbose);
void CalculateSymmetryPointGroupKLattice(xstructure& str,bool ossverbose,ostream& oss);
void CalculateSymmetryPointGroupKLattice(xstructure& str,bool ossverbose);
void CalculateSymmetryPointGroupKLattice(xstructure& str);
void CalculateSymmetryPointGroupKCrystal(xstructure& str,bool ossverbose,ostream& oss,bool fffverbose);  //ME20200114
void CalculateSymmetryPointGroupKCrystal(xstructure& str,bool ossverbose,ostream& oss);  //ME20200114
void CalculateSymmetryPointGroupKCrystal(xstructure& str,bool ossverbose);  //ME20200114
void CalculateSymmetryPointGroupKCrystal(xstructure& str);  //ME20200114
void CalculateSymmetryPointGroupKPatterson(xstructure& str,bool ossverbose,ostream& oss,bool fffverbose);  //ME20200129
void CalculateSymmetryPointGroupKPatterson(xstructure& str,bool ossverbose,ostream& oss);  //ME20200129
void CalculateSymmetryPointGroupKPatterson(xstructure& str,bool ossverbose);  //ME20200129
void CalculateSymmetryPointGroupKPatterson(xstructure& str);  //ME20200129
xstructure Rotate(const xstructure& a,const xmatrix<double>& rm);
xstructure GetLTCell(const xmatrix<double>& lt,const xstructure& str);
xstructure GetLTFVCell(const xvector<double>& nvec,const double phi,const xstructure& str);
xstructure ShiftPos(const xstructure& a,const xvector<double>& shift,const int& flag);
xstructure ShiftCPos(const xstructure& a,const xvector<double>& shift);
xstructure ShiftFPos(const xstructure& a,const xvector<double>& shift);
double MaxStructureLattice(const xstructure& str);
double MinStructureLattice(const xstructure& str);
double AtomDist(const xstructure& str,const _atom& atom1,const _atom& atom2);
bool SameAtom(const xstructure& str,const _atom& atom1,const _atom& atom2);
bool SameAtom(const _atom& atom1,const _atom& atom2);
bool DifferentAtom(const xstructure& str,const _atom& atom1,const _atom& atom2);
xmatrix<double> GetDistMatrix(const xstructure& a); //CO20171025
vector<double> GetNBONDXX(const xstructure& a);
int GenerateGridAtoms(xstructure& str,int i1,int i2,int j1,int j2,int k1,int k2);
int GenerateGridAtoms_20190520(xstructure& str,int i1,int i2,int j1,int j2,int k1,int k2); //DX20191218 [ORIG]
int GenerateGridAtoms_20191218(xstructure& str,int i1,int i2,int j1,int j2,int k1,int k2); //DX20191218 [NEW]
int GenerateGridAtoms(xstructure& str,int d1,int d2,int d3);
int GenerateGridAtoms(xstructure& str,int d);
int GenerateGridAtoms(xstructure& str,const xvector<int>& dims);
int GenerateGridAtoms(xstructure& str);
int GenerateGridAtoms(xstructure& str,const double& radius);

void l2ijk(const xstructure& str,const int &l,int &i,int &j,int &k);
void l2ijk(const xstructure& str,const int &l,xvector<int>& ijk);
xvector<int> l2ijk(const xstructure& str,const int &l);
void ijk2l(const xstructure& str,int &l,const int &i,const int &j,const int &k);
void ijk2l(const xstructure& str,int &l,const xvector<int>& ijk);
int ijk2l(const xstructure& str,const int &i,const int &j,const int &k);
int ijk2l(const xstructure& str,const xvector<int>& ijk);
xvector<double> r_lattice(const xstructure& str,const int &l);
xvector<double> r_lattice(const xstructure& str,const int &i,const int &j,const int &k);
xvector<double> r_lattice(const xstructure& str,const xvector<int>& ijk);
xstructure input2AIMSxstr(istream& input);
xstructure input2ABINITxstr(istream& input);
xstructure input2QExstr(istream& input);
xstructure input2VASPxstr(istream& input);
xstructure input2ELKxstr(istream& input); //DX20200313

// ----------------------------------------------------------------------------
// functions related to AtomEnvironment - DX20191122
vector<AtomEnvironment> getAtomEnvironments(const xstructure& xstr, uint mode=ATOM_ENVIRONMENT_MODE_1);
vector<AtomEnvironment> getLFAAtomEnvironments(const xstructure& xstr, const string& lfa, const vector<string>& LFAs, uint mode=ATOM_ENVIRONMENT_MODE_1);
void minimumCoordinationShellLatticeOnly(const xmatrix<double>& lattice,
    double& min_dist, uint& frequency, vector<xvector<double> >& coordinates); //DX20191122
void minimumCoordinationShellLatticeOnly(const xmatrix<double>& lattice,
    double& min_dist, uint& frequency, vector<xvector<double> >& coordinates, double radius); //DX20191122
void minimumCoordinationShellLatticeOnly(const xmatrix<double>& lattice, xvector<int>& dims,
    vector<xvector<double> >& l1, vector<xvector<double> >& l2, vector<xvector<double> >& l3, 
    vector<int>& a_index, vector<int>& b_index, vector<int>& c_index, 
    double& min_dist, uint& frequency, vector<xvector<double> >& coordinates,
    double radius); //DX20191122
void minimumCoordinationShell(const xstructure& xstr, uint center_index, 
    double& min_dist, uint& frequency, vector<xvector<double> >& coordinates); //DX20191122
void minimumCoordinationShell(const xstructure& xstr, uint center_index, 
    double& min_dist, uint& frequency, vector<xvector<double> >& coordinates, const string& type); //DX20191122

//makefile tests
bool CeramGenTest(ostream& oss=cout);
bool CeramGenTest(ofstream& FileMESSAGE,ostream& oss=cout);
bool EgapTest(ostream& oss=cout);
bool EgapTest(ofstream& FileMESSAGE,ostream& oss=cout);
bool gcdTest(ostream& oss=cout);
bool gcdTest(ofstream& FileMESSAGE,ostream& oss=cout);
bool smithTest(ostream& oss=cout);
bool smithTest(ofstream& FileMESSAGE,ostream& oss=cout);

// ----------------------------------------------------------------------------
// Structure Prototypes
// aflow_xproto.cpp
#define _HTQC_PROJECT_STRING_ "HTQC Project"
#define _TERNARY_PROJECT_STRING_ "HTQC^3 Project"
#define _ICSD_STRING_ "(icsd library)"
#define _ICSD_PROJECT_STRING_ "ICSD Project"
#define _ICSD_AFLOWLIB_STRING_ "(icsd_aflowlib library)"

// aflow_xproto.cpp
namespace aflowlib {
  string PrototypeCleanLatticeString(const string& latticeIN);
}
double NearestNeighbour(const xstructure &str_in);

// for HTQC
#define STRUCTURE_MODE_NONE             0
#define STRUCTURE_MODE_RAW              1
#define STRUCTURE_MODE_ABC              2
#define STRUCTURE_MODE_WYC              3
#define STRUCTURE_MODE_ICSD             4
#define STRUCTURE_MODE_HTQC_ICSD        5
#define STRUCTURE_MODE_USE              6
#define STRUCTURE_MODE_REMOVE           7
#define STRUCTURE_MODE_SPECIES          8
#define STRUCTURE_MODE_SWAP_AB          9
#define STRUCTURE_MODE_SWAP_BC         10
#define STRUCTURE_MODE_SWAP_AC         11
#define STRUCTURE_MODE_SWAP_XY         12
#define STRUCTURE_MODE_PRIM            13
#define STRUCTURE_MODE_CONVENTIONAL    14
#define STRUCTURE_MODE_VOLUME          15
#define LIBRARY_MODE_ICSD               0
#define LIBRARY_MODE_ICSD_AFLOWLIB      1
#define LIBRARY_MODE_HTQC               2
#define LIBRARY_MODE_HTQC_ICSD          3
#define LIBRARY_MODE_HTQC_ICSD_AFLOWLIB 4
#define LIBRARY_MODE_LIB0               5
#define LIBRARY_MODE_LIB3               6
#define LIBRARY_MODE_LIB4               7
#define LIBRARY_MODE_LIB5               8
#define LIBRARY_MODE_LIB6               9
#define LIBRARY_MODE_LIB7               10
#define LIBRARY_MODE_LIB8               11
#define LIBRARY_MODE_LIB9               12
#define LIBRARY_MODE_PROTOTYPE          13
#define LIBRARY_MODE_XSTRUCTURE         14
#define LIBRARY_MODE_AUID               15
#define LIBRARY_MODE_ARUN               16  //ME20181226
string* LOAD_Library_ICSD(string file);

namespace aflowlib {
  struct _PROTO_PARAMS{
    string label;
    string parameters;
    deque<string> vatomX;
    deque<double> vvolumeX;
    double volume_in;
    int mode;
    bool flip_option;
  };

  xstructure PrototypePure(ostream &FileMESSAGE,string label,string parameters,string atA,double volA);
  xstructure PrototypePure(ostream &FileMESSAGE,string label,string parameters,string atA);
  xstructure PrototypePure(ostream &FileMESSAGE,string label,string parameters);
  xstructure PrototypePureHTQC(ostream &FileMESSAGE,string label,string parameters,string atA,double volA);
  xstructure PrototypePureHTQC(ostream &FileMESSAGE,string label,string parameters,string atA);
  xstructure PrototypePureHTQC(ostream &FileMESSAGE,string label,string parameters);
  uint PrototypeLibrariesSpeciesNumber(const string& label,ostream& oss=cerr); //CO20181226
  // xstructure PrototypeLibraries(ostream &oss,string label,string parameters,int mode=LIBRARY_MODE_HTQC);
  // xstructure PrototypeLibraries(ostream &oss,string label,string parameters,deque<string> &vatomX,int mode=LIBRARY_MODE_HTQC);
  uint GetAllPrototypeLabels(vector<string>& prototype_labels, string library="all"); //DX20181009
  uint GetAllPrototypeLabels(vector<string>& prototype_labels, vector<string>& compositions, string library="all"); //DX20181009
  uint GetAllPrototypeLabels(vector<string>& prototype_labels, vector<string>& compositions, 
      vector<uint>& space_group_numbers, vector<vector<vector<string> > >& Wyckoff_letter_strings, 
      string library="all");
  vector<string> GetPrototypesBySpeciesNumber(uint number_of_species, string library="all"); //DX20181009
  vector<string> GetPrototypesByStoichiometry(vector<uint> stoichiometry, string library="all"); //DX20181009
  vector<string> GetPrototypesByStoichiometry(vector<uint> stoichiometry, vector<string>& protototype_composition, vector<uint>& prototype_space_group_numbers, 
      vector<vector<vector<string> > >& prototype_grouped_Wyckoff_letters, string library="all");
  vector<string> GetPrototypesBySymmetry(vector<uint>& stoichiometry, uint& space_group_number, vector<GroupedWyckoffPosition>& grouped_Wyckoff_positions, uint setting, string library="all"); //DX20181010
  //vector<string> GetPrototypesBySymmetry(vector<uint>& stoichiometry, uint& space_group_number, vector<vector<vector<string> > >& grouped_possible_Wyckoff_letters, uint setting, string library="all"); //DX20181010
  vector<string> GetPrototypesBySymmetry(vector<uint>& stoichiometry, uint& space_group_number, vector<GroupedWyckoffPosition>& grouped_Wyckoff_positions, vector<uint>& prototype_space_groups, uint setting, string library="all"); //DX20181010
  //vector<string> GetPrototypesBySymmetry(vector<uint>& stoichiometry, uint& space_group_number, vector<vector<vector<string> > >& grouped_possible_Wyckoff_letters, vector<uint>& prototype_space_groups, uint setting, string library="all"); //DX20181010
  xstructure PrototypeLibraries(ostream &oss,string label,string parameters,int mode);
  xstructure PrototypeLibraries(ostream &oss,string label,string parameters,deque<string> &vatomX,int mode);
  xstructure PrototypeLibraries(ostream &oss,string label,string parameters,deque<string> &vatomX,deque<double> &vvolumeX,double volume_in,int mode);//=LIBRARY_MODE_HTQC);
  xstructure PrototypeLibraries(ostream &oss,string label,string parameters,deque<string> &vatomX,deque<double> &vvolumeX,double volume_in,int mode,bool flip_option);
  xstructure PrototypeLibraries(ostream &oss,_PROTO_PARAMS *PARAMS);

  string PrototypesHelp(void);
  string PrototypesIcsdHelp(string options);
  string CALCULATED(string options);
  string CALCULATED_ICSD_RANDOM(void);
  // aflow_xproto_gus.cpp
  xstructure PrototypeBinaryGUS(ostream &oss,string label);
  xstructure PrototypeBinaryGUS(ostream &oss,string label,string atA,string atB);
  xstructure PrototypeBinaryGUS(ostream &oss,string label,string atA,double volA,string atB,double volB,double vol_in);
}

extern string PrototypeBinaryGUS_Cache_Library[];

// ----------------------------------------------------------------------------
// aflow_anrl.cpp
//DX20180710 - updated - #define DOI_ANRL " [ANRL doi: arXiv:1607.02532]"
#define DOI_ANRL " [ANRL doi: 10.1016/j.commatsci.2017.01.017 (part 1), doi: 10.1016/j.commatsci.2018.10.043 (part 2)]" //DX20180710 - updated //DX20190214 updated part 2 doi
#define DOI_POCC " [POCC doi: 10.1021/acs.chemmater.6b01449]"

namespace anrl {
  xstructure PrototypeANRL(ostream &oss,string label,string parameters,deque<string> &vatomX,deque<double> &vvolumeX,double volume_in,int mode,bool flip_option);
  uint PrototypeANRL_LoadList(vector<string>& vproto,
      vector<string>& vproto_label,
      vector<uint>& vproto_nspecies,
      vector<uint>& vproto_natoms,
      vector<uint>& vproto_spacegroup,
      vector<uint>& vproto_nunderscores,
      vector<uint>& vproto_nparameters,
      vector<string>& vproto_Pearson_symbol,
      vector<string>& vproto_params,
      vector<string>& vproto_Strukturbericht,
      vector<string>& vproto_prototype,
      vector<string>& vproto_dialect);
  vector<string> getANRLParameters(string anrl_label, string library="", int choice=-1, bool keep_original_lattice_parameter=false); //DX20181009 //DX20190227 - added keep_original_lattice_parameter
  bool vproto2tokens(string proto,
      string& label,
      uint& nspecies,
      uint& natoms,
      uint& spacegroup,
      uint& nunderscores,
      uint& nparameters,
      string& Pearson_symbol,
      string& params,
      string& Strukturbericht,
      string& prototype,
      string& dialect);
  vector<uint> extractStoichiometry(string& anrl_label);
  bool PrototypeANRL_Consistency(uint vparameters_size,uint proto_nparameters,string proto_prototype,
      string proto_label,string proto_Strukturbericht,string proto_Pearson_symbol,
      uint proto_spacegroup, string proto_params, uint print_mode); //DX20180710 - added print_mode //DX20200207 - oss no longer needed
  string groupedWyckoffPosition2ANRLString(const vector<GroupedWyckoffPosition>& grouped_positions, bool alphabetize);
  vector<string> getANRLLatticeParameterString(char& lattice_type);
  vector<double> getANRLLatticeParameterValuesFromWyccar(const vector<string>& wyccar_ITC, char lattice_type, char lattice_centering, uint setting); //DX20191031
  vector<double> getANRLLatticeParameterValuesFromABCAngles(const xstructure& xstr, char lattice_type, char lattice_centering, uint setting); //DX20191031
  vector<double> getANRLLatticeParameterValues(const vector<double>& all_lattice_parameters, char lattice_type, char lattice_centering, uint setting); //DX20191031
  uint getANRLSettingChoice(int spacegroup); //DX20191031 - removed reference
  string structure2anrl(istream& input, aurostd::xoption& vpflow);           // xoption
  string structure2anrl(xstructure& xstr, bool recalculate_symmetry=true);   // use default options //DX20191031 - added recalculate_symmetry
  string structure2anrl(xstructure& xstr, double tolerance);                 // specify symmetry tolerance //CO20190520 - removed pointers for bools and doubles, added const where possible
  string structure2anrl(xstructure& xstr, uint setting);                     // specify setting
  string structure2anrl(xstructure& xstr, double tolerance, uint setting, bool recalculate_symmetry=true);  // main function //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190829 - added recalculate_symmetry //DX20191031 - removed reference
  xstructure rhl2hex(xstructure& str, double& a, double& c); 
}

// ----------------------------------------------------------------------------
// Various prototypes to be moved somewhere sometime
// PROTOTYPES
// uint argsprint(vector<string> argv);
// ----------------------------------------------------------------------------
// aflow.cpp
//[CO20200502 - DUPLICATE?]string aflow_get_time_string(void);
//[CO20200502 - DUPLICATE?]string aflow_get_time_string_short(void);
//[CO20200502 - DUPLICATE?]string strPID(void);
//[CO20200502 - DUPLICATE?]string strTID(void);  //CO20200502 - threadID
int AFLOW_main(vector<string> &argv);
namespace aflow {
  string License_Preamble_aflow(void);
  string Intro_aflow(string x);
  string Intro_sflow(string x);
  string Intro_HELP(string x);
  string Banner(string type);
}
int VASP_Main(vector<string> argv);
int GRND_Main(vector<string> argv);
namespace KBIN {
  int KBIN_Main(vector<string> argv);
}
string MessageTime(void);
string MessageHostTime(const _aflags& aflags);
string MessageDir(const _aflags& aflags);
string MessageDirTime(const _aflags& aflags);
string MessageDirHostTime(const _aflags& aflags);
//[CO20200624 - REDUNDANT]bool AFLOW_BlackList(string hostname);

// ----------------------------------------------------------------------------
// aflow_pthreads.cpp
namespace AFLOW_PTHREADS {
  int GetTotalCPUs(void);
  bool Check_Threads(vector<string> argv,const bool& VERBOSE);
  void Clean_Threads(void);
  void No_Threads(void);
  bool Available_Free_Threads(int &fthread);
  bool Wait_Available_Free_Threads(int &fthread,const double& pthread_wait,const bool& VERBOSE);
  bool Wait_Available_Free_Threads(int &fthread,const bool& VERBOSE);
}
// interfaces
namespace KBIN {
  void RUN_Directory_PTHREADS(_aflags &aflags);
  void *_threaded_interface_RUN_Directory(void *ptr);
} // namespace KBIN
namespace aurostd { // Multithreaded add on to aurostd
  bool multithread_execute(deque<string> vcommand,int NUM_THREADS,bool VERBOSE);
  bool multithread_execute(deque<string> vcommand,int NUM_THREADS);
  bool multithread_execute(deque<string> vcommand);
  bool multithread_execute(vector<string> vcommand,int NUM_THREADS,bool VERBOSE);
  bool multithread_execute(vector<string> vcommand,int NUM_THREADS);
  bool multithread_execute(vector<string> vcommand);
} // namespace aurostd
namespace AFLOW_PTHREADS {
  bool MULTI_sh(vector<string> argv);
  bool MULTI_compress(string cmd,vector<string> argv);
  bool MULTI_zip(vector<string> argv);
  bool MULTI_bz2xz(vector<string> argv);bool MULTI_xz2bz2(vector<string> argv);
  bool MULTI_gz2xz(vector<string> argv);
}
namespace sflow {
  void KILL(string options);
  void JUST(string options,istream& input,string mode);
  void QSUB(string options);
  void QSUB(string options,string cmd);
  void QDEL(string options);
  void QDEL(string options,string cmd);
}
vector<vector<int> > getThreadDistribution(const int&, const int&);  //ME20190218

// ----------------------------------------------------------------------------
// aflow_kbin.cpp
//int KbinCheckInputFiles(string Directory,ofstream& FileERROR);
namespace KBIN {
  void MPI_Extract(string AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags);
  void RUN_Directory(_aflags& aflags);
  void AFLOW_RUN_Directory(const _aflags& aflags);
  void RUN_DirectoryScript(const _aflags& aflags,const string& script,const string& output);
  bool CompressDirectory(const _aflags& aflags,const _kflags& kflags);
  bool CompressDirectory(const _aflags& aflags);
  void Clean(const _aflags& aflags);
  void Clean(const string directory);
  void XClean(string options);
  void GenerateAflowinFromVASPDirectory(_aflags& aflags);
  void StartStopCheck(const string &AflowIn,string str1,string str2,bool &flag,bool &flagS);
  void StartStopCheck(const string &AflowIn,string str1,bool &flag,bool &flagS);
  bool Legitimate_aflowin(string aflowindir,const bool& osswrite,ostringstream& oss);
  bool Legitimate_aflowin(string aflowindir);
  void getAflowInFromAFlags(const _aflags& aflags,string& AflowIn_file,string& AflowIn,ostream& oss=cout); //CO20191110
  void getAflowInFromAFlags(const _aflags& aflags,string& AflowIn_file,string& AflowIn,ofstream& FileMESSAGE,ostream& oss=cout); //CO20191110
  void getAflowInFromDirectory(const string& directory,string& AflowIn_file,string& AflowIn,ostream& oss=cout); //CO20191110
  void getAflowInFromDirectory(const string& directory,string& AflowIn_file,string& AflowIn,ofstream& FileMESSAGE,ostream& oss=cout); //CO20191110
  int get_NCPUS();  //ME20200219
  int get_NCPUS(const _kflags&);  //ME20200219
}

// ----------------------------------------------------------------------------
// aflow_modules.cpp
//ME20181027
namespace KBIN {
  void setModules(_xvasp&);
  void setModules(_xinput&);
  void readModulesFromAflowIn(const string&, _kflags&, _xvasp&);
  void readModulesFromAflowIn(const string&, _kflags&, _xinput&);
  vector<aurostd::xoption> loadDefaultsAPL();
  bool writeFlagAPL(const string& key,const xoption& xopt); //CO20181226  //ME20190113
  void readParametersAPL(const string&, _moduleOptions&, _xinput&);
  vector<aurostd::xoption> loadDefaultsQHA();  //AS20200302
  void readParametersQHA(const string&, _moduleOptions&, _xinput&); //AS20200302
  vector<aurostd::xoption> loadDefaultsAAPL();
  bool writeFlagAAPL(const string& key,const xoption& xopt);  //CO20181226  //ME20190113
  void readParametersAAPL(const string&, _moduleOptions&, _xinput&);
  vector<aurostd::xoption> loadDefaultsAEL();
  bool writeFlagAEL(const string& key,const xoption& xopt); 
  vector<aurostd::xoption> loadDefaultsAGL();
  bool writeFlagAGL(const string& key,const xoption& xopt); 

}

// ----------------------------------------------------------------------------
// aflow_qsub.cpp
namespace KBIN {
  bool QSUB_Extract(_xqsub& xqsub,string AflowIn,ifstream &FileAFLOWIN,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags);
  bool QSUB_RunFinished(_aflags &aflags,ofstream &FileMESSAGE,bool=FALSE);
  void QSUB_WaitFinished(_aflags &aflags,ofstream &FileMESSAGE,bool=FALSE);
  bool QSUB_Extract_Mode1(_xqsub& xqsub,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags);
  bool QSUB_Extract_Mode2(_xqsub& xqsub,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags);
  bool QSUB_Extract_Mode3(_xqsub& xqsub,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags);
}

// ----------------------------------------------------------------------------
// aflow_ialien.cpp
namespace ALIEN {
  bool Produce_INPUT(_xalien& xalien,string AflowIn,ifstream &FileAFLOWIN,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_alienflags &alienflags);
  bool Modify_INPUT(_xalien& xalien,ofstream &FileMESSAGE,_aflags &aflags,_alienflags &alienflags);
  bool Write_INPUT(_xalien& xalien);
  bool Produce_INPUT_FILE(_xalien& xalien,string AflowIn,ifstream &FileAFLOWIN,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_alienflags &alienflags);
  bool Modify_INPUT_FILE(_xalien& xalien,ofstream &FileMESSAGE,_aflags &aflags,_alienflags &alienflags);
}

// ----------------------------------------------------------------------------
// aflow_kalien.cpp
namespace ALIEN {
  _alienflags Get_Alienflags_from_AflowIN(string &AflowIn);
  bool Run_Directory(ofstream& FileERROR,_aflags& aflags,_kflags& kflags);
}

// ----------------------------------------------------------------------------
// aflow_matlab.cpp aflow_matlab_funcs.cpp
bool KBIN_MATLAB_Extract(string AflowIn,ifstream &FileAFLOWIN,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags);
bool KBIN_MATLAB_Run(_kflags &kflags);
_kflags KBIN_MATLAB_Get_Matlabflags_from_AflowIN(string &AflowIn);
bool KBIN_MATLAB_Directory(ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags);
string MATLAB_FUNCS_param(void);
//string MATLAB_FUNCS_plotband(string DIRECTORY,string OPTION1);
string MATLAB_FUNCS_plotband(void);

// ----------------------------------------------------------------------------
// aflow_gnuplot_plotbz.cpp
//string GNUPLOT_FUNCS_plotbz(string DIRECTORY,string OPTION1);
string GNUPLOT_FUNCS_plotbz(void);

// ----------------------------------------------------------------------------
// aflow_ifrozsl.cpp

namespace KBIN {
  void VASP_RunPhonons_FROZSL(_xvasp &xvasp,string AflowIn,_aflags &aflags,_kflags &kflags,_vflags &vflags,ofstream &FileMESSAGE);
}

namespace FROZSL {
  bool Extract_INPUT(const string& AflowIn,ofstream &FileMESSAGE,stringstream &input_file,_aflags &aflags,_kflags &kflags);
  bool Setup_frozsl_init_input(const string& AflowIn,ofstream &FileMESSAGE,stringstream &input_file,_aflags &aflags,_kflags &kflags);
  bool Already_Calculated_Input(const string& AflowIn);
  bool WGET_INPUT(ofstream &FileMESSAGE,string AflowIn,_aflags &aflags,_kflags &kflags);
  bool WGET_OUTPUT(ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags);
  bool input_TO_poscar(ofstream &FileMESSAGE,stringstream &input_file,_aflags &aflags,_kflags &kflags);
  string Generate_Input_file(ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags);
  bool File_INPUT(const string& AflowIn,ofstream &FileMESSAGE,stringstream &input_file,_aflags &aflags,_kflags &kflags);
  bool Write(string data,string directory);
  bool Delete(string data,string directory);
}
namespace FINDSYM {
  bool Write(string data,string directory);
}

// ----------------------------------------------------------------------------
// aflow_kvasp.cpp

namespace KBIN {
  _kflags VASP_Get_Kflags_from_AflowIN(const string &AflowIn,_aflags &aflags,ostream& oss=cout);
  _kflags VASP_Get_Kflags_from_AflowIN(const string &AflowIn,ofstream &FileMESSAGE,_aflags &aflags,ostream& oss=cout);
  _vflags VASP_Get_Vflags_from_AflowIN(const string &AflowIn,_aflags &aflags,_kflags& kflags,ostream& oss=cout);
  _vflags VASP_Get_Vflags_from_AflowIN(const string &AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags& kflags,ostream& oss=cout);
  bool VASP_Fix_Machine_Kflags_from_AflowIN(ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_vflags &vflags);
  bool VASP_Directory(ofstream& FileERROR,_aflags& aflags,_kflags& kflags);
  void VASP_BackupOriginal(_aflags aflags);
  // [OBSOLETE] G++6 not needed  void VASP_Wait(_xvasp& xvasp,_aflags &aflags,_kflags &kflags,_vflags &vflags,ofstream &FileMESSAGE);
  bool VASP_Run(_xvasp &xvasp,_aflags &aflags,_kflags &kflags,_vflags &vflags,ofstream &FileMESSAGE);
  bool VASP_Run(_xvasp &xvasp,_aflags &aflags,_kflags &kflags,_vflags &vflags,string relaxA,string relaxB,bool qmwrite,ofstream &FileMESSAGE);
  bool VASP_Run(_xvasp &xvasp,_aflags &aflags,_kflags &kflags,_vflags &vflags,string relaxA,bool qmwrite,ofstream &FileMESSAGE);
  bool VASP_RunFinished(_xvasp &xvasp,_aflags &aflags,ofstream &FileMESSAGE,bool=FALSE);
  void WaitFinished(_xvasp &xvasp,_aflags &aflags,ofstream &FileMESSAGE,bool=FALSE);
  void VASP_Error(_xvasp &xvasp,string="",string="",string="");
  void VASP_Error(_xvasp &xvasp,ofstream &FileMESSAGE,string="",string="",string="");
  string VASP_Analyze(_xvasp &xvasp,bool qmwrite);
  void VASP_CompressDirectory(_xvasp xvasp,_kflags &kflags);
  void VASP_Backup(_xvasp& xvasp,bool qmwrite,string relax);
  void VASP_CONTCAR_Save(_xvasp xvasp,string relax);
  void VASP_Recycle(_xvasp xvasp,string relax);
  void VASP_Recycle(_xvasp xvasp,int relax_number);
  void VASP_RecycleExtraFile(_xvasp xvasp,string xfile,string relax);
  void VASP_RecycleExtraFile(_xvasp xvasp,string xfile,int relax_number);
  bool VASP_CheckUnconvergedOSZICAR(string dir);
  void GetStatDiel(string& outcar, xvector<double>& eigr, xvector<double>& eigi); // CAMILO
  void GetDynaDiel(string& outcar, xvector<double>& eigr, xvector<double>& eigi); // CAMILO
  string getVASPVersionString(const string&);  //ME20190219
  string getVASPVersionNumber(const string&);  //CO20200610
  double getVASPVersion(const string&);  //CO20200610
}

// ----------------------------------------------------------------------------
// aflow_ivasp.cpp
namespace KBIN {
  bool VASP_Produce_INPUT(_xvasp& xvasp,const string& AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_vflags &vflags,bool load_POSCAR_from_xvasp=false);
  bool VASP_Modify_INPUT(_xvasp& xvasp,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_vflags &vflags);
  bool VASP_Produce_and_Modify_INPUT(_xvasp& xvasp,const string& AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_vflags &vflags,bool load_POSCAR_from_xvasp=false); //CO20180418
  bool VASP_Write_ppAUID_FILE(const string& directory,const vector<string>& vppAUIDs,const vector<string>& species);
  bool VASP_Write_ppAUID_FILE(const string& directory,const deque<string>& vppAUIDs,const deque<string>& species);
  bool VASP_Write_ppAUID_AFLOWIN(const string& directory,const vector<string>& vppAUIDs,const vector<string>& species);
  bool VASP_Write_ppAUID_AFLOWIN(const string& directory,const deque<string>& vppAUIDs,const deque<string>& species);
  bool VASP_Write_INPUT(_xvasp& xvasp,_vflags &vflags);
  bool VASP_Produce_INCAR(_xvasp& xvasp,const string& AflowIn,ofstream& FileERROR,_aflags& aflags,_kflags& kflags,_vflags& vflags);
  bool VASP_Modify_INCAR(_xvasp& xvasp,ofstream& FileERROR,_aflags& aflags,_kflags& kflags,_vflags& vflags);
  bool VASP_Reread_INCAR(_xvasp& xvasp,ofstream &FileMESSAGE,_aflags &aflags);
  bool VASP_Produce_POSCAR(_xvasp& xvasp,const string& AflowIn,ofstream& FileERROR,_aflags& aflags,_kflags& kflags,_vflags& vflags);
  bool VASP_Produce_POSCAR(_xvasp& xvasp);
  bool VASP_Modify_POSCAR(_xvasp& xvasp,const string& AflowIn,ofstream& FileERROR,_aflags& aflags,_vflags& vflags);
  void convertPOSCARFormat(_xvasp&, const _kflags&);  //ME20190220
  bool VASP_Convert_Unit_Cell(_xvasp&, _vflags&, _aflags&, ofstream&, ostringstream&); //ME20181216
  bool VASP_Reread_POSCAR(_xvasp& xvasp,ofstream &FileMESSAGE,_aflags &aflags);
  bool VASP_Produce_KPOINTS(_xvasp& xvasp,const string& AflowIn,ofstream& FileERROR,_aflags& aflags,_kflags& kflags,_vflags& vflags);
  bool VASP_Modify_KPOINTS(_xvasp& xvasp,ofstream& FileERROR,_aflags& aflags,_vflags& vflags);
  bool VASP_Reread_KPOINTS(_xvasp& xvasp,ofstream &FileMESSAGE,_aflags &aflags);
  bool VASP_Find_DATA_POTCAR(const string& species_pp,string &FilePotcar,string &DataPotcar,string &AUIDPotcar);
  bool VASP_Find_FILE_POTCAR(const string& species_pp,string &FilePotcar,string &DataPotcar,string &AUIDPotcar);
  bool VASP_Produce_POTCAR(_xvasp& xvasp,const string& AflowIn,ofstream& FileERROR,_aflags& aflags,_kflags& kflags,_vflags& vflags);
  bool VASP_Modify_POTCAR(_xvasp& xvasp,ofstream& FileERROR,_aflags& aflags,_vflags& vflags);
  bool VASP_Reread_POTCAR(_xvasp& xvasp,ofstream &FileMESSAGE,_aflags &aflags);
  string VASP_PseudoPotential_CleanName(const string& specieIN);
  string VASP_PseudoPotential_CleanName_20190712(const string& specieIN); //CO20190712
  string VASP_PseudoPotential_CleanName_20190101(const string& specieIN); //CO20190712
  void VASP_PseudoPotential_CleanName_InPlace(string& species,bool capital_letters_only=false); //CO20190712
  bool VASP_PseudoPotential_CleanName_TEST(void); //CO20190712
  uint VASP_SplitAlloySpecies(string alloy_in, vector<string> &speciesX);
  uint VASP_SplitAlloySpecies(string alloy_in, vector<string> &speciesX, vector<double> &natomsX);
  bool VASP_SplitAlloySpecies(string alloy_in, string &specieA, string &specieB);
  bool VASP_SplitAlloySpecies(string alloy_in, string &specieA, string &specieB, string &specieC);
  bool VASP_SplitAlloySpecies(vector<string> alloy, vector<string> &speciesA, vector<string> &speciesB);
  bool VASP_SplitAlloySpecies(vector<string> alloy, vector<string> &speciesA, vector<string> &speciesB, vector<string> &speciesC);
  uint VASP_SplitAlloyPseudoPotentials(string alloy_in, vector<string> &species_ppX);
  uint VASP_SplitAlloyPseudoPotentials(string alloy_in, vector<string> &species_ppX, vector<double> &natomsX);
  bool VASP_SplitAlloyPseudoPotentials(string alloy, string &species_ppA, string &species_ppB);
  bool VASP_SplitAlloyPseudoPotentials(string alloy, string &species_ppA, string &species_ppB, string &species_ppC);
  bool VASP_SplitAlloyPseudoPotentials(vector<string> alloy, vector<string> &species_ppsA, vector<string> &species_ppsB);
  bool VASP_SplitAlloyPseudoPotentials(vector<string> alloy, vector<string> &species_ppsA, vector<string> &species_ppsB, vector<string> &species_ppsC);
  void VASP_MPI_Autotune(_xvasp& xvasp,_aflags &aflags,bool VERBOSE);
  void XVASP_INCAR_System_Auto(_xvasp& xvasp,bool VERBOSE);
  void XVASP_INCAR_Relax_ON(_xvasp& xvasp,bool VERBOSE);
  void XVASP_INCAR_Relax_ON(_xvasp& xvasp,_vflags& vflags,int number); // for steps
  void XVASP_INCAR_Static_ON(_xvasp& xvasp,_vflags& vflags);
  void XVASP_INCAR_Relax_Static_ON(_xvasp& xvasp,_vflags& vflags);
  void XVASP_INCAR_Relax_Static_Bands_ON(_xvasp& xvasp,_vflags& vflags);
  void XVASP_INCAR_RWIGS_Static(_xvasp& xvasp,_vflags& vflags,ofstream &FileMESSAGE,bool OPERATION);
  void XVASP_INCAR_Precision(_xvasp& xvasp,_vflags& vflags);
  void XVASP_INCAR_Metagga(_xvasp& xvasp,_vflags& vflags);
  void XVASP_INCAR_Ivdw(_xvasp& xvasp,_vflags& vflags);
  void XVASP_INCAR_ABMIX(_xvasp& xvasp,_vflags& vflags);
  int XVASP_INCAR_GetNBANDS(_xvasp& xvasp,bool ispin);
  bool XVASP_INCAR_PREPARE_GENERIC(string command,_xvasp& xvasp,_vflags& vflags,string svalue,int ivalue,double dvalue,bool bvalue);
  //  bool XVASP_INCAR_PREPARE_GENERIC(string command,_xvasp& xvasp,_kflags &kflags,_vflags& vflags,string svalue,int ivalue,double dvalue,bool bvalue);
  // ALGO, ENMAX_MULTIPLY, IMIX, IALGO, TYPE, PAW_CORRECTIONS, NBANDS, PSTRESS, EDIFFG, POTIM, SPIN, LS_COUPLING, AUTO_MAGMOM, NWS, SYM, WAVECAR, CHGCAR
  void XVASP_INCAR_ADJUST_ICHARG(_xvasp&, _vflags&, _aflags&, int, ofstream&);  //ME20191028
  void XVASP_INCAR_SPIN_REMOVE_RELAX(_xvasp& xvasp,_aflags &aflags,_vflags& vflags,int step,ofstream &FileMESSAGE);
  void XVASP_KPOINTS_IBZKPT_UPDATE(_xvasp& xvasp,_aflags &aflags,_vflags& vflags,int step,ofstream &FileMESSAGE);
  void XVASP_INCAR_LDAU_OFF(_xvasp& xvasp,bool VERBOSE);
  void XVASP_INCAR_LDAU_ON(_xvasp& xvasp,_vflags& vflags,uint type);
  void XVASP_INCAR_LDAU_ADIABATIC(_xvasp& xvasp,int step);
  void XVASP_INCAR_LDAU_CUTOFF(_xvasp& xvasp,bool VERBOSE);
  void XVASP_INCAR_KPOINTS_Dielectric_SET(_xvasp& xvasp,_kflags &kflags,_vflags& vflags,string mode_dielectric);
  void XVASP_INCAR_REMOVE_ENTRY(_xvasp& xvasp,string ENTRY,string COMMENT,bool VERBOSE);

  bool XVASP_KPOINTS_KPOINTS(_xvasp &xvasp,ofstream &FileMESSAGE,bool VERBOSE);
  bool XVASP_KPOINTS_KPOINTS(_xvasp &xvasp);
  // bool XVASP_KPOINTS_EVEN(_xvasp& xvasp); TO REMOVE
  // bool XVASP_KPOINTS_ODD(_xvasp& xvasp); TO REMOVE
  bool XVASP_KPOINTS_OPERATION(_xvasp& xvasp,string operation);
  // bool XVASP_KPOINTS_Kshift_Gamma_EVEN(_xvasp& xvasp); TO REMOVE
  // bool XVASP_KPOINTS_Kshift_Gamma_ODD(_xvasp& xvasp); TO REMOVE
  // bool XVASP_KPOINTS_Kscheme(_xvasp& xvasp,string kscheme);
  bool XVASP_KPOINTS_Fix_KPPRA(_xvasp &xvasp,int NK,ofstream &FileMESSAGE,bool VERBOSE);
  bool XVASP_KPOINTS_Fix_KSHIFT(_xvasp &xvasp,_xvasp &rxvasp,bool KAUTOSHIFT,bool VERBOSE);
  bool XVASP_KPOINTS_Fix_KPOINTS(_xvasp &xvasp,int NK,ofstream &FileMESSAGE,bool VERBOSE);
  void XVASP_string2numbers(_xvasp& xvasp);
  void XVASP_numbers2string(_xvasp& xvasp);
  void XVASP_Afix_Clean(_xvasp& xvasp,string preserve_name);
  void XVASP_Afix_ROTMAT(_xvasp& xvasp,int mode,bool verbose,_aflags &aflags,ofstream &FileMESSAGE);
  void XVASP_Afix_NBANDS(_xvasp& xvasp,int& nbands,bool VERBOSE);
  void XVASP_Afix_POTIM(_xvasp& xvasp,double& potim,bool VERBOSE);
  double XVASP_Afix_GENERIC(string mode,_xvasp& xvasp,_kflags& kflags,_vflags& vflags,double=0.0,int=0);

  string ExtractSystemName(const string& directory);  //ME20200217
  string ExtractSystemNameFromAFLOWIN(const string& directory);  //ME20200217
  string ExtractSystemNameFromVASP(const string& directory);  //ME20200217
  double ExtractEfermiOUTCAR(string directory);
  xstructure GetMostRelaxedStructure(string directory); //CO20180627
  vector<string> ExtractAtomicSpecies(const string& directory,ostream& oss=cout);
  vector<string> ExtractAtomicSpecies(const string& directory,ofstream& FileMESSAGE,ostream& oss=cout);

}

// ----------------------------------------------------------------------------
// aflow_avasp.cpp
#define _AVASP_PSEUDOPOTENTIAL_AUTO_ string("AUTO")
#define _AVASP_PSEUDOPOTENTIAL_DELIMITER_ string(":")
#define _AVASP_PSEUDOPOTENTIAL_POTENTIAL_TYPE_ string("TYPE") //CO20191020
#define _AVASP_PSEUDOPOTENTIAL_POTENTIAL_COMPLETE_ string("COMPLETE")

struct _AVASP_PROTO{
  vector<string> ucell;
  deque<int> vkppra;
  vector<double> vpressure;
  aurostd::xoption vparams;
};

bool AVASP_MakePrototype_AFLOWIN(_AVASP_PROTO *PARAMS);
bool AVASP_MakePrototype_AFLOWIN_20181226(_AVASP_PROTO *PARAMS);
bool AVASP_MakePrototype_AFLOWIN_20180101(_AVASP_PROTO *PARAMS);
bool AVASP_MakePrototypeICSD_AFLOWIN(_AVASP_PROTO *PARAMS,bool flag_AFLOW_IN_ONLY_IF_MISSING);
void AVASP_Get_LDAU_Parameters(string species,bool &LDAU,vector<string>& vLDAUspecies,
    vector<uint>& vLDAUtype,vector<int>& vLDAUL, vector<double>& vLDAUU, vector<double> &vLDAUJ);
string AVASP_Get_PseudoPotential_PAW_PBE_KIN(string species);
string AVASP_Get_PseudoPotential_PAW_PBE(string species);
string AVASP_Get_PseudoPotential_PAW_GGA(string species);
string AVASP_Get_PseudoPotential_PAW_LDA_KIN(string species);
string AVASP_Get_PseudoPotential_PAW_LDA(string species);
string AVASP_Get_PseudoPotential_PBE(string species);
string AVASP_Get_PseudoPotential_GGA(string species);
string AVASP_Get_PseudoPotential_LDA(string species);
bool AVASP_populateXVASP(const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,_xvasp& xvasp);
void AVASP_populateXVASP_ARUN(const _aflags&,const _kflags&, const _vflags&,_xvasp&);  //ME20181030
void setStatic(_xvasp&);  //ME20181102
void setPreserveUnitCell(_xvasp&);   //ME20181102
void AVASP_fix_volumes_masses_XVASP(_xvasp&,bool skip_volume=false); //ME20181103 //CO20181226
bool AVASP_MakeSingleAFLOWIN(_xvasp& xvaspin,stringstream &_aflowin,bool flag_WRITE,int=-1,bool flag_PRINT=TRUE);   // last is pthread number, if <0 then serial
bool AVASP_MakeSingleAFLOWIN_20181226(_xvasp& xvaspin,stringstream &_aflowin,bool flag_WRITE,int=-1,bool flag_PRINT=TRUE);   // last is pthread number, if <0 then serial
bool AVASP_MakeSingleAFLOWIN_20180101(_xvasp& xvaspin,stringstream &_aflowin,bool flag_WRITE,int=-1,bool flag_PRINT=TRUE);   // last is pthread number, if <0 then serial
bool AVASP_MakeSingleAFLOWIN(_xvasp& xvasp_in,bool flag_WRITE,int=-1,bool flag_PRINT=TRUE);  // last is pthread number, if <0 then serial
bool AVASP_MakeSingleAFLOWIN(_xvasp& xvasp_in,int=-1,bool flag_PRINT=TRUE);  // last is pthread number, if <0 then serial
bool AVASP_DefaultValuesBinary_AFLOWIN(_xvasp &xvasp);
bool AVASP_MakeSinglePOSCAR(_xvasp& xvaspin);
bool Alloys_LibraryU(vector<string> &alloy,vector<string> &pseudosA,vector<string> &pseudosB);
bool Alloys_LibraryG(vector<string> &alloy,vector<string> &pseudosA,vector<string> &pseudosB);
bool Alloys_LibraryX(vector<string> &alloy,vector<string> &pseudosA,vector<string> &pseudosB);
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// aflow_ovasp.cpp
class xOUTCAR;
class xDOSCAR;
class xEIGENVAL;
class xPOTCAR;
class xVASPRUNXML;
class xIBZKPT;
class xKPOINTS;
class xCHGCAR;
class xVASPOUT;
class xQMVASP;  //CO20190803
namespace aflowlib { class _aflowlib_entry;}

// -------------------------------------------------------------------------------------------------
class xOUTCAR : public xStream { //CO20200404 - xStream integration for logging
  public:
    xOUTCAR(ostream& oss=cout);                         // default, just allocate  //CO20200404 - xStream integration for logging
    xOUTCAR(ofstream& FileMESSAGE,ostream& oss=cout);   // default, just allocate  //CO20200404 - xStream integration for logging
    xOUTCAR(const string& fileIN,bool=TRUE,ostream& oss=cout);                        // constructor from filename, QUIET  //CO20200404 - xStream integration for logging
    xOUTCAR(const string& fileIN,ofstream& FileMESSAGE,bool=TRUE,ostream& oss=cout);  // constructor from filename, QUIET  //CO20200404 - xStream integration for logging
    bool initialize(const string& fileIN,ostream& oss,bool=TRUE);  //ME20200427  //CO20200508
    bool initialize(const string& fileIN,ofstream& FileMESSAGE,ostream& oss,bool=TRUE);  //ME20200427  //CO20200508
    bool initialize(const string& fileIN,bool=TRUE);  //ME20200427  //CO20200508

    xOUTCAR(const xOUTCAR& b);                                    // constructor copy
    ~xOUTCAR();                                                   // kill everything
    const xOUTCAR& operator=(const xOUTCAR &b);                   // copy
    void clear(void);                                             // clear

    bool m_initialized;  //CO20200404 - xStream integration for logging

    // CONTENT
    string content;vector<string> vcontent;string filename;       // the content, and lines of it
    string SYSTEM;
    int NIONS;
    double Efermi;
    bool isLSCOUPLING;
    int nelectrons; //AS20200528
    double natoms;                                                // for aflowlib_libraries.cpp
    double energy_cell,energy_atom;                               // for aflowlib_libraries.cpp
    double enthalpy_cell,enthalpy_atom;                           // for aflowlib_libraries.cpp
    double eentropy_cell,eentropy_atom;                           // for aflowlib_libraries.cpp
    double PV_cell,PV_atom;                                       // for aflowlib_libraries.cpp
    xmatrix<double> stress;                                       // for aflowlib_libraries.cpp
    double mag_cell,mag_atom;                                     // for aflowlib_libraries.cpp
    vector<double> vmag;                                          // for aflowlib_libraries.cpp
    vector<xvector<double> > vmag_noncoll;                        //DX20171205 - non-collinear
    double volume_cell,volume_atom;                               // for aflowlib_libraries.cpp
    double pressure;                                              // for aflowlib_libraries.cpp // SAME AS PSTRESS
    double pressure_residual;                                     // for aflowlib_libraries.cpp
    double Pulay_stress;                                          // for aflowlib_libraries.cpp
    vector<aurostd::xvector<double> > vforces;                    // for aflowlib_libraries.cpp
    vector<aurostd::xvector<double> > vpositions_cartesian;       // for aflowlib_libraries.cpp
    double ENCUT,EDIFF,EDIFFG,POTIM,TEIN,TEBEG,TEEND,SMASS,NPACO,APACO,PSTRESS;     // 
    int NBANDS,NKPTS,NSW,NBLOCK,KBLOCK,IBRION,NFREE,ISIF,IWAVPR,ISYM,ISPIN;   // for aflowlib_libraries.cpp
    double total_energy_change;                                   // for aflowlib_libraries.cpp
    // DOS related values
    double EMIN,EMAX,SIGMA;                                       // eV - energy-range for DOS
    int ISMEAR;                                                   // broadening in eV -4-tet -1-fermi 0-gaus
    //  Electronic relaxation
    int IALGO;              //  algorithm                         // for aflowlib_libraries.cpp
    string LDIAG;           //   sub-space diagonalisation        // for aflowlib_libraries.cpp
    int IMIX,INIMIX,MIXPRE; //     mixing-type and parameters     // for aflowlib_libraries.cpp
    double AMIX,BMIX,AMIX_MAG,BMIX_MAG,AMIN,WC; // parameters     // for aflowlib_libraries.cpp
    // Intra band minimization
    double WEIMIN,EBREAK,DEPER,TIME;  // for aflowlib_libraries.cpp
    // begin shared xPOTCAR
    double ENMAX;vector<double> vENMAX;                            // eV
    double ENMIN;vector<double> vENMIN;                            // eV
    double POMASS_sum,POMASS_min,POMASS_max;vector<double> vPOMASS;// mass
    double ZVAL_sum,ZVAL_min,ZVAL_max;vector<double> vZVAL;   // valence
    double EATOM_min,EATOM_max;vector<double> vEATOM;        // eV
    double RCORE_min,RCORE_max;vector<double> vRCORE;        // outmost cutoff radius
    double RWIGS_min,RWIGS_max;vector<double> vRWIGS;        // wigner-seitz radius (au A)
    double EAUG_min,EAUG_max;vector<double> vEAUG;            // augmentation
    double RAUG_min,RAUG_max;vector<double> vRAUG;            // augmentation
    double RMAX_min,RMAX_max;vector<double> vRMAX;            // unicity
    vector<string> vTITEL;                                         // unicity
    vector<string> vLEXCH;                                         // unicity
    // end shared xPOTCAR
    string pp_type;
    vector<string> species;                                        // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
    vector<int>    species_Z;                                      // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
    vector<string> species_pp;                                     // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
    vector<string> species_pp_type;                                // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
    vector<string> species_pp_version;                             // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
    vector<string> species_pp_AUID;                                // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
    vector<string> species_pp_AUID_collisions;                     // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
    vector<double> species_pp_groundstate_energy;                  // meV/atom
    vector<string> species_pp_groundstate_structure;               // name that we have, maybe ANRL
    deque<deque<double> > species_pp_vLDAU;                       // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
    bool isKIN;                                                   // METAGGA
    bool isMETAGGA;string METAGGA;                                // METAGGA
    string string_LDAU;                                           // for aflowlib_libraries.cpp
    uint nweights,nkpoints_irreducible;                           // kpoints reading
    vector<aurostd::xvector<double> > vkpoint_reciprocal;         // kpoints reading
    vector<aurostd::xvector<double> > vkpoint_cartesian;          // kpoints reading
    vector<double> vweights;                                      // kpoints reading
    double calculation_time;                                      // for aflowlib_libraries.cpp - calculation_time
    double calculation_memory;                                    // for aflowlib_libraries.cpp - calculation_memory
    uint calculation_cores;                                       // for aflowlib_libraries.cpp - calculation_cores
    xstructure xstr;                                              // for GetBandGap()
    vector<string> GetCorrectPositions(string line,uint expected_count);                //CO20170725 - vasp issues with lattice spacing (negative sign) 
    bool GetProperties(const stringstream& stringstreamIN,bool=TRUE);          // get everything QUIET
    bool GetProperties(const string& stringIN,bool=TRUE);                      // get everything QUIET
    bool GetPropertiesFile(const string& fileIN,bool=TRUE);                    // get everything QUIET
    bool GetPropertiesFile(const string& fileIN,uint natoms_check,bool=TRUE);       // get everything QUIET  //CO20200404 - added default for bool
    bool GetPropertiesUrlFile(const string& url,const string& file,bool=TRUE); // get everything from an aflowlib entry
    vector<int> band_index;
    vector<int> carrier_spin;
    vector<string> carrier_type;
    vector<vector<double> > extrema_cart_coord;
    vector<vector<double> > effective_mass_axes;
    vector<int> equivalent_valley;
    vector<double> effective_mass_DOS;
    vector<double> effective_mass_COND;
    vector<double> mass_elec_dos;
    vector<double> mass_hole_dos;
    vector<double> mass_elec_conduction;
    vector<double> mass_hole_conduction;
    // BAND GAPS
    bool GetXStructure();
    int isKPointLine(uint iline,xvector<double>& kpoint); //if returns 0 if not KPointLine, -1 means it gave *** for kpoint
    int isKPointLine(uint iline);                         //if returns 0 if not KPointLine, -1 means it gave *** for kpoint
    bool GetStartingKPointLines(vector<uint>& ilines);
    bool GetNextKPointLine(uint& iline);
    bool ProcessKPoint(uint iline,double EFERMI,vector<double>& b_energies,vector<double>& b_occs);
    bool GetBandEdge(vector<double>& b_energies,vector<double>& b_occs,double EFERMI,uint& iedge,double efermi_tol=AUROSTD_NAN,double energy_tol=1e-4,double occ_tol=1e-5);
    bool identicalKPoints(vector<xvector<double> >& vkpoints,uint kpt1,uint kpt2,double tol=1e-12);
    bool identicalKPoints(xvector<double>& kpoint1,xvector<double>& kpoint2,double tol=1e-12);
    bool removeDuplicateKPoints(vector<xvector<double> >& vkpoints,vector<uint>& vikpt);
    bool removeDuplicateKPoints(vector<vector<xvector<double> > >& vkpoints,vector<uint>& vikpt,vector<uint>& vispin);
    double minimumDistanceKPoints(vector<xvector<double> >& vkpoints,uint ikp1,uint ikp2);
    double minimumDistanceKPoints(xvector<double>& kpoint1,xvector<double>& kpoint2);
    struct bandEnergyOcc{
      double energy;
      double occ;
      //bool operator<(const bandEnergyOcc& other) const {return energy<other.energy;}
    };
    struct bandEnergyOccCompare{
      bandEnergyOccCompare(double _energy_tol) : energy_tol(_energy_tol) {};
      double energy_tol;
      bool operator()(const bandEnergyOcc& a,const bandEnergyOcc b) const;
    };
    bool orderBands(vector<double>& b_energies,vector<double>& b_occs,double energy_tol=1e-4);
    enum BROAD_TYPES {empty,metal,insulator};                   //bandgap types
    enum EMPTY_TYPES {empty_all,empty_partial};                 //bandgap types
    enum INSULATOR_TYPES {insulator_direct,insulator_indirect}; //bandgap types
    enum GAP_TYPES {zero_gap,non_zero_gap};                     //bandgap types
    //bool GetBandGap(void); //CO20171002 - need POSCAR for kpt_tol
    bool GetBandGap(double EFERMI=AUROSTD_NAN,double efermi_tol=AUROSTD_NAN,double energy_tol=1e-4,double occ_tol=1e-5);
    bool GetBandGap_Camilo(double kpt_tol);
    vector<double> conduction_band_min;
    double         conduction_band_min_net;
    vector<double> valence_band_max;
    double         valence_band_max_net;
    vector<double> Egap;
    double         Egap_net;
    vector<double> Egap_fit;
    double         Egap_fit_net;
    vector<string> Egap_type;
    string         Egap_type_net;
    //[CO20200404 - OBSOLETE]string ERROR;
    //int number_bands,number_kpoints; //CO20171006 - camilo garbage
    //int ISPIN; // turn this into spin = 0 if ISPIN = 1 //CO20171006 - camilo garbage
    //int spin;  //CO20171006 - camilo garbage
    friend ostream& operator<<(ostream&, const xOUTCAR&);  //ME20190623
  private:                       //
    void free();                 // free space
    void copy(const xOUTCAR& b); //
};

// EFFECTIVE MASSES //CO20200404 - moved from "friend" of xOUTCAR
bool GetEffectiveMass(xOUTCAR& outcar, xDOSCAR& doscar, xEIGENVAL& eigenval, xstructure xstr,ostream& oss=cout); //CO20200404
bool GetEffectiveMass(xOUTCAR& outcar, xDOSCAR& doscar, xEIGENVAL& eigenval, xstructure xstr,ofstream& FileMeSSAGE,ostream& oss=cout); //CO20200404

//-------------------------------------------------------------------------------------------------
class xDOSCAR : public xStream { //CO20200404 - xStream integration for logging
  public:
    xDOSCAR(ostream& oss=cout);                           // default, just allocate  //CO20200404 - xStream integration for logging
    xDOSCAR(ofstream& FileMESSAGE,ostream& oss=cout);     // constructor from filename QUIET //CO20200404 - xStream integration for logging
    xDOSCAR(const string& fileIN,bool=TRUE,ostream& oss=cout);                          // constructor from filename QUIET //CO20200404 - xStream integration for logging
    xDOSCAR(const string& fileIN,ofstream& FileMESSAGE,bool=TRUE,ostream& oss=cout);    // constructor from filename QUIET //CO20200404 - xStream integration for logging
    bool initialize(const string& fileIN,ostream& oss,bool=TRUE);  //ME20200427  //CO20200508
    bool initialize(const string& fileIN,ofstream& FileMESSAGE,ostream& oss,bool=TRUE);  //ME20200427  //CO20200508
    bool initialize(const string& fileIN, bool=TRUE); //ME20200427

    xDOSCAR(const xDOSCAR& b);                                    // constructor copy
    ~xDOSCAR();                                                   // kill everything
    const xDOSCAR& operator=(const xDOSCAR &b);                   // copy
    void clear(void);                                             // clear

    bool m_initialized;  //CO20200404 - xStream integration for logging

    // CONTENT
    string content;vector<string> vcontent;string filename;       // the content, and lines of it
    string title;
    uint spin;
    double Vol,POTIM;
    xvector<double> lattice;
    double temperature;
    bool RWIGS;
    double Efermi;
    double spinF;
    double energy_max;
    double energy_min;
    uint number_energies;
    uint number_atoms;  //ME20190614
    bool partial;  //ME20190614
    double denergy;
    deque<double> venergy;                                        // venergy.at(energy_number) 
    deque<double> venergyEf;                                      // venergyEf.at(energy_number) 
    //ME20190614 BEGIN
    //[OBSOLETE]  deque<deque<double> > vDOS;                                   // vDOS.at(energy_number).at(spin)
    deque<deque<double> > viDOS;                                  // viDOS.at(spin).at(energy_number)
    //[OBSOLETE]  deque<deque<double> > vDOSs;                                  // vDOSs.at(energy_number).at(spin)
    //[OBSOLETE]  deque<deque<double> > vDOSp;                                  // vDOSp.at(energy_number).at(spin)
    //[OBSOLETE]  deque<deque<double> > vDOSd;                                  // vDOSd.at(energy_number).at(spin)
    deque<deque<deque<deque<double> > > > vDOS;                   // vDOS.at(atom).at(orbital).at(spin).at(energy_number); 0 = total for atoms and orbitals
    //ME20190614 END
    //ME20190620 BEGIN
    bool isLSCOUPLING;  // Contains spin-orbit coupling
    bool lmResolved;  // Is it lm-resolved?
    string carstring;  // The fourth line of the DOSCAR
    //ME20190620 END
    //[CO20200404 - OBSOLETE]string ERROR; //CO20191004
    vector<double> conduction_band_min;     //CO20191004
    double         conduction_band_min_net; //CO20191004
    vector<double> valence_band_max;        //CO20191004
    double         valence_band_max_net;    //CO20191004
    vector<double> Egap;                    //CO20191004
    double         Egap_net;                //CO20191004
    vector<double> Egap_fit;                //CO20191004
    double         Egap_fit_net;            //CO20191004
    vector<string> Egap_type;               //CO20191004
    string         Egap_type_net;           //CO20191004
    bool GetProperties(const stringstream& stringstreamIN,bool=TRUE);       // get everything QUIET
    bool GetProperties(const string& stringIN,bool=TRUE);                   // get everything QUIET
    bool GetPropertiesFile(const string& fileIN,bool=TRUE);                 // get everything QUIET
    bool GetPropertiesUrlFile(const string& url,const string& file,bool=TRUE); // get everything from an aflowlib entry
    void convertSpinOFF2ON(); //CO20191217 - copies everything from spin channel 1 to spin channel 2
    bool checkDOS(string& ERROR_out) const;  //CO20191010
    bool GetBandGap(double EFERMI=AUROSTD_NAN,double efermi_tol=AUROSTD_NAN,double energy_tol=1e-3,double occ_tol=1e-4); //CO20191110
    deque<deque<deque<deque<double> > > > GetVDOSSpecies(const xstructure& xstr) const; //vDOS.at(species).at(spin).at(energy_number)  //CO20191110
    deque<deque<deque<deque<double> > > > GetVDOSSpecies(deque<int> num_each_type) const; //vDOS.at(species).at(spin).at(energy_number)  //CO20191110
    friend ostream& operator<<(ostream&, const xDOSCAR&);  //ME20190623
  private:                                                        //
    void free();                                                  // free space
    void copy(const xDOSCAR& b);                                  //
};
//-------------------------------------------------------------------------------------------------
class xEIGENVAL : public xStream { //CO20200404 - xStream integration for logging
  public:
    xEIGENVAL(ostream& oss=cout);                           // default, just allocate  //CO20200404 - xStream integration for logging
    xEIGENVAL(ofstream& FileMESSAGE,ostream& oss=cout);     // constructor from filename QUIET //CO20200404 - xStream integration for logging
    xEIGENVAL(const string& fileIN,bool=TRUE,ostream& oss=cout);                          // constructor from filename QUIET //CO20200404 - xStream integration for logging
    xEIGENVAL(const string& fileIN,ofstream& FileMESSAGE,bool=TRUE,ostream& oss=cout);    // constructor from filename QUIET //CO20200404 - xStream integration for logging
    bool initialize(const string& fileIN,ostream& oss,bool=TRUE);  //ME20200427  //CO20200508
    bool initialize(const string& fileIN,ofstream& FileMESSAGE,ostream& oss,bool=TRUE);  //ME20200427  //CO20200508
    bool initialize(const string& fileIN, bool=TRUE); //ME20200427

    xEIGENVAL(const xEIGENVAL& b);                                // constructor copy
    ~xEIGENVAL();                                                 // kill everything
    const xEIGENVAL& operator=(const xEIGENVAL &b);               // copy
    void clear(void);                                             // clear

    bool m_initialized;  //CO20200404 - xStream integration for logging

    // CONTENT
    string content;vector<string> vcontent;string filename;       // the content, and lines of it
    string title;
    uint number_atoms;  //ME20190623
    uint number_loops;  //ME20190623
    uint spin;
    double Vol,POTIM;
    xvector<double> lattice;
    double temperature;
    uint number_electrons,number_kpoints,number_bands;
    deque<double> vweight;                                        // vweight.at(kpoint number)
    deque<xvector<double> > vkpoint;                              // vkpoint.at(kpoint number)[1,2,3]=xyz.
    deque<deque<deque<double> > > venergy;                        // venergy.at(kpoint number).at(band number).at(spin number)
    string carstring;  //ME20190620 - the fourth line of the EIGENVAL file
    bool GetProperties(const stringstream& stringstreamIN,bool=TRUE);       // get everything QUIET
    bool GetProperties(const string& stringIN,bool=TRUE);                   // get everything QUIET
    bool GetPropertiesFile(const string& fileIN,bool=TRUE);                 // get everything QUIET
    bool GetPropertiesUrlFile(const string& url,const string& file,bool=TRUE); // get everything from an aflowlib entry
    double energy_max;  //ME20190614
    double energy_min;  //ME20190614
    friend ostream& operator<<(ostream&, const xEIGENVAL&);  //ME20190623
  private:                                                        //
    void free();                                                  // free space
    void copy(const xEIGENVAL& b);                                //
};
//-------------------------------------------------------------------------------------------------
class xPOTCAR : public xStream { //CO20200404 - xStream integration for logging
  public:
    xPOTCAR(ostream& oss=cout);                           // default, just allocate  //CO20200404 - xStream integration for logging
    xPOTCAR(ofstream& FileMESSAGE,ostream& oss=cout);     // constructor from filename QUIET //CO20200404 - xStream integration for logging
    xPOTCAR(const string& fileIN,bool=TRUE,ostream& oss=cout);                          // constructor from filename QUIET //CO20200404 - xStream integration for logging
    xPOTCAR(const string& fileIN,ofstream& FileMESSAGE,bool=TRUE,ostream& oss=cout);    // constructor from filename QUIET //CO20200404 - xStream integration for logging
    bool initialize(const string& fileIN,ostream& oss,bool=TRUE);  //ME20200427  //CO20200508
    bool initialize(const string& fileIN,ofstream& FileMESSAGE,ostream& oss,bool=TRUE);  //ME20200427  //CO20200508
    bool initialize(const string& fileIN, bool=TRUE); //ME20200427

    ~xPOTCAR();                                                   // kill everything
    xPOTCAR(const xPOTCAR& b);                                    // constructor copy
    const xPOTCAR& operator=(const xPOTCAR &b);                   // copy
    void clear(void);                                             // clear

    bool m_initialized;  //CO20200404 - xStream integration for logging

    // CONTENT
    string content;vector<string> vcontent;                       // the content and the lines
    string filename;                                              // the filename - THIS IS A GLOBAL PROPERTY OF THE WHOLE POTCAR
    string title;
    bool   POTCAR_PAW;
    string POTCAR_TYPE;
    bool   POTCAR_KINETIC;
    bool   POTCAR_GW;
    bool   POTCAR_AE;
    double ENMAX;vector<double> vENMAX;                            // eV
    double ENMIN;vector<double> vENMIN;                            // eV
    double POMASS_sum,POMASS_min,POMASS_max;vector<double> vPOMASS;// mass
    double ZVAL_sum,ZVAL_min,ZVAL_max;vector<double> vZVAL;        // valence
    double EATOM_min,EATOM_max;vector<double> vEATOM;              // eV
    double RCORE_min,RCORE_max;vector<double> vRCORE;              // outmost cutoff radius
    double RWIGS_min,RWIGS_max;vector<double> vRWIGS;              // wigner-seitz radius (au A)
    double EAUG_min,EAUG_max;vector<double> vEAUG;                 // augmentation
    double RAUG_min,RAUG_max;vector<double> vRAUG;                 // augmentation
    double RMAX_min,RMAX_max;vector<double> vRMAX;                 // unicity
    vector<string> vTITEL;                                         // unicity
    vector<string> vLEXCH;                                         // unicity
    string pp_type;
    vector<string> species;                                        // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
    vector<int>    species_Z;                                      // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
    vector<string> species_pp;                                     // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
    vector<string> species_pp_type;                                // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
    vector<string> species_pp_version;                             // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
    vector<string> species_pp_AUID;                                // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
    vector<string> species_pp_AUID_collisions;                     // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
    vector<double> species_pp_groundstate_energy;                  // meV/atom
    vector<string> species_pp_groundstate_structure;               // name that we have, maybe ANRL
    bool GetProperties(const stringstream& stringstreamIN,bool=TRUE);       // get everything QUIET
    bool GetProperties(const string& stringIN,bool=TRUE);                   // get everything QUIET
    bool GetPropertiesFile(const string& fileIN,bool=TRUE);                 // get everything QUIET
    bool GetPropertiesUrlFile(const string& url,const string& file,bool=TRUE); // get everything from an aflowlib entry
    // objects/functions for references energies defined only with one specie
    // [OBSOLETE] vector<string> vsymbol;                                        // Ta 
    // [OBSOLETE] vector<string> vname;                                          // Ta_pv // FIX COPY CONSTRUCTOR   
    // [OBSOLETE] vector<string> vdate;                            // 07Sep2000  // FIX COPY CONSTRUCTOR
    string AUID;                                                   // crc32 - THIS IS A GLOBAL PROPERTY OF THE WHOLE POTCAR
    friend ostream& operator<<(ostream &,const xPOTCAR&);          // print // FIX COPY CONSTRUCTOR
    // xPOTCAR xPOTCAR_initialize(uint Z);                         // function to clean up the name // FIX COPY CONSTRUCTOR
  private:                                                         //
    void free();                                                   // free space
    void copy(const xPOTCAR& b);                                   //
};

extern std::vector<xPOTCAR> vxpseudopotential;        // store starting from ONE
uint xPOTCAR_Initialize(void);
bool xPOTCAR_PURE_Printer(xPOTCAR& xPOT,ostream& oss,bool LVERBOSE=FALSE);
xPOTCAR xPOTCAR_Finder(vector<string>& species_pp_AUID,vector<string>& species_pp_AUID_collisions,const string& TITEL,const string& LEXCH,const double& EATOM,const double& RMAX,bool LVERBOSE=FALSE);
xPOTCAR xPOTCAR_Finder(const string& AUID,bool LVERBOSE=FALSE);
bool xPOTCAR_EnthalpyReference_AUID(string AUID,string METAGGA=""); // returns if available
bool xPOTCAR_EnthalpyReference_AUID(string AUID,string METAGGA,string& gs,double& enthalpy_atom,double& volume_atom,double& spin_atom);

// -------------------------------------------------------------------------------------------------
class xVASPRUNXML : public xStream { //CO20200404 - xStream integration for logging
  public:
    xVASPRUNXML(ostream& oss=cout);                         // default, just allocate  //CO20200404 - xStream integration for logging
    xVASPRUNXML(ofstream& FileMESSAGE,ostream& oss=cout);   // default, just allocate  //CO20200404 - xStream integration for logging
    xVASPRUNXML(const string& fileIN,bool=TRUE,ostream& oss=cout);                        // constructor from filename, QUIET  //CO20200404 - xStream integration for logging
    xVASPRUNXML(const string& fileIN,ofstream& FileMESSAGE,bool=TRUE,ostream& oss=cout);  // constructor from filename, QUIET  //CO20200404 - xStream integration for logging
    bool initialize(const string& fileIN,ostream& oss,bool=TRUE);  //ME20200427  //CO20200508
    bool initialize(const string& fileIN,ofstream& FileMESSAGE,ostream& oss,bool=TRUE);  //ME20200427  //CO20200508
    bool initialize(const string& fileIN,bool=TRUE);  //ME20200427  //CO20200508

    xVASPRUNXML(const xVASPRUNXML& b);                            // constructor copy
    ~xVASPRUNXML();                                               // kill everything
    const xVASPRUNXML& operator=(const xVASPRUNXML &b);           // copy
    void clear(void);                                             // clear

    bool m_initialized;  //CO20200404 - xStream integration for logging

    // CONTENT
    string content;vector<string> vcontent;string filename;       // the content, and lines of it
    double natoms;                                                // for aflowlib_libraries.cpp
    xmatrix<double> stress;                                       // for aflowlib_libraries.cpp
    vector<aurostd::xvector<double> > vkpoint;                    // for aflowlib_libraries.cpp
    vector<aurostd::xvector<double> > vweights;                   // for aflowlib_libraries.cpp
    vector<aurostd::xvector<double> > vforces;                    // for aflowlib_libraries.cpp
    bool GetProperties(const stringstream& stringstreamIN,bool=TRUE);       // get everything QUIET
    bool GetProperties(const string& stringIN,bool=TRUE);                   // get everything QUIET
    bool GetPropertiesFile(const string& fileIN,bool=TRUE);                 // get everything QUIET
    bool GetPropertiesUrlFile(const string& url,const string& file,bool=TRUE); // get everything from an aflowlib entry
    bool GetForces(const string&, bool=true);  //ME20190204
    bool GetForcesFile(const string&, bool=true);  //ME20190204
    bool GetForces(stringstream&, bool=true);  //ME20190204
  private:                       //
    void free();                 // free space
    void copy(const xVASPRUNXML& b); //
};
// -------------------------------------------------------------------------------------------------
class xIBZKPT : public xStream { //CO20200404 - xStream integration for logging
  public:
    xIBZKPT(ostream& oss=cout);                         // default, just allocate  //CO20200404 - xStream integration for logging
    xIBZKPT(ofstream& FileMESSAGE,ostream& oss=cout);   // default, just allocate  //CO20200404 - xStream integration for logging
    xIBZKPT(const string& fileIN,bool=TRUE,ostream& oss=cout);                        // constructor from filename, QUIET  //CO20200404 - xStream integration for logging
    xIBZKPT(const string& fileIN,ofstream& FileMESSAGE,bool=TRUE,ostream& oss=cout);  // constructor from filename, QUIET  //CO20200404 - xStream integration for logging
    bool initialize(const string& fileIN,ostream& oss,bool=TRUE);  //ME20200427  //CO20200508
    bool initialize(const string& fileIN,ofstream& FileMESSAGE,ostream& oss,bool=TRUE);  //ME20200427  //CO20200508
    bool initialize(const string& fileIN,bool=TRUE);  //ME20200427  //CO20200508

    xIBZKPT(const xIBZKPT& b);                                    // constructor copy
    ~xIBZKPT();                                                   // kill everything
    const xIBZKPT& operator=(const xIBZKPT &b);                   // copy
    void clear(void);                                             // clear

    bool m_initialized;  //CO20200404 - xStream integration for logging

    // CONTENT
    string content;vector<string> vcontent;string filename;       // the content, and lines of it
    uint nweights;                                                // for aflowlib_libraries.cpp
    uint nkpoints_irreducible;                                    // for aflowlib_libraries.cpp
    vector<aurostd::xvector<double> > vkpoint;                    // for aflowlib_libraries.cpp
    vector<uint> vweights;                                        // for aflowlib_libraries.cpp
    uint ntetrahedra;                                             // for aflowlib_libraries.cpp
    double wtetrahedra;                                           // for aflowlib_libraries.cpp
    vector<aurostd::xvector<int> > vtetrahedra;                   // for aflowlib_libraries.cpp
    bool GetProperties(const stringstream& stringstreamIN,bool=TRUE);       // get everything QUIET
    bool GetProperties(const string& stringIN,bool=TRUE);                   // get everything QUIET
    bool GetPropertiesFile(const string& fileIN,bool=TRUE);                 // get everything QUIET
    bool GetPropertiesUrlFile(const string& url,const string& file,bool=TRUE); // get everything from an aflowlib entry
  private:                       //
    void free();                 // free space
    void copy(const xIBZKPT& b); //
};
// -------------------------------------------------------------------------------------------------
class xKPOINTS : public xStream { //CO20200404 - xStream integration for logging
  public:
    xKPOINTS(ostream& oss=cout);                         // default, just allocate  //CO20200404 - xStream integration for logging
    xKPOINTS(ofstream& FileMESSAGE,ostream& oss=cout);   // default, just allocate  //CO20200404 - xStream integration for logging
    xKPOINTS(const string& fileIN,bool=TRUE,ostream& oss=cout);                        // constructor from filename, QUIET  //CO20200404 - xStream integration for logging
    xKPOINTS(const string& fileIN,ofstream& FileMESSAGE,bool=TRUE,ostream& oss=cout);  // constructor from filename, QUIET  //CO20200404 - xStream integration for logging
    bool initialize(const string& fileIN,ostream& oss,bool=TRUE);  //ME20200427  //CO20200508
    bool initialize(const string& fileIN,ofstream& FileMESSAGE,ostream& oss,bool=TRUE);  //ME20200427  //CO20200508
    bool initialize(const string& fileIN,bool=TRUE);  //ME20200427  //CO20200508

    xKPOINTS(const xKPOINTS& b);                                   // constructor copy
    ~xKPOINTS();                                                    // kill everything
    const xKPOINTS& operator=(const xKPOINTS &b);                  // copy
    void clear(void);                                              // clear

    bool m_initialized;  //CO20200404 - xStream integration for logging

    // CONTENT
    string content;vector<string> vcontent;string filename;        // the content, and lines of it
    string title;                                                  // first line
    int mode;                                                      // sort of mode
    string grid_type;                                              // if grid specified
    bool is_KPOINTS_NNN,is_KPOINTS_PATH;                           // control parameters
    xvector<int>    nnn_kpoints; // N*N*N                          // triplet of kpoints
    xvector<double> ooo_kpoints; // ORIGIN                         // triplet of origin
    int nkpoints;                                                  // total kpoints
    string path_mode,path;vector<string> vpath;int path_grid;      // path if any
    vector<xvector<double> > vkpoints;                             //ME20190614 - k-point coordinates of the path
    bool GetProperties(const stringstream& stringstreamIN,bool=TRUE);       // get everything QUIET
    bool GetProperties(const string& stringIN,bool=TRUE);                   // get everything QUIET
    bool GetPropertiesFile(const string& fileIN,bool=TRUE);                 // get everything QUIET
    bool GetPropertiesUrlFile(const string& url,const string& file,bool=TRUE); // get everything from an aflowlib entry
    friend ostream& operator<<(ostream&, const xKPOINTS&);  //ME20190623
    string createStandardTitlePath(const xstructure&);  //ME20190623
  private:                       //
    void free();                 // free space
    void copy(const xKPOINTS& b); //
};
// -------------------------------------------------------------------------------------------------
class xCHGCAR : public xStream { //CO20200404 - xStream integration for logging
  public:
    xCHGCAR(ostream& oss=cout);                         // default, just allocate  //CO20200404 - xStream integration for logging
    xCHGCAR(ofstream& FileMESSAGE,ostream& oss=cout);   // default, just allocate  //CO20200404 - xStream integration for logging
    xCHGCAR(const string& fileIN,bool=TRUE,ostream& oss=cout);                        // constructor from filename, QUIET  //CO20200404 - xStream integration for logging
    xCHGCAR(const string& fileIN,ofstream& FileMESSAGE,bool=TRUE,ostream& oss=cout);  // constructor from filename, QUIET  //CO20200404 - xStream integration for logging
    bool initialize(const string& fileIN,ostream& oss,bool=TRUE);  //ME20200427  //CO20200508
    bool initialize(const string& fileIN,ofstream& FileMESSAGE,ostream& oss,bool=TRUE);  //ME20200427  //CO20200508
    bool initialize(const string& fileIN,bool=TRUE);  //ME20200427  //CO20200508

    xCHGCAR(const xCHGCAR& b);                                     // constructor copy
    ~xCHGCAR();                                                    // kill everything
    const xCHGCAR& operator=(const xCHGCAR &b);                    // copy
    void clear(void);                                              // clear

    bool m_initialized;  //CO20200404 - xStream integration for logging

    // CONTENT
    string content;vector<string> vcontent;string filename;        // the content, and lines of it
    xvector<int>     grid; // N*N*N                                // triplet of grid
    vector<string>   vstring; // ORIGIN                            // string of values
    xvector<double>  vvalues; // ORIGIN                            // xvector of values
    //[OBSOLETE ME20180705]xtensor3<double> tvalues; // ORIGIN       // xtensor of values
    xtensor<double> tvalues; // ORIGIN                             // xtensor of values ME20180705
    bool GetProperties(const stringstream& stringstreamIN,bool=TRUE);          // get everything QUIET
    bool GetProperties(const string& stringIN,bool=TRUE);                      // get everything QUIET
    bool GetPropertiesFile(const string& fileIN,bool=TRUE);                    // get everything QUIET
    bool GetPropertiesUrlFile(const string& url,const string& file,bool=TRUE); // get everything from an aflowlib entry
  private:                       //
    void free();                 // free space
    void copy(const xCHGCAR& b); //
};

class xQMVASP : public xStream {  //CO20191110 //CO20200404 - xStream integration for logging
  public:
    xQMVASP(ostream& oss=cout);                           // default, just allocate  //CO20200404 - xStream integration for logging
    xQMVASP(ofstream& FileMESSAGE,ostream& oss=cout);     // constructor from filename QUIET //CO20200404 - xStream integration for logging
    xQMVASP(const string& fileIN,bool=TRUE,ostream& oss=cout);                          // constructor from filename QUIET //CO20200404 - xStream integration for logging
    xQMVASP(const string& fileIN,ofstream& FileMESSAGE,bool=TRUE,ostream& oss=cout);    // constructor from filename QUIET //CO20200404 - xStream integration for logging
    bool initialize(const string& fileIN,ostream& oss,bool=TRUE);  //ME20200427  //CO20200508
    bool initialize(const string& fileIN,ofstream& FileMESSAGE,ostream& oss,bool=TRUE);  //ME20200427  //CO20200508
    bool initialize(const string& fileIN, bool=TRUE); //ME20200427

    ~xQMVASP();                                                    // kill everything
    xQMVASP(const xQMVASP& b);                                     // constructor copy
    const xQMVASP& operator=(const xQMVASP &b);                    // copy
    void clear(void);                                              // clear

    bool m_initialized;  //CO20200404 - xStream integration for logging

    // CONTENT
    string content;vector<string> vcontent;string filename;        // the content, and lines of it
    double H_atom_relax;
    double H_atom_static;
    vector<aurostd::xvector<double> > vforces;                     // for APL - only one (no relax vs. static), get most relaxed forces  //CO20191112
    bool GetProperties(const stringstream& stringstreamIN,bool=TRUE);          // get everything QUIET
    bool GetProperties(const string& stringIN,bool=TRUE);                      // get everything QUIET
    bool GetPropertiesFile(const string& fileIN,bool=TRUE);                    // get everything QUIET
    bool GetPropertiesUrlFile(const string& url,const string& file,bool=TRUE); // get everything from an aflowlib entry
  private:                       //
    void free();                 // free space
    void copy(const xQMVASP& b); //
};

// -------------------------------------------------------------------------------------------------
// aflow_kaims.cpp
namespace KBIN {
  _aimsflags AIMS_Get_AIMSflags_from_AflowIN(string& AflowIn,_aflags& aflags,_kflags& kflags);
  _aimsflags AIMS_Get_AIMSflags_from_AflowIN(string& AflowIn,ofstream& FileMESSAGE,_aflags& aflags,_kflags& kflags);
  bool AIMS_Directory(ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags);
}
// -------------------------------------------------------------------------------------------------
// aflow_iaims.cpp
namespace KBIN {
  bool AIMS_Produce_INPUT(_xaims& xaims,string AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_aimsflags &aimsflags);
  bool AIMS_Modify_INPUT(_xaims& xaims,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_aimsflags &aimsflags);
  bool AIMS_Write_INPUT(_xaims& xaims,_aimsflags &aimsflags);
  bool AIMS_Write_CONTROL(_xaims& xaims,_aimsflags &aimsflags);
  bool AIMS_Write_GEOM(_xaims& xaims,_aimsflags &aimsflags);
  bool AIMS_Produce_CONTROL(_xaims& xaims,string AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_aimsflags &aimsflags);
  bool AIMS_Modify_CONTROL(_xaims& xaims,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_aimsflags &aimsflags);
  bool AIMS_Reread_CONTROL(_xaims& xaims,ofstream &FileMESSAGE,_aflags &aflags);
  bool AIMS_Produce_GEOM(_xaims& xaims,string AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_aimsflags &aimsflags);
  bool AIMS_Produce_GEOM(_xaims& xaims);
  bool AIMS_Modify_GEOM(_xaims& xaims,string AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_aimsflags &aimsflags);
  bool AIMS_Reread_GEOM(_xaims& xaims,ofstream &FileMESSAGE,_aflags &aflags);
  bool XAIMS_CONTROL_PREPARE_GENERIC(string command,_xaims& xaims,_aimsflags& aimsflags,string svalue,int ivalue,double dvalue,bool OPTION);
  void XAIMS_CONTROL_REMOVE_ENTRY(_xaims& xaims,string ENTRY,string COMMENT,bool VERBOSE);
}
// -------------------------------------------------------------------------------------------------
// aflow_oaims.cpp
class xAIMSOUT;
class xAIMSOUT {
  public:
    xAIMSOUT();                                                    // default, just allocate
    ~xAIMSOUT();                                                   // kill everything
    xAIMSOUT(const string& fileIN,bool=TRUE);                      // constructor from filename, QUIET
    xAIMSOUT(const xAIMSOUT& b);                                    // constructor copy
    const xAIMSOUT& operator=(const xAIMSOUT &b);                   // copy
    void clear(void);                                             // clear
    // CONTENT
    string content;vector<string> vcontent;string filename;       // the content, and lines of it
    vector<aurostd::xvector<double> > vforces;                    // for aflowlib_libraries.cpp
    double natoms;
    //[CO20200404 - OBSOLETE]string ERROR;
    bool GetProperties(const stringstream& stringstreamIN,bool=TRUE);          // get everything QUIET
    bool GetProperties(const string& stringIN,bool=TRUE);                      // get everything QUIET
    bool GetPropertiesFile(const string& fileIN,bool=TRUE);                    // get everything QUIET
    bool GetPropertiesFile(const string& fileIN,uint natoms_check,bool);       // get everything QUIET
    bool GetPropertiesUrlFile(const string& url,const string& file,bool=TRUE); // get everything from an aflowlib entry
  private:                       //
    void free();                 // free space
    void copy(const xAIMSOUT& b); //
};
// -----------------------------------------------------------------------------------------------
bool PrintBandGap       (string& WorkDir, ostream &oss);
bool PrintBandGap_DOS   (string& WorkDir, ostream &oss); //CO20191110
bool PrintEffectiveMass (string& WorkDir, ostream &oss);
bool PrintEigCurv       (string& WorkDir, ostream &oss);
// -----------------------------------------------------------------------------------------------
bool ParseKPOINTS(stringstream& File_Kpoints, int& GRIDS, vector<xvector<double> >& special_kpts, vector<xvector<double> >& unique_kpts, vector<int>& repeat_kpts_num);
bool AdjacencyList_KPT(vector<xvector<double> >& special_kpts,vector<xvector<double> >& unique_kpts,vector<xvector<int> >& connect_kpts,vector<int>& connect_kpts_num);
bool AdjacencyList_EIG(vector<xvector<double> >& unique_kpts,vector<xvector<int> >& connect_kpts,vector<int>& connect_kpts_num,xEIGENVAL& xeigenval,vector<xvector<double> >& unique_kpts_EIG,vector<xvector<int> >& connect_kpts_EIG,vector<xvector<double> >& vkpoint_eig);
bool RepeatsList(vector<xvector<double> >& unique_kpts_EIG,vector<int>& repeat_kpts_num,vector<xvector<double> >& vkpoint_eig,vector<xvector<int> >& repeat_kpts_EIG);
bool VertexPaths(vector<xvector<int> >& repeat_kpts_EIG,vector<xvector<int> >& connect_kpts_EIG,vector<int>& repeat_kpts_num,int& GRIDS,vector<xvector<int> >& vrtx_path);
bool RepeatedEdges(vector<xvector<int> >& vrtx_path,vector<xvector<int> >& repeat_kpts_EIG,vector<int>& repeat_kpts_num,vector<xvector<int> >& ndx_edges);
bool VertexBranches(vector<xvector<int> >& ndx_edges,vector<int>& repeat_kpts_num,vector<xvector<int> >& repeat_kpts_EIG,vector<vector<xvector<int> > >& branches);
bool PathDataStuct(xEIGENVAL& xeigenval,vector<xvector<double> >& vkpoint_eig,vector<vector<xvector<int> > >& branches,vector<vector< vector<int> > >& branches_indx,vector<vector< vector<xvector<double> > > >& branches_kpts,vector<vector< vector<vector<vector<double> > > > >& branches_bnds);
bool IBZextrema(xEIGENVAL& xeigenval, vector<xvector<double> >& vkpoint_eig, vector<vector<xvector<int> > >& branches);
void CompareDoublesChar(bool& MATCH, double& number1, double& number2);
void CompareEdges(vector<vector<xvector<int> > >& branches, vector<xvector<int> >& vertex_edges, xvector<int>& test_edge, bool& MATCH);
void NaiveCurvatures(xvector<double>& eigvec, vector<xvector<double> >& posvec, vector<double>& curvature);
double StencilLinear1D(vector<xvector<double> >& positions, xvector<double>& eigenvals);
//-------------------------------------------------------------------------------------------------
struct kEn_st {
  xvector<double> kpoint;
  double energy[2];
  int band_index;
  int band_type; // 0 -- valence band; 1 -- conduction band
};
#define _SIGMA 1.0 // default standard deviation of input data
// range of energy point to fit the ellipse curve
const double _FIT_ENERGY_RANGE = 0.026; // eV range of band
const int _FIT_POINTS_NUMBER = 8; // minimum fit points in Irreducible BZ
//range of band extremes to determine the number of bands for effective mass calculations
const double _BANDS_ENERGY_RANGE = 0.026; // eV
// used to determine cluster of points can be changed to other values
const double _BANDS_PARAMETER_MIN_RATIO = 0.2;
// factor unit
// mass is in unit of electron mass
const double _MASS_FACTOR = 3.80998; // hbar^2*10^{20}/(2.0*me*eV)
bool comparison_kEn_str_up          (const kEn_st& k1, const kEn_st& k2);
bool comparison_kEn_str_dn          (const kEn_st& k1, const kEn_st& k2);
bool comparison_kEn_str_position    (const kEn_st& k1, const kEn_st& k2);
bool comparison_kEn_str_band_type_up(const kEn_st& k1, const kEn_st& k2);
bool comparison_kEn_str_band_type_dn(const kEn_st& k1, const kEn_st& k2);
bool is_equal_position_kEn_str      (const kEn_st& k1, const kEn_st& k2);
bool near_to                        (const xvector<double> & k1, const xvector<double> & k2, const vector<double> & max_distance);
// [OBSOLETE] bool GetEffectiveMass(xOUTCAR& outcar,xDOSCAR& doscar,xEIGENVAL& eigenval,xstructure xstr,ostream& oss,const bool& osswrite);
//-------------------------------------------------------------------------------------------------
//ME20190614 - plotter functions
namespace plotter {
  // Plot setup --------------------------------------------------------------
  // Plot options
  aurostd::xoption getPlotOptions(const aurostd::xoption&, const string&, bool=false);
  aurostd::xoption getPlotOptionsEStructure(const aurostd::xoption&, const string&, bool=false);
  aurostd::xoption getPlotOptionsPhonons(const aurostd::xoption&, const string&);

  // Plot functions
  void generateHeader(stringstream&, const aurostd::xoption&, bool=false);
  void savePlotGNUPLOT(const aurostd::xoption&, const stringstream&);
  void setFileName(aurostd::xoption&, string="");
  void setTitle(aurostd::xoption&,ostream& oss=cout); //CO20200404
  void setTitle(aurostd::xoption&,ofstream& FileMESSAGE,ostream& oss=cout); //CO20200404
  string formatDefaultPlotTitle(const aurostd::xoption&,ostream& oss=cout); //CO20200404
  string formatDefaultPlotTitle(const aurostd::xoption&,ofstream& FileMESSAGE,ostream& oss=cout); //CO20200404
  vector<double> getCompositionFromHTQCPrototype(const string&, const string&);  //ME20190813
  vector<double> getCompositionFromANRLPrototype(const string&);
  string formatDefaultTitlePOCC(const aurostd::xoption&,ostream& oss=cout); //CO20200404
  string formatDefaultTitlePOCC(const aurostd::xoption&,ofstream& FileMESSAGE,ostream& oss=cout); //CO20200404
  string formatDefaultTitlePOCC_20191004(const aurostd::xoption&,ostream& oss=cout); //CO20191110 //CO20200404
  string formatDefaultTitlePOCC_20191004(const aurostd::xoption&,ofstream& FileMESSAGE,ostream& oss=cout); //CO20191110 //CO20200404
  string formatDefaultTitlePOCC_20190101(const aurostd::xoption&,ostream& oss=cout);  //CO20200404
  vector<double> getCompositionFromPoccString(const string&, bool&);

  // Electronic structure ----------------------------------------------------
  void patchDefaultTitleAFLOWIN(xoption& plotoptions);  //CO20191110
  // Plot functions
  void PLOT_DOS(aurostd::xoption&,ostream& oss=cout); //CO20200404
  void PLOT_DOS(aurostd::xoption&,ofstream& FileMESSAGE,ostream& oss=cout); //CO20200404
  void PLOT_DOS(aurostd::xoption&,const xDOSCAR&,ostream& oss=cout); //CO20191110 //CO20200404
  void PLOT_DOS(aurostd::xoption&,const xDOSCAR&,ofstream& FileMESSAGE,ostream& oss=cout); //CO20191110 //CO20200404
  void PLOT_DOS(aurostd::xoption&,stringstream&,ostream& oss=cout); //CO20200404
  void PLOT_DOS(aurostd::xoption&,stringstream&,ofstream& FileMESSAGE,ostream& oss=cout); //CO20200404
  void PLOT_DOS(aurostd::xoption&,stringstream&,const xDOSCAR&,ostream& oss=cout);  //CO20191110  //CO20200404
  void PLOT_DOS(aurostd::xoption&,stringstream&,const xDOSCAR&,ofstream& FileMESSAGE,ostream& oss=cout);  //CO20191110  //CO20200404

  void PLOT_PDOS(aurostd::xoption&,ostream& oss=cout);  //CO20200404
  void PLOT_PDOS(aurostd::xoption&,ofstream& FileMESSAGE,ostream& oss=cout);  //CO20200404
  void PLOT_PDOS(aurostd::xoption&, const xDOSCAR&,ostream& oss=cout); //CO20191110 //CO20200404
  void PLOT_PDOS(aurostd::xoption&, const xDOSCAR&,ofstream& FileMESSAGE,ostream& oss=cout); //CO20191110 //CO20200404
  void PLOT_PDOS(aurostd::xoption&, stringstream&,ostream& oss=cout); //CO20200404
  void PLOT_PDOS(aurostd::xoption&, stringstream&,ofstream& FileMESSAGE,ostream& oss=cout); //CO20200404
  void PLOT_PDOS(aurostd::xoption&, stringstream&, const xDOSCAR&,ostream& oss=cout);  //CO20191110 //CO20200404
  void PLOT_PDOS(aurostd::xoption&, stringstream&, const xDOSCAR&,ofstream& FileMESSAGE,ostream& oss=cout);  //CO20191110 //CO20200404

  void PLOT_BAND(aurostd::xoption&,ostream& oss=cout);  //CO20200404
  void PLOT_BAND(aurostd::xoption&,ofstream& FileMESSAGE,ostream& oss=cout);  //CO20200404
  void PLOT_BAND(aurostd::xoption&, stringstream&,ostream& oss=cout); //CO20200404
  void PLOT_BAND(aurostd::xoption&, stringstream&,ofstream& FileMESSAGE,ostream& oss=cout); //CO20200404
  void BANDDOS2JSON(ostream&, string);
  void PLOT_BANDDOS(aurostd::xoption&,ostream& oss=cout); //CO20200404
  void PLOT_BANDDOS(aurostd::xoption&,ofstream& FileMESSAGE,ostream& oss=cout); //CO20200404
  void PLOT_BANDDOS(aurostd::xoption&, stringstream&,ostream& oss=cout);  //CO20200404
  void PLOT_BANDDOS(aurostd::xoption&, stringstream&,ofstream& FileMESSAGE,ostream& oss=cout);  //CO20200404

  // Helper functions
  xstructure getStructureWithNames(const aurostd::xoption&,const string& carstring="CAR",ostream& oss=cout);  //CO20191110 //CO20200404
  xstructure getStructureWithNames(const aurostd::xoption&,ofstream& FileMESSAGE,const string& carstring="CAR",ostream& oss=cout);  //CO20191110 //CO20200404
  string getLatticeFromKpointsTitle(const string&);
  void shiftEfermiToZero(xEIGENVAL&, double);
  void setEMinMax(aurostd::xoption&, double, double);

  // DOS
  void generateDosPlot(stringstream&, const xDOSCAR&, const aurostd::xoption&,ostream& oss=cout);  //CO20200404
  void generateDosPlot(stringstream&, const xDOSCAR&, const aurostd::xoption&,ofstream& FileMESSAGE,ostream& oss=cout);  //CO20200404

  // Bands
  void generateBandPlot(stringstream&, const xEIGENVAL&, const xKPOINTS&, const xstructure&, const aurostd::xoption&);
  string convertKPointLabel(const string&, const string&);
  string convertKPointLetter(string, const string&);

  // Gnuplot
  void generateDosPlotGNUPLOT(stringstream&, const xDOSCAR&, const deque<double>&,
      const deque<deque<deque<double> > >&, const vector<string>&, const aurostd::xoption&);
  double getDosLimits(const aurostd::xoption&, const xDOSCAR&, const deque<deque<deque<double> > >&, const deque<double>&);
  void generateBandPlotGNUPLOT(stringstream&, const xEIGENVAL&, const vector<double>&,
      const vector<double>&, const vector<string>&, const aurostd::xoption&);
  string getFormattedUnit(const string&);

  // Phonons -----------------------------------------------------------------
  void PLOT_PHDOS(aurostd::xoption&,ostream& oss=cout); //CO20200404
  void PLOT_PHDOS(aurostd::xoption&,ofstream& FileMESSAGE,ostream& oss=cout); //CO20200404
  void PLOT_PHDOS(aurostd::xoption&, stringstream&,ostream& oss=cout);  //CO20200404
  void PLOT_PHDOS(aurostd::xoption&, stringstream&,ofstream& FileMESSAGE,ostream& oss=cout);  //CO20200404

  void PLOT_PHDISP(aurostd::xoption&,ostream& oss=cout);  //CO20200404
  void PLOT_PHDISP(aurostd::xoption&,ofstream& FileMESSAGE,ostream& oss=cout);  //CO20200404
  void PLOT_PHDISP(aurostd::xoption&, stringstream&,ofstream& FileMESSAGE,ostream& oss=cout); //CO20200404
  void PLOT_PHDISPDOS(aurostd::xoption&,ostream& oss=cout); //CO20200404
  void PLOT_PHDISPDOS(aurostd::xoption&,ofstream& FileMESSAGE,ostream& oss=cout); //CO20200404
  void PLOT_PHDISPDOS(aurostd::xoption&, stringstream&,ostream& oss=cout);  //CO20204004
  void PLOT_PHDISPDOS(aurostd::xoption&, stringstream&,ofstream& FileMESSAGE,ostream& oss=cout);  //CO20204004

  void convertEnergies(xEIGENVAL&, const string&);
  void convertEnergies(xDOSCAR&, const string&);
  double getEnergyConversionFactor(const string&);

  // Properties plotter ------------------------------------------------------
  void PLOT_THERMO(aurostd::xoption&,ostream& oss=cout);  //CO20200404
  void PLOT_THERMO(aurostd::xoption&,ofstream& FileMESSAGE,ostream& oss=cout);  //CO20200404
  void PLOT_THERMO(aurostd::xoption&, stringstream&,ostream& oss=cout); //CO20200404
  void PLOT_THERMO(aurostd::xoption&, stringstream&,ofstream& FileMESSAGE,ostream& oss=cout); //CO20200404
  void PLOT_TCOND(aurostd::xoption&,ostream& oss=cout); //CO20200404
  void PLOT_TCOND(aurostd::xoption&,ofstream& FileMESSAGE,ostream& oss=cout); //CO20200404
  void PLOT_TCOND(aurostd::xoption&, stringstream&,ostream& oss=cout);  //CO20200404
  void PLOT_TCOND(aurostd::xoption&, stringstream&,ofstream& FileMESSAGE,ostream& oss=cout);  //CO20200404

  // General plots -----------------------------------------------------------
  void plotSingleFromSet(xoption&, stringstream&, const vector<vector<double> >&, int,ostream& oss=cout); //CO20200404
  void plotSingleFromSet(xoption&, stringstream&, const vector<vector<double> >&, int,ofstream& FileMESSAGE,ostream& oss=cout); //CO20200404
  void plotMatrix(xoption& plotoptions, stringstream&,ostream& oss=cout); //CO20200404
  void plotMatrix(xoption& plotoptions, stringstream&,ofstream& FileMESSAGE,ostream& oss=cout); //CO20200404
  void setPlotLabels(aurostd::xoption&, const string&, const string&, const string&, const string&);
  vector<vector<double> > readAflowDataFile(aurostd::xoption&);
  void generatePlotGNUPLOT(stringstream&, const xoption&, const vector<vector<double> >&);
}

//-------------------------------------------------------------------------------------------------
// aflow_estructure_dos.cpp

namespace aurostd {
  int CountWordsinString(string& input);  // put aurostd
  int CountWordsinString_web(string input); // put aurostd
}

namespace estructure {
  string PEDOS_GENERATE_GNUPLOTSCRIPT(const string&,const string&,const double&,const double&,const double&,const double&,const int&,const vector<vector<vector<double> > >&,const string&);
  bool isSpecialKPOINT(string kpoint);  //CO20170830
  string fixSpecialKPOINT_GNUPLOT(string kpoint,bool json=false);  //CO20170830
  string fixSpecialKPOINT_HTML(string kpoint);  //CO20170830
  string fixSpecialKPOINT_LATEX(string kpoint);  //CO20170830
  string fixKPOINT_GNUPLOT(string kpoint,bool json=false);  //CO20170830
  string fixKPOINT_HTML(string kpoint);  //CO20170830
  string fixKPOINT_LATEX(string kpoint);  //CO20170830
  string fixKPOINT_SPECIALONLY(string kpoint);  //CO20170830
  void PLOT_BANDDOS(string options);
  void PLOT_BAND(string options);
  void PLOT_DOS(string options);
  void PLOT_PEDOS(string options);
  void PLOT_PEDOSALL(string options);
  void PLOT_PEDOSALL_AFLOWLIB(string options,_aflags& aflags);
  void PLOT_BAND2(string options);
  // [OBSOLETE]  void PLOT_BAND3(string options);
  void PLOT_BAND_SPINSPLIT(string options);
  void PLOT_DOSWEB(string options);
  // manipulation
  string changeICSDNameGunplot(string ICSDName);
  void CombineTDOSAndTOTALPDOS(const vector<vector<double> >& TDOS, const vector<vector<double> >& TOTALPDOS, vector<vector<double> >& vvDOS);
  double GET_TDOSDATA(const string& str_dir, vector<vector<double> >& TDOS);
  double GET_TDOSDATA(stringstream& ss_dosfile, stringstream& ss_outcarfile, vector<vector<double> >& TDOS);
  double GET_TOTALPDOSDATA(const string& str_dir, vector<vector<double> >& TOTALPDOS);
  double GET_TOTALPDOSDATA(stringstream& ss_dosfile, stringstream& ss_outfile, vector<vector<double> >& TOTALPDOS);
  double GET_PDOSDATA(const string& str_dir, vector<vector<vector<double> > >& PDOS);
  double GET_PDOSDATA(stringstream& ss_dosfile, stringstream& ss_outfile, vector<vector<vector<double> > >& PDOS);
  // [OBSOLETE]  void GET_DOS_DATA(vector<string>& argv);
  bool GET_DOS_DATA(stringstream& ss_dosfile, stringstream& ss_outfile, double& Efermi, vector<vector<double> >& TDOS, vector<vector<double> >& TOTALPDOS); //CO20180216
  bool GET_DOS_DATA(const string& str_dir,   double& Efermi, vector<vector<double> >& TDOS, vector<vector<double> >& TOTALPDOS, vector<vector<vector<double> > >& PDOS);  //CO20180216
  bool GET_DOS_DATA(stringstream& ss_dosfile, stringstream& ss_outfile,   double& Efermi, vector<vector<double> >& TDOS, vector<vector<double> >& TOTALPDOS, vector<vector<vector<double> > >& PDOS); //CO20180216
  void FormatSpinofPDOS(vector<vector<vector<double> > >& vvva);

  // Functions for serializing bands data to JSON
  // Added by EG
  bool DOSDATA_JSON(aurostd::xoption& vpflow,ostream& oss=cout);
  bool DOSDATA_JSON(aurostd::xoption& vpflow,string directory,stringstream& json,bool wrapping_brackets=true);
  bool BANDSDATA_JSON(aurostd::xoption& vpflow,ostream& oss=cout);
  bool BANDSDATA_JSON(aurostd::xoption& vpflow,string directory,stringstream& json,bool wrapping_brackets=true);
  //uint DOSDATA_JSON(string options);
  //uint DOSDATA_JSON(string options, ostream& json);
  //uint BANDSDATA_JSON(string options);
  //uint BANDSDATA_JSON(string options, string json_dir);
  //uint BANDSDATA_JSON(string options, ostream& json);
  string linelabel2HTML(string linelabel);
  uint inequivalentAtomsJSON( vector<vector<vector<double> > >& PDOS, vector<int>& iatoms, vector<double>& numbers, vector<string>& vspecies, ostream& json);
  uint constructInequivalentAtomPDOSJSON(vector<vector<vector<double> > >& PDOS, int iatom, ostream& json); 
  // End of bands data JSON serializers

}


// ----------------------------------------------------------------------------
// aflow_poccupation_*.cpp
//  #include "aflow_pocc.h"

// aflow_poccupation_params.cpp
namespace pocc {
  bool poccInput(); //CO20170805

  string ReturnAtomSpecies(string atom);
  string ReturnAtomSpeciesPotential(string atom);
  string ReturnUFFParameters(string atom);
  class UFFPara {
    public:
      UFFPara(); // constructor
      ~UFFPara(); // destructor
      string symbol;
      double r1,theta0,x1,D1,zeta,Z1,Vi,Uj,Xi,hard,radius; 
      void GetUFFParameters(string);
    private:
      void free(); //free space
  };
  string ReturnAtomProperties(string atom);
  //Atomic Properties Database
  class Atom {
    public:
      Atom();
      ~Atom();
      string name,symbol;
      int number; //atomic number
      double mass,radius,Xi; //atomic, weight radius /pauling electronegativity
      void GetAtomicProperties(string);
    private:
      void free();
  };
} // namespace pocc

// aflow_poccupation_forcefield.cpp
namespace pocc {
  class Bond{
    public:
      Bond();
      Bond(const Bond& b);
      ~Bond();
      _atom bgn,end;
      double length;
      void Set(xstructure , _atom , _atom );
      const Bond & operator=(const Bond &other);
      bool operator==(const Bond &other) const;
      bool operator!=(const Bond &other) const;
      friend ostream& operator<<(ostream&,const Bond&);
    private:
      void free();
      void copy(const Bond& b);
  };
  void SetUFFPara(_atom atomi, _atom atomj, double& R0, double& Kij, double& Xij, double& Dij);
  double CalculateBondEnergy(xstructure xstr, _atom atomi, _atom atomj);
  double CalculateNonBondEnergy(xstructure xstr, _atom atomi, _atom atomj);
  double CalculateUFFEnergy(xstructure xstr);
  void RemoveSameBond(vector<Bond>& Bonds_orig, vector<Bond>& Bonds_new);
  void ExtractBonds(const xstructure& xstr, deque<deque<_atom> >& neigh_mat_bonded, deque<deque<_atom> >& neigh_mat_nonbonded);
  void AnalyzeBonds(const xstructure& xstr, vector<Bond>& Bonds, vector<Bond>& NonBonds);
  void UFFENERGY(istream& input);
}

// ----------------------------------------------------------------------------
// aflow_mix.cpp  aflow_nomix.cpp   aflow_mix_pauling.cpp
#define MISCIBILITY_SYSTEM_NOT_STUDIED  3
#define MISCIBILITY_SYSTEM_SOLUTION     2
#define MISCIBILITY_SYSTEM_MISCIBLE     1
#define MISCIBILITY_SYSTEM_NOMIX        0
#define MISCIBILITY_SYSTEM_UNKNOWN     -1
#define MISCIBILITY_SYSTEM_CUTOFF     200
#define MIEDEMA_MIX_SLOPE 3.069              // Miedema Rule Table 1a Physica 100B (1980) 1-28

int MiscibilityCheck(int speciesA,int speciesB);                  // aflow_mix.cpp
int MiscibilityCheck(string speciesA,string speciesB);            // aflow_mix.cpp
int MiscibilityExperimentsCheck(int speciesA,int speciesB);       // aflow_mix.cpp
int MiscibilityExperimentsCheck(string speciesA,string speciesB); // aflow_mix.cpp
int MiscibilityMiedemaCheck(int speciesA,int speciesB);           // aflow_mix.cpp
int MiscibilityMiedemaCheck(string speciesA,string speciesB);     // aflow_mix.cpp
int MiscibilityMiedemaCheck(string system_in);                    // aflow_mix.cpp
int MiscibilityHumeRotheryCheck(int speciesA,int speciesB);       // aflow_mix.cpp
int MiscibilityHumeRotheryCheck(string speciesA,string speciesB); // aflow_mix.cpp
int MiscibilityHumeRotheryCheck(string system_in);                // aflow_mix.cpp
int MiscibilityCheck(string system_in);                           // aflow_nomix.cpp
int MiscibilityExperimentsCheck(string system_in);                    // aflow_mix_pauling.cpp

// ----------------------------------------------------------------------------
// symmetry prototypes
// aflow_symmetry.cpp
namespace SYM {
  //DX+CO START
  bool ApplyAtomValidate(const _atom &atom_in,_atom& atom_out,const _sym_op &symop,const xstructure& a); //DX+CO
  bool ApplyAtomValidate(const _atom &atom_in,_atom& atom_out,const _sym_op &symop,const xstructure& a, bool _incell_,bool roff); //DX+CO
  bool ApplyAtomValidate(const _atom &atom_in,_atom& atom_out,const _sym_op &symop,const xstructure& a,bool skew, bool _incell_,bool roff,double _eps_); //DX+CO
  bool ApplyAtomValidate(const _atom &atom_in,_atom& atom_out,const _sym_op &symop,const xmatrix<double>& lattice,const xmatrix<double>& c2f, const xmatrix<double>& f2c,bool skew, bool _incell_,bool roff,double _eps_); //DX+CO
  _atom ApplyAtom(const _atom& atom_in,const _sym_op& symop,const xstructure& str); //DX
  //DX+CO END
  _atom ApplyAtom(const _atom& atom_in,const _sym_op& symop,const xstructure& str,bool _incell_);
  //DX+CO START
  _atom ApplyAtom(const _atom& atom_in,const _sym_op& symop,const xstructure& str,bool _incell_,bool roff); //DX
  _atom ApplyAtom(const _atom& atom_in,const _sym_op& symop,const xstructure& str,bool _incell_,bool roff,bool validatePosition); //DX
  _atom ApplyAtom(const _atom& atom_in,const _sym_op& symop,const xmatrix<double>& lattice,const xmatrix<double>& c2f, const xmatrix<double>& f2c,bool skew,bool _incell_,bool roff,bool validatePosition,double eps); //DX
  _atom ApplyAtom_20161115(const _atom& atom_in,const _sym_op& symop,const xmatrix<double>& lattice,const xmatrix<double>& c2f, const xmatrix<double>& f2c,bool skew,bool _incell_,bool roff,bool validatePosition,double eps); //DX
  _atom ApplyAtom_20160101(const _atom& atom_in,const _sym_op& symop,const xstructure& str,bool _incell_); //DX
  //DX+CO END
  xvector<double> ApplyCpos(const xvector<double> &cpos_in,const _sym_op &symop,const xstructure& str,bool _incell_);
  xvector<double> ApplyCpos(const xvector<double> &cpos_in,const _sym_op &symop,const xstructure& str);
  xvector<double> ApplyFpos(const xvector<double> &fpos_in,const _sym_op &symop,const xstructure& str,bool _incell_);
  xvector<double> ApplyFpos(const xvector<double> &fpos_in,const _sym_op &symop,const xstructure& str);
  xvector<int> ApplyIJK(const xvector<int> &ijk_in,const _sym_op &symop,const xstructure& str);
  int  ApplyL(const int &l_in,const _sym_op &symop,const xstructure& str);
  xstructure ApplyXstructure(const _sym_op &symop,const xstructure& str);
  xstructure ApplyXstructure(const _sym_op &symop,const xstructure& str,bool _incell_);
  //DX+CO START
  bool AtomsEquivalent(xstructure& str,_atom& atom1,_atom& atom2); //DX
  bool AtomsEquivalent_20161115(xstructure& str,_atom& atom1,_atom& atom2); //DX
  bool AtomsEquivalent_20160101(xstructure& str,_atom& atom1,_atom& atom2); //DX
  bool AtomsEquivalent(xstructure& str, _atom& atom1, _atom& atom2,double& eps); //DX
  bool AtomsEquivalent(xstructure& str, _atom& a, _atom& b, bool skew, double tol); //DX //CO20190520 - removed pointers for bools and doubles, added const where possible
  bool AtomsEquivalent_20161115(xstructure& str,_atom& atom1,_atom& atom2,double& eps); //DX
  bool AtomsEquivalent_20160101(xstructure& str,_atom& atom1,_atom& atom2,double& eps); //DX
  bool AtomsEquivalent_Basis(xstructure& str, int atom1_indx,int atom2_indx);
  //DX+CO END
  bool CposEquivalent(const xstructure& str,const xvector<double>& cpos1,const xvector<double>& cpos2,const double& eps);
  bool FposEquivalent(const xstructure& str,const xvector<double>& fpos1,const xvector<double>& fpos2,const double& eps);
  bool CposEquivalent(const xstructure& str,const xvector<double>& cpos1,const xvector<double>& cpos2);
  bool FposEquivalent(const xstructure& str,const xvector<double>& fpos1,const xvector<double>& fpos2);
  bool TypePointGroupOperation(const xmatrix<double>& Uc,const xmatrix<double>& Uf,string& _string,bool& _inversion,double& _angle,
      xvector<double>& _axis,xmatrix<double>& _generator, xvector<double>& _generator_coefficients, 
      xmatrix<xcomplex<double> >& _SU2_matrix, xvector<xcomplex<double> >& _su2_coefficients, double _eps_);  // calculate the symmetry inversion,type,axis,generator //DX20171206 - Added generator coefficients //DX20171207 - Added Uf //DX20180117 - Added SU2 and su2 coefficients
  bool TypePointGroupOperationInternational(const xmatrix<double>& Uc,string& _stringHM,string& _stringSC,
      const bool& _inversion,const double& _angle,
      const xvector<double>& _axis,const xmatrix<double>& _generator, xvector<double>& _generator_coefficients, 
      xmatrix<xcomplex<double> >& _SU2_matrix, xvector<xcomplex<double> >& _su2_coefficients, double _eps_);  // International symbol = Hermann-Mauguin notation & Schonflies notation //DX20171206 - Added generator coefficients //DX20180117 - Added SU2 and su2 coefficients
  //DX+CO START
  uint AddSymmetryToStructure(xstructure &a,const uint& iat,
      const xmatrix<double> &Uc,const xmatrix<double> &Uf,const xvector<double> &ctau,const xvector<double> &ftau,
      const xvector<double> &ctrasl,const xvector<double> &ftrasl,
      const std::vector<int> &basis_atoms_map,const std::vector<int> &basis_types_map,bool basis_map_calculated,char group);
  uint AddSymmetryToStructure(xstructure &a,const uint& iat,
      const xmatrix<double> &Uc,const xmatrix<double> &Uf,const xvector<double> &ctau,const xvector<double> &ftau,
      const xvector<double> &ctrasl,const xvector<double> &ftrasl,
      const std::vector<int> &basis_atoms_map,const std::vector<int> &basis_types_map,bool basis_map_calculated,char group,bool roff); //DX
  uint AddSymmetryToStructure(xstructure& a,const xmatrix<double>& Uc,const xmatrix<double>& Uf,
      const xvector<double>& ctau,const xvector<double>& ftau,const xvector<double>& ctrasl,
      const xvector<double>& ftrasl,
      const std::vector<int>& basis_atoms_map,const std::vector<int>& basis_types_map,bool basis_map_calculated,char group);
  uint AddSymmetryToStructure(xstructure& a,const xmatrix<double>& Uc,const xmatrix<double>& Uf,
      const xvector<double>& ctau,const xvector<double>& ftau,const xvector<double>& ctrasl,
      const xvector<double>& ftrasl,
      const std::vector<int>& basis_atoms_map,const std::vector<int>& basis_types_map,bool basis_map_calculated,char group,bool roff); //DX
  uint AddSymmetryToStructure(xstructure &a,const uint& iat,const xmatrix<double> &Uc,const xmatrix<double> &Uf,
      const std::vector<int> &basis_atoms_map,const std::vector<int> &basis_types_map,bool basis_map_calculated,char group);
  uint AddSymmetryToStructure(xstructure &a,const uint& iat,const xmatrix<double> &Uc,const xmatrix<double> &Uf,
      const std::vector<int> &basis_atoms_map,const std::vector<int> &basis_types_map,bool basis_map_calculated,char group,bool roff); //DX
  uint AddSymmetryToStructure(xstructure &a,const xmatrix<double> &Uc,const xmatrix<double> &Uf,
      const std::vector<int> &basis_atoms_map,const std::vector<int> &basis_types_map,bool basis_map_calculated,char group); //DX
  uint AddSymmetryToStructure(xstructure &a,const xmatrix<double> &Uc,const xmatrix<double> &Uf,
      const std::vector<int> &basis_atoms_map,const std::vector<int> &basis_types_map,bool basis_map_calculated,char group,bool roff); //DX
  bool PointGroupsIdentical(const vector<_sym_op>& vpg1,const vector<_sym_op>& vpg2, double eps, bool is_same_lattice=false); //DX20171207 - added is_same_lattice
  //GG START
  bool CalculateQuaternion(_sym_op& a);
  //GG STOP
  bool ComplexSU2Rotations(xmatrix<xcomplex<double> > & _SU2_matrix, xvector<xcomplex<double> >& _su2_coefficients, double& theta, xvector<double>& _axis); //DX20180117 - add SU(2) and su(2) coefficients
  //DX+CO END
  bool CalculatePointGroup(ofstream& FileMESSAGE,xstructure& a,_aflags& aflags,bool _write_,const bool& osswrite,ostream& oss,string format="txt");      // POINT GROUP      _PGROUP_
  uint CalculatePointGroup(const xmatrix<double>& lattice, vector<_sym_op > pgroup, ofstream &FileMESSAGE,bool _write_,const bool& osswrite,ostream& oss,double _eps_);
  uint CalculatePointGroup(const xmatrix<double>& lattice, vector<_sym_op > pgroup, bool _write_,const bool& osswrite,ostream& oss,double _eps_);     // POINT GROUP      _PGROUP_
  uint CalculatePointGroup(const xmatrix<double>& lattice,double _eps_);     // POINT GROUP      _PGROUP_
  uint CalculatePointGroup(const xmatrix<double>& lattice);     // POINT GROUP      _PGROUP_
  bool CalculatePointGroup(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format="txt");      // POINT GROUP      _PGROUP_
  //DX+CO START
  bool CalculatePointGroup_20160101(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_); //DX
  bool CalculatePointGroup_20160801(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format="txt"); //DX
  //DX+CO END
  bool CalculatePointGroupKLattice(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,string format="txt");  // POINT GROUP KLATTICE     _PGROUPK_
  bool CalculatePointGroupKCrystal(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,string format="txt");  // POINT GROUP KCRYSTAL     _PGROUPK_XTAL_ //DX20171205 - New group: reciprocal space counterpart of pgroup_xtal
  bool TransformSymmetryFromRealToReciprocal(ofstream &FileMESSAGE, xstructure& real_space_crystal, xstructure& reciprocal_space,
      _aflags& aflags, const bool& osswrite, ostream& oss, string& pgroup_type); //DX20170808 - New klattice routine //DX20171205 - Added pgroup_type option to account for pgroupk_xtal
  bool CalculateSitePointGroup(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,string format="txt");  // SITE POINT GROUP _AGROUP_
  //DX+CO START
  bool CalculateSitePointGroup(ofstream &FileMESSAGE,xstructure &a,int CALCULATION_MODE,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,string format="txt"); // SITE POINT GROUP _AGROUP_
  bool CalculateSitePointGroup_20160801(ofstream &FileMESSAGE,xstructure &a,int CALCULATION_MODE,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format="txt"); // SITE POINT GROUP _AGROUP_ //DX
  bool CalculateSitePointGroup_20160101(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_); // SITE POINT GROUP _AGROUP_ //DX
  bool CalculateSitePointGroup_EquivalentSites(xstructure &a,double _eps_); //DX
  bool CalculateSitePointGroup_EquivalentSites(xstructure &a,bool get_full_basis,double _eps_); //DX
  //DX+CO END
  bool CalculatePointGroupCrystal(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,string format="txt");     // POINT GROUP      _PGROUP_
  bool CalculatePointGroupCrystal(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format="txt");      // POINT GROUP      _PGROUP_
  //DX+CO START
  bool CalculatePointGroupCrystal_20170814(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format="txt");      // POINT GROUP      _PGROUP_ //DX
  bool CalculatePointGroupCrystal_20160801(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format="txt");      // POINT GROUP      _PGROUP_ //DX
  bool CalculatePointGroupCrystal_20160101(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_);      // POINT GROUP      _PGROUP_ //DX
  bool CalculatePointGroupKPatterson(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,string format="txt");     // POINT GROUP PATTERSON     _PGROUPK_PATTERSON_ //DX20200129 
  bool CalculatePointGroupKPatterson(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format="txt");      // POINT GROUP PATTERSON     _PGROUPK_PATTERSON_ //DX20200129
  bool PointGroupMap(xstructure& a, string& pgname, string& operations, char group); //DX20170906
  bool PointGroupLookUpTable(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,string format);
  //DX+CO END
  void CalculateSitePointGroup2(xstructure &a,bool ComMidss); // for --agroup2 and --agroup2m
  //DX START
  //xstructure and _sym_op
  bool getFullSymBasis(const xstructure& a, _sym_op& symOp,bool map_types,vector<int>& basis_atoms_map,vector<int>& basis_types_map);
  bool getFullSymBasis(const xstructure& a, _sym_op& symOp,bool map_types,double tolerance,vector<int>& basis_atoms_map,vector<int>& basis_types_map); //CO20190520 - removed pointers for bools and doubles, added const where possible
  bool getFullSymBasis(const xstructure& a, _sym_op& symOp,bool map_types,bool skew,double tolerance,vector<int>& basis_atoms_map,vector<int>& basis_types_map); //CO20190520 - removed pointers for bools and doubles, added const where possible
  bool getFullSymBasis(const xstructure& a, _sym_op& symOp,bool map_types,bool skew,double tolerance,vector<int>& basis_atoms_map,vector<int>& basis_types_map); //CO20190520 - removed pointers for bools and doubles, added const where possible
  //atoms, c2f, f2c and _sym_op
  bool getFullSymBasis(const deque<_atom>& atoms,const xmatrix<double>& lattice,const xmatrix<double>& c2f, const xmatrix<double>& f2c, _sym_op& symOp,bool map_types,bool skew,double tolerance,vector<int>& basis_atoms_map,vector<int>& basis_types_map); //CO20190520 - removed pointers for bools and doubles, added const where possible
  bool getFullSymBasis(const deque<_atom>& atoms,const xmatrix<double>& lattice,const xmatrix<double>& c2f, const xmatrix<double>& f2c, _sym_op& symOp,bool map_types,bool skew,double tolerance,vector<int>& basis_atoms_map,vector<int>& basis_types_map); //CO20190520 - removed pointers for bools and doubles, added const where possible
  bool getFullSymBasis_20170729(const deque<_atom>& atoms,const xmatrix<double>& lattice,const xmatrix<double>& c2f, const xmatrix<double>& f2c, _sym_op& symOp,bool map_types,bool skew,double tolerance,vector<int>& basis_atoms_map,vector<int>& basis_types_map); //CO20190520 - removed pointers for bools and doubles, added const where possible
  //bool getFullSymBasis(deque<_atom>& atoms,xmatrix<double>& Uf, xmatrix<double>& c2f, xmatrix<double>& f2c, bool skew, double tolerance, vector<int>& basis_atoms_map,vector<int>& basis_types_map); //CO20190520 - removed pointers for bools and doubles, added const where possible
  //bool getFullSymBasis(deque<_atom>& atoms,xmatrix<double>& Uf, xmatrix<double>& c2f, xmatrix<double>& f2c, string& str_Hermann_Mauguin, bool skew, double tolerance, vector<int>& basis_atoms_map,vector<int>& basis_types_map); //CO20190520 - removed pointers for bools and doubles, added const where possible
  //bool getFullSymBasis(deque<_atom>& atoms,xmatrix<double>& Uf, xmatrix<double>& c2f, xmatrix<double>& f2c, xvector<double>& ftau, bool skew, double tolerance, vector<int>& basis_atoms_map,vector<int>& basis_types_map); //CO20190520 - removed pointers for bools and doubles, added const where possible
  //bool getFullSymBasis(deque<_atom>& atoms,xmatrix<double>& Uf, xmatrix<double>& c2f, xmatrix<double>& f2c, string& str_Hermann_Mauguin, xvector<double>& ftau, bool skew, double tolerance, vector<int>& basis_atoms_map,vector<int>& basis_types_map); //CO20190520 - removed pointers for bools and doubles, added const where possible
  bool CalculateFactorGroup_20160801(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format="txt");
  //DX END
  bool CalculateFactorGroup(ofstream& FileMESSAGE,xstructure& a,_aflags& aflags,bool _write_,const bool& osswrite,ostream& oss,string format="txt");     // FACTOR GROUP     _FGROUP_
  bool CalculateFactorGroup(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format="txt");      // FACTOR GROUP      _FGROUP_
  //DX START
  bool AtomsMapped(const _atom& a, const _atom& b, const xmatrix<double>& lattice, bool skew, double tol); //DX20190620
  bool AtomsMapped(const _atom& a, const _atom& b, const xmatrix<double>& lattice, const xmatrix<double>& f2c, bool skew, double tol);   //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - lattice and f2c as input
  xvector<double> minimizeDistanceCartesianMethod(const xvector<double>& cpos1, const xvector<double>& cpos2, const xmatrix<double>& lattice); //DX20190613
  xvector<double> minimizeDistanceCartesianMethod(const xvector<double>& cpos1, const xvector<double>& cpos2, const xmatrix<double>& lattice, xvector<int>& ijk); //DX20190613
  xvector<double> minimizeDistanceFractionalMethod(const xvector<double>& fpos1, const xvector<double>& fpos2); //DX20190613
  xvector<double> minimizeDistanceFractionalMethod(const xvector<double>& fdiff); //DX20190613
  xvector<double> minimizeDistanceFractionalMethod(const xvector<double>& fdiff, xvector<int>& ijk); //DX20190613
  //DX20190613 [OBOSLETE] bool minimizeCartesianDistance(const xvector<double>& coord1, const xvector<double>& coord2, xvector<double>& out, const xmatrix<double>& c2f, const xmatrix<double>& f2c, double tol); //CO20190520 - removed pointers for bools and doubles, added const where possible
  //DX20190613 [OBOSLETE] bool minimizeCartesianDistance(const xvector<double>& coord1, const xvector<double>& coord2, xvector<double>& out, const xmatrix<double>& c2f, const xmatrix<double>& f2c, xvector<int>& ijk, bool& restriction, double tol); //CO20190520 - removed pointers for bools and doubles, added const where possible
  //DX20190613 [OBOSLETE] double minimumCartesianDistance(const xvector<double>& coord1, const xvector<double>& coord2, const xmatrix<double>& lattice);
  //DX20190613 [OBOSLETE] double minimumCartesianDistance(const xvector<double>& coord1, const xvector<double>& coord2, const xmatrix<double>& lattice,xvector<double>& min_vec,xvector<int>& ijk);
  //DX20190613 [OBOSLETE] xvector<double> minimumCartesianVector(const xvector<double>&, const xvector<double>&,
  //DX20190613 [OBOSLETE]                                        const xmatrix<double>&);  //ME20180730
  //DX20190613 [OBOSLETE] xvector<double> minimumCartesianVector(const xvector<double>&, const xvector<double>&,
  //DX20190613 [OBOSLETE]                                        const xmatrix<double>&, xvector<int>&);  //ME20180730
  //DX20190613 [OBOSLETE] bool PBC(xvector<double>& v_in, xvector<int>& ijk, bool& restriction);
  //DX20190613 [OBOSLETE] bool PBC(xvector<double>& v_in);
  xvector<double> FPOSDistFromFPOS(const xvector<double>& fpos1,const xvector<double>& fpos2,const xmatrix<double>& lattice,bool skew=false);  //CO20190525
  xvector<double> FPOSDistFromFPOS(const xvector<double>& fpos1,const xvector<double>& fpos2,const xmatrix<double>& lattice,const xmatrix<double>& c2f,const xmatrix<double>& f2c,bool skew=false);  //CO20190525
  xvector<double> CPOSDistFromFPOS(const xvector<double>& fpos1,const xvector<double>& fpos2,const xmatrix<double>& lattice,bool skew=false);  //DX20190620
  xvector<double> CPOSDistFromFPOS(const xvector<double>& fpos1,const xvector<double>& fpos2,const xmatrix<double>& lattice,const xmatrix<double>& f2c,bool skew=false); //DX20190620

  bool FPOSMatch(const deque<_atom>& atom_set, const _atom& atom2, uint& match_type, const xmatrix<double>& lattice, bool skew, double tol); //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - lattice and f2c as input
  bool FPOSMatch(const deque<_atom>& atom_set, const _atom& atom2, uint& match_type, const xmatrix<double>& lattice, const xmatrix<double>& f2c, bool skew, double tol); //DX20190620 - overload
  bool FPOSMatch(const _atom& atom1, const _atom& atom2, const xmatrix<double>& lattice, bool skew, double tol); //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - lattice and f2c as input
  bool FPOSMatch(const _atom& atom1, const _atom& atom2, const xmatrix<double>& lattice, const xmatrix<double>& f2c, bool skew, double tol); //DX20190620 - overload 
  bool FPOSMatch(const xvector<double>& atom1, const xvector<double>& atom2, const xmatrix<double>& lattice, bool skew, double tol); //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - lattice and f2c as input
  bool FPOSMatch(const xvector<double>& atom1, const xvector<double>& atom2, const xmatrix<double>& lattice, const xmatrix<double>& f2c, bool skew, double tol); //DX20190620 - overload 
  bool validateAtomPosition(const _atom& atom,const xmatrix<double>& c2f,const xmatrix<double>& f2c,bool skew,double& _eps_); //CO20190520 - removed pointers for bools and doubles, added const where possible
  bool validateAtomPosition(const xvector<double>& cpos,const xvector<double>& fpos,const xmatrix<double>& c2f,const xmatrix<double>& f2c,bool skew,double& _eps_); //CO20190520 - removed pointers for bools and doubles, added const where possible
  bool MapAtom(const deque<_atom>& a_vec, const _atom& b, bool map_types, const xmatrix<double>& lattice, bool skew, double tol); //DX20190620
  bool MapAtom(const deque<_atom>& a_deq, const _atom& b, bool map_types, const xmatrix<double>& lattice, const xmatrix<double>& f2c, bool skew, double tol); //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - lattice and f2c as input
  bool MapAtom(const vector<_atom> a_vec, const _atom b, bool map_types, const xmatrix<double>& lattice, bool skew, double tol); //DX20190620
  bool MapAtom(const vector<_atom> a_vec, const _atom b, bool map_types, const xmatrix<double>& lattice, const xmatrix<double>& f2c, bool skew, double tol); //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - lattice and f2c as input
  bool MapAtom(const xvector<double>& a, const xvector<double>& b, const xmatrix<double>& lattice, bool skew, double tol); //DX20190620
  bool MapAtom(const xvector<double>& a, const xvector<double>& b, const xmatrix<double>& lattice, const xmatrix<double>& f2c, bool skew, double tol); //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - lattice and f2c as input
  bool MapAtom(const _atom& a, const _atom& b, bool map_types, const xmatrix<double>& lattice, bool skew, double tol); //DX20190620 
  bool MapAtom(const _atom& a, const _atom& b, bool map_types, const xmatrix<double>& lattice, const xmatrix<double>& f2c, bool skew, double tol); //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - lattice and f2c as input, remove "Atom" prefix from name
  bool BringInCellTolABC(xstructure& a, xvector<double> tol_abc_res);
  //[CO20190515 - not needed and is ambiguous with overload]bool MapAtomWithBasis(vector<_atom>& vec, _atom& a, bool map_types, deque<uint>& index_to_check, xmatrix<double>& c2f, xmatrix<double>& f2c, bool skew, double tol,bool fast=true);    //CO20190520 - removed pointers for bools and doubles, added const where possible
  bool MapAtomWithBasis(const vector<_atom>& vec, const _atom& a, bool map_types, deque<uint>& index_to_check, const xmatrix<double>& lattice, bool skew, double tol, uint& mapped_index,bool fast=true); //DX20190620
  bool MapAtomWithBasis(const vector<_atom>& vec, const _atom& a, bool map_types, deque<uint>& index_to_check, const xmatrix<double>& lattice, const xmatrix<double>& f2c, bool skew, double tol, uint& mapped_index,bool fast=true);   //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - lattice and f2c as input
  //[CO20190515 - not needed and is ambiguous with overload]bool MapAtomWithBasis(deque<_atom>& vec, _atom& a, bool map_types, deque<uint>& index_to_check, xmatrix<double>& c2f, xmatrix<double>& f2c, bool skew, double tol,bool fast=true); //CO20190520 - removed pointers for bools and doubles, added const where possible
  bool MapAtomWithBasis(const deque<_atom>& vec, const _atom& a, bool map_types, deque<uint>& index_to_check, const xmatrix<double>& lattice, bool skew, double tol, uint& mapped_index,bool fast=true); //DX20190620
  bool MapAtomWithBasis(const deque<_atom>& vec, const _atom& a, bool map_types, deque<uint>& index_to_check, const xmatrix<double>& lattice, const xmatrix<double>& f2c, bool skew, double tol, uint& mapped_index,bool fast=true); //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190619 - lattice and f2c as input, remove "Atom" prefix from name
  bool isLatticeSkewed(const xmatrix<double>& lattice, double& min_dist, double tol); //CO20190520 - removed pointers for bools and doubles, added const where possible
  double minimumDistance(const xstructure& xstr);
  double minimumDistance(const deque<_atom>& atoms);  //CO20190808 - for NON periodic systems
  double minimumDistance(const deque<_atom>& atoms, const xmatrix<double>& lattice,double scale=1.0);
  double defaultTolerance(const xstructure& xstr);
  bool checkAngle(xvector<double>& v1, xvector<double>& v2, double input_angle, double tolerance); //CO20190520 - removed pointers for bools and doubles, added const where possible
  bool checkAngle(xvector<double>& v1, xvector<double>& v2, double input_angle, bool& is_deg, double tolerance); //CO20190520 - removed pointers for bools and doubles, added const where possible
  bool checkAngle(double& mod_v1, double& mod_v2, double angle1, double angle2, double tolerance); //CO20190520 - removed pointers for bools and doubles, added const where possible
  bool checkAngle(double& mod_v1, double& mod_v2, double angle1, double angle2, bool& is_deg, double tolerance); //CO20190520 - removed pointers for bools and doubles, added const where possible
  //DX20170905 [OBSOLETE] bool change_tolerance(xstructure& xstr, double tolerance, double& orig_tolerance, int& count , double& min_dist, bool& no_scan); //CO20190520 - removed pointers for bools and doubles, added const where possible
  bool change_tolerance(xstructure& xstr, double& tolerance, double& min_dist, bool& no_scan); //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190524 - need pointer for tolerance, otherwise it will not update
  deque<deque<_atom> > break_up_by_type(deque<_atom>& expanded_crystal);
  vector<vector<_atom> > break_up_by_type(vector<_atom> expanded_crystal);
  //DX20190905 [OBSOLETE] double mod_one(double d); //DX 
  //DX20190905 [OBSOLETE] _atom mod_one_atom(const _atom& atom_in); //CO
  //DX20190905 [OBSOLETE] xvector<double> mod_one_xvec(xvector<double> a); //DX
  bool CheckForIdentity(const xstructure& xstr); //DX
  bool checkSuperCellLatticePoints(xstructure& xstr, int& num_lattice_points, char& centering, uint& expand_size); //DX
  bool ComparePointGroupAndSpaceGroupString(xstructure& xstr, int& multiplicity_of_primitive, bool& derivative_structure); //DX
  bool CalculateFactorGroup_20160101(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_);      // FACTOR GROUP      _FGROUP_
  bool CalculateSpaceGroup_20160101(ofstream& FileMESSAGE,xstructure& a,_aflags& aflags,bool _write_,const bool& osswrite,ostream& oss);      // SPACE GROUP      _SGROUP_
  bool CalculateSpaceGroup_20160801(ofstream& FileMESSAGE,xstructure& a,_aflags& aflags,bool _write_,const bool& osswrite,ostream& oss,string format="txt");      // SPACE GROUP      _SGROUP_
  //DX END
  bool CalculateSpaceGroup(ofstream& FileMESSAGE,xstructure& a,_aflags& aflags,bool _write_,const bool& osswrite,ostream& oss,string format="txt");      // SPACE GROUP      _SGROUP_

  bool CalculateInequivalentAtoms(xstructure&);
  bool CalculateInequivalentAtoms(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,string format="txt"); // EQUIVALENT ATOMS _IATOMS_
  //DX+CO START
  bool CalculateInequivalentAtoms(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format="txt"); // EQUIVALENT ATOMS _IATOMS_ //DX
  bool CalculateInequivalentAtoms(ofstream &FileMESSAGE,xstructure &a,bool rely_on_basis,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,string format="txt"); // EQUIVALENT ATOMS _IATOMS_ //DX
  bool CalculateInequivalentAtoms(ofstream &FileMESSAGE,xstructure &a,bool rely_on_basis,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format="txt"); // EQUIVALENT ATOMS _IATOMS_ //DX
  bool CalculateInequivalentAtoms_20160801(ofstream &FileMESSAGE,xstructure &a,bool rely_on_basis,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss,double _eps_,string format="txt"); // EQUIVALENT ATOMS _IATOMS_ //DX  
  bool CalculateInequivalentAtoms_20160101(ofstream &FileMESSAGE,xstructure &a,_aflags &aflags,bool _write_,const bool& osswrite,ostream& oss); // EQUIVALENT ATOMS _IATOMS_ //DX
  //DX+CO END
}
string AgroupSymmetryToJson(vector<vector<_sym_op> >& group, char& mode); //DX20170803 - For Python wrapper
string EquivalentAtomsToJson(vector<vector<int> >& iatoms); //DX20170803 - For Python wrapper
string SymmetryToJson(vector<_sym_op>& group, char& mode); //DX20170803 - For Python wrapper
bool KBIN_SymmetryWrite(ofstream& FileMESSAGE,xstructure& a,_aflags& aflags,char group,const bool& osswrite,ostream& oss,const string& format="txt");
//bool KBIN_SymmetryToScreen(xstructure& a, string& format, ostream& oss); //DX20170803 - For Python wrapper
bool KBIN_SymmetryToScreen(xstructure& a, string& format, ostream& oss, char mode='\0'); //DX20170822 - For Python wrapper
bool KBIN_StepSymmetryPerform(xstructure& a,string AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,const bool& osswrite,ostream& oss);
//DX+CO START
bool KBIN_StepSymmetryPerform_20161205(xstructure& a,string AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,const bool& osswrite,ostream& oss); //DX
bool KBIN_StepSymmetryPerform_20160101(xstructure& a,string AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,const bool& osswrite,ostream& oss); //DX
//DX+CO END
vector<double> PointGroupHistogramCheck(xstructure& a);
//bool SYM_CalculatePointGroup(ofstream& FileMESSAGE,xstructure& a,_aflags& aflags,bool _write_,const bool& osswrite,ostream& oss);      // MIKNOWSKI BASIS REDUCTION

// ----------------------------------------------------------------------------
// aflow_spacegroup.cpp
// spagegroup tables and routines
extern string LibrarySPACEGROUP;
namespace spacegroup{
  class _spacegroup {
    public:
      // constructors/destructors
      _spacegroup(void);    // do nothing
      _spacegroup(const _spacegroup& b);    // do nothing
      ~_spacegroup();        // do nothing
      // OPERATORS                                                  // --------------------------------------
      const _spacegroup& operator=(const _spacegroup& b);             // some operators
      // CONTENT
      uint number;
      uint option;
      string stroption;
      string name;
      string sginfo;
      std::vector<_sym_op> fgroup;                                  // rotations/inversions + incell_translations operations
      std::vector<_sym_op> pgroup;                                  // rotations/inversions
    private:                                                       // ---------------------------------------
      void free();                                                  // to free everything
      void copy(const _spacegroup& b);                               // the flag is necessary because sometimes you need to allocate the space.
  };

  extern std::vector<_spacegroup> vspacegroups;
  uint SpaceGroupInitialize(void);
  bool SpaceGroupNumberStructure(xstructure &str);
  bool SpaceGroupOptionRequired(uint spacegroup);
}

// ----------------------------------------------------------------------------
// aflow_lattice.cpp
// lattice and brillouin zones
namespace LATTICE {
  bool lattice_is_working(string lat);
  string Lattice2TypeAndCentering(const string& lattice_type); //DX20191031
  string SpaceGroup2Lattice(uint sg);
  string SpaceGroup2LatticeTypeAndCentering(uint sg); //DX20191031
  uint Lattice2SpaceGroup(string lattice,vector<uint>& vsg);
  string SpaceGroup2LatticeVariation(uint sg,const xstructure& str);
  string ConventionalLattice_SpaceGroup(uint sg,double a,double b,double c);
  string ConventionalLattice_SpaceGroup(uint sg,const xstructure& str);
  xvector<double> Getabc_angles_Conventional(const xmatrix<double>& rlattice, string lattice,const int& mode);
  bool fix_sts_sp(xstructure& str_sp,xmatrix<double> &rlattice,xmatrix<double> &plattice);
  bool Standard_Lattice_Structure(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc,bool full_sym=true);
  bool Standard_Lattice_StructureDefault(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc,bool full_sym=true);
  bool Standard_Lattice_StructureCoarse(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc);
  bool Standard_Lattice_StructureNormal(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc);
  bool Standard_Lattice_StructureMedium(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc);
  bool Standard_Lattice_StructurePrecise(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc);
  bool Standard_Lattice_StructureUltra(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc);
  //bool Standard_Lattice_Structure(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc,double eps,double epsang); //SC OLD VERSION
  bool Standard_Lattice_Structure(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc,double eps,double epsang,int& time,double symeps);
  bool Standard_Lattice_Structure(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc,double eps,double epsang,int& time,double symeps,bool histogram);
  bool Standard_Lattice_Structure_20170101(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc,double eps,double epsang,int& time,double symeps,bool histogram);
  bool Standard_Lattice_Structure_20170718(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc,bool full_sym=true);
  //bool Standard_Lattice_Structure(const xstructure& str_in,xstructure& str_sp,xstructure& str_sc,double eps,double epsang,int& time,int mode);
  bool Bravais_Lattice_Structure(xstructure& str_in,xstructure& str_sp,xstructure& str_sc,double eps,double epsang); // calculate everything
  bool Bravais_Lattice_StructureDefault(xstructure& str_in,xstructure& str_sp,xstructure& str_sc,bool full_sym=true); // calculate everything
  //DX START
  bool Bravais_Lattice_StructureDefault_20170401(xstructure& str_in,xstructure& str_sp,xstructure& str_sc,bool full_sym=true); // calculate everything
  bool Bravais_Lattice_StructureDefault_20160101(xstructure& str_in,xstructure& str_sp,xstructure& str_sc); // calculate everything
  //DX END
  bool Lattice(const xmatrix<double>& lattice,xmatrix<double>& lattice_sp,xmatrix<double>& lattice_sc,string& bravais_lattice_type,string& bravais_lattice_variation_type,string& bravais_lattice_system,double eps,double epsang);
  string Bravais_Lattice_Type(const xmatrix<double>& lattice,xmatrix<double>& lattice_sp,xmatrix<double>& lattice_sc,double eps,double epsang);
  string Bravais_Lattice_Type(const xmatrix<double>& lattice,xmatrix<double>& lattice_sp,xmatrix<double>& lattice_sc);
  string Bravais_Lattice_Type(const xmatrix<double>& lattice,double eps,double epsang);
  string Bravais_Lattice_Type(const xmatrix<double>& lattice);
  string Bravais_Lattice_Variation_Type(const xmatrix<double>& lattice,xmatrix<double>& lattice_sp,xmatrix<double>& lattice_sc,double eps,double epsang);
  string Bravais_Lattice_Variation_Type(const xmatrix<double>& lattice,xmatrix<double>& lattice_sp,xmatrix<double>& lattice_sc);
  string Bravais_Lattice_Variation_Type(const xmatrix<double>& lattice,double eps,double epsang);
  string Bravais_Lattice_Variation_Type(const xmatrix<double>& lattice);
  string Bravais_Lattice_System(const xmatrix<double>& lattice,xmatrix<double>& lattice_sp,xmatrix<double>& lattice_sc,double eps,double epsang);
  string Bravais_Lattice_System(const xmatrix<double>& lattice,xmatrix<double>& lattice_sp,xmatrix<double>& lattice_sc);
  string Bravais_Lattice_System(const xmatrix<double>& lattice,double eps,double epsang);
  string Bravais_Lattice_System(const xmatrix<double>& lattice);
  string Primitive_Lattice_Type(const xstructure& str);
  string Bravais_Lattice_System(const xstructure& str);
  string Conventional_Lattice_Type(const xstructure& str);
  xstructure Standard_Primitive_Lattice_Structure(const xstructure& str);
  xstructure Standard_Conventional_Lattice_Structure(const xstructure& str);
  string Get_Primitive_Lattice_Structure(const xstructure& str);
  xmatrix<double> sc2sp(const xmatrix<double>& rlattice, string lattice,bool inverseflag);
  xmatrix<double> sp2sc(const xmatrix<double>& rlattice, string lattice,bool inverseflag);
  void BZPLOTDATA(string options,istream& poscar, int mode);
}
xvector<double> Vrotate(xvector<double> v, xvector<double> axisrot, double theta);
void CheckLatticeHistogram();
namespace LATTICE {
  // kpoints and brillouin zones
  string KPOINTS_Directions(xstructure str_in,double grid,bool &foundBZ);
  string KPOINTS_Directions(string lattice_type,xmatrix<double> sp, double _grid,int iomode,bool &foundBZ);
  string KPOINTS_Directions(string lattice_type,xmatrix<double> sp, xmatrix<double> transformation_matrix, double _grid,int iomode,bool &foundBZ); //DX20181101
}

// ----------------------------------------------------------------------------
// neighbours prototypes
// aflow_neighbours.cpp
bool StepNeighboursPerform(xstructure& a,string AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags);

// ----------------------------------------------------------------------------
// surface prototypes
// aflow_surface.cpp
namespace surface {
  double PointInTriangleContribution(const xvector<double>& _point,const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3);
  double PointInRhombusContribution(const xvector<double>& _point,const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3,const xvector<double>& v4);
  double TriangleArea(const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3);
  bool PlaneGetABCD(double& a,double& b,double& c,double& d,const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3);
  double PlaneDistance(const xvector<double>& r,const double& a,const double& b,const double& c,const double& d);
  double PlaneDistance(const xvector<double>& r,const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3);
  xvector<double> PlaneGetProjection(const xvector<double>& r,const double& a,const double& b,const double& c,const double& d);
  xvector<double> PlaneGetProjection(const xvector<double>& r,const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3);
  xvector<double> PlaneGetHKL(const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3,const xvector<double>& a1,const xvector<double>& a2,const xvector<double>& a3);
  bool PlaneGetVVV(const xvector<double>& hkl,double& area,xvector<double>& v1,xvector<double>& v2,xvector<double>& v3,xvector<double>& v4,const xvector<double>& a1,const xvector<double>& a2,const xvector<double>& a3);
  double GetPlaneDensityAtoms(const xstructure& _str,const xvector<double>& hkl,const double& roughness,const int& type_at);
  double GetPlaneDensityAtoms(const xstructure& _str,const xvector<double>& hkl,const double& roughness);
  double GetPlaneDensityBBonds(const xstructure& _str,const xvector<double>& hkl,const double& roughness,const double& bbdistance,const int& type_at1,const int& type_at2);
  double GetPlaneDensityBBonds(const xstructure& _str,const xvector<double>& hkl,const double& roughness,const double& bbdistance);
  double GetNNeighbours(const xstructure& _str,const int& type_at1,const int& type_at2);
  double GetNNeighbours(const xstructure& _str);
  string PrintHKLSigma(int num_types,int num_types_combinations);
  string PrintHKLSigmaBB(int num_types,int num_types_combinations,const double& bbfrac,const double& bbdistance,const xmatrix<double>& bbdistances);
  bool GetSurfaceHKL(const xstructure& _str,_aflags& aflags,const xvector<double>& _hkl,vector<vector<double> >& planesreducible,vector<vector<double> >& planesirreducible,ostream& oss);
  bool GetSurfaceHKLSearch(const xstructure& _str,_aflags& aflags,const xvector<double>& iparams,vector<vector<double> >& planesreducible,vector<vector<double> >& planesirreducible,vector<vector<uint> >& planesirreducible_images,ostream& oss,const string& smode);
}

namespace slab { // ROMAN CHEPULSKYY
  xstructure MAKE_SLAB(string options, istream& cin);
  xstructure MAKE_SLAB(string options, xstructure& str_in);
  void POSCAR_reading(istream& cin);
  double VectorAbsValue(int Layer, int NinLayer,const xmatrix<double>& UnitCellVector,const vector<vector<vector<int> > >& LayerSitesDirCoords);
  double VectorScalarMult(int Layer1, int NinLayer1, int Layer2, int NinLayer2,const xmatrix<double>& UnitCellVector,const vector<vector<vector<int> > >& LayerSitesDirCoords);
  double CosAngle(int Layer1, int NinLayer1, int Layer2, int NinLayer2,const xmatrix<double>& UnitCellVector,const vector<vector<vector<int> > >& LayerSitesDirCoords);
  double hkl_CartCoord_Length(xvector<double>& hkl_CartCoord,const xvector<double>& hkl);
} // namespace slab

#define DEFAULT_TOTAL_LAYERS 10 //CO20190601 - move to .aflow.rc eventually
#define DEFAULT_V3_ANGLE_DEVIATION 5.0  //CO20190803 - move to .aflow.rc eventually
namespace slab { //CO20190601
  xvector<double> HKLPlane2Normal(const xstructure& xstr_in,int h,int k,int l);  //CO20190321
  xvector<double> HKLPlane2Normal(const xmatrix<double>& lattice,int h,int k,int l); //CO20190321
  xvector<double> HKLPlane2Normal(const xstructure& xstr_in,const xvector<int>& hkl);  //CO20190321
  xvector<double> HKLPlane2Normal(const xmatrix<double>& lattice,const xvector<int>& hkl); //CO20190321
  bool Normal2HKLPlane(const xstructure& xstr_in,const xvector<double>& n,xvector<int>& hkl);  //CO20190321
  bool Normal2HKLPlane(const xmatrix<double>& lattice,const xvector<double>& n,xvector<int>& hkl); //CO20190321
  vector<xvector<double> > getHKLPlaneIntercepts(const xstructure& xstr_in,int h,int k,int l); //CO20190321
  vector<xvector<double> > getHKLPlaneIntercepts(const xmatrix<double>& lattice,int h,int k,int l);  //CO20190321
  vector<xvector<double> > getHKLPlaneIntercepts(const xstructure& xstr_in,const xvector<int>& hkl); //CO20190321
  vector<xvector<double> > getHKLPlaneIntercepts(const xmatrix<double>& lattice,const xvector<int>& hkl);  //CO20190321
  double getSpacingHKLPlane(const xstructure& xstr_in,int h,int k,int l);  //CO20190321
  double getSpacingHKLPlane(const xmatrix<double>& lattice,int h,int k,int l); //CO20190321
  double getSpacingHKLPlane(const xstructure& xstr_in,const xvector<int>& hkl);  //CO20190321
  double getSpacingHKLPlane(const xmatrix<double>& lattice,const xvector<int>& hkl); //CO20190321
  double getAngleHKLPlanes(const xstructure& xstr_in,int h1,int k1,int l1,int h2,int k2,int l2);  //CO20190321
  double getAngleHKLPlanes(const xmatrix<double>& lattice,int h1,int k1,int l1,int h2,int k2,int l2); //CO20190321
  double getAngleHKLPlanes(const xstructure& xstr_in,const xvector<int>& hkl1,const xvector<int>& hkl2);  //CO20190321
  double getAngleHKLPlanes(const xmatrix<double>& lattice,const xvector<int>& hkl1,const xvector<int>& hkl2); //CO20190321
  void BringInBoundary(xvector<double>& vec,double padding=0.0);
  //[CO20190520 - plugged into BringInBoundary() with no padding]void Bring2OppositeBoundary(xvector<double>& vec);
  xvector<double> getNextAtomInPath(const xstructure& xstr_in,const xvector<double>& l_cpos,const xvector<double>& cpos_starting,vector<uint>& atoms2skip,uint& loop_iteration,bool outside_current_cell=false); //CO20190321
  double getDistanceBetweenImages(const xstructure& xstr_in,const xvector<double>& n,bool outside_cell=false);  //CO20190321
  bool distanceBetweenImages_HKL(const xstructure& xstr_in,const xvector<double>& n,double& distance_between_images,bool outside_cell=false);  //CO20190321
  bool distanceBetweenImages_Tracing(const xstructure& xstr_in,const xvector<double>& n,double& distance_between_images,bool outside_cell=false);  //CO20190321

  //[CO20190808 - OBSOLETE: does not work, need to redefine lattice vectors, use CreateSlab_SurfaceLattice() instead]  //easy inputs
  //[CO20190808 - OBSOLETE: does not work, need to redefine lattice vectors, use CreateSlab_SurfaceLattice() instead]  xstructure CreateSlab_RigidRotation(const xstructure& xstr_in,const xvector<int>& hkl,int total_layers,double vacuum,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,ostream& oss=cout); //CO20190321
  //[CO20190808 - OBSOLETE: does not work, need to redefine lattice vectors, use CreateSlab_SurfaceLattice() instead]  xstructure CreateSlab_RigidRotation(const xstructure& xstr_in,const xvector<int>& hkl,int total_layers,double vacuum,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,const _aflags& aflags,ostream& oss=cout); //CO20190321
  //[CO20190808 - OBSOLETE: does not work, need to redefine lattice vectors, use CreateSlab_SurfaceLattice() instead]  xstructure CreateSlab_RigidRotation(const xstructure& xstr_in,const xvector<int>& hkl,int total_layers,double vacuum,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,ofstream& FileMESSAGE,ostream& oss=cout); //CO20190321
  //[CO20190808 - OBSOLETE: does not work, need to redefine lattice vectors, use CreateSlab_SurfaceLattice() instead]  xstructure CreateSlab_RigidRotation(const xstructure& xstr_in,const xvector<int>& hkl,int total_layers,double vacuum,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss=cout); //CO20190321
  //[CO20190808 - OBSOLETE: does not work, need to redefine lattice vectors, use CreateSlab_SurfaceLattice() instead]  xstructure CreateSlab_RigidRotation(const xstructure& xstr_in,const xvector<int>& hkl,int total_layers,double vacuum,ostream& oss=cout); //CO20190321
  //[CO20190808 - OBSOLETE: does not work, need to redefine lattice vectors, use CreateSlab_SurfaceLattice() instead]  xstructure CreateSlab_RigidRotation(const xstructure& xstr_in,const xvector<int>& hkl,int total_layers,double vacuum,const _aflags& aflags,ostream& oss=cout); //CO20190321
  //[CO20190808 - OBSOLETE: does not work, need to redefine lattice vectors, use CreateSlab_SurfaceLattice() instead]  xstructure CreateSlab_RigidRotation(const xstructure& xstr_in,const xvector<int>& hkl,int total_layers,double vacuum,ofstream& FileMESSAGE,ostream& oss=cout); //CO20190321
  //[CO20190808 - OBSOLETE: does not work, need to redefine lattice vectors, use CreateSlab_SurfaceLattice() instead]  xstructure CreateSlab_RigidRotation(const xstructure& xstr_in,const xvector<int>& hkl,int total_layers,double vacuum,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss=cout); //CO20190321
  //[CO20190808 - OBSOLETE: does not work, need to redefine lattice vectors, use CreateSlab_SurfaceLattice() instead]  //load from xoptions
  //[CO20190808 - OBSOLETE: does not work, need to redefine lattice vectors, use CreateSlab_SurfaceLattice() instead]  xstructure CreateSlab_RigidRotation(const aurostd::xoption& vpflow,istream& input,xvector<int>& hkl,int& total_layers,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,ostream& oss=cout); //CO20190321
  //[CO20190808 - OBSOLETE: does not work, need to redefine lattice vectors, use CreateSlab_SurfaceLattice() instead]  xstructure CreateSlab_RigidRotation(const aurostd::xoption& vpflow,const xstructure& xstr_in,xvector<int>& hkl,int& total_layers,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,ostream& oss=cout);  //CO20190321
  //[CO20190808 - OBSOLETE: does not work, need to redefine lattice vectors, use CreateSlab_SurfaceLattice() instead]  xstructure CreateSlab_RigidRotation(const aurostd::xoption& vpflow,istream& input,xvector<int>& hkl,int& total_layers,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,const _aflags& aflags,ostream& oss=cout); //CO20190321
  //[CO20190808 - OBSOLETE: does not work, need to redefine lattice vectors, use CreateSlab_SurfaceLattice() instead]  xstructure CreateSlab_RigidRotation(const aurostd::xoption& vpflow,const xstructure& xstr_in,xvector<int>& hkl,int& total_layers,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,const _aflags& aflags,ostream& oss=cout);  //CO20190321
  //[CO20190808 - OBSOLETE: does not work, need to redefine lattice vectors, use CreateSlab_SurfaceLattice() instead]  xstructure CreateSlab_RigidRotation(const aurostd::xoption& vpflow,istream& input,xvector<int>& hkl,int& total_layers,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,ofstream& FileMESSAGE,ostream& oss=cout); //CO20190321  
  //[CO20190808 - OBSOLETE: does not work, need to redefine lattice vectors, use CreateSlab_SurfaceLattice() instead]  xstructure CreateSlab_RigidRotation(const aurostd::xoption& vpflow,const xstructure& xstr_in,xvector<int>& hkl,int& total_layers,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,ofstream& FileMESSAGE,ostream& oss=cout);  //CO20190321
  //[CO20190808 - OBSOLETE: does not work, need to redefine lattice vectors, use CreateSlab_SurfaceLattice() instead]  xstructure CreateSlab_RigidRotation(const aurostd::xoption& vpflow,istream& input,xvector<int>& hkl,int& total_layers,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss=cout); //CO20190321
  //[CO20190808 - OBSOLETE: does not work, need to redefine lattice vectors, use CreateSlab_SurfaceLattice() instead]  xstructure CreateSlab_RigidRotation(const aurostd::xoption& vpflow,const xstructure& xstr_in,xvector<int>& hkl,int& total_layers,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss=cout); //CO20190321
  //[CO20190808 - OBSOLETE: does not work, need to redefine lattice vectors, use CreateSlab_SurfaceLattice() instead]  //input directly
  //[CO20190808 - OBSOLETE: does not work, need to redefine lattice vectors, use CreateSlab_SurfaceLattice() instead]  xstructure CreateSlab_RigidRotation(istream& input,const xvector<int>& hkl,int total_layers,double vacuum,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,ofstream& FileMESSAGE,ostream& oss=cout); //CO20190321  
  //[CO20190808 - OBSOLETE: does not work, need to redefine lattice vectors, use CreateSlab_SurfaceLattice() instead]  xstructure CreateSlab_RigidRotation(const xstructure& xstr_in,const xvector<int>& hkl,int total_layers,double vacuum,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,ofstream& FileMESSAGE,ostream& oss=cout);  //CO20190321
  //[CO20190808 - OBSOLETE: does not work, need to redefine lattice vectors, use CreateSlab_SurfaceLattice() instead]  xstructure CreateSlab_RigidRotation(istream& input,const xvector<int>& hkl,int total_layers,double vacuum,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss=cout); //CO20190321
  //[CO20190808 - OBSOLETE: does not work, need to redefine lattice vectors, use CreateSlab_SurfaceLattice() instead]  xstructure CreateSlab_RigidRotation(const xstructure& xstr_in,const xvector<int>& hkl,int total_layers,double vacuum,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss=cout); //CO20190321
  //[CO20190808 - OBSOLETE: does not work, need to redefine lattice vectors, use CreateSlab_SurfaceLattice() instead]  xstructure CreateSlab_RigidRotation(istream& input,const xvector<int>& hkl,int total_layers,double vacuum,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,ostream& oss=cout); //CO20190321
  //[CO20190808 - OBSOLETE: does not work, need to redefine lattice vectors, use CreateSlab_SurfaceLattice() instead]  xstructure CreateSlab_RigidRotation(const xstructure& xstr_in,const xvector<int>& hkl,int total_layers,double vacuum,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,ostream& oss=cout);  //CO20190321
  //[CO20190808 - OBSOLETE: does not work, need to redefine lattice vectors, use CreateSlab_SurfaceLattice() instead]  xstructure CreateSlab_RigidRotation(istream& input,const xvector<int>& hkl,int total_layers,double vacuum,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,const _aflags& aflags,ostream& oss=cout); //CO20190321
  //[CO20190808 - OBSOLETE: does not work, need to redefine lattice vectors, use CreateSlab_SurfaceLattice() instead]  xstructure CreateSlab_RigidRotation(const xstructure& xstr_in,const xvector<int>& hkl,int total_layers,double vacuum,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,const _aflags& aflags,ostream& oss=cout);  //CO20190321

  xmatrix<double> getSlabLattice(istream& input,const xvector<int>& hkl,xmatrix<double>& lattice_slab_origbasis,double ang_dev=DEFAULT_V3_ANGLE_DEVIATION,double vlen_max_strict=AUROSTD_MAX_DOUBLE); //CO20190321
  xmatrix<double> getSlabLattice(const xstructure& _xstr_in,const xvector<int>& hkl,xmatrix<double>& lattice_slab_origbasis,double ang_dev=DEFAULT_V3_ANGLE_DEVIATION,double vlen_max_strict=AUROSTD_MAX_DOUBLE); //CO20190321
  //easy inputs
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in,const xvector<int>& hkl,int total_layers,double vacuum,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,double v3len_max_strict=AUROSTD_MAX_DOUBLE,ostream& oss=cout); //CO20190321
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in,const xvector<int>& hkl,int total_layers,double vacuum,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,const _aflags& aflags,double v3len_max_strict=AUROSTD_MAX_DOUBLE,ostream& oss=cout); //CO20190321
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in,const xvector<int>& hkl,int total_layers,double vacuum,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,ofstream& FileMESSAGE,double v3len_max_strict=AUROSTD_MAX_DOUBLE,ostream& oss=cout); //CO20190321
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in,const xvector<int>& hkl,int total_layers,double vacuum,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,const _aflags& aflags,ofstream& FileMESSAGE,double v3len_max_strict=AUROSTD_MAX_DOUBLE,ostream& oss=cout); //CO20190321
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in,const xvector<int>& hkl,int total_layers,double vacuum,double v3len_max_strict=AUROSTD_MAX_DOUBLE,ostream& oss=cout); //CO20190321
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in,const xvector<int>& hkl,int total_layers,double vacuum,const _aflags& aflags,double v3len_max_strict=AUROSTD_MAX_DOUBLE,ostream& oss=cout); //CO20190321
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in,const xvector<int>& hkl,int total_layers,double vacuum,ofstream& FileMESSAGE,double v3len_max_strict=AUROSTD_MAX_DOUBLE,ostream& oss=cout); //CO20190321
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in,const xvector<int>& hkl,int total_layers,double vacuum,const _aflags& aflags,ofstream& FileMESSAGE,double v3len_max_strict=AUROSTD_MAX_DOUBLE,ostream& oss=cout); //CO20190321
  //load from xoptions
  xstructure CreateSlab_SurfaceLattice(const aurostd::xoption& vpflow,istream& input,xvector<int>& hkl,int& total_layers,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,double v3len_max_strict=AUROSTD_MAX_DOUBLE,ostream& oss=cout); //CO20190321
  xstructure CreateSlab_SurfaceLattice(const aurostd::xoption& vpflow,const xstructure& xstr_in,xvector<int>& hkl,int& total_layers,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,double v3len_max_strict=AUROSTD_MAX_DOUBLE,ostream& oss=cout);  //CO20190321
  xstructure CreateSlab_SurfaceLattice(const aurostd::xoption& vpflow,istream& input,xvector<int>& hkl,int& total_layers,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,const _aflags& aflags,double v3len_max_strict=AUROSTD_MAX_DOUBLE,ostream& oss=cout); //CO20190321
  xstructure CreateSlab_SurfaceLattice(const aurostd::xoption& vpflow,const xstructure& xstr_in,xvector<int>& hkl,int& total_layers,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,const _aflags& aflags,double v3len_max_strict=AUROSTD_MAX_DOUBLE,ostream& oss=cout);  //CO20190321
  xstructure CreateSlab_SurfaceLattice(const aurostd::xoption& vpflow,istream& input,xvector<int>& hkl,int& total_layers,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,ofstream& FileMESSAGE,double v3len_max_strict=AUROSTD_MAX_DOUBLE,ostream& oss=cout); //CO20190321  
  xstructure CreateSlab_SurfaceLattice(const aurostd::xoption& vpflow,const xstructure& xstr_in,xvector<int>& hkl,int& total_layers,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,ofstream& FileMESSAGE,double v3len_max_strict=AUROSTD_MAX_DOUBLE,ostream& oss=cout);  //CO20190321
  xstructure CreateSlab_SurfaceLattice(const aurostd::xoption& vpflow,istream& input,xvector<int>& hkl,int& total_layers,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,const _aflags& aflags,ofstream& FileMESSAGE,double v3len_max_strict=AUROSTD_MAX_DOUBLE,ostream& oss=cout); //CO20190321
  xstructure CreateSlab_SurfaceLattice(const aurostd::xoption& vpflow,const xstructure& xstr_in,xvector<int>& hkl,int& total_layers,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,const _aflags& aflags,ofstream& FileMESSAGE,double v3len_max_strict=AUROSTD_MAX_DOUBLE,ostream& oss=cout); //CO20190321
  //input directly
  xstructure CreateSlab_SurfaceLattice(istream& input,const xvector<int>& hkl,int total_layers,double vacuum,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,ofstream& FileMESSAGE,double v3len_max_strict=AUROSTD_MAX_DOUBLE,ostream& oss=cout); //CO20190321  
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in,const xvector<int>& hkl,int total_layers,double vacuum,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,ofstream& FileMESSAGE,double v3len_max_strict=AUROSTD_MAX_DOUBLE,ostream& oss=cout);  //CO20190321
  xstructure CreateSlab_SurfaceLattice(istream& input,const xvector<int>& hkl,int total_layers,double vacuum,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,const _aflags& aflags,ofstream& FileMESSAGE,double v3len_max_strict=AUROSTD_MAX_DOUBLE,ostream& oss=cout); //CO20190321
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in,const xvector<int>& hkl,int total_layers,double vacuum,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,const _aflags& aflags,ofstream& FileMESSAGE,double v3len_max_strict=AUROSTD_MAX_DOUBLE,ostream& oss=cout); //CO20190321
  xstructure CreateSlab_SurfaceLattice(istream& input,const xvector<int>& hkl,int total_layers,double vacuum,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,double v3len_max_strict=AUROSTD_MAX_DOUBLE,ostream& oss=cout); //CO20190321
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in,const xvector<int>& hkl,int total_layers,double vacuum,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,double v3len_max_strict=AUROSTD_MAX_DOUBLE,ostream& oss=cout);  //CO20190321
  xstructure CreateSlab_SurfaceLattice(istream& input,const xvector<int>& hkl,int total_layers,double vacuum,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,const _aflags& aflags,double v3len_max_strict=AUROSTD_MAX_DOUBLE,ostream& oss=cout); //CO20190321
  xstructure CreateSlab_SurfaceLattice(const xstructure& xstr_in,const xvector<int>& hkl,int total_layers,double vacuum,xmatrix<double>& rotation,xstructure& xstr_slab_newbasis,vector<int>& sc2pcMap_slab,vector<int>& pc2scMap_slab,const _aflags& aflags,double v3len_max_strict=AUROSTD_MAX_DOUBLE,ostream& oss=cout);  //CO20190321

  bool slabTest(ostream& oss=cout);  //CO20190520
  bool slabTest(ofstream& FileMESSAGE,ostream& oss=cout);  //CO20190520
} // namespace slab

// ----------------------------------------------------------------------------
// aflow_defects.cpp
class acage {
  public:
    // constructors/destructors                                   // --------------------------------------
    acage();                                                      // constructor default
    acage(const acage& b);                                        // constructor copy
    ~acage();                                                     // destructor
    // CONTENT                                                    // --------------------------------------
    void clear(void);                                             // clear everything
    xvector<double> origin_fpos;// default 3D
    xvector<double> origin_cpos;// default 3D
    double radius;
    uint coordination_position;
    int cages_position;
    int cages_irrtype;
    deque<_atom> atoms;
  private:                                                        // ---------------------------------------
    void free();                                                  // to free everything
    void copy(const acage& b);                                    // the flag is necessary because sometimes you need to allocate the space.
};
class _isort_acage_radius {                   // sorting through reference
  public:
    bool operator()(const acage& a1, const acage& a2) const {
      return (bool) (a1.radius>a2.radius);}
};

bool GetSphereFromFourPoints(xvector<double>& orig,double& radius,const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3,const xvector<double>& v4);
bool GetCircumCircleeFromThreePoints(xvector<double>& orig,double& radius,const xvector<double>& v1,const xvector<double>& v2,const xvector<double>& v3);
bool GetCircleFromTwoPoints(xvector<double>& orig,double& radius,const xvector<double>& v1,const xvector<double>& v2);
bool FPositionInsideCell(const xvector<double>& r);
bool EmptySphere(const deque<_atom>& grid_atoms,const xvector<double>& origin_cpos,const double& radius);
bool EmptySphere(const xstructure& str,const xvector<double>& origin_cpos,const double& radius);
uint CoordinationPoint(const deque<_atom>& atoms,deque<_atom>& ratoms,const xvector<double>& point,const double& rmin,const double& rmax);
uint CoordinationPoint(const deque<_atom>& atoms,deque<_atom>& ratoms,const xvector<double>& point,const double& rmin);
uint CoordinationPoint(const xstructure& str,deque<_atom>& ratoms,const xvector<double>& point,const double& rmin,const double& rmax);
uint CoordinationPoint(const xstructure& str,deque<_atom>& ratoms,const xvector<double>& point,const double& rmin);
bool AddCageToCages(const xstructure& str,const xvector<double>& origin_cpos,const xvector<double>& origin_fpos,const double& radius,
    const int& cage_points_type,const double& roughness,
    vector<acage>& cages,vector<acage>& cagesX,
    const bool& osswrite1,ostream& oss1, const bool& osswrite2,ostream& oss2,int ithread);
uint GetCages4(const xstructure& str,const double& roughness,vector<acage>& cages,vector<acage>& cages4,
    const bool& osswrite1,ostream& oss1, const bool& osswrite2,ostream& oss2);
uint GetCages3(const xstructure& str,const double& roughness,vector<acage>& cages,vector<acage>& cages3,
    const bool& osswrite1,ostream& oss1, const bool& osswrite2,ostream& oss2);
uint GetCages2(const xstructure& str,const double& roughness,vector<acage>& cages,vector<acage>& cages2,
    const bool& osswrite1,ostream& oss1, const bool& osswrite2,ostream& oss2);
bool GetCages(const xstructure& _str,_aflags& aflags,
    vector<acage>& cagesirreducible,vector<acage>& cagesreducible,vector<acage>& cages4,
    vector<acage>& cages3,vector<acage>& cages2,const double& _roughness,const bool& osswrite,ostream& oss);
// ----------------------------------------------------------------------------
// aflow_pocc //CO20180502
namespace KBIN {
  void VASP_RunPOCC(const string& directory,ostream& oss=std::cout);  //CO20200624
  void VASP_RunPOCC(const string& directory,ofstream& FileMESSAGE,ostream& oss=std::cout);  //CO20200624
  void VASP_RunPOCC(const _xvasp& xvasp,const string& AflowIn,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags,ofstream& FileMESSAGE,ostream& oss=std::cout);
}
// ----------------------------------------------------------------------------
// aflow_phonons.cpp
namespace KBIN {
  bool relaxStructureAPL_VASP(int, const string&, aurostd::xoption&, const aurostd::xvector<int>&, bool, _xvasp&, _aflags&, _kflags&, _vflags&, ofstream&,ostream& oss=std::cout);  //ME20181107
  bool runRelaxationsAPL_VASP(int, const string&, _xvasp&, _aflags&, _kflags&, _vflags&, ofstream&);  //ME20200427
  void VASP_RunPhonons_APL(_xvasp &xvasp,string AflowIn,_aflags &aflags,_kflags &kflags,_vflags &vflags,ofstream &FileMESSAGE, ostream& oss=std::cout);
  void RunPhonons_APL(_xinput &xinput,string AflowIn,_aflags &aflags,_kflags &kflags,_xflags &xflags,ofstream &FileMESSAGE, ostream& oss=std::cout);  //now it's general
  void RunPhonons_APL_20181216(_xinput &xinput,string AflowIn,_aflags &aflags,_kflags &kflags,_xflags &xflags,ofstream &FileMESSAGE, ostream& oss=std::cout);  //now it's general //CO20181216
  // ----------------------------------------------------------------------------
  // aflow_agl_debye.cpp
  uint relaxStructureAGL_VASP(const string& AflowIn, _xvasp& xvasp, _aflags& aflags, _kflags& kflags, _vflags& vflags, ofstream& FileMessage);  //CT20200501
  void VASP_RunPhonons_AGL(_xvasp &xvasp,string AflowIn,_aflags &aflags,_kflags &kflags,_vflags &vflags,ofstream &FileMESSAGE);
  void VASP_RunPhonons_AGL_postprocess(const string& directory_LIB, string& AflowInName, string& FileLockName);  //CT20200624
  // ----------------------------------------------------------------------------
  // aflow_ael_elasticity.cpp
  uint relaxStructureAEL_VASP(const string& AflowIn, _xvasp& xvasp, _aflags& aflags, _kflags& kflags, _vflags& vflags, ofstream& FileMessage);  //CT20200501
  void VASP_RunPhonons_AEL(_xvasp &xvasp,string AflowIn,_aflags &aflags,_kflags &kflags,_vflags &vflags,ofstream &FileMESSAGE);
  void VASP_RunPhonons_AEL_postprocess(const string& directory_LIB, string& AflowInName, string& FileLockName);  //CT20200624
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------
// aconvasp_aflow.cpp

xstructure PutInCell(const xstructure& a);     // Bring all atoms in the cell (will be moved to external function)
xstructure PutInCompact(const xstructure& a);  // Bring all atoms in a compact shape (will be moved to external function)
xstructure GetPrim(const xstructure& a);
bool IsTranslationFVector(const xstructure& a,const xvector<double>& ftvec);
bool IsTranslationCVector(const xstructure& a,const xvector<double>& ctvec);
xvector<double> GetMom1(const xstructure& a);          // get moment_1 position of the atoms
xstructure SetMom1(const xstructure& a,const xvector<double>& mom1_in);     // set moment_1 position of atoms
xvector<double> AtomCDisp(const _atom& at1,const _atom& at2);
double AtomDist(const xstructure& str,const _atom& atom1,const _atom& atom2);  // with structure
double AtomDist(const _atom& at1, const _atom& at2);  // without structure
xvector<double> GetCDispFromOrigin(const _atom& atom);
double GetDistFromOrigin(const _atom& atom);
//void GetUnitCellRep(const xvector<double>& ppos,xvector<double>& p_cell0,xvector<int>& ijk,const xmatrix<double>& lattice,const bool coord_flag);
_atom ConvertAtomToLat(const _atom& in_at,const xmatrix<double>& lattice);
double GetXrayScattFactor(const string& name,double lambda=XRAY_RADIATION_COPPER_Kalpha,bool clean=true); //CO20190322
xmatrix<double> RecipLat(const xmatrix<double>& lat);
double Normal(const double& x,const double& mu,const double& sigma);
xstructure SetLat(const xstructure& a,const xmatrix<double>& in_lat);
xmatrix<double> GetLat(const xstructure& a);

namespace pflow {
  // [OBSOLETE] void GetNeighData(const deque<_atom>& in_atom_vec,const xstructure& in_str,const double& rmin,const double& rmax,deque<deque<_atom> >& neigh_mat);
  // [OBSOLETE] void GetStrNeighData(const xstructure& str,const double cutoff,deque<deque<_atom> >& neigh_mat);
  double GetVol(const xmatrix<double>& lat);
  double GetVol(const aurostd::matrix<double>& lat); //CO20200404 pflow::matrix()->aurostd::matrix()
  double GetSignedVol(const xmatrix<double>& lat);
  double GetSignedVol(const aurostd::matrix<double>& lat); //CO20200404 pflow::matrix()->aurostd::matrix()
  xmatrix<double> RecipLat(const xmatrix<double>& lat);
  aurostd::matrix<double> RecipLat(const aurostd::matrix<double>& lat);  //CO20200404 pflow::matrix()->aurostd::matrix()
  _atom SetCpos(const _atom& a,const vector<double>& in_cpos);
  _atom SetFpos(const _atom& a,const vector<double>& in_fpos);
  vector<double> vecF2C(const aurostd::matrix<double>& lat,const vector<double>& vf);  //CO20200404 pflow::matrix()->aurostd::matrix()
  vector<double> vecC2F(const aurostd::matrix<double>& lat,const vector<double>& vc);  //CO20200404 pflow::matrix()->aurostd::matrix()
  _atom SetName(const _atom& a,const string& in_name);
  _atom SetType(const _atom& a,const int in_type);
  _atom SetNum(const _atom& a,const int in_num);
  // [RF20200415 - duplicate from xatom]vector<int> GetTypes(const xstructure& a);
  // [RF20200415 - duplicate from xatom]vector<string> GetNames(const xstructure& a);
  // [RF20200415 - duplicate from xatom]vector<string> GetCleanNames(const xstructure& a);
  // [RF20200415 - duplicate from xatom]vector<double> GetSpins(const xstructure& a);
  aurostd::matrix<double> GetFpos(const xstructure& str);  //CO20200404 pflow::matrix()->aurostd::matrix()
  aurostd::matrix<double> GetCpos(const xstructure& str);  //CO20200404 pflow::matrix()->aurostd::matrix()
  xstructure SetNumEachType(const xstructure& a,const deque<int>& in_num_each_type);
  deque<int> GetNumEachType(const xstructure& a);
  xstructure SetLat(const xstructure& a,const aurostd::matrix<double>& in_lat);  //CO20200404 pflow::matrix()->aurostd::matrix()
  aurostd::matrix<double> GetLat(const xstructure& a); //CO20200404 pflow::matrix()->aurostd::matrix()
  double GetScale(const xstructure& a);
  aurostd::matrix<double> GetScaledLat(const xstructure& a); //CO20200404 pflow::matrix()->aurostd::matrix()
  xstructure AddAllAtomPos(const xstructure& a,const aurostd::matrix<double>& in_pos,const int in_coord_flag); //CO20200404 pflow::matrix()->aurostd::matrix()
  xstructure SetAllAtomPos(const xstructure& a,const aurostd::matrix<double>& in_pos,const int in_coord_flag); //CO20200404 pflow::matrix()->aurostd::matrix()
  xstructure SetAllAtomNames(const xstructure& a,const vector<string>& in_names);
  xstructure SetNamesWereGiven(const xstructure& a,const vector<int>& in_names_were_given);
  xstructure SetOrigin(const xstructure& a,const vector<double>& in_origin);
  xstructure SetOrigin(const xstructure& a,const xvector<double>& in_origin);
  bool VVequal(const vector<double>& a,const vector<double>& b);
  bool VVequal(const vector<int>& a,const vector<int>& b);
  bool VVequal(const deque<double>& a,const deque<double>& b);
  bool VVequal(const deque<int>& a,const deque<int>& b);
  vector<double> SmoothFunc(const vector<double>& func,const double& sigma);
  void VVset(aurostd::matrix<double>& mat,const double& value);  //CO20200404 pflow::matrix()->aurostd::matrix()
  void VVset(vector<vector< int> >& mat,const int& value);
  double norm(const vector<double>& v);
  double getcos(const vector<double>& a,const vector<double>& b);
  //  vector<double> Getabc_angles(const aurostd::matrix<double>& lat);   // confuses namespace  //CO20200404 pflow::matrix()->aurostd::matrix()
  vector<double> Sort_abc_angles(const vector<double>& abc_angles);
  void Vout(const vector<double>& a,ostream& out);
  void Vout(const vector<int>& a,ostream& out);
  void Vout(const vector<string>& a,ostream& out);
  void Mout(const aurostd::matrix<double>& m,ostream& out);  //CO20200404 pflow::matrix()->aurostd::matrix()
  void Mout(const vector<vector<double> >& m,ostream& out);
  vector<double> SVprod(const double& s,const vector<double>& b);
  vector<int> SVprod(const int& s,const vector<int>& b);
  vector<double> VVsum(const vector<double>& a,const vector<double>& b);
  vector<double> VVsum(const vector<double>& a,const vector<int>& b);
  vector<double> VVdiff(const vector<double>& a,const vector<double>& b);
  double VVprod(const vector<double>& a,const vector<double>& b);
  double VVprod(const vector<double>& a,const vector<int>& b);
  aurostd::matrix<double> MMmult(const aurostd::matrix<double>& a,const aurostd::matrix<double>& b); //CO20200404 pflow::matrix()->aurostd::matrix()
  vector<double> MVmult(const aurostd::matrix<double>& A,const vector<double>& v); //CO20200404 pflow::matrix()->aurostd::matrix()
  vector<double> VMmult(const vector<double>& v,const aurostd::matrix<double>& A); //CO20200404 pflow::matrix()->aurostd::matrix()
  vector<double> VMmult(const vector<int>& v,const aurostd::matrix<double>& A);  //CO20200404 pflow::matrix()->aurostd::matrix()
  vector<double> VVcross(const vector<double>& a,const vector<double>& b);
  double VVdot(const vector<double>& a,const vector<double>& b);
  int GetNumAtoms(const xstructure& a);
  void SetSpline(const vector<double>& x,const vector<double>& y,const double& yp1,const double& ypn,vector<double>& y2);
  void GetSplineInt(const vector<double>& xa,const vector<double>& ya,vector<double>& y2a,const double& x,double& y);
  void PrintSpline(const vector<double>& x,const vector<double>& y,const int& npts,ostream& outf);
}

// ----------------------------------------------------------------------------
// aflow_wyckoff.cpp

xvector<double> wv(const double &x,const double &y,const double &z);
void wa(_atom& a,xstructure &str);
xstructure WyckoffPOSITIONS(uint spacegroup, xstructure strin);
xstructure WyckoffPOSITIONS(uint spacegroup, uint option, xstructure strin);


// ----------------------------------------------------------------------------
// aflow_makefile.cpp

namespace makefile {
  void getDependencies(const string& filename,vector<string>& files_already_explored,vector<string>& dfiles);
  void getDependencies(const string& filename,vector<string>& files_already_explored,vector<string>& dfiles,bool& mt_required);
  //[CO20200508 - OBSOLETE]void replaceMakefileDefinitions(const vector<string>& vvariables,vector<string>& vdefinitions);
  //[CO20200508 - OBSOLETE]void replaceMakefileDefinitions(const vector<string>& vvariables,vector<vector<string> >& vvdefinitions);
  //[CO20200508 - OBSOLETE]void splitMakefileDefinitions(const string& definitions,vector<string>& vdefinitions);
  //[CO20200508 - OBSOLETE]void splitMakefileDefinitions(const vector<string>& vdefinitions,vector<vector<string> >& vvdefinitions);
  //[CO20200508 - OBSOLETE]void readMakefileVariables(const string& directory,vector<string>& vvariables,vector<vector<string> >& vvdefinitions);
  //[CO20200508 - OBSOLETE]void readMakefileVariables(const vector<string>& vlines,vector<string>& vvariables,vector<vector<string> >& vvdefinitions);
  void updateDependenciesAUROSTD(vector<string>& vdependencies);
  void updateDependenciesVariable(const vector<string>& vdeps_var,const string& var,vector<string>& vdeps_replace);
  void createMakefileAFLOW(const string& directory=".");
}

// ----------------------------------------------------------------------------
// aflow_apennsy stuff
#include "aflow_apennsy.h"

// ----------------------------------------------------------------------------
// aflowlib.h stuff
#include "aflowlib.h"

// ----------------------------------------------------------------------------
// aflowlib.h stuff
#include "aflowlib.h"

// ----------------------------------------------------------------------------
// aflow_pflow.h stuff
#include "aflow_pflow.h"

// ----------------------------------------------------------------------------
// aflow_xelement.h stuff
#include "aflow_xelement.h"




#endif
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *                                                                         *
// ***************************************************************************

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2021           *
// *                                                                         *
// ***************************************************************************
// this file contains the routines to prepare VASP input files
// Stefano Curtarolo - 2007 Duke
// fixed for xz - 2008 (SC)

#ifndef _AFLOW_IVASP_CPP
#define _AFLOW_IVASP_CPP

#include "aflow.h"
#define _incarpad_ 48
#define _IVASP_DOUBLE2STRING_PRECISION_ 7
#define DIELECTRIC_DK 0.1
#define DEFAULT_EFIELD_PEAD 0.001

// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// INPUT
namespace KBIN {
  bool VASP_Produce_INPUT(_xvasp& xvasp,const string& AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_vflags &vflags,bool load_POSCAR_from_xvasp) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(AflowIn.length()==0) {throw aurostd::xerror(_AFLOW_FILE_NAME_,"KBIN::VASP_Produce_INPUT():","empty AflowIn",_FILE_CORRUPT_);} //CO20200624
    bool Krun=TRUE;
    if(load_POSCAR_from_xvasp){
      if(Krun) Krun=(Krun && KBIN::VASP_Produce_POSCAR(xvasp));     // produce POSCAR before KPOINTS  //CO20180420 - good for POCC
    } else {
      if(Krun) Krun=(Krun && KBIN::VASP_Produce_POSCAR(xvasp,AflowIn,FileMESSAGE,aflags,kflags,vflags));     // produce POSCAR before KPOINTS
    }
    if(Krun) Krun=(Krun && KBIN::VASP_Produce_INCAR(xvasp,AflowIn,FileMESSAGE,aflags,kflags,vflags));
    if(LDEBUG) cerr << "VASP_Produce_INPUT: Checking if generate aflow.in only before calling VASP_Produce_POTCAR" << endl;  //CT20180719
    if (!XHOST.GENERATE_AFLOWIN_ONLY) { //CT20180719
      if(LDEBUG) cerr << "VASP_Produce_INPUT: Calling VASP_Produce_POTCAR" << endl;  //CT20180719
      if(Krun) Krun=(Krun && KBIN::VASP_Produce_POTCAR(xvasp,AflowIn,FileMESSAGE,aflags,kflags,vflags));
      if(LDEBUG) {  //CO20181226 check species_pp patching (may cause problems with auto aflow.in generation)
        cerr << "VASP_Produce_INPUT: Checking species_pp patch" << endl;
        for(uint i=0;i<xvasp.str.species_pp.size();i++){
          cerr << "xstr.species[" << i << "]=" << xvasp.str.species[i] << endl;
          cerr << "xstr.species_pp[" << i << "]=" << xvasp.str.species_pp[i] << endl;
        }
      }
      convertPOSCARFormat(xvasp, kflags); //ME20190220
    } //CT20180719
    if(Krun) Krun=(Krun && KBIN::VASP_Produce_KPOINTS(xvasp,AflowIn,FileMESSAGE,aflags,kflags,vflags));
    return Krun;
  }
}

namespace KBIN {
  bool VASP_Modify_INPUT(_xvasp& xvasp,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_vflags &vflags) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    bool Krun=TRUE;
    // if(Krun) Krun=(Krun && KBIN::VASP_Modify_POSCAR(xvasp,FileMESSAGE,aflags,vflags));  // moved up
    if(Krun) Krun=(Krun && KBIN::VASP_Modify_INCAR(xvasp,FileMESSAGE,aflags,kflags,vflags));
    if(LDEBUG) cerr << "VASP_Modify_INPUT: Checking if generate aflow.in only before calling VASP_Produce_POTCAR" << endl; //CT20180719
    if (!XHOST.GENERATE_AFLOWIN_ONLY) { //CT20180719
      if(LDEBUG) cerr << "VASP_Modify_INPUT: Calling VASP_Modify_POTCAR" << endl;  //CT20180719 
      if(Krun) Krun=(Krun && KBIN::VASP_Modify_POTCAR(xvasp,FileMESSAGE,aflags,vflags));
    } //CT20180719
    if(Krun) Krun=(Krun && KBIN::VASP_Modify_KPOINTS(xvasp,FileMESSAGE,aflags,vflags));
    return Krun;
  }
}

namespace KBIN {
  bool VASP_Produce_and_Modify_INPUT(_xvasp& xvasp,const string& AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_vflags &vflags,bool load_POSCAR_from_xvasp) {  //CO20180418
    bool Krun=TRUE;
    if(Krun) Krun=(Krun && KBIN::VASP_Produce_INPUT(xvasp,AflowIn,FileMESSAGE,aflags,kflags,vflags,load_POSCAR_from_xvasp));
    if(Krun) Krun=(Krun && KBIN::VASP_Modify_INPUT(xvasp,FileMESSAGE,aflags,kflags,vflags));
    return Krun;
  }
}


namespace KBIN {
  bool VASP_Write_ppAUID_FILE(const string& directory,const vector<string>& vppAUIDs,const vector<string>& species) {
    vector<string> vWRITE;
    vWRITE.push_back("AFLOW: AUID of pseudopotentials: NAUID="+aurostd::utype2string<uint>(vppAUIDs.size()));
    if(vppAUIDs.size()!=species.size()) {
      cout << "WARNING: KBIN::VASP_Write_ppAUID_FILE: vppAUIDs.size()=" << vppAUIDs.size() << "!=species.size()=" << species.size() << endl;
      cerr << "WARNING: KBIN::VASP_Write_ppAUID_FILE: vppAUIDs.size()=" << vppAUIDs.size() << "!=species.size()=" << species.size() << endl;
      return FALSE;
    }
    for(uint i=0;i<vppAUIDs.size();i++)
      vWRITE.push_back(vppAUIDs.at(i)+" "+aurostd::PaddedPOST(KBIN::VASP_PseudoPotential_CleanName(species.at(i)),3));
    aurostd::vectorstring2file(vWRITE,directory+"/"+DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT);
    return TRUE;
  }
  bool VASP_Write_ppAUID_FILE(const string& directory,const deque<string>& vppAUIDs,const deque<string>& species) {
    return VASP_Write_ppAUID_FILE(directory,aurostd::deque2vector(vppAUIDs),aurostd::deque2vector(species));
  }
  bool VASP_Write_ppAUID_AFLOWIN(const string& directory,const vector<string>& vppAUIDs,const vector<string>& species) {
    if(vppAUIDs.size()!=species.size()) {
      cout << "WARNING: KBIN::VASP_Write_ppAUID_AFLOWIN: vppAUIDs.size()=" << vppAUIDs.size() << "!=species.size()=" << species.size() << endl;
      cerr << "WARNING: KBIN::VASP_Write_ppAUID_AFLOWIN: vppAUIDs.size()=" << vppAUIDs.size() << "!=species.size()=" << species.size() << endl;
      return FALSE;
    }
    struct stat fileInfo;
    stat(string(directory+"/"+_AFLOWIN_).c_str(), &fileInfo);
    string date=std::ctime(&fileInfo.st_mtime);
    if (!date.empty() && date[date.length()-1] == '\n') date.erase(date.length()-1); // remove last newline

    string WRITE="[VASP_POTCAR_AUID]";
    for(uint i=0;i<vppAUIDs.size();i++) {
      WRITE+=vppAUIDs.at(i);
      if(i<vppAUIDs.size()-1) WRITE+=",";
    }
    aurostd::execute("echo \""+WRITE+"\" >> "+directory+"/"+_AFLOWIN_);
    aurostd::execute("touch -m --date=\""+date+"\" "+directory+"/"+_AFLOWIN_);
    return TRUE;
  }
  bool VASP_Write_ppAUID_AFLOWIN(const string& directory,const deque<string>& vppAUIDs,const deque<string>& species) {
    return VASP_Write_ppAUID_AFLOWIN(directory,aurostd::deque2vector(vppAUIDs),aurostd::deque2vector(species));
  }
}

namespace KBIN {
  bool VASP_Write_INPUT(_xvasp& xvasp,_vflags &vflags,const string& ext_module) {        // AFLOW_FUNCTION_IMPLEMENTATION
    string soliloquy=XPID+"KBIN::VASP_Write_INPUT():";
    ifstream DirectoryStream;
    DirectoryStream.open(xvasp.Directory.c_str(),std::ios::in);
    if(!DirectoryStream) {
      ostringstream aus;
      aus << "XXXXX  MAKING DIRECTORY = " << xvasp.Directory << endl;
      aurostd::PrintMessageStream(aus,XHOST.QUIET); // return FALSE;
      string str="mkdir "+xvasp.Directory;
      system(str.c_str());
    }
    DirectoryStream.close();
    if(XHOST.AVOID_RUNNING_VASP){  //CO20200624
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"VASP should NOT be running",_INPUT_ILLEGAL_);  //better to throw to avoid VASP_Backup(), etc.
      //return false;
    }
    bool Krun=TRUE;
    // VASP VASP WRITE
    if(Krun) Krun=(Krun && aurostd::stringstream2file(xvasp.POSCAR,string(xvasp.Directory+"/POSCAR")));
    if(Krun) Krun=(Krun && aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR")));
    if(Krun) Krun=(Krun && aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS")));
    if(Krun) Krun=(Krun && aurostd::stringstream2file(xvasp.POTCAR,string(xvasp.Directory+"/POTCAR")));
    if(Krun) Krun=(Krun && VASP_Write_ppAUID_FILE(xvasp.Directory,xvasp.POTCAR_AUID,xvasp.str.species));
    if(Krun) Krun=(Krun && VASP_Write_ppAUID_AFLOWIN(xvasp.Directory,xvasp.POTCAR_AUID,xvasp.str.species));

    // VASP BACKUP VASP WRITE
    if(Krun && xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed"))  Krun=(Krun && aurostd::stringstream2file(xvasp.POSCAR_orig,string(xvasp.Directory+"/POSCAR.orig"+ext_module)));
    if(Krun && xvasp.aopts.flag("FLAG::XVASP_INCAR_changed"))   Krun=(Krun && aurostd::stringstream2file(xvasp.INCAR_orig,string(xvasp.Directory+"/INCAR.orig"+ext_module)));
    if(Krun && xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed")) Krun=(Krun && aurostd::stringstream2file(xvasp.KPOINTS_orig,string(xvasp.Directory+"/KPOINTS.orig"+ext_module)));
    if(Krun && xvasp.aopts.flag("FLAG::XVASP_POTCAR_changed"))  Krun=(Krun && aurostd::stringstream2file(xvasp.POTCAR_orig,string(xvasp.Directory+"/POTCAR.orig"+ext_module)));
    // AS20210302 adding ext_module to be able to distinguish between different orig files (for example, we do not want APL to overwrite an existing POSCAR.orig file)

    if(vflags.KBIN_VASP_INCAR_VERBOSE) {;} // DUMMY

    return Krun;
  }
}

// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// PseudoPotential_CleanName
// gets rid of all junk in the name
namespace KBIN {
  string VASP_PseudoPotential_CleanName(const string& speciesIN) {return VASP_PseudoPotential_CleanName_20190712(speciesIN);} //CO20190712
  string VASP_PseudoPotential_CleanName_20190712(const string& speciesIN) { //CO20190712
    //the old function assumed only a single species input
    //now, the species input can be a full species string: Mn_pvPt, etc.
    //need to remove ALL instances of pp info
    //no longer need for loops
    //also, this new function leverages InPlace substitution, instead of creating new strings every time (FASTER)

    string species=speciesIN;
    aurostd::VASP_PseudoPotential_CleanName_InPlace(species);
    return species;
  }
  string VASP_PseudoPotential_CleanName_20190101(const string& speciesIN) {
    string species=speciesIN;
    uint i=1,imax=2;
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_GW");  //CO20190712
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_new");  //CO20190712
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_AE");  //CO20190712
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_sv");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_pv");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_d2");  //CO20190712 - BEFORE _d //potpaw_LDA/potpaw_LDA.05May2010/As_d2_GW
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_d");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_soft");  //CO20190712 - BEFORE _s
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_s");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_200eV");
    //[CO20190712 - moved up before _s]for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_soft");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_2_n");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_h");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_1");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_2");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_3");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"1.75"); //CO20190712 - potpaw_LDA.52/potpaw_LDA.52.19Apr2012/H1.75
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"1.66"); //CO20190712 - potpaw_LDA.52/potpaw_LDA.52.19Apr2012/H1.66
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"1.33"); //CO20190712 - potpaw_PBE.54/potpaw_PBE.54.04Sep2015/H1.33
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"1.25"); //CO20190712 - before all other decimal numbers
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"1.5"); //CO20190712 - potpaw_PBE/potpaw_PBE.06May2010/H1.5
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,".75");  //CO20190712 - before 0.5
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,".25");  //CO20190712 - potpaw_LDA.52/potpaw_LDA.52.19Apr2012/H.25
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,".66"); //CO20190712 - potpaw_LDA.52/potpaw_LDA.52.19Apr2012/H.66
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,".33"); //CO20190712 - potpaw_PBE.54/potpaw_PBE.54.04Sep2015/H.33
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,".42"); //CO20190712 - potpaw_PBE.54/potpaw_PBE.54.04Sep2015/H.42
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,".58"); //CO20190712 - before 0.5 //potpaw_LDA.52/potpaw_LDA.52.19Apr2012/H.58
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,".5");
    //[CO20190712 - moved up 0.5]for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,".75");
    //[CO20190712 - moved up before all other decimal numbers]for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"1.25");

    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"+1");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"+3");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"+5");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"+7");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"-1");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"-3");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"-5");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"-7");
    //  for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"1");
    //CO START
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,DEFAULT_VASP_POTCAR_DIR_POT_LDA+"/");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,DEFAULT_VASP_POTCAR_DIR_POT_GGA+"/");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,DEFAULT_VASP_POTCAR_DIR_POT_PBE+"/");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA+"/");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA+"/");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE+"/");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN+"/");
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN+"/");
    //CO END
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"__"); //CO20190712 - BEFORE _ - potpaw_LDA/potpaw_LDA.05May2010/Si_sv_GW__
    for(i=1;i<=imax;i++) species=aurostd::RemoveSubStringFirst(species,"_");  //CO20190712  //potpaw_LDA/potpaw_LDA.05May2010/Si_sv_GW_
    return species;
  }
}

namespace KBIN {
  bool VASP_PseudoPotential_CleanName_TEST(void){ //CO20190712
    bool LDEBUG=(TRUE || XHOST.DEBUG);
    string soliloquy = XPID + "KBIN::VASP_PseudoPotential_CleanName_TEST():";
    stringstream message;

    ostream& oss=cout;
    ofstream FileMESSAGE;
    _aflags aflags; aflags.Directory=XHOST.vflag_control.getattachedscheme("DIRECTORY_CLEAN");

    string compound_to_test="";
    vector<string> velements;
    vector<double> vcomposition;

    compound_to_test="W_sv_GWCs_sv_GWCr_sv_GWPt_ZORALi_AE_GW2";
    bool keep_pp=false;
    velements=aurostd::getElements(compound_to_test,vcomposition,pp_string,FileMESSAGE,false,true,keep_pp,oss);
    if(LDEBUG){cerr << soliloquy << " compound_to_test=\"" << compound_to_test << "\", velements=" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(velements,"\""),",") << endl;}
    if(aurostd::joinWDelimiter(velements,"")!="CrCsLiPtW"){
      message << "getElements() failed [" << compound_to_test << " with keep_pp=" << keep_pp << "]";
      pflow::logger(_AFLOW_FILE_NAME_, soliloquy, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
      return false;
    }

    compound_to_test="W_sv_GWCs_sv_GWCr_sv_GWPt_ZORALi_AE_GW2";
    keep_pp=true;
    velements=aurostd::getElements(compound_to_test,vcomposition,pp_string,FileMESSAGE,false,true,keep_pp,oss);
    if(LDEBUG){cerr << soliloquy << " compound_to_test=\"" << compound_to_test << "\", velements=" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(velements,"\""),",") << endl;}
    if(aurostd::joinWDelimiter(velements,"")!="Cr_sv_GWCs_sv_GWLi_AE_GW2Pt_ZORAW_sv_GW"){
      message << "getElements() failed [" << compound_to_test << " with keep_pp=" << keep_pp << "]";
      pflow::logger(_AFLOW_FILE_NAME_, soliloquy, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
      return false;
    }

    compound_to_test="H1.25";
    keep_pp=false;
    velements=aurostd::getElements(compound_to_test,vcomposition,pp_string,FileMESSAGE,false,true,keep_pp,oss);
    if(LDEBUG){cerr << soliloquy << " compound_to_test=\"" << compound_to_test << "\", velements=" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(velements,"\""),",") << endl;}
    if(aurostd::joinWDelimiter(velements,"")!="H"){
      message << "getElements() failed [" << compound_to_test << " with keep_pp=" << keep_pp << "]";
      pflow::logger(_AFLOW_FILE_NAME_, soliloquy, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
      return false;
    }

    compound_to_test="H1.25";
    keep_pp=true;
    velements=aurostd::getElements(compound_to_test,vcomposition,pp_string,FileMESSAGE,false,true,keep_pp,oss);
    if(LDEBUG){cerr << soliloquy << " compound_to_test=\"" << compound_to_test << "\", velements=" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(velements,"\""),",") << endl;}
    if(aurostd::joinWDelimiter(velements,"")!="H1.25"){
      message << "getElements() failed [" << compound_to_test << " with keep_pp=" << keep_pp << "]";
      pflow::logger(_AFLOW_FILE_NAME_, soliloquy, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
      return false;
    }

    compound_to_test="Mn.5W5H1.25";
    keep_pp=false;
    velements=aurostd::getElements(compound_to_test,vcomposition,composition_string,FileMESSAGE,false,true,keep_pp,oss);
    if(LDEBUG){cerr << soliloquy << " compound_to_test=\"" << compound_to_test << "\", velements=" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(velements,"\""),",") << ", vcomposition=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vcomposition,5),",") << endl;}
    if(!(aurostd::joinWDelimiter(velements,"")=="HMnW" && vcomposition.size()==3 && aurostd::isequal(vcomposition[0],1.25) && aurostd::isequal(vcomposition[1],0.5) && aurostd::isequal(vcomposition[2],5.0) )){
      message << "getElements() failed [" << compound_to_test << " with keep_pp=" << keep_pp << "]";
      pflow::logger(_AFLOW_FILE_NAME_, soliloquy, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
      return false;
    }

    //return true;

    //full test of pp available
    vector<string> directories_2_search;
    vector<string> potcar_paths;
    vector<string> path_parts;
    vector<string> element_parts,_element_parts;
    string element_raw="",element_clean="";
    string element="",suffix="";
    aurostd::string2tokens(DEFAULT_VASP_POTCAR_DIRECTORIES,directories_2_search,",");
    for(uint i=0;i<directories_2_search.size();i++){
      if(aurostd::IsDirectory(directories_2_search[i])){
        const string& dir=directories_2_search[i];
        if(LDEBUG){cerr << soliloquy << " searching dir=" << dir << endl;}
        aurostd::string2vectorstring(aurostd::execute2string("find "+dir+" -name POTCAR*"),potcar_paths);
        if(LDEBUG){cerr << soliloquy << " potcar_paths=" << aurostd::joinWDelimiter(potcar_paths," ") << endl;}
        for(uint j=0;j<potcar_paths.size();j++){
          if(!aurostd::RemoveWhiteSpacesFromTheBack(potcar_paths[j]).empty()){
            const string& potcar_path=potcar_paths[j];
            aurostd::string2tokens(potcar_path,path_parts,"/");
            if(aurostd::WithinList(path_parts,"TESTS")){continue;}  //skip VASP/TESTS
            if(aurostd::WithinList(path_parts,"BENCHS")){continue;}  //skip VASP/BENCHS
            if(path_parts.size()<2){continue;}
            element_raw=path_parts[path_parts.size()-2];
            if(aurostd::substring2bool(element_raw,"runelements")){continue;} //garbage pp
            if(aurostd::substring2bool(element_raw,"Free")){continue;} //garbage pp
            if(LDEBUG){cerr << soliloquy << " element_raw=" << element_raw << endl;}
            aurostd::string2tokens(element_raw,element_parts,"_");
            if(element_parts.size()==0){continue;}
            element=suffix="";
            if(element_parts.size()>1){
              _element_parts.clear();
              for(uint k=1;k<element_parts.size();k++){_element_parts.push_back(element_parts[k]);} //skip element
              suffix="_"+aurostd::joinWDelimiter(_element_parts,"_");
            }
            element=element_parts[0];
            element=KBIN::VASP_PseudoPotential_CleanName(element);  //clean H1.25
            if(LDEBUG){cerr << soliloquy << " element=\"" << element << "\", suffix=\"" << suffix << "\"" << endl;}
            element_clean=KBIN::VASP_PseudoPotential_CleanName(element_raw);
            if(element!=element_clean){ //idempotent
              message << "\"" << element << "\" != \"" << element_clean << "\" == KBIN::VASP_PseudoPotential_CleanName(\"" << element_raw << "\")";
              pflow::logger(_AFLOW_FILE_NAME_, soliloquy, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
              return false;
            }
            if(GetAtomNumber(element)==0){
              message << "(GetAtomNumber(\"" << element << "\")==0 [element_raw=\"" << element_raw << "\",element_clean=\"" << element_clean << "\"" << "]";
              pflow::logger(_AFLOW_FILE_NAME_, soliloquy, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
              return false;
            }
          }
        }
      }
    }
    return true;
  }
}

namespace KBIN {
  uint VASP_SplitAlloySpecies(string alloy_in, vector<string> &speciesX) {
    return XATOM_SplitAlloySpecies(alloy_in,speciesX);
  }
}

namespace KBIN {
  uint VASP_SplitAlloySpecies(string alloy_in, vector<string> &speciesX, vector<double> &natomsX) {
    return XATOM_SplitAlloySpecies(alloy_in,speciesX,natomsX);
  }
}

namespace KBIN {
  bool VASP_SplitAlloySpecies(string alloy_in, string &speciesA, string &speciesB) {
    string alloy=KBIN::VASP_PseudoPotential_CleanName(KBIN::VASP_PseudoPotential_CleanName(alloy_in));
    if(0) {
      int i=0;
      speciesA=alloy[i++];
      if(alloy[i]>='a' && alloy[i]<='z') speciesA=speciesA+alloy[i++];
      speciesB=alloy[i++];
      if(alloy[i]>='a' && alloy[i]<='z') speciesB=speciesB+alloy[i++];
    }
    if(1) {
      speciesA="";speciesB="";
      int speciesN=0;
      for(uint i=0;i<alloy.length();i++) {
        if(alloy[i]>='A' && alloy[i]<='Z') speciesN++;
        if(speciesN==1) speciesA+=alloy[i];
        if(speciesN==2) speciesB+=alloy[i];
      }
    }
    speciesA=aurostd::CleanStringASCII(speciesA);
    speciesB=aurostd::CleanStringASCII(speciesB);
    return TRUE;
  }
}

namespace KBIN {
  bool VASP_SplitAlloySpecies(string alloy_in, string &speciesA, string &speciesB, string &speciesC) {
    string alloy=KBIN::VASP_PseudoPotential_CleanName(KBIN::VASP_PseudoPotential_CleanName(alloy_in));
    speciesA="";speciesB="";speciesC="";
    int speciesN=0;
    for(uint i=0;i<alloy.length();i++) {
      if(alloy[i]>='A' && alloy[i]<='Z') speciesN++;
      if(speciesN==1) speciesA+=alloy[i];
      if(speciesN==2) speciesB+=alloy[i];
      if(speciesN==3) speciesC+=alloy[i];
    }
    speciesA=aurostd::CleanStringASCII(speciesA);
    speciesB=aurostd::CleanStringASCII(speciesB);
    speciesC=aurostd::CleanStringASCII(speciesC);
    return TRUE;
  }
}

namespace KBIN {
  bool VASP_SplitAlloySpecies(vector<string> alloy, vector<string> &speciesA, vector<string> &speciesB) {
    for (uint k=0;k<alloy.size();k++) {
      speciesA.push_back("");speciesB.push_back("");
      KBIN::VASP_SplitAlloySpecies(alloy.at(k),speciesA.at(k),speciesB.at(k));
    }
    return TRUE;
  }
}

namespace KBIN {
  bool VASP_SplitAlloySpecies(vector<string> alloy, vector<string> &speciesA, vector<string> &speciesB, vector<string> &speciesC) {
    for (uint k=0;k<alloy.size();k++) {
      speciesA.push_back("");speciesB.push_back("");
      KBIN::VASP_SplitAlloySpecies(alloy.at(k),speciesA.at(k),speciesB.at(k),speciesC.at(k));
    }
    return TRUE;
  }
}

namespace KBIN {
  uint VASP_SplitAlloyPseudoPotentials(string alloy_in, vector<string> &species_ppX) {
    return XATOM_SplitAlloyPseudoPotentials(alloy_in,species_ppX);
  }
}

namespace KBIN {
  uint VASP_SplitAlloyPseudoPotentials(string alloy_in, vector<string> &species_ppX, vector<double> &natomsX) {
    return XATOM_SplitAlloyPseudoPotentials(alloy_in,species_ppX,natomsX);
  }
}

namespace KBIN {
  bool VASP_SplitAlloyPseudoPotentials(string alloy, string &species_ppA, string &species_ppB) {
    if(0) {
      uint i=0;
      species_ppA=alloy[i++];
      while((alloy[i]>='Z' || alloy[i]<='A') && i<=alloy.length()) species_ppA=species_ppA+alloy[i++];
      species_ppB=alloy[i++];
      while((alloy[i]>='Z' || alloy[i]<='A') && i<=alloy.length()) species_ppB=species_ppB+alloy[i++];
    }
    if(1) {
      species_ppA="";species_ppB="";
      int speciesN=0;
      for(uint i=0;i<alloy.length();i++) {
        if(alloy[i]>='A' && alloy[i]<='Z') speciesN++;
        if(speciesN==1) species_ppA+=alloy[i];
        if(speciesN==2) species_ppB+=alloy[i];
      }
    }
    species_ppA=aurostd::CleanStringASCII(species_ppA);
    species_ppB=aurostd::CleanStringASCII(species_ppB);
    return TRUE;
  }
}

namespace KBIN {
  bool VASP_SplitAlloyPseudoPotentials(string alloy, string &species_ppA, string &species_ppB, string &species_ppC) {
    species_ppA="";species_ppB="";species_ppC="";
    int speciesN=0;
    for(uint i=0;i<alloy.length();i++) {
      if(alloy[i]>='A' && alloy[i]<='Z') speciesN++;
      if(speciesN==1) species_ppA+=alloy[i];
      if(speciesN==2) species_ppB+=alloy[i];
      if(speciesN==3) species_ppC+=alloy[i];
    }
    species_ppA=aurostd::CleanStringASCII(species_ppA);
    species_ppB=aurostd::CleanStringASCII(species_ppB);
    species_ppC=aurostd::CleanStringASCII(species_ppC);
    return TRUE;
  }
}

namespace KBIN {
  bool VASP_SplitAlloyPseudoPotentials(vector<string> alloy, vector<string> &pseudosA, vector<string> &pseudosB) {
    for (uint k=0;k<alloy.size();k++) {
      pseudosA.push_back("");pseudosB.push_back("");
      KBIN::VASP_SplitAlloyPseudoPotentials(alloy.at(k),pseudosA.at(k),pseudosB.at(k));
    }
    return TRUE;
  }
}

namespace KBIN {
  bool VASP_SplitAlloyPseudoPotentials(vector<string> alloy, vector<string> &pseudosA, vector<string> &pseudosB, vector<string> &pseudosC) {
    for (uint k=0;k<alloy.size();k++) {
      pseudosA.push_back("");pseudosB.push_back("");
      KBIN::VASP_SplitAlloyPseudoPotentials(alloy.at(k),pseudosA.at(k),pseudosB.at(k),pseudosC.at(k));
    }
    return TRUE;
  }
}


// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// INCAR
namespace KBIN {
  bool VASP_Produce_INCAR(_xvasp& xvasp,const string& AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_vflags &vflags) { // AFLOW_FUNCTION_IMPLEMENTATION
    if(AflowIn.length()==0) {throw aurostd::xerror(_AFLOW_FILE_NAME_,"KBIN::VASP_Produce_INCAR():","empty AflowIn",_FILE_CORRUPT_);} //CO20200624
    if(!kflags.AFLOW_MODE_VASP) {cerr << XPID << "KBIN::VASP_Produce_INCAR: should kflags.AFLOW_MODE_VASP be set ??" << endl;}
    ostringstream aus;
    bool Krun=TRUE;
    xvasp.INCAR.str(std::string());
    xvasp.NCPUS=0;
    xvasp.aopts.flag("FLAG::XVASP_INCAR_generated",FALSE);
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",FALSE);

    aus << "00000  MESSAGE INCAR   generation in " << xvasp.Directory << "  " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    bool KBIN_VASP_INCAR_MODE_EMPTY=!vflags.KBIN_VASP_INCAR_MODE.flag("IMPLICIT") && !vflags.KBIN_VASP_INCAR_MODE.flag("EXPLICIT") && !vflags.KBIN_VASP_INCAR_MODE.flag("EXTERNAL");

    // IMPLICIT or EXPLICIT or EXTERNAL for INCAR
    Krun=(Krun && (vflags.KBIN_VASP_INCAR_MODE.flag("IMPLICIT") ||
          vflags.KBIN_VASP_INCAR_MODE.flag("EXPLICIT") ||
          vflags.KBIN_VASP_INCAR_MODE.flag("EXTERNAL") || KBIN_VASP_INCAR_MODE_EMPTY));
    if(!Krun) {
      aurostd::StringstreamClean(aus);
      aus << "EEEEE  [VASP_INCAR_MODE_IMPLICIT] or [VASP_INCAR_MODE_EXPLICIT] or [VASP_INCAR_MODE_EXTERNAL] must be specified " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl; //ME20181113
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      Krun=FALSE;
      return Krun;
    }
    // EMPTY ************************************************** INCAR
    if(Krun && KBIN_VASP_INCAR_MODE_EMPTY) {  // [VASP_INCAR_MODE_EMPTY] construction
      aus << "00000  MESSAGE INCAR   generation EMPTY file from " << _AFLOWIN_ << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.INCAR << "#AFLOW INCAR automatically generated" << endl;
    }
    // IMPLICIT ************************************************** INCAR
    if(Krun && vflags.KBIN_VASP_INCAR_MODE.flag("IMPLICIT")) {  // [VASP_INCAR_MODE_IMPLICIT] construction
      aus << "00000  MESSAGE INCAR   generation IMPLICIT file from " << _AFLOWIN_ << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(vflags.KBIN_VASP_INCAR_FILE.flag("SYSTEM_AUTO")) {
        xvasp.INCAR << "#AFLOW INCAR automatically generated" << endl;
        xvasp.INCAR << "SYSTEM=" << xvasp.str.title << endl;
        xvasp.INCAR << "#PROTOTYPE=" << xvasp.str.prototype << endl;
        xvasp.INCAR << "#INFO=" << xvasp.str.info << endl;
        // if(LDEBUG) cerr << xvasp.INCAR.str() << endl;
      } else if (vflags.AFLOW_SYSTEM.isentry) {  //ME20181121
        xvasp.INCAR << "SYSTEM=" << vflags.AFLOW_SYSTEM.content_string << std::endl;
      }
    }
    // EXPLICIT ************************************************** INCAR
    if(Krun && vflags.KBIN_VASP_INCAR_MODE.flag("EXPLICIT")) {  // [VASP_INCAR_MODE_EXPLICIT] construction
      if(vflags.KBIN_VASP_INCAR_FILE.flag("KEYWORD") && !vflags.KBIN_VASP_INCAR_MODE.flag("EXPLICIT_START_STOP")) {
        aus << "00000  MESSAGE INCAR   generation EXPLICIT file from " << _AFLOWIN_ << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        xvasp.INCAR << vflags.KBIN_VASP_INCAR_EXPLICIT.str(); //ME20181113
        // [OBSOLETE]    aurostd::ExtractToStringstreamEXPLICIT(FileAFLOWIN,xvasp.INCAR,"[VASP_INCAR_FILE]");
        // [ME20181226 OBSOLETE]    aurostd::ExtractToStringstreamEXPLICIT(AflowIn,xvasp.INCAR,"[VASP_INCAR_FILE]");
      } else if(!vflags.KBIN_VASP_INCAR_FILE.flag("KEYWORD") && vflags.KBIN_VASP_INCAR_MODE.flag("EXPLICIT_START_STOP")) {
        aus << "00000  MESSAGE INCAR   generation EXPLICIT file from " << _AFLOWIN_ << " with START/STOP  " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        xvasp.INCAR << vflags.KBIN_VASP_INCAR_EXPLICIT_START_STOP.str();
        // [ME20181226 OBSOLETE] if(aurostd::substring2bool(AflowIn,"[VASP_INCAR_MODE_EXPLICIT]START") && aurostd::substring2bool(AflowIn,"[VASP_INCAR_MODE_EXPLICIT]STOP"))
        // [OBSOLETE]	aurostd::ExtractToStringstreamEXPLICIT(FileAFLOWIN,xvasp.INCAR,"[VASP_INCAR_MODE_EXPLICIT]START","[VASP_INCAR_MODE_EXPLICIT]STOP");
        // [ME20181226 OBSOLETE] aurostd::ExtractToStringstreamEXPLICIT(AflowIn,xvasp.INCAR,"[VASP_INCAR_MODE_EXPLICIT]START","[VASP_INCAR_MODE_EXPLICIT]STOP");
      } else {
        aus << "EEEEE  [VASP_INCAR_MODE_EXPLICIT] do not confuse aflow !!" << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aus << "EEEEE  [VASP_INCAR_MODE_EXPLICIT] Possible modes " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aus << "----------------------------------------------------------------------------------------------------" << endl;
        aus << "[AFLOW] INCAR EXPLICIT MODE without START/STOP (default)" << endl;
        aus << "[VASP_INCAR_MODE_EXPLICIT]" << endl;
        aus << "[VASP_INCAR_FILE]SYSTEM=AuTi.274_LDA" << endl;
        aus << "[VASP_INCAR_FILE]#Prototype Ni2In" << endl;
        aus << "[VASP_INCAR_FILE]PREC=med" << endl;
        aus << "[VASP_INCAR_FILE]ISMEAR=1" << endl;
        aus << "[VASP_INCAR_FILE]SIGMA=0.2" << endl;
        aus << "[VASP_INCAR_FILE]IBRION=2" << endl;
        aus << "[VASP_INCAR_FILE]NSW=51" << endl;
        aus << "[VASP_INCAR_FILE]ISIF=3" << endl;
        aus << "[VASP_INCAR_FILE]ENMAX=333.666" << endl;
        aus << "[VASP_INCAR_FILE]NBANDS=47" << endl;
        aus << "[VASP_INCAR_FILE]MAGMOM=  5 5 5 5 5 5" << endl;
        aus << "[VASP_INCAR_FILE]ISPIND=2" << endl;
        aus << "[VASP_INCAR_FILE]ISPIN=2" << endl;
        aus << "[AFLOW]" << endl;
        aus << "----------------------------------------------------------------------------------------------------" << endl;
        aus << "[AFLOW] INCAR EXPLICIT MODE with START/STOP" << endl;
        aus << "[VASP_INCAR_MODE_EXPLICIT]" << endl;
        aus << "[VASP_INCAR_MODE_EXPLICIT]START" << endl;
        aus << "SYSTEM=AuTi.274_LDA" << endl;
        aus << "#Prototype Ni2In" << endl;
        aus << "PREC=med" << endl;
        aus << "ISMEAR=1" << endl;
        aus << "SIGMA=0.2" << endl;
        aus << "IBRION=2" << endl;
        aus << "NSW=51" << endl;
        aus << "ISIF=3" << endl;
        aus << "ENMAX=333.666" << endl;
        aus << "NBANDS=47" << endl;
        aus << "MAGMOM=  5 5 5 5 5 5" << endl;
        aus << "ISPIND=2" << endl;
        aus << "ISPIN=2" << endl;
        aus << "[VASP_INCAR_MODE_EXPLICIT]STOP" << endl;
        aus << "[AFLOW]" << endl;
        aus << "----------------------------------------------------------------------------------------------------" << endl;
        aus << "EEEEE  [VASP_INCAR_MODE_EXPLICIT] Note " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aus << "EEEEE  [VASP_INCAR_MODE_EXPLICIT]START must be present and no [VASP_INCAR_FILE]" << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aus << "EEEEE  [VASP_INCAR_MODE_EXPLICIT]STOP  must be present and no [VASP_INCAR_FILE]" << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aus << "EEEEE  or [VASP_INCAR_FILE] present and NO START/STOP" << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
    }
    // EXTERNAL **************************************************
    if(Krun && vflags.KBIN_VASP_INCAR_MODE.flag("EXTERNAL")) {  // [VASP_INCAR_MODE_EXTERNAL] construction
      string file;
      aus << "00000  MESSAGE INCAR   generation EXTERNAL file from " << _AFLOWIN_ << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(vflags.KBIN_VASP_INCAR_FILE.flag("COMMAND") && vflags.KBIN_VASP_INCAR_FILE.flag("FILE")) {
        aus << "EEEEE   [VASP_INCAR_MODE]FILE=  and  [VASP_INCAR_MODE]COMMAND=  can not be used together " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
      if(!vflags.KBIN_VASP_INCAR_FILE.flag("COMMAND") && (vflags.KBIN_VASP_INCAR_FILE.flag("FILE") || !vflags.KBIN_VASP_INCAR_FILE.flag("FILE"))) {
        if(vflags.KBIN_VASP_INCAR_FILE.flag("FILE")) {
          // [ME20181226 OBSOLETE] file=aurostd::substring2string(AflowIn,"[VASP_INCAR_FILE]FILE=",TRUE);
          file = vflags.KBIN_VASP_INCAR_FILE.getattachedscheme("FILE"); //ME20181113
          aus << "00000  MESSAGE INCAR   generation from file=" << file << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        } else {
          file=DEFAULT_VASP_EXTERNAL_INCAR;
          aus << "00000  MESSAGE INCAR   generation from DEFAULT file=" << DEFAULT_VASP_EXTERNAL_INCAR << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if(!aurostd::FileExist(file)) {
          aus << "EEEEE  ERROR INCAR file=" << file << " does not exist! " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        if(aurostd::FileEmpty(file)) {
          aus << "EEEEE  ERROR INCAR file=" << file << " is empty! " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        xvasp.INCAR << aurostd::file2string(file);
      }
      if(vflags.KBIN_VASP_INCAR_FILE.flag("COMMAND") && !vflags.KBIN_VASP_INCAR_FILE.flag("FILE")) {
        // [ME20181226 OBSOLETE] file=aurostd::substring2string(AflowIn,"[VASP_INCAR_FILE]COMMAND=",FALSE);
        file = vflags.KBIN_VASP_INCAR_FILE.getattachedscheme("COMMAND"); //ME20181113
        aus << "00000  MESSAGE INCAR   generation from command= '" << file << "' " << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        file=file+" > ./_aflow_INCAR."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";    // create temp  //CO20200502 - threadID
        aurostd::execute(file);                           // create temp
        file="./_aflow_INCAR."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";            // file name  //CO20200502 - threadID
        if(!aurostd::FileExist(file)) {  // could not write (directory protected)
          aus << "EEEEE  ERROR INCAR file=" << file << " does not exist! " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        if(aurostd::FileEmpty(file)) {  // contains nothing good
          aus << "EEEEE  ERROR INCAR file=" << file << " is empty! " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        xvasp.INCAR << aurostd::file2string(file);       // load INCAR
        aurostd::RemoveFile("./_aflow_INCAR."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp");     // remove temp //CO20200502 - threadID
      }
    }
    // INCAR DONE **************************************************
    xvasp.INCAR << "#incar" << endl;
    xvasp.INCAR_orig << xvasp.INCAR.str();
    xvasp.aopts.flag("FLAG::XVASP_INCAR_generated",TRUE);
    return Krun;
  };  // KBIN::VASP_Produce_INCAR
}



namespace KBIN {
  bool VASP_Modify_INCAR(_xvasp& xvasp,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_vflags &vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    ostringstream aus;
    bool Krun=TRUE;

    // bool vflags.KBIN_VASP_INCAR_VERBOSE=TRUE;
    if(Krun && kflags.KBIN_MPI) {
      xvasp.NCPUS=kflags.KBIN_MPI_NCPUS;
      if(kflags.KBIN_MPI_AUTOTUNE) {
        aus << "00000  MESSAGE INCAR-MPI: found AUTOTUNE option " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aus << "00000  MESSAGE INCAR-MPI: input files WILL be auto-tuned for PARALLEL execution with " << kflags.KBIN_MPI_NCPUS << " CPUs " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        KBIN::VASP_MPI_Autotune(xvasp,aflags,vflags.KBIN_VASP_INCAR_VERBOSE);
        xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
      } else {
        aus << "00000  MESSAGE INCAR-MPI: AUTOTUNE option NOT found! (aflow_ivasp.cpp) " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aus << "00000  MESSAGE INCAR-MPI: input files MUST be appropriate for PARALLEL execution with " << kflags.KBIN_MPI_NCPUS << " CPUs " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    }
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_SYSTEM_AUTO.isentry) {                                                        /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]SYSTEM_AUTO - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_System_Auto(xvasp,vflags.KBIN_VASP_INCAR_VERBOSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    } else if (vflags.AFLOW_SYSTEM.isentry) {  //ME20181121
      xvasp.INCAR << "SYSTEM=" << vflags.AFLOW_SYSTEM.content_string << std::endl;
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // if(Krun && !vflags.KBIN_VASP_RUN.flag("STATIC") && !vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("STATIC") && ! vflags.KBIN_VASP_RUN_RELAX_STATIC_PATCH_STATIC) //o change for relax if STATIC !!!
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.isentry && 
        (vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("ALL")               || vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS") ||
         vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_SHAPE")        || vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_VOLUME") ||
         vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS_CELL_VOLUME"))) {      /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("ALL"))                 aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]RELAX_ALL - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS"))                aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]RELAX_IONS - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_SHAPE"))          aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]RELAX_CELL_SHAPE - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_VOLUME"))         aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]RELAX_CELL_VOLUME - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS_CELL_VOLUME"))    aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]RELAX_IONS_CELL_VOLUME - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.isentry) {
        if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.xscheme=="ENERGY")         aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]RELAX_MODE=ENERGY (default) - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.xscheme=="FORCES")         aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]RELAX_MODE=FORCES - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.xscheme=="ENERGY_FORCES")  aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]RELAX_MODE=ENERGY_FORCES - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.xscheme=="FORCES_ENERGY")  aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]RELAX_MODE=FORCES_ENERGY - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      } else {
        aus << "00000  MESSAGE-DEFAULT RELAX_MODE=" << "ENERGY" << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_Relax_ON(xvasp,vflags,1); // relax number (start)
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // with KBIN_VASP_FORCE_OPTION_EDIFFG_EQUAL.isentry && KBIN_VASP_FORCE_OPTION_EDIFFG_EQUAL.content_double
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_EDIFFG_EQUAL.isentry) {      /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]EDIFFG=" << vflags.KBIN_VASP_FORCE_OPTION_EDIFFG_EQUAL.content_double << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("EDIFFG",xvasp,vflags,"",0,0.0,FALSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // with KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry) {             /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]NBANDS - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("NBANDS",xvasp,vflags,"",0,0.0,FALSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // with KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.isentry && KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.content_int
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.isentry && vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.content_int>0) {      /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]NBANDS=" << vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.content_int << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("NBANDS",xvasp,vflags,"",0,0.0,FALSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // with KBIN_VASP_FORCE_OPTION_PSTRESS_EQUAL.isentry && KBIN_VASP_FORCE_OPTION_PSTRESS_EQUAL.content_double
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_PSTRESS_EQUAL.isentry) {      /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]PSTRESS=" << vflags.KBIN_VASP_FORCE_OPTION_PSTRESS_EQUAL.content_double << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("PSTRESS",xvasp,vflags,"",0,0.0,FALSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // with KBIN_VASP_FORCE_OPTION_POTIM_EQUAL.isentry && KBIN_VASP_FORCE_OPTION_POTIM_EQUAL.content_double
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_POTIM_EQUAL.isentry && vflags.KBIN_VASP_FORCE_OPTION_POTIM_EQUAL.content_double>0) {      /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]POTIM=" << vflags.KBIN_VASP_FORCE_OPTION_POTIM_EQUAL.content_double << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("POTIM",xvasp,vflags,"",0,0.0,FALSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // SPIN
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                     /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_SPIN.isentry) {
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]SPIN=" << (vflags.KBIN_VASP_FORCE_OPTION_SPIN.option?"ON":"OFF") << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      } else {
        //      aus << "00000  MESSAGE-DEFAULT SPIN=" << (DEFAULT_VASP_FORCE_OPTION_SPIN?"ON":"OFF") << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aus << "00000  MESSAGE-DEFAULT SPIN=" << "NEGLECT" << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(vflags.KBIN_VASP_FORCE_OPTION_SPIN.isentry) 
        KBIN::XVASP_INCAR_PREPARE_GENERIC("SPIN",xvasp,vflags,"",0,0.0,vflags.KBIN_VASP_FORCE_OPTION_SPIN.option);  // CHANGE ONLY IF SPIN IS MENTIONED
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
      if(vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1) {
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]SPIN=REMOVE_RELAX_1 - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      if(vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2) {
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]SPIN=REMOVE_RELAX_2 - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    }

    // LSCOUPLING must be after spin
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                     /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_LSCOUPLING.isentry) {
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]LSCOUPLING=" << (vflags.KBIN_VASP_FORCE_OPTION_LSCOUPLING.option?"ON":"OFF") << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      } else {
        aus << "00000  MESSAGE-DEFAULT LSCOUPLING=" << (DEFAULT_VASP_FORCE_OPTION_LSCOUPLING?"ON":"OFF") << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(vflags.KBIN_VASP_FORCE_OPTION_LSCOUPLING.isentry || DEFAULT_VASP_FORCE_OPTION_LSCOUPLING) KBIN::XVASP_INCAR_PREPARE_GENERIC("LS_COUPLING",xvasp,vflags,"",0,0.0,vflags.KBIN_VASP_FORCE_OPTION_LSCOUPLING.option);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // AUTO_MAGMOM must be after LSCOUPLING, AFTER SPIN
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                     /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_SPIN.isentry&&vflags.KBIN_VASP_FORCE_OPTION_SPIN.option) {   //CO, MAGMOM should only be on if SPIN is SPECIFIED and ON (as default is to neglect)
        if(vflags.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.isentry) {
          aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]AUTO_MAGMOM=" << (vflags.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.option?"ON":"OFF") << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        } else {
          aus << "00000  MESSAGE-DEFAULT AUTO_MAGMOM=" << (DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM?"ON":"OFF") << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        }
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        if(vflags.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.isentry || DEFAULT_VASP_FORCE_OPTION_AUTO_MAGMOM) KBIN::XVASP_INCAR_PREPARE_GENERIC("AUTO_MAGMOM",xvasp,vflags,"",0,0.0,vflags.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.option);
        xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
      }
    }

    // BADER 
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                     /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry) {
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]BADER=" << (vflags.KBIN_VASP_FORCE_OPTION_BADER.option?"ON":"OFF") << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      } else {
        aus << "00000  MESSAGE-DEFAULT BADER=" << (DEFAULT_VASP_FORCE_OPTION_BADER?"ON":"OFF") << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      //  if(vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry || DEFAULT_VASP_FORCE_OPTION_BADER) KBIN::XVASP_INCAR_BADER(xvasp,vflags,vflags.KBIN_VASP_FORCE_OPTION_BADER.option); WILL BE DONE WHEN NEEDED
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // ELF 
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                     /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry) {
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]ELF=" << (vflags.KBIN_VASP_FORCE_OPTION_ELF.option?"ON":"OFF") << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      } else {
        aus << "00000  MESSAGE-DEFAULT ELF=" << (DEFAULT_VASP_FORCE_OPTION_ELF?"ON":"OFF") << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      //  if(vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry || DEFAULT_VASP_FORCE_OPTION_ELF) KBIN::XVASP_INCAR_ELF(xvasp,vflags,vflags.KBIN_VASP_FORCE_OPTION_ELF.option); WILL BE DONE WHEN NEEDED
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // NSW_EQUAL
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_NSW_EQUAL) {           /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_NSW_EQUAL_VALUE>0) {
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]NSW=" << vflags.KBIN_VASP_FORCE_OPTION_NSW_EQUAL_VALUE << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        KBIN::XVASP_INCAR_PREPARE_GENERIC("NSW",xvasp,vflags,"",0,0.0,FALSE);	
        xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
      }
    }

    // SYM
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                     /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_SYM.isentry) {
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]SYM=" << (vflags.KBIN_VASP_FORCE_OPTION_SYM.option?"ON":"OFF") << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      } else {
        aus << "00000  MESSAGE-DEFAULT SYM=" << (DEFAULT_VASP_FORCE_OPTION_SYM?"ON":"OFF") << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("SYM",xvasp,vflags,"",0,0.0,vflags.KBIN_VASP_FORCE_OPTION_SYM.option);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // WAVECAR
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                     /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_WAVECAR.isentry) {
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]WAVECAR=" << (vflags.KBIN_VASP_FORCE_OPTION_WAVECAR.option?"ON":"OFF") << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      } else {
        aus << "00000  MESSAGE-DEFAULT WAVECAR=" << (DEFAULT_VASP_FORCE_OPTION_WAVECAR?"ON":"OFF") << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("WAVECAR",xvasp,vflags,"",0,0.0,vflags.KBIN_VASP_FORCE_OPTION_WAVECAR.option);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // CHGCAR
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                     /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_CHGCAR.isentry) {
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]CHGCAR=" << (vflags.KBIN_VASP_FORCE_OPTION_CHGCAR.option?"ON":"OFF") << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      } else {
        aus << "00000  MESSAGE-DEFAULT CHGCAR=" << (DEFAULT_VASP_FORCE_OPTION_CHGCAR?"ON":"OFF") << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("CHGCAR",xvasp,vflags,"",0,0.0,vflags.KBIN_VASP_FORCE_OPTION_CHGCAR.option);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    //ME20191028
    if (Krun && vflags.KBIN_VASP_FORCE_OPTION_CHGCAR_FILE.isentry) {
      string chgcar = vflags.KBIN_VASP_FORCE_OPTION_CHGCAR_FILE.content_string;
      if (chgcar[0] != '/') chgcar = aurostd::CleanFileName(aflags.Directory + "/" + chgcar);  // relative path
      if (aurostd::FileExist(chgcar)) {
        if (aurostd::IsCompressed(chgcar)) {
          string ext = aurostd::GetCompressionExtension(chgcar);
          aurostd::CopyFile(chgcar, aflags.Directory + "/CHGCAR" + ext);  // relative path
          aurostd::UncompressFile(aflags.Directory + "/CHGCAR" + ext);
        } else {
          aurostd::CopyFile(chgcar, aflags.Directory + "/CHGCAR");
        }
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        KBIN::XVASP_INCAR_PREPARE_GENERIC("ICHARG", xvasp, vflags, "", 1, 0.0, true);
        xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
      } else {
        string function = "KBIN::VASP_Modify_INCAR()";
        string message = "Cannot use CHGCAR file " + chgcar + ". File not found.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_NOT_FOUND_);
      }
    }

    // LDAU0
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_LDAU0.isentry) {          /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]LDAU=OFF - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_LDAU_OFF(xvasp,vflags.KBIN_VASP_INCAR_VERBOSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // LDAU1
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_LDAU1.isentry) {          /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]LDAU1=ON - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(vflags.KBIN_VASP_LDAU_SPECIES!="") aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]LDAU_SPECIES=\"" << vflags.KBIN_VASP_LDAU_SPECIES << "\" - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      if(vflags.KBIN_VASP_LDAU_PARAMETERS!="") aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]LDAU_PARAMETERS=\"" << vflags.KBIN_VASP_LDAU_PARAMETERS << "\" - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]LDAU_AFLOW_AUTO_flag=\"" << vflags.KBIN_VASP_LDAU_AFLOW_AUTO_flag << "\" - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_LDAU_ON(xvasp,vflags,1);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // LDAU2
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_LDAU2.isentry) {          /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]LDAU2=ON - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(vflags.KBIN_VASP_LDAU_SPECIES!="") aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]LDAU_SPECIES=\"" << vflags.KBIN_VASP_LDAU_SPECIES << "\" - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      if(vflags.KBIN_VASP_LDAU_PARAMETERS!="") aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]LDAU_PARAMETERS=\"" << vflags.KBIN_VASP_LDAU_PARAMETERS << "\" - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]LDAU_AFLOW_AUTO_flag=\"" << vflags.KBIN_VASP_LDAU_AFLOW_AUTO_flag << "\" - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_LDAU_ON(xvasp,vflags,2);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // LDAU_ADIABATIC
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.isentry) {          /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]LDAU=ADIABATIC (steps=" << vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.content_int << ") - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      //  KBIN::XVASP_INCAR_LDAU_ADIABATIC_OFF(xvasp,vflags.KBIN_VASP_INCAR_VERBOSE); // NOTHING TO CHANGE, IT IS ON THE FLY
      //  xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);                                                  // NOTHING TO CHANGE, IT IS ON THE FLY
    }

    // LDAU_CUTOFF
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_LDAU_CUTOFF.isentry) {          /*************** INCAR **************/
      //    cerr << "vflags.KBIN_VASP_FORCE_OPTION_LDAU_CUTOFF.isentry=" << vflags.KBIN_VASP_FORCE_OPTION_LDAU_CUTOFF.isentry << endl;
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]LDAU=CUTOFF - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      //  KBIN::XVASP_INCAR_LDAU_CUTOFF_OFF(xvasp,vflags.KBIN_VASP_INCAR_VERBOSE); // NOTHING TO CHANGE, IT IS ON THE FLY
      //  xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);                                                  // NOTHING TO CHANGE, IT IS ON THE FLY
    }

    // PREC
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                     /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_PREC.isentry==TRUE) {                                                         /*************** INCAR **************/
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]PREC=" << vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_PREC.preserved)   aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]PREC_preserved - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      } else {
        aus << "00000  MESSAGE-DEFAULT PREC=" << DEFAULT_VASP_FORCE_OPTION_PREC_SCHEME << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_Precision(xvasp,vflags);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // with KBIN_VASP_FORCE_OPTION_ENMAX_MULTIPLY_EQUAL.isentry
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_ENMAX_MULTIPLY_EQUAL.isentry) {             /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]ENMAX_MULTIPLY="<< vflags.KBIN_VASP_FORCE_OPTION_ENMAX_MULTIPLY_EQUAL.content_double << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("ENMAX_MULTIPLY",xvasp,vflags,"",0,0.0,FALSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // ALGO
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                    /*************** INCAR **************/
      //ME20181108 - PREC=PHONONS overrides ALGO
      if (vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme == "PHONONS") {
        aus << "00000  MESSAGE-OPTION ALGO = NORMAL (override by PREC=PHONONS) - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        vflags.KBIN_VASP_FORCE_OPTION_ALGO.scheme2scheme(vflags.KBIN_VASP_FORCE_OPTION_ALGO.xscheme, "NORMAL");
      } else if(vflags.KBIN_VASP_FORCE_OPTION_ALGO.isentry) {                                                              /*************** INCAR **************/
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]ALGO=" << vflags.KBIN_VASP_FORCE_OPTION_ALGO.xscheme << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved)   aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]ALGO_PRESERVED - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      } else {
        aus << "00000  MESSAGE-DEFAULT ALGO=" << DEFAULT_VASP_FORCE_OPTION_ALGO_SCHEME << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("ALGO",xvasp,vflags,"",0,0.0,FALSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // METAGGA
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                    /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.isentry) {                                                              /*************** INCAR **************/
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]METAGGA=" << vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;;
      } else {
        aus << "00000  MESSAGE-DEFAULT METAGGA=" << DEFAULT_VASP_FORCE_OPTION_METAGGA_SCHEME << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_Metagga(xvasp,vflags);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // IVDW
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                    /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_IVDW.isentry) {                                                              /*************** INCAR **************/
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]IVDW=" << vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;;
      } else {
        aus << "00000  MESSAGE-DEFAULT IVDW=" << DEFAULT_VASP_FORCE_OPTION_IVDW_SCHEME << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_Ivdw(xvasp,vflags);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // ABMIX
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                    /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_ABMIX.isentry) {                                                              /*************** INCAR **************/
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]ABMIX=" << vflags.KBIN_VASP_FORCE_OPTION_ABMIX.xscheme << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        // the rest is neglected... no AMIX BMIX AMIX_MAG BMIX_MAG
      } else {
        // neglected     aus << "00000  MESSAGE-DEFAULT ABMIX=" << DEFAULT_VASP_FORCE_OPTION_ABMIX_SCHEME << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_ABMIX(xvasp,vflags);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // TYPE
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                    /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_TYPE.isentry) {                                                              /*************** INCAR **************/
        if(vflags.KBIN_VASP_FORCE_OPTION_TYPE.xscheme.at(0)=='D') aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]TYPE=" << vflags.KBIN_VASP_FORCE_OPTION_TYPE.xscheme << "  - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_TYPE.xscheme.at(0)=='M') aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]TYPE=" << vflags.KBIN_VASP_FORCE_OPTION_TYPE.xscheme << "  - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_TYPE.xscheme.at(0)=='S') aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]TYPE=" << vflags.KBIN_VASP_FORCE_OPTION_TYPE.xscheme << "  - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_TYPE.xscheme.at(0)=='I') aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]TYPE=" << vflags.KBIN_VASP_FORCE_OPTION_TYPE.xscheme << "  - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      } else {
        aus << "00000  MESSAGE-DEFAULT TYPE=" << DEFAULT_VASP_FORCE_OPTION_TYPE_SCHEME << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("TYPE",xvasp,vflags,vflags.KBIN_VASP_FORCE_OPTION_TYPE.xscheme,0,0.0,FALSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    //CO20181128
    //this must come AFTER TYPE
    // with KBIN_VASP_FORCE_OPTION_ISMEAR_EQUAL.isentry && KBIN_VASP_FORCE_OPTION_ISMEAR_EQUAL.content_int
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_EQUAL.isentry) {      /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]ISMEAR=" << vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_EQUAL.content_int << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("ISMEAR",xvasp,vflags,"",0,0.0,FALSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    //CO20181128
    //this must come AFTER TYPE
    // with KBIN_VASP_FORCE_OPTION_SIGMA_EQUAL.isentry && KBIN_VASP_FORCE_OPTION_SIGMA_EQUAL.content_double
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_SIGMA_EQUAL.isentry) {      /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]SIGMA=" << vflags.KBIN_VASP_FORCE_OPTION_SIGMA_EQUAL.content_double << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("SIGMA",xvasp,vflags,"",0,0.0,FALSE);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // CONVERT_UNIT_CELL
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {                                                    /*************** INCAR **************/
      if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.isentry) {                                                              /*************** INCAR **************/
        if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_PRIMITIVE"))  aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]CONVERT_UNIT_CELL=STANDARD_PRIMITIVE - "<< Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_CONVENTIONAL"))  aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]CONVERT_UNIT_CELL=STANDARD_CONVENTIONAL - "<< Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("NIGGLI"))  aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]CONVERT_UNIT_CELL=NIGGLI - "<< Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("MINKOWSKI"))  aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]CONVERT_UNIT_CELL=MINKOWSKI - "<< Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("INCELL"))  aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]CONVERT_UNIT_CELL=INCELL - "<< Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("COMPACT"))  aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]CONVERT_UNIT_CELL=COMPACT - "<< Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("WIGNERSEITZ"))  aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]CONVERT_UNIT_CELL=WIGNERSEITZ - "<< Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("CARTESIAN"))  aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]CONVERT_UNIT_CELL=CARTESIAN - "<< Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("FRACTIONAL"))  aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]CONVERT_UNIT_CELL=FRACTIONAL - "<< Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        if(vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("PRESERVE"))  aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]CONVERT_UNIT_CELL=PRESERVE - "<< Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl; //CO
      }
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }

    // print the AFIX
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.isentry) {                                                    /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]IGNORE_AFIX=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.content_string << "  - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      for(uint i=0;i<vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.vxscheme.size();i++) 
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]IGNORE_AFIX=" << vflags.KBIN_VASP_FORCE_OPTION_IGNORE_AFIX.vxscheme.at(i) << "  - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }

    // DO THE PAW CORRECTIONS
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && 0) {                                                    /*************** INCAR **************/
      if(xvasp.POTCAR_PAW==TRUE) {
        aus << "00000  MESSAGE-DEFAULT PAW_CORRECTIONS" << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        KBIN::XVASP_INCAR_PREPARE_GENERIC("PAW_CORRECTIONS",xvasp,vflags,"",0,0.0,FALSE);
        xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
      }
    }

    // STATIC at the end after all the relax stuff has been added
    if(Krun && vflags.KBIN_VASP_RUN.flag("STATIC")) {                                                                      /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_RUN_STATIC] - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_Static_ON(xvasp,vflags);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("STATIC")) {                                                             /*************** INCAR **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]STATIC - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_Static_ON(xvasp,vflags);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }

    // these functions are forced in aflow_kvasp.cpp inside the RELAX_STATIC_BANDS run
    // if(Krun && vflags.KBIN_VASP_RUN_RELAX_STATIC_PATCH_STATIC) {                                                   /*************** INCAR **************/
    // aus << "00000  MESSAGE-OPTION  [VASP_RUN_RELAX_STATIC] - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
    // aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // KBIN::XVASP_INCAR_Relax_Static_ON(xvasp,vflags.KBIN_VASP_INCAR_VERBOSE);
    // xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    // }
    // if(Krun && vflags.KBIN_VASP_RUN_RELAX_STATIC_BANDS_PATCH_STATIC_BANDS) {                                       /*************** INCAR **************/
    // aus << "00000  MESSAGE-OPTION  [VASP_RUN_RELAX_STATIC_BANDS] - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
    // aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // KBIN::XVASP_INCAR_Relax_Static_Bands_ON(xvasp,vflags.KBIN_VASP_INCAR_VERBOSE);     // FIX
    // xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    // }


    // INTERCEPT ERRORS AND WRAP UP after all the incar has been prepared
    if(Krun && aurostd::substring2bool(xvasp.INCAR,"LEPSILON") && !aurostd::substring2bool(xvasp.INCAR,"#LEPSILON")) {  /*************** INCAR **************/
      //ME20191205 - also remove NCORE
      aus << "00000  MESSAGE REMOVE ENTRIES NPAR and NCORE because of LEPSILON - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;  //ME20191205
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_REMOVE_ENTRY(xvasp,"NPAR","LEPSILON",vflags.KBIN_VASP_INCAR_VERBOSE);
      KBIN::XVASP_INCAR_REMOVE_ENTRY(xvasp,"NCORE","LEPSILON",vflags.KBIN_VASP_INCAR_VERBOSE);  //ME20191205
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    }
    if(Krun && aurostd::substring2bool(xvasp.INCAR,"LCALCEPS") && !aurostd::substring2bool(xvasp.INCAR,"#LCALCEPS")) {  /*************** INCAR **************/
      //ME20191205 - also remove NCORE
      aus << "00000  MESSAGE REMOVE ENTRIES NPAR and NCORE because of LCALCEPS - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;  //ME20191205
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_REMOVE_ENTRY(xvasp,"NPAR","LCALCEPS",vflags.KBIN_VASP_INCAR_VERBOSE);
      KBIN::XVASP_INCAR_REMOVE_ENTRY(xvasp,"NCORE","LCALCEPS",vflags.KBIN_VASP_INCAR_VERBOSE);  //ME20191205
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    } // KEVIN
    if(Krun && aurostd::substring2bool(xvasp.INCAR,"IBRION") && !aurostd::substring2bool(xvasp.INCAR,"#IBRION")) {  /*************** INCAR **************/
      uint IBRION=aurostd::substring2utype<uint>(xvasp.INCAR.str(),"IBRION=");
      if(IBRION==8) {
        aus << "00000  MESSAGE REMOVE ENTRIES NPAR and NCORE because of IBRION=" << IBRION << "  - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;  //ME20191205
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        KBIN::XVASP_INCAR_REMOVE_ENTRY(xvasp,"NPAR","IBRION=8",vflags.KBIN_VASP_INCAR_VERBOSE);
        KBIN::XVASP_INCAR_REMOVE_ENTRY(xvasp,"NCORE","IBRION=8",vflags.KBIN_VASP_INCAR_VERBOSE);  //ME20191205
        xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
      }
    }

    // ------------------------------------
    // end
    xvasp.aopts.flag("FLAG::XVASP_INCAR_generated",TRUE);
    return Krun;
  };  // KBIN::VASP_Produce_INCAR
}


namespace KBIN {
  bool VASP_Reread_INCAR(_xvasp& xvasp,ofstream &FileMESSAGE,_aflags &aflags) { // AFLOW_FUNCTION_IMPLEMENTATION
    ostringstream aus;
    bool Krun=TRUE;
    if(!aurostd::FileExist(xvasp.Directory+"/INCAR")) {
      aus << "EEEEE  KBIN::VASP_Reread_INCAR: INCAR not present in directory: " << xvasp.Directory << " - "  << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      Krun=FALSE;
      return Krun;
    }
    xvasp.INCAR_orig.str(std::string()); xvasp.INCAR_orig << xvasp.INCAR.str();
    xvasp.INCAR.str(std::string()); xvasp.INCAR << aurostd::file2string(xvasp.Directory+"/INCAR"); // DID REREAD
    // xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    return Krun;
  }
}

// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// POSCAR
namespace KBIN {
  bool VASP_Produce_POSCAR(_xvasp& xvasp,const string& AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_vflags &vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(AflowIn.length()==0) {throw aurostd::xerror(_AFLOW_FILE_NAME_,"KBIN::VASP_Produce_POSCAR():","empty AflowIn",_FILE_CORRUPT_);} //CO20200624
    if(!kflags.AFLOW_MODE_VASP) {cerr << XPID << "KBIN::VASP_Produce_POSCAR: should kflags.AFLOW_MODE_VASP be set ??" << endl;}
    ostringstream aus;
    bool Krun=TRUE;
    xvasp.POSCAR.str(std::string());xvasp.POSCAR.clear();
    xvasp.POSCAR_orig.str(std::string());xvasp.POSCAR_orig.clear();
    xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",FALSE);
    xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",FALSE);
    //
    // IMPLICIT or EXPLICIT or EXTERNAL for POSCAR
    Krun=(Krun && (vflags.KBIN_VASP_POSCAR_MODE.flag("IMPLICIT") ||
          vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT") ||
          vflags.KBIN_VASP_POSCAR_MODE.flag("EXTERNAL")));
    if(!Krun) {
      aurostd::StringstreamClean(aus);
      aus << "EEEEE  [VASP_POSCAR_MODE_IMPLICIT] or [VASP_POSCAR_MODE_EXPLICIT] or [VASP_POSCAR_MODE_EXPLICIT] must be specified " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      Krun=FALSE;
      return Krun;
    }

    // IMPLICIT **************************************************
    if(Krun && vflags.KBIN_VASP_POSCAR_MODE.flag("IMPLICIT")) {  // [VASP_POSCAR_MODE_IMPLICIT] construction
      aus << "00000  MESSAGE POSCAR  generation IMPLICIT file from " << _AFLOWIN_ << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(!vflags.KBIN_VASP_POSCAR_FILE.flag("PROTOTYPE")) {
        aus << "EEEEE  [VASP_POSCAR_FILE] In POSCAR_MODE_IMPLICIT you must specify PROTOTYPE " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      } else {
        std::vector<string> tokens,tokens2,atomABC;
        std::string structure,label,parameters="";  // FIX NRL
        vector<double> volumeABC;
        structure=aurostd::substring2string(AflowIn,"[VASP_POSCAR_FILE]PROTOTYPE=",TRUE);
        aurostd::string2tokens(structure,tokens,";");
        label=tokens[0];
        for(uint i=1;i<tokens.size();i++) {
          // find SPECIES
          if(aurostd::substring2bool(tokens[i],"SPECIES=",TRUE) || aurostd::substring2bool(tokens[i],"SPECIE=",TRUE)) {
            aurostd::string2tokens(aurostd::substring2string(tokens[i],"=",TRUE),tokens2,",");
            for(uint j=0;j<tokens2.size();j++)
              atomABC.push_back(tokens2[j]);
          }
          // find VOLUMES
          if(aurostd::substring2bool(tokens[i],"VOLUMES=",TRUE) || aurostd::substring2bool(tokens[i],"VOLUME=",TRUE)) {
            aurostd::string2tokens(aurostd::substring2string(tokens[i],"=",TRUE),tokens2,",");
            for(uint j=0;j<tokens2.size();j++)
              volumeABC.push_back(aurostd::string2utype<double>(tokens2[j]));
          }
        }
        // for(uint j=0;j<atomABC.size();j++) cerr << atomABC.at(j) << endl;
        // for(uint j=0;j<volumeABC.size();j++) cerr << volumeABC.at(j) << endl;
        bool done=FALSE;
        if(atomABC.size()==2 && volumeABC.size()==0) {
          done=TRUE;
          deque<string> atomX;deque<double> volumeX;
          for(uint isp=0;isp<=1;isp++) {atomX.push_back(atomABC[isp]);volumeX.push_back(GetAtomVolume(atomABC[isp]));} //CO20181129
          //[CO20181106]for(uint isp=0;isp<=1;isp++) {atomX.push_back(atomABC[isp]);volumeX.push_back(GetAtomVolume(KBIN::VASP_PseudoPotential_CleanName(atomABC[isp])));}
          xvasp.str=aflowlib::PrototypeLibraries(FileMESSAGE,label,parameters,atomX,volumeX,-1.0,LIBRARY_MODE_HTQC);
          // xvasp.str=aflowlib::PrototypeLibraries(FileMESSAGE,label,parameters,atomABC[0],GetAtomVolume(KBIN::VASP_PseudoPotential_CleanName(atomABC[0])),atomABC[1],GetAtomVolume(KBIN::VASP_PseudoPotential_CleanName(atomABC[1])),-1.0); // OLD WAY
        }
        if(atomABC.size()==2 && volumeABC.size()==1) {
          done=TRUE;
          deque<string> atomX;deque<double> volumeX;
          for(uint isp=0;isp<=1;isp++) {atomX.push_back(atomABC[isp]);volumeX.push_back(0.0);}
          // xvasp.str=aflowlib::PrototypeLibraries(FileMESSAGE,label,parameters,atomABC[0],0.0,atomABC[1],0.0,volumeABC[0]); // OLD WAY
          xvasp.str=aflowlib::PrototypeLibraries(FileMESSAGE,label,parameters,atomX,volumeX,volumeABC[0],LIBRARY_MODE_HTQC);
        }
        if(atomABC.size()==2 && volumeABC.size()==2) {
          done=TRUE;
          deque<string> atomX;deque<double> volumeX;
          for(uint isp=0;isp<=1;isp++) {atomX.push_back(atomABC[isp]);volumeX.push_back(volumeABC[isp]);};
          // xvasp.str=aflowlib::PrototypeLibraries(FileMESSAGE,label,parameters,atomABC[0],volumeABC[0],atomABC[1],volumeABC[1],-1.0); // OLD WAY
          xvasp.str=aflowlib::PrototypeLibraries(FileMESSAGE,label,parameters,atomX,volumeX,-1.0,LIBRARY_MODE_HTQC);
        }
        if(done==FALSE) {
          aus << "EEEEE  POSCAR_MODE_IMPLICIT error in the PROTOTYPE definition" << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aus << "EEEEE  [VASP_POSCAR_FILE]PROTOTYPE=" << structure << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        // done

        xvasp.POSCAR << xvasp.str;
        // cerr << xvasp.POSCAR.str() << endl;
      }
    }
    // EXPLICIT **************************************************
    if(Krun && vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT")) {  // [VASP_POSCAR_MODE_EXPLICIT] construction
      if(vflags.KBIN_VASP_POSCAR_FILE.flag("KEYWORD") && !vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP") && !vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT")) {
        aus << "00000  MESSAGE POSCAR  generation EXPLICIT file from " << _AFLOWIN_ << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // [OBSOLETE]    aurostd::ExtractToStringstreamEXPLICIT(FileAFLOWIN,xvasp.POSCAR,"[VASP_POSCAR_FILE]");
        aurostd::ExtractToStringstreamEXPLICIT(AflowIn,xvasp.POSCAR,"[VASP_POSCAR_FILE]");
        xvasp.str=xstructure(xvasp.POSCAR,IOVASP_AUTO);  // load structure
        xvasp.POSCAR.str(std::string());xvasp.POSCAR.clear();
        xvasp.str.iomode=IOVASP_POSCAR;
        xvasp.POSCAR << xvasp.str;
      } else if(!vflags.KBIN_VASP_POSCAR_FILE.flag("KEYWORD") && (vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP") || vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT"))) {
        aus << "00000  MESSAGE POSCAR  generation EXPLICIT file from " << _AFLOWIN_ << " with START/STOP  " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        // normal get ONE of ONE
        if(vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP")) {
          if(aurostd::substring2bool(AflowIn,"[VASP_POSCAR_MODE_EXPLICIT]START") &&
              aurostd::substring2bool(AflowIn,"[VASP_POSCAR_MODE_EXPLICIT]STOP"))
            // [OBSOLETE]	  aurostd::ExtractLastToStringstreamEXPLICIT(FileAFLOWIN,xvasp.POSCAR,"[VASP_POSCAR_MODE_EXPLICIT]START","[VASP_POSCAR_MODE_EXPLICIT]STOP");
            aurostd::ExtractLastToStringstreamEXPLICIT(AflowIn,xvasp.POSCAR,"[VASP_POSCAR_MODE_EXPLICIT]START","[VASP_POSCAR_MODE_EXPLICIT]STOP");
          xvasp.str=xstructure(xvasp.POSCAR,IOVASP_AUTO);   // load structure
        }
        // get ONE of MANY
        if(vflags.KBIN_VASP_POSCAR_MODE.flag("EXPLICIT_START_STOP_POINT")) {
          xvasp.str=vflags.KBIN_VASP_POSCAR_MODE_EXPLICIT_VSTRUCTURE.at(xvasp.POSCAR_index);
        }
        // GOT IT
        xvasp.POSCAR.str(std::string());xvasp.POSCAR.clear();
        xvasp.str.iomode=IOVASP_POSCAR;
        xvasp.POSCAR << xvasp.str;
      } else {
        aus << "EEEEE  [VASP_POSCAR_MODE_EXPLICIT] do not confuse aflow !!" << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aus << "EEEEE  [VASP_POSCAR_MODE_EXPLICIT] Possible modes " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aus << "----------------------------------------------------------------------------------------------------" << endl;
        aus << "[AFLOW] POSCAR EXPLICIT MODE without START/STOP (default)" << endl;
        aus << "[VASP_POSCAR_MODE_EXPLICIT]" << endl;
        aus << "[VASP_POSCAR_FILE]POSCAR of the structure example" << endl;
        aus << "[VASP_POSCAR_FILE]-98.5397" << endl;
        aus << "[VASP_POSCAR_FILE]   4.18890 0.00000 0.00000" << endl;
        aus << "[VASP_POSCAR_FILE]  -2.09445 3.62769 0.00000" << endl;
        aus << "[VASP_POSCAR_FILE]   0.00000 0.00000 5.12300" << endl;
        aus << "[VASP_POSCAR_FILE]2 4" << endl;
        aus << "[VASP_POSCAR_FILE]Direct" << endl;
        aus << "[VASP_POSCAR_FILE]0.33333 0.66666 0.25000 Au" << endl;
        aus << "[VASP_POSCAR_FILE]0.66666 0.33333 0.75000 Au" << endl;
        aus << "[VASP_POSCAR_FILE]0.00000 0.00000 0.00000 Ti" << endl;
        aus << "[VASP_POSCAR_FILE]0.00000 0.00000 0.50000 Ti" << endl;
        aus << "[VASP_POSCAR_FILE]0.33333 0.66666 0.75000 Ti" << endl;
        aus << "[VASP_POSCAR_FILE]0.66666 0.33333 0.25000 Ti" << endl;
        aus << "[AFLOW]" << endl;
        aus << "----------------------------------------------------------------------------------------------------" << endl;
        aus << "[AFLOW] POSCAR EXPLICIT MODE with START/STOP" << endl;
        aus << "[VASP_POSCAR_MODE_EXPLICIT]" << endl;
        aus << "[VASP_POSCAR_MODE_EXPLICIT]START" << endl;
        aus << "POSCAR of the structure example with START/STOP" << endl;
        aus << "-98.5397" << endl;
        aus << "   4.18890 0.00000 0.00000" << endl;
        aus << "  -2.09445 3.62769 0.00000" << endl;
        aus << "   0.00000 0.00000 5.12300" << endl;
        aus << "2 4" << endl;
        aus << "Direct" << endl;
        aus << "0.33333 0.66666 0.25000 Au" << endl;
        aus << "0.66666 0.33333 0.75000 Au" << endl;
        aus << "0.00000 0.00000 0.00000 Ti" << endl;
        aus << "0.00000 0.00000 0.50000 Ti" << endl;
        aus << "0.33333 0.66666 0.75000 Ti" << endl;
        aus << "0.66666 0.33333 0.25000 Ti" << endl;
        aus << "[VASP_POSCAR_MODE_EXPLICIT]STOP" << endl;
        aus << "[AFLOW]" << endl;
        aus << "----------------------------------------------------------------------------------------------------" << endl;
        aus << "EEEEE  [VASP_POSCAR_MODE_EXPLICIT] Note " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aus << "EEEEE  [VASP_POSCAR_MODE_EXPLICIT]START must be present and no [VASP_POSCAR_FILE]" << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aus << "EEEEE  [VASP_POSCAR_MODE_EXPLICIT]STOP  must be present and no [VASP_POSCAR_FILE]" << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aus << "EEEEE  or [VASP_POSCAR_FILE] present and NO START/STOP" << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
    }
    // EXTERNAL **************************************************
    if(Krun && vflags.KBIN_VASP_POSCAR_MODE.flag("EXTERNAL")) {  // [VASP_POSCAR_MODE_EXTERNAL] construction
      string file;
      aus << "00000  MESSAGE POSCAR  generation EXTERNAL file from " << _AFLOWIN_ << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(vflags.KBIN_VASP_POSCAR_FILE.flag("COMMAND") && vflags.KBIN_VASP_POSCAR_FILE.flag("FILE")) {
        aus << "EEEEE   [VASP_POSCAR_MODE]FILE=  and  [VASP_POSCAR_MODE]COMMAND=  can not be used together " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
      if(!vflags.KBIN_VASP_POSCAR_FILE.flag("COMMAND") && (vflags.KBIN_VASP_POSCAR_FILE.flag("FILE") || !vflags.KBIN_VASP_POSCAR_FILE.flag("FILE"))) {
        if(vflags.KBIN_VASP_POSCAR_FILE.flag("FILE")) {
          file=aurostd::substring2string(AflowIn,"[VASP_POSCAR_FILE]FILE=",TRUE);
          aus << "00000  MESSAGE POSCAR  generation from file=" << file << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        } else {
          file=DEFAULT_VASP_EXTERNAL_POSCAR;
          aus << "00000  MESSAGE POSCAR  generation from DEFAULT file=" << DEFAULT_VASP_EXTERNAL_POSCAR << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if(!aurostd::FileExist(file)) {
          aus << "EEEEE  ERROR POSCAR file=" << file << " does not exist! " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        if(aurostd::FileEmpty(file)) {
          aus << "EEEEE  ERROR POSCAR file=" << file << " is empty! " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        xvasp.POSCAR << aurostd::file2string(file);
        xvasp.str=xstructure(xvasp.POSCAR,IOVASP_POSCAR);  // load structure
      }
      if(vflags.KBIN_VASP_POSCAR_FILE.flag("COMMAND") && !vflags.KBIN_VASP_POSCAR_FILE.flag("FILE")) {
        file=aurostd::substring2string(AflowIn,"[VASP_POSCAR_FILE]COMMAND=",FALSE);
        aus << "00000  MESSAGE POSCAR  generation from command= '" << file << "' " << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        file=file+" > ./_aflow_POSCAR."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";    // create temp //CO20200502 - threadID
        aurostd::execute(file);                           // create temp
        file="./_aflow_POSCAR."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";            // file name //CO20200502 - threadID
        if(!aurostd::FileExist(file)) {  // could not write (directory protected)
          aus << "EEEEE  ERROR POSCAR file=" << file << " does not exist! " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        if(aurostd::FileEmpty(file)) {  // contains nothing good
          aus << "EEEEE  ERROR POSCAR file=" << file << " is empty! " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        xvasp.POSCAR << aurostd::file2string(file);       // load POSCAR
        xvasp.str=xstructure(xvasp.POSCAR,IOVASP_POSCAR);              // load structure
        aurostd::RemoveFile("./_aflow_POSCAR."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp");     // remove temp  //CO20200502 - threadID
      }
    }
    // POSCAR DONE **************************************************
    xvasp.POSCAR_orig << xvasp.POSCAR.str();
    xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
    xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",FALSE);

    // must modify POSCAR before calculating everything else
    if(!kflags.KBIN_POCC){  //CO20180420 - VERY SLOW AND UNNECESSARY - let's do these mods when we read-in/run in the ARUNS
      if(Krun) Krun=(Krun && KBIN::VASP_Modify_POSCAR(xvasp,AflowIn,FileMESSAGE,aflags,vflags));
    }

    // CHECK for negative determinant
    if(det(xvasp.str.scale*xvasp.str.lattice)<0.0) {
      aus << "EEEEE  POSCAR ERROR: the triple product of the basis vectors is negative                      " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aus << "EEEEE  POSCAR ERROR: exchange two basis vectors and adjust the atomic positions accordingly   " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      Krun=FALSE;
      return Krun;
    }

    // some useful LDEBUG
    if(Krun && LDEBUG) {
      // STRUCTURE IS GENERATED
      FileMESSAGE <<  endl;
      FileMESSAGE <<  "******** STRUCTURE IN CARTESIAN *****************************" << endl;
      xvasp.str.SetCoordinates(_COORDS_CARTESIAN_);
      FileMESSAGE <<  xvasp.str << endl;
      FileMESSAGE <<  "******** STRUCTURE IN FRACTIONAL ****************************" << endl;
      xvasp.str.SetCoordinates(_COORDS_FRACTIONAL_);
      FileMESSAGE <<  xvasp.str << endl;
      FileMESSAGE <<  "*************************************************************" << endl;
      // xvasp.str.write_klattice_flag=TRUE;
      FileMESSAGE <<  "SCALE" << endl;
      FileMESSAGE <<  xvasp.str.scale << endl;
      FileMESSAGE <<  "DIRECT LATTICE (with scale)" << endl;
      FileMESSAGE <<  xvasp.str.scale*(xvasp.str.lattice) << endl;
      FileMESSAGE <<  "RECIPROCAL LATTICE" << endl;
      FileMESSAGE <<  (xvasp.str.klattice) << endl;
      FileMESSAGE <<  "ORTOGONALITY (a*b')/2pi=I" << endl;
      FileMESSAGE <<  (xvasp.str.scale*(xvasp.str.lattice))*trasp(xvasp.str.klattice)/(2.0*pi) << endl;
      FileMESSAGE <<  "*************************************************************" << endl;
    }
    // done produced and modified
    return Krun;
  };  // KBIN::VASP_Produce_POSCAR
}

namespace KBIN {
  bool VASP_Produce_POSCAR(_xvasp& xvasp) {        // AFLOW_FUNCTION_IMPLEMENTATION
    bool Krun=TRUE;
    xvasp.POSCAR.str(std::string());xvasp.POSCAR.clear();
    xvasp.POSCAR_orig.str(std::string());xvasp.POSCAR_orig.clear();
    xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",FALSE);
    xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",FALSE);
    xvasp.POSCAR << xvasp.str;
    // POSCAR done
    xvasp.POSCAR_orig << xvasp.POSCAR.str();
    xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
    return Krun;
  };  // KBIN::VASP_Produce_POSCAR
}

namespace KBIN {
  bool VASP_Modify_POSCAR(_xvasp& xvasp,const string& AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_vflags &vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    ostringstream aus;
    bool Krun=TRUE;
    if(xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated")==FALSE) {
      aus << "EEEEE  KBIN::VASP_Modify_POSCAR: can`t modify POSCAR if it does not exist ! " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      Krun=FALSE;
      return Krun;
    }
    // return Krun;
    // xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",FALSE);
    //LDEBUG=TRUE;

    // POSCAR must be modified before doing the KPOINTS
    //ME20181121 - Outsourced to be used by other functions
    Krun = VASP_Convert_Unit_Cell(xvasp, vflags, aflags, FileMESSAGE, aus);

    if(LDEBUG) cerr << "vflags.KBIN_VASP_POSCAR_FILE_VOLUME.flag(\"EQUAL_EQUAL\")=" << vflags.KBIN_VASP_POSCAR_FILE_VOLUME.flag("EQUAL_EQUAL") << endl;
    if(Krun && vflags.KBIN_VASP_POSCAR_MODE.flag("IMPLICIT") && vflags.KBIN_VASP_POSCAR_FILE_VOLUME.flag("EQUAL_EQUAL")) {  // [VASP_POSCAR_FILE]VOLUME=                       /*************** POSCAR **************/
      double factor=aurostd::substring2utype<double>(AflowIn,"[VASP_POSCAR_FILE]VOLUME=",FALSE);
      aus << "00000  MESSAGE POSCAR  IMPLICIT Volume = " << factor << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(factor<=0.0) {
        aus << "EEEEE  KBIN::VASP_Modify_POSCAR Volume can not be <=0 (factor=" << factor << ")" << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
      aus << "00000  MESSAGE POSCAR  IMPLITIC Old Volume= " << xvasp.str.Volume() << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.SetVolume(factor);
      aus << "00000  MESSAGE POSCAR  IMPLITIC New Volume= " << xvasp.str.Volume() << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.title=xvasp.str.title+" [Forced Volume =]";
      xvasp.POSCAR.str(std::string());xvasp.POSCAR.clear();
      xvasp.POSCAR << xvasp.str;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag(\"EQUAL_EQUAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag("EQUAL_EQUAL") << endl;
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag("EQUAL_EQUAL")) {  // [VASP_FORCE_OPTION]VOLUME=                                                    /*************** POSCAR **************/
      double factor1=aurostd::substring2utype<double>(AflowIn,"[VASP_FORCE_OPTION]VOLUME=",FALSE);
      double factor=vflags.KBIN_VASP_FORCE_OPTION_VOLUME.getattachedutype<double>("EQUAL_EQUAL");
      aus << "00000  MESSAGE POSCAR  FORCE Volume = " << factor << " (factor1=" << factor1 << ")  " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(factor<=0.0) {
        aus << "EEEEE  KBIN::VASP_Modify_POSCAR Volume can not be <=0 (factor=" << factor << ")" << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
      aus << "00000  MESSAGE POSCAR  FORCE Old Volume= " << xvasp.str.Volume() << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.SetVolume(factor);
      aus << "00000  MESSAGE POSCAR  FORCE New Volume= " << xvasp.str.Volume() << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.title=xvasp.str.title+" [Forced Volume =]";
      xvasp.POSCAR.str(std::string());xvasp.POSCAR.clear();
      xvasp.POSCAR << xvasp.str;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_POSCAR_FILE_VOLUME.flag(\"MULTIPLY_EQUAL\")=" << vflags.KBIN_VASP_POSCAR_FILE_VOLUME.flag("MULTIPLY_EQUAL") << endl;
    if(Krun && vflags.KBIN_VASP_POSCAR_MODE.flag("IMPLICIT") && vflags.KBIN_VASP_POSCAR_FILE_VOLUME.flag("MULTIPLY_EQUAL")) {  // [VASP_POSCAR_FILE]VOLUME*=             /*************** POSCAR **************/
      double factor=aurostd::substring2utype<double>(AflowIn,"[VASP_POSCAR_FILE]VOLUME*=",FALSE);
      //     double factor=aurostd::string2utype<double>(vflags.KBIN_VASP_POSCAR_FILE_VOLUME.getattachedscheme("MULTIPLY_EQUAL"));
      //      cerr << "CORMAC MULTIPLY_EQUAL=" << factor << endl;
      aus << "00000  MESSAGE POSCAR  IMPLICIT Volume *= " << factor << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(factor<=0.0) {
        aus << "EEEEE  KBIN::VASP_Modify_POSCAR Volume can not be <=0 (factor=" << factor << ")" << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
      aus << "00000  MESSAGE POSCAR  IMPLITIC Old Volume= " << xvasp.str.Volume() << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.SetVolume(xvasp.str.Volume()*factor);
      aus << "00000  MESSAGE POSCAR  IMPLITIC New Volume= " << xvasp.str.Volume() << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.title=xvasp.str.title+" [Forced Volume *=]";
      xvasp.POSCAR.str(std::string());xvasp.POSCAR.clear();
      xvasp.POSCAR << xvasp.str;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag(\"MULTIPLY_EQUAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag("MULTIPLY_EQUAL") << endl;
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag("MULTIPLY_EQUAL")) {  // [VASP_FORCE_OPTION]VOLUME*=                                                    /*************** POSCAR **************/
      double factor1=aurostd::substring2utype<double>(AflowIn,"[VASP_FORCE_OPTION]VOLUME*=",FALSE);
      double factor=vflags.KBIN_VASP_FORCE_OPTION_VOLUME.getattachedutype<double>("MULTIPLY_EQUAL");
      aus << "00000  MESSAGE POSCAR  FORCE Volume *= " << factor << " (factor1=" << factor1 << ")  " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(factor<=0.0) {
        aus << "EEEEE  KBIN::VASP_Modify_POSCAR Volume can not be <=0 (factor=" << factor << ")" << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
      aus << "00000  MESSAGE POSCAR  FORCE Old Volume= " << xvasp.str.Volume() << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.SetVolume(xvasp.str.Volume()*factor);
      aus << "00000  MESSAGE POSCAR  FORCE New Volume= " << xvasp.str.Volume() << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.title=xvasp.str.title+" [Forced Volume *=]";
      xvasp.POSCAR.str(std::string());xvasp.POSCAR.clear();
      xvasp.POSCAR << xvasp.str;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_POSCAR_FILE_VOLUME.flag(PLUS_EQUAL)=" << vflags.KBIN_VASP_POSCAR_FILE_VOLUME.flag("PLUS_EQUAL") << endl;
    if(Krun && vflags.KBIN_VASP_POSCAR_MODE.flag("IMPLICIT") && vflags.KBIN_VASP_POSCAR_FILE_VOLUME.flag("PLUS_EQUAL")) {  // [VASP_POSCAR_FILE]VOLUME+=               /*************** POSCAR **************/
      double factor=aurostd::substring2utype<double>(AflowIn,"[VASP_POSCAR_FILE]VOLUME+=",FALSE);
      aus << "00000  MESSAGE POSCAR  IMPLICIT Volume += " << factor << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      aus << "00000  MESSAGE POSCAR  IMPLITIC Old Volume= " << xvasp.str.Volume() << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.SetVolume(xvasp.str.Volume()+factor);
      aus << "00000  MESSAGE POSCAR  IMPLITIC New Volume= " << xvasp.str.Volume() << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.title=xvasp.str.title+" [Forced Volume +=]";
      xvasp.POSCAR.str(std::string());xvasp.POSCAR.clear();
      xvasp.POSCAR << xvasp.str;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag(PLUS_EQUAL)=" << vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag("PLUS_EQUAL") << endl;
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_VOLUME.flag("PLUS_EQUAL")) {  // [VASP_FORCE_OPTION]VOLUME+=                                                    /*************** POSCAR **************/
      double factor1=aurostd::substring2utype<double>(AflowIn,"[VASP_FORCE_OPTION]VOLUME+=",FALSE);
      double factor=vflags.KBIN_VASP_FORCE_OPTION_VOLUME.getattachedutype<double>("PLUS_EQUAL");
      aus << "00000  MESSAGE POSCAR  FORCE Volume += " << factor << " (factor1=" << factor1 << ")  " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      aus << "00000  MESSAGE POSCAR  FORCE Old Volume= " << xvasp.str.Volume() << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.SetVolume(xvasp.str.Volume()+factor);
      aus << "00000  MESSAGE POSCAR  FORCE New Volume= " << xvasp.str.Volume() << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.title=xvasp.str.title+" [Forced Volume +=]";
      xvasp.POSCAR.str(std::string());xvasp.POSCAR.clear();
      xvasp.POSCAR << xvasp.str;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
    }

    // POSCAR done
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {
      if(0) {
        aus << "00000  MESSAGE-OPTION  XXXXX" << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    }
    // ------------------------------------
    // end
    // xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
    return Krun;
  }; // KBIN::VASP_Produce_POSCAR
}

namespace KBIN {
  //ME20190220
  // Convert POSCAR into the format that is consistent with the binary version   
  //ME20200114 - Throw a warning when the version could not be determined and
  // leave as is. Aborting is not desirable for instances where VASP does not
  // need to be run (e.g. post-processing).
  void convertPOSCARFormat(_xvasp& xvasp, const _kflags& kflags) {
    string vaspVersion = getVASPVersionNumber(kflags.KBIN_BIN); //CO20200610
    if (vaspVersion.empty()) {
      string function = "KBIN::convertPOSCARFormat()";
      stringstream message;
      message << "Could not determine VASP version " << kflags.KBIN_BIN << "."
        << " POSCAR may have the wrong format.";
      pflow::logger(_AFLOW_FILE_NAME_, function, message, std::cout, _LOGGER_WARNING_);
      return;
    }
    int vaspv = vaspVersion[0] - '0';
    if ((vaspv == 4) && !xvasp.str.is_vasp4_poscar_format) {
      xvasp.str.is_vasp4_poscar_format = true;
      xvasp.str.is_vasp5_poscar_format = false;
      xvasp.POSCAR.str(std::string());
      xvasp.POSCAR.clear();
      xvasp.POSCAR << xvasp.str;
    } else if ((vaspv == 5) && !xvasp.str.is_vasp5_poscar_format) {
      xvasp.str.is_vasp4_poscar_format = false;
      xvasp.str.is_vasp5_poscar_format = true;
      xvasp.POSCAR.str(std::string());
      xvasp.POSCAR.clear();
      xvasp.POSCAR << xvasp.str;
    }
  }

  // CONVERT_UNIT_CELL STUFF
  bool VASP_Convert_Unit_Cell(_xvasp& xvasp, _vflags& vflags, _aflags& aflags, ofstream& FileMESSAGE, ostringstream& aus) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    bool Krun = true;

    if(LDEBUG) aus << "00000  MESSAGE-OPTION  KBIN::VASP_Convert_Unit_Cell: BEGIN " << endl;  //SC20200410
    if(LDEBUG) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);                      //SC20200410

    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.isentry) {        /*************** POSCAR **************/
      //    aus << "00000  MESSAGE POSCAR  [VASP_FORCE_OPTION]CONVERT_UNIT_CELL=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.content_string << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]CONVERT_UNIT_CELL=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.content_string << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"PRESERVE\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("PRESERVE") << endl;
    // POSCAR must be modified before doing the KPOINTS
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("PRESERVE")) {  // [VASP_FORCE_OPTION]CONVERT_UNIT_CELL_PRESERVE construction /*************** POSCAR **************/
      aus << "00000  MESSAGE POSCAR  PRESERVE Unit Cell " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"STANDARD_PRIMITIVE\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_PRIMITIVE") << endl;
    // POSCAR must be modified before doing the KPOINTS
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_PRIMITIVE")) {  // [VASP_FORCE_OPTION]CONVERT_UNIT_CELL_STANDARD_PRIMITIVE construction /*************** POSCAR **************/
      aus << "00000  MESSAGE POSCAR  STANDARD_PRIMITIVE Unit Cell " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.Standard_Primitive_UnitCellForm();
      xvasp.str.sortAtomsEquivalent();  //CO20190216 - critical for prospective APL calculations, better to sort before relaxations and not have to sort again in the future
      //CO START
      //CO, fix issue with iatoms tag becoming atom names
      bool write_inequivalent_flag=xvasp.str.write_inequivalent_flag;
      //CO END
      string bravais_lattice_type=xvasp.str.bravais_lattice_type,bravais_lattice_variation_type=xvasp.str.bravais_lattice_variation_type,pearson_symbol=xvasp.str.pearson_symbol;
      xvasp.str.title=xvasp.str.title+" [Standard_Primitive Unit Cell Form]";
      xvasp.POSCAR.str(std::string());xvasp.POSCAR.clear();
      //CO START
      //CO, fix if write_inequivalent_flag is present
      xvasp.str.write_inequivalent_flag=FALSE;
      //CO, fix if write_inequivalent_flag is present
      xvasp.POSCAR << xvasp.str;
      xvasp.str.clear();xvasp.POSCAR >> xvasp.str;  //CO, this is important, clear all symmetry stuff as the whole lattice has changed //DX20191220 - uppercase to lowercase clear
      //CO add these flags to prevent recalculation and wasted effort
      xvasp.str.Standard_Lattice_calculated=TRUE;
      xvasp.str.Standard_Lattice_primitive=TRUE;
      //CO add these flags to prevent recalculation and wasted effort
      //CO, fix issue with iatoms tag becoming atom names
      xvasp.str.write_inequivalent_flag=write_inequivalent_flag;
      //CO END
      xvasp.str.bravais_lattice_type=bravais_lattice_type;xvasp.str.bravais_lattice_variation_type=bravais_lattice_variation_type;xvasp.str.pearson_symbol=pearson_symbol;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
      aus << "00000  MESSAGE POSCAR  STANDARD_PRIMITIVE Unit Cell Lattice = ["+xvasp.str.bravais_lattice_type << "," << xvasp.str.bravais_lattice_variation_type << "," << xvasp.str.pearson_symbol << "]" << "  " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"STANDARD_CONVENTIONAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_CONVENTIONAL") << endl;
    // POSCAR must be modified before doing the KPOINTS
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("STANDARD_CONVENTIONAL")) {  // [VASP_FORCE_OPTION]CONVERT_UNIT_CELL_STANDARD_CONVENTIONAL construction /*************** POSCAR **************/
      aus << "00000  MESSAGE POSCAR  STANDARD_CONVENTIONAL Unit Cell " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.Standard_Conventional_UnitCellForm();
      xvasp.str.sortAtomsEquivalent();  //CO20190216 - critical for prospective APL calculations, better to sort before relaxations and not have to sort again in the future
      //CO START
      //CO, fix issue with iatoms tag becoming atom names
      bool write_inequivalent_flag=xvasp.str.write_inequivalent_flag;
      //CO END
      string bravais_lattice_type=xvasp.str.bravais_lattice_type,bravais_lattice_variation_type=xvasp.str.bravais_lattice_variation_type,pearson_symbol=xvasp.str.pearson_symbol;
      xvasp.str.title=xvasp.str.title+" [Standard_Conventional Unit Cell Form]";
      xvasp.POSCAR.str(std::string());xvasp.POSCAR.clear();
      //CO START
      //CO, fix if write_inequivalent_flag is present
      xvasp.str.write_inequivalent_flag=FALSE;
      //CO, fix if write_inequivalent_flag is present
      xvasp.POSCAR << xvasp.str;
      xvasp.str.clear();xvasp.POSCAR >> xvasp.str;  //CO, this is important, clear all symmetry stuff as the whole lattice has change //DX20191220 - uppercase to lowercase clear
      //CO add these flags to prevent recalculation and wasted effort
      xvasp.str.Standard_Lattice_calculated=TRUE;
      xvasp.str.Standard_Lattice_conventional=TRUE;
      //CO add these flags to prevent recalculation and wasted effort
      //CO, fix issue with iatoms tag becoming atom names
      xvasp.str.write_inequivalent_flag=write_inequivalent_flag;
      //CO END
      xvasp.str.bravais_lattice_type=bravais_lattice_type;xvasp.str.bravais_lattice_variation_type=bravais_lattice_variation_type;xvasp.str.pearson_symbol=pearson_symbol;
      // xvasp.str.clear();xvasp.POSCAR >> xvasp.str; //DX20191220 - uppercase to lowercase clear
      // cout << xvasp.str << endl;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
      aus << "00000  MESSAGE POSCAR  STANDARD_CONVENTIONAL Unit Cell Lattice = ["+xvasp.str.bravais_lattice_type << "," << xvasp.str.bravais_lattice_variation_type << "," << xvasp.str.pearson_symbol << "]" << "  " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl; //CO
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET); //CO
      // cerr << det(xvasp.str.lattice) << endl;
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"NIGGLI\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("NIGGLI") << endl;
    // POSCAR must be modified before doing the KPOINTS
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("NIGGLI")) {  // [VASP_FORCE_OPTION]CONVERT_UNIT_CELL_NIGGLI construction                     /*************** POSCAR **************/
      aus << "00000  MESSAGE POSCAR  NIGGLI Unit Cell Reduction " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.NiggliUnitCellForm();
      xvasp.str.sortAtomsEquivalent();  //CO20190216 - critical for prospective APL calculations, better to sort before relaxations and not have to sort again in the future
      xvasp.str.title=xvasp.str.title+" [Niggli Unit Cell Form]";
      xvasp.POSCAR.str(std::string());xvasp.POSCAR.clear();
      xvasp.POSCAR << xvasp.str;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
      // cerr << det(xvasp.str.lattice) << endl;
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(MINKOWSKI)=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("MINKOWSKI") << endl;
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("MINKOWSKI")) {  // [VASP_FORCE_OPTION]CONVERT_UNIT_CELL_MINKOWSKI construction               /*************** POSCAR **************/
      aus << "00000  MESSAGE POSCAR  MINKOWSKI Basis Reduction " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.MinkowskiBasisReduction();
      xvasp.str.sortAtomsEquivalent();  //CO20190216 - critical for prospective APL calculations, better to sort before relaxations and not have to sort again in the future
      xvasp.str.title=xvasp.str.title+" [Minkowski Basis Reduction]";
      xvasp.POSCAR.str(std::string());xvasp.POSCAR.clear();
      xvasp.POSCAR << xvasp.str;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
      // cerr << det(xvasp.str.lattice) << endl;
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"INCELL\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("INCELL") << endl;
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("INCELL")) {  // [VASP_FORCE_OPTION]CONVERT_UNIT_CELL_INCELL construction                     /*************** POSCAR **************/
      aus << "00000  MESSAGE POSCAR  INCELL Unit Cell Basis " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.BringInCell();
      xvasp.str.title=xvasp.str.title+" [Bring In Cell Basis]";
      xvasp.POSCAR.str(std::string());xvasp.POSCAR.clear();
      xvasp.POSCAR << xvasp.str;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"COMPACT\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("COMPACT") << endl;
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("COMPACT")) {  // [VASP_FORCE_OPTION]CONVERT_UNIT_CELL_COMPACT construction                   /*************** POSCAR **************/
      aus << "00000  MESSAGE POSCAR  COMPACT Unit Cell Basis " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.BringInCompact();
      xvasp.str.title=xvasp.str.title+" [Bring In Compact Basis]";
      xvasp.POSCAR.str(std::string());xvasp.POSCAR.clear();
      xvasp.POSCAR << xvasp.str;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"WIGNERSEITZ\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("WIGNERSEITZ") << endl;
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("WIGNERSEITZ")) {  // [VASP_FORCE_OPTION]CONVERT_UNIT_CELL_WIGNERSEITZ construction           /*************** POSCAR **************/
      aus << "00000  MESSAGE POSCAR  WIGNERSEITZ Unit Cell Basis " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.BringInWignerSeitz();
      xvasp.str.title=xvasp.str.title+" [WignerSeitz Basis]";
      xvasp.POSCAR.str(std::string());xvasp.POSCAR.clear();
      xvasp.POSCAR << xvasp.str;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"CARTESIAN\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("CARTESIAN") << endl;
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("CARTESIAN")) {  // [VASP_FORCE_OPTION]CONVERT_UNIT_CELL_CARTESIAN construction               /*************** POSCAR **************/
      aus << "00000  MESSAGE POSCAR  CARTESIAN Basis Coordinates" << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.SetCoordinates(_COORDS_CARTESIAN_);
      //[CO20190216 - seems mistake]xvasp.str.title=xvasp.str.title+" [WignerSeitz Basis]";
      xvasp.POSCAR.str(std::string());xvasp.POSCAR.clear();
      xvasp.POSCAR << xvasp.str;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"FRACTIONAL\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("FRACTIONAL") << endl;
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("FRACTIONAL")) {  // [VASP_FORCE_OPTION]CONVERT_UNIT_CELL_FRACTIONAL construction            /*************** POSCAR **************/
      aus << "00000  MESSAGE POSCAR  FRACTIONAL Basis Coordinate" << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.SetCoordinates(_COORDS_FRACTIONAL_);
      xvasp.POSCAR.str(std::string());xvasp.POSCAR.clear();
      xvasp.POSCAR << xvasp.str;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
    }

    if(LDEBUG) cerr << "vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag(\"DIRECT\")=" << vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("DIRECT") << endl;
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("DIRECT")) {  // [VASP_FORCE_OPTION]CONVERT_UNIT_CELL_DIRECT construction                    /*************** POSCAR **************/
      aus << "00000  MESSAGE POSCAR  DIRECT Basis Coordinate" << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.SetCoordinates(_COORDS_FRACTIONAL_);
      xvasp.POSCAR.str(std::string());xvasp.POSCAR.clear();
      xvasp.POSCAR << xvasp.str;
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated",TRUE);
      xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
    }

    if(LDEBUG) aus << "00000  MESSAGE-OPTION  KBIN::VASP_Convert_Unit_Cell: END " << endl;  //SC20200410
    if(LDEBUG) aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);                    //SC20200410

    return Krun;
  }
}

namespace KBIN {
  bool VASP_Reread_POSCAR(_xvasp& xvasp,ofstream &FileMESSAGE,_aflags &aflags) { // AFLOW_FUNCTION_IMPLEMENTATION
    ostringstream aus;
    bool Krun=TRUE;
    if(!aurostd::FileExist(xvasp.Directory+"/POSCAR")) {
      aus << "EEEEE  KBIN::VASP_Reread_POSCAR: POSCAR not present in directory: " << xvasp.Directory << " - "  << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      Krun=FALSE;
      return Krun;
    }
    xvasp.POSCAR_orig.str(std::string()); xvasp.POSCAR_orig << xvasp.POSCAR.str();
    xvasp.POSCAR.str(std::string()); xvasp.POSCAR << aurostd::file2string(xvasp.Directory+"/POSCAR"); // DID REREAD
    xvasp.aopts.flag("FLAG::XVASP_POSCAR_changed",TRUE);
    return Krun;
  }
}


// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// KPOINTS
namespace KBIN {
  bool VASP_Produce_KPOINTS(_xvasp& xvasp,const string& AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_vflags &vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    if(AflowIn.length()==0) {throw aurostd::xerror(_AFLOW_FILE_NAME_,"KBIN::VASP_Produce_KPOINTS():","empty AflowIn",_FILE_CORRUPT_);} //CO20200624
    if(!kflags.AFLOW_MODE_VASP) {cerr << XPID << "KBIN::VASP_Produce_KPOINTS: should kflags.AFLOW_MODE_VASP be set ??" << endl;}
    ostringstream aus;
    bool Krun=TRUE;
    bool IMPLICIT=FALSE;
    xvasp.KPOINTS.str(std::string());
    xvasp.str.kpoints_k1=0;xvasp.str.kpoints_k2=0;xvasp.str.kpoints_k3=0;        // RESET KPOINTS
    xvasp.str.kpoints_s1=0.0;xvasp.str.kpoints_s2=0.0;xvasp.str.kpoints_s3=0.0;  // RESET SHIFTS

    xvasp.str.kpoints_kmax=0;
    xvasp.str.kpoints_kscheme.clear();
    xvasp.aopts.flag("FLAG::XVASP_KPOINTS_generated",FALSE);
    xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",FALSE);

    // IMPLICIT or EXPLICIT or EXTERNAL for KPOINTS
    Krun=(Krun && (vflags.KBIN_VASP_KPOINTS_MODE.flag("IMPLICIT") ||
          vflags.KBIN_VASP_KPOINTS_MODE.flag("EXPLICIT") ||
          vflags.KBIN_VASP_KPOINTS_MODE.flag("EXTERNAL")));
    if(!Krun) {
      aurostd::StringstreamClean(aus);
      aus << "EEEEE  [VASP_KPOINTS_MODE_IMPLICIT] or [VASP_KPOINTS_MODE_EXPLICIT] or [VASP_KPOINTS_MODE_EXTERNAL] must be specified " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl; //ME20181113
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      Krun=FALSE;
      return Krun;
    }
    if(Krun && xvasp.aopts.flag("FLAG::XVASP_POSCAR_generated")==FALSE) {
      aus      << "EEEEE  INCAR   generation: POSCAR must be called before " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      Krun=FALSE;
      return Krun;
    }
    // EXPLICIT **************************************************
    // kpoints explicit through string xvasp.KPOINTS
    if(Krun && vflags.KBIN_VASP_KPOINTS_MODE.flag("EXPLICIT")) {  // [VASP_KPOINTS_MODE_EXPLICIT] construction
      IMPLICIT=FALSE;
      if(vflags.KBIN_VASP_KPOINTS_FILE.flag("KEYWORD") && !vflags.KBIN_VASP_KPOINTS_MODE.flag("EXPLICIT_START_STOP")) {
        aus << "00000  MESSAGE KPOINTS generation EXPLICIT file from " << _AFLOWIN_ << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        xvasp.KPOINTS << vflags.KBIN_VASP_KPOINTS_EXPLICIT.str(); //ME20181113
        // [OBSOLETE]      aurostd::ExtractToStringstreamEXPLICIT(FileAFLOWIN,xvasp.KPOINTS,"[VASP_KPOINTS_FILE]");
        // [ME20181226 OBSOLETE]      aurostd::ExtractToStringstreamEXPLICIT(AflowIn,xvasp.KPOINTS,"[VASP_KPOINTS_FILE]");
        KBIN::XVASP_string2numbers(xvasp);
      } else if(!vflags.KBIN_VASP_KPOINTS_FILE.flag("KEYWORD") && vflags.KBIN_VASP_KPOINTS_MODE.flag("EXPLICIT_START_STOP")) {
        string START="[VASP_KPOINTS_MODE_EXPLICIT]START";
        string STOP="[VASP_KPOINTS_MODE_EXPLICIT]STOP";
        aus << "00000  MESSAGE KPOINTS generation EXPLICIT file from " << _AFLOWIN_ << " with START/STOP  " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        xvasp.KPOINTS << vflags.KBIN_VASP_KPOINTS_EXPLICIT_START_STOP.str(); //ME20181113
        // [ME20181226 OBSOLETE] if(aurostd::substring2bool(AflowIn,START) && aurostd::substring2bool(AflowIn,STOP))
        // [OBSOLETE]	aurostd::ExtractToStringstreamEXPLICIT(FileAFLOWIN,xvasp.KPOINTS,START,STOP);
        // [ME20181226 OBSOLETE] aurostd::ExtractToStringstreamEXPLICIT(AflowIn,xvasp.KPOINTS,START,STOP);
        KBIN::XVASP_string2numbers(xvasp);
      } else {
        aus << "EEEEE  [VASP_KPOINTS_MODE_EXPLICIT] do not confuse aflow !!" << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aus << "EEEEE  [VASP_KPOINTS_MODE_EXPLICIT] Possible modes " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aus << "----------------------------------------------------------------------------------------------------" << endl;
        aus << "[AFLOW] KPOINTS EXPLICIT MODE without START/STOP (default)" << endl;
        aus << "[VASP_KPOINTS_MODE_EXPLICIT]" << endl;
        aus << "[VASP_KPOINTS_FILE]KPOINTS of the structure example with START/STOP" << endl;
        aus << "[VASP_KPOINTS_FILE]0" << endl;
        aus << "[VASP_KPOINTS_FILE]Monkhorst-Pack" << endl;
        aus << "[VASP_KPOINTS_FILE]7 7 5" << endl;
        aus << "[VASP_KPOINTS_FILE]0 0 0" << endl;
        aus << "[AFLOW]" << endl;
        aus << "----------------------------------------------------------------------------------------------------" << endl;
        aus << "[AFLOW] KPOINTS EXPLICIT MODE with START/STOP" << endl;
        aus << "[VASP_KPOINTS_MODE_EXPLICIT]" << endl;
        aus << "[VASP_KPOINTS_MODE_EXPLICIT]START" << endl;
        aus << "KPOINTS of the structure example with START/STOP" << endl;
        aus << "0" << endl;
        aus << "Monkhorst-Pack" << endl;
        aus << "7 7 5" << endl;
        aus << "0 0 0" << endl;
        aus << "[VASP_KPOINTS_MODE_EXPLICIT]STOP" << endl;
        aus << "[AFLOW]" << endl;
        aus << "----------------------------------------------------------------------------------------------------" << endl;
        aus << "EEEEE  [VASP_KPOINTS_MODE_EXPLICIT] Note " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aus << "EEEEE  [VASP_KPOINTS_MODE_EXPLICIT]START must be present and no [VASP_KPOINTS_FILE]" << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aus << "EEEEE  [VASP_KPOINTS_MODE_EXPLICIT]STOP  must be present and no [VASP_KPOINTS_FILE]" << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aus << "EEEEE  or [VASP_KPOINTS_FILE] present and NO START/STOP" << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
    }
    // IMPLICIT **************************************************
    // kpoints implicit through string xvasp.KPOINTS
    if(Krun && vflags.KBIN_VASP_KPOINTS_MODE.flag("IMPLICIT")) // IT MIGHT NOT CONTAIN OPTION SO IT IS ALL DEFAULT && vflags.KBIN_VASP_KPOINTS_FILE)  // [VASP_KPOINTS_MODE_IMPLICIT] construction
    { //CO20200106 - patching for auto-indenting
      IMPLICIT=TRUE;
      aus << "00000  MESSAGE KPOINTS generation IMPLICIT file from " << _AFLOWIN_ << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      string stringKPPRA;
      int NK=1;
      xvasp.str.kpoints_k1=1;xvasp.str.kpoints_k2=1;xvasp.str.kpoints_k3=1;
      xvasp.str.kpoints_mode=0;
      xvasp.str.kpoints_kscheme.clear();
      bool LOCAL_KBIN_VASP_KPOINTS_KMODE_isentry=vflags.KBIN_VASP_KPOINTS_KMODE.isentry;  // DEFAULT
      bool LOCAL_KBIN_VASP_KPOINTS_KPPRA_isentry=vflags.KBIN_VASP_KPOINTS_KPPRA.isentry;  // DEFAULT
      bool LOCAL_KBIN_VASP_KPOINTS_KSCHEME_isentry=vflags.KBIN_VASP_KPOINTS_KSCHEME.isentry;  // DEFAULT
      bool LOCAL_KBIN_VASP_KPOINTS_KSHIFT_isentry=vflags.KBIN_VASP_KPOINTS_KSHIFT.isentry;  // DEFAULT
      int LOCAL_KBIN_VASP_KPOINTS_KMODE_content_int=vflags.KBIN_VASP_KPOINTS_KMODE.content_int;  // DEFAULT
      int LOCAL_KBIN_VASP_KPOINTS_KPPRA_content_int=vflags.KBIN_VASP_KPOINTS_KPPRA.content_int;  // DEFAULT
      string LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string=vflags.KBIN_VASP_KPOINTS_KSCHEME.content_string;  // DEFAULT
      string LOCAL_KBIN_VASP_KPOINTS_KSHIFT_content_string=vflags.KBIN_VASP_KPOINTS_KSHIFT.content_string;  // DEFAULT

      string STRING_KPOINTS_TO_SHOW="KPOINTS";

      if(vflags.KBIN_VASP_RUN.flag("STATIC")) { // do the switching
        STRING_KPOINTS_TO_SHOW="KPOINTS_STATIC";
        LOCAL_KBIN_VASP_KPOINTS_KSCHEME_isentry=vflags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.isentry;  // overrides kscheme
        LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string=vflags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.content_string;  // overrides kscheme
        if(vflags.KBIN_VASP_KPOINTS_STATIC_KMODE.isentry) {
          LOCAL_KBIN_VASP_KPOINTS_KMODE_isentry=vflags.KBIN_VASP_KPOINTS_STATIC_KMODE.isentry;
          LOCAL_KBIN_VASP_KPOINTS_KMODE_content_int=vflags.KBIN_VASP_KPOINTS_STATIC_KMODE.content_int;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << "_KMODE overrides KPOINTS_KMODE generation IMPLICIT file from " << _AFLOWIN_ << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << "_KMODE: LOCAL_KBIN_VASP_KPOINTS_KMODE_isentry=" << vflags.KBIN_VASP_KPOINTS_STATIC_KMODE.isentry << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << "_KMODE: LOCAL_KBIN_VASP_KPOINTS_KMODE_content_int=" << vflags.KBIN_VASP_KPOINTS_STATIC_KMODE.content_int << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if(vflags.KBIN_VASP_KPOINTS_STATIC_KPPRA.isentry) {
          LOCAL_KBIN_VASP_KPOINTS_KPPRA_isentry=vflags.KBIN_VASP_KPOINTS_STATIC_KPPRA.isentry;
          LOCAL_KBIN_VASP_KPOINTS_KPPRA_content_int=vflags.KBIN_VASP_KPOINTS_STATIC_KPPRA.content_uint;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << "_KPPRA overrides KPOINTS_KPPRA generation IMPLICIT file from " << _AFLOWIN_ << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << "_KPPRA: LOCAL_KBIN_VASP_KPOINTS_KPPRA_isentry=" << vflags.KBIN_VASP_KPOINTS_STATIC_KPPRA.isentry << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << "_KPPRA: LOCAL_KBIN_VASP_KPOINTS_KPPRA_content_int=" << vflags.KBIN_VASP_KPOINTS_STATIC_KPPRA.content_uint << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if(vflags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.isentry) {
          LOCAL_KBIN_VASP_KPOINTS_KSCHEME_isentry=vflags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.isentry;
          LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string=vflags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.content_string;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << "_KSCHEME overrides KPOINTS_KSCHEME generation IMPLICIT file from " << _AFLOWIN_ << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << "_KSCHEME: LOCAL_KBIN_VASP_KPOINTS_KSCHEME_isentry=" << vflags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.isentry << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << "_KSCHEME: LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string=" << vflags.KBIN_VASP_KPOINTS_STATIC_KSCHEME.content_string << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if(vflags.KBIN_VASP_KPOINTS_STATIC_KSHIFT.isentry) {
          LOCAL_KBIN_VASP_KPOINTS_KSHIFT_isentry=vflags.KBIN_VASP_KPOINTS_STATIC_KSHIFT.isentry;
          LOCAL_KBIN_VASP_KPOINTS_KSHIFT_content_string=vflags.KBIN_VASP_KPOINTS_STATIC_KSHIFT.content_string;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << "_KSHIFT overrides KPOINTS_KSHIFT generation IMPLICIT file from " << _AFLOWIN_ << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << "_KSHIFT: LOCAL_KBIN_VASP_KPOINTS_KSHIFT_isentry=" << vflags.KBIN_VASP_KPOINTS_STATIC_KSHIFT.isentry << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << "_KSHIFT: LOCAL_KBIN_VASP_KPOINTS_KSHIFT_content_string=" << vflags.KBIN_VASP_KPOINTS_STATIC_KSHIFT.content_string << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
      }

      if(LOCAL_KBIN_VASP_KPOINTS_KMODE_content_int==0 || LOCAL_KBIN_VASP_KPOINTS_KMODE_isentry==FALSE) {
        // KSCHEME ******************************
        if(LOCAL_KBIN_VASP_KPOINTS_KSCHEME_isentry==FALSE) {
          xvasp.str.kpoints_kscheme="Monkhorst-Pack";  // DEFAULT FIX
          xvasp.str.kpoints_kscheme=DEFAULT_KSCHEME;  // DEFAULT FIX
          aus << "00000  MESSAGE-DEFAULT " << STRING_KPOINTS_TO_SHOW << "_KSCHEME=" << xvasp.str.kpoints_kscheme << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        } else {
          xvasp.str.kpoints_kscheme="Monkhorst-Pack";  // DEFAULT FIX
          xvasp.str.kpoints_kscheme=DEFAULT_KSCHEME;  // DEFAULT FIX
          if(LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string.at(0)=='M' || LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string.at(0)=='m') xvasp.str.kpoints_kscheme="Monkhorst-Pack";
          if(LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string.at(0)=='G' || LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string.at(0)=='g') xvasp.str.kpoints_kscheme="Gamma";
          //	if(LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string.at(0)=='A' || LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string.at(0)=='a') xvasp.str.kpoints_kscheme="Auto";
          //	if(LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string.at(0)=='L' || LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string.at(0)=='l') xvasp.str.kpoints_kscheme=DEFAULT_KSCHEME;
          if(LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string.at(0)=='A' || LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string.at(0)=='a') xvasp.str.kpoints_kscheme=DEFAULT_KSCHEME;
          if(LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string.at(0)=='C' || LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string.at(0)=='c') xvasp.str.kpoints_kscheme="Cartesian";
          if(LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string.at(0)=='K' || LOCAL_KBIN_VASP_KPOINTS_KSCHEME_content_string.at(0)=='k') xvasp.str.kpoints_kscheme="Cartesian";
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << "_KSCHEME=" << xvasp.str.kpoints_kscheme << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if(xvasp.str.kpoints_kscheme==DEFAULT_KSCHEME) {
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << " Calculating structure lattice" << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
          xvasp.str.GetLatticeType();
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << " Found Lattice=" << xvasp.str.bravais_lattice_variation_type << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          xvasp.str.kpoints_kscheme="Monkhorst-Pack";
          //          if(xvasp.str.bravais_lattice_variation_type=="HEX" || xvasp.str.bravais_lattice_variation_type=="FCC") xvasp.str.kpoints_kscheme="Gamma";  // also add RHL
          if(xvasp.str.bravais_lattice_variation_type=="HEX" || xvasp.str.bravais_lattice_variation_type=="FCC" || xvasp.str.bravais_lattice_variation_type=="RHL") xvasp.str.kpoints_kscheme="Gamma";  //CO+DX20200404 - 3-fold symmetries REQUIRE Gamma-centered
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << "_KSCHEME=\"" << xvasp.str.kpoints_kscheme << "\" " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        // 
        // KPPRA ******************************
        if(LOCAL_KBIN_VASP_KPOINTS_KPPRA_isentry==FALSE) {
          NK=1; // DEFAULT FIX
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << " KPPRA=NNNN is missing, taking NNNN=" << NK << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        } else {
          NK=LOCAL_KBIN_VASP_KPOINTS_KPPRA_content_int;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << " KPPRA=" << NK << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        stringKPPRA=KPPRA(xvasp.str,NK);
        // VERBOSE on LOCK /SCREEN
        if(TRUE) FileMESSAGE << stringKPPRA;
        if(TRUE) {
          aus.precision(5);
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << " KPPRA routine ["
            << xvasp.str.kpoints_k1 << "," << xvasp.str.kpoints_k2 << "," << xvasp.str.kpoints_k3 << "]=" << xvasp.str.kpoints_k1*xvasp.str.kpoints_k2*xvasp.str.kpoints_k3 << "=["
            << modulus(xvasp.str.klattice(1)/((double) xvasp.str.kpoints_k1)) << ","
            << modulus(xvasp.str.klattice(2)/((double) xvasp.str.kpoints_k2)) << ","
            << modulus(xvasp.str.klattice(3)/((double) xvasp.str.kpoints_k3)) << "]   " << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        // KSHIFT ******************************
        if(LOCAL_KBIN_VASP_KPOINTS_KSHIFT_isentry==FALSE) {
          xvasp.str.kpoints_s1=xvasp.str.kpoints_s2=xvasp.str.kpoints_s3=0.0;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << " KSHIFT= X X X  is missing, taking X X X = 0.0 0.0 0.0 " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        } else {
          xvasp.str.kpoints_s1=xvasp.str.kpoints_s2=xvasp.str.kpoints_s3=0.0;
          stringstream oaus;
          oaus.str(LOCAL_KBIN_VASP_KPOINTS_KSHIFT_content_string);
          oaus >> xvasp.str.kpoints_s1 >> xvasp.str.kpoints_s2 >> xvasp.str.kpoints_s3;
          aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << " KSHIFT= " << xvasp.str.kpoints_s1 << " " << xvasp.str.kpoints_s2 << " " << xvasp.str.kpoints_s3 << " " << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        // create KPOINTS ******************************
        xvasp.KPOINTS.str(std::string());xvasp.KPOINTS.clear();
        xvasp.KPOINTS_orig.str(std::string());xvasp.KPOINTS_orig.clear();
        xvasp.KPOINTS << "KPOINTS File, automatically generated by KPPRA - aflow with Kpoints=" << xvasp.str.kpoints_kppra  << "  [KPPRA=" << NK << "]" << endl;
        xvasp.KPOINTS << xvasp.str.kpoints_mode << endl;  // MODE AUTO 0
        xvasp.KPOINTS << xvasp.str.kpoints_kscheme << endl;
        xvasp.KPOINTS << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << endl;
        xvasp.KPOINTS.precision(3);
        xvasp.KPOINTS << xvasp.str.kpoints_s1 << " " << xvasp.str.kpoints_s2 << " " << xvasp.str.kpoints_s3 << endl;

        aurostd::StringstreamClean(aus);
        aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << " K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aus << "00000  MESSAGE " << STRING_KPOINTS_TO_SHOW << " Kpoints=" << xvasp.str.kpoints_kppra  << "  [KPPRA=" << NK  << "]    with " << xvasp.str.kpoints_kscheme << "   " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);			
      } else {
        aurostd::StringstreamClean(aus);
        aus << "EEEEE             Only [VASP_KPOINTS_FILE]KMODE=0 is supported ! " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
    }
    // EXTERNAL **************************************************
    if(Krun && vflags.KBIN_VASP_KPOINTS_MODE.flag("EXTERNAL")) {  // [VASP_KPOINTS_MODE_EXTERNAL] construction
      string file;
      aus << "00000  MESSAGE KPOINTS  generation EXTERNAL file from " << _AFLOWIN_ << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(vflags.KBIN_VASP_KPOINTS_FILE.flag("COMMAND") && vflags.KBIN_VASP_KPOINTS_FILE.flag("FILE")) {
        aus << "EEEEE   [VASP_KPOINTS_MODE]FILE=  and  [VASP_KPOINTS_MODE]COMMAND=  can not be used together " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
      if(!vflags.KBIN_VASP_KPOINTS_FILE.flag("COMMAND") && (vflags.KBIN_VASP_KPOINTS_FILE.flag("FILE") || !vflags.KBIN_VASP_KPOINTS_FILE.flag("FILE"))) {
        if(vflags.KBIN_VASP_KPOINTS_FILE.flag("FILE")) {
          // [ME20181226 OBSOLETE] file=aurostd::substring2string(AflowIn,"[VASP_KPOINTS_FILE]FILE=",TRUE);
          file = vflags.KBIN_VASP_KPOINTS_FILE.getattachedscheme("FILE"); //ME20181113
          aus << "00000  MESSAGE KPOINTS  generation from file=" << file << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        } else {
          file=DEFAULT_VASP_EXTERNAL_KPOINTS;
          aus << "00000  MESSAGE KPOINTS  generation from DEFAULT file=" << DEFAULT_VASP_EXTERNAL_KPOINTS << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if(!aurostd::FileExist(file)) {
          aus << "EEEEE  ERROR KPOINTS file=" << file << " does not exist! " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        if(aurostd::FileEmpty(file)) {
          aus << "EEEEE  ERROR KPOINTS file=" << file << " is empty! " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        xvasp.KPOINTS << aurostd::file2string(file);
      }
      if(vflags.KBIN_VASP_KPOINTS_FILE.flag("COMMAND") && !vflags.KBIN_VASP_KPOINTS_FILE.flag("FILE")) {
        // [ME20181226 OBSOLETE] file=aurostd::substring2string(AflowIn,"[VASP_KPOINTS_FILE]COMMAND=",FALSE);
        file = vflags.KBIN_VASP_KPOINTS_FILE.getattachedscheme("COMMAND"); //ME20181113
        aus << "00000  MESSAGE KPOINTS  generation from command= '" << file << "' " << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        file=file+" > ./_aflow_KPOINTS."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";    // create temp  //CO20200502 - threadID
        aurostd::execute(file);                           // create temp
        file="./_aflow_KPOINTS."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";            // file name  //CO20200502 - threadID
        if(!aurostd::FileExist(file)) {  // could not write (directory protected)
          aus << "EEEEE  ERROR KPOINTS file=" << file << " does not exist! " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        if(aurostd::FileEmpty(file)) {  // contains nothing good
          aus << "EEEEE  ERROR KPOINTS file=" << file << " is empty! " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        xvasp.KPOINTS << aurostd::file2string(file);       // load KPOINTS
        aurostd::RemoveFile("./_aflow_KPOINTS."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp");     // remove temp //CO20200502 - threadID
      }
      KBIN::XVASP_string2numbers(xvasp);                          // create KPOINTS numbers
    }
    // KPOINTS DONE **************************************************
    if(IMPLICIT==FALSE) KBIN::XVASP_string2numbers(xvasp);
    xvasp.str.kpoints_kmax=max(xvasp.str.kpoints_k1,xvasp.str.kpoints_k2,xvasp.str.kpoints_k3);
    xvasp.str.kpoints_kppra=xvasp.str.kpoints_k1*xvasp.str.kpoints_k2*xvasp.str.kpoints_k3*xvasp.str.atoms.size();
    // KPOINTS done
    xvasp.KPOINTS_orig << xvasp.KPOINTS.str();
    xvasp.aopts.flag("FLAG::XVASP_KPOINTS_generated",TRUE);
    // cerr << xvasp.str.kpoints_kscheme << endl;
    return Krun;
  };  // KBIN::VASP_Produce_KPOINTS
}

namespace KBIN {
  bool VASP_Modify_KPOINTS(_xvasp& xvasp,ofstream &FileMESSAGE,_aflags &aflags,_vflags &vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    ostringstream aus;
    bool Krun=TRUE;

    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("KEEPK")==FALSE && 0) {             /*************** KPOINTS **************/
      vflags.KBIN_VASP_WRITE_KPOINTS=TRUE;
      if(_isodd(xvasp.str.kpoints_k1)) xvasp.str.kpoints_k1++;
      if(_isodd(xvasp.str.kpoints_k2)) xvasp.str.kpoints_k2++;
      if(_isodd(xvasp.str.kpoints_k3)) xvasp.str.kpoints_k3++;
      xvasp.str.kpoints_kmax=max(xvasp.str.kpoints_k1,xvasp.str.kpoints_k2,xvasp.str.kpoints_k3);
      xvasp.str.kpoints_kppra=xvasp.str.kpoints_k1*xvasp.str.kpoints_k2*xvasp.str.kpoints_k3*xvasp.str.atoms.size();
      KBIN::XVASP_numbers2string(xvasp);
      xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
      aus << "00000  MESSAGE KPOINTS Option Tune - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << " " << xvasp.str.kpoints_kmax << "] - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }

    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.isentry) {                        /*************** KPOINTS **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]KPOINTS=" << vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.content_string << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("KEEPK")) {                        /*************** KPOINTS **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]KPOINTS=KEEPK - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      //  KBIN::XVASP_KPOINTS_OPERATION(xvasp,""); //    KBIN::XVASP_KPOINTS_KEEPK(xvasp);
      xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
    }
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("EVEN")) {                        /*************** KPOINTS **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]KPOINTS=EVEN - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_KPOINTS_OPERATION(xvasp,"Xeven,Yeven,Zeven"); //    KBIN::XVASP_KPOINTS_EVEN(xvasp);
      xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
    }
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("ODD")) {                         /*************** KPOINTS **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]KPOINTS=ODD - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_KPOINTS_OPERATION(xvasp,"Xodd,Yodd,Zodd");  // KBIN::XVASP_KPOINTS_ODD(xvasp);
      xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
    }
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("KSHIFT_GAMMA_EVEN")) {           /*************** KPOINTS **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]KPOINTS=KSHIFT_GAMMA_EVEN - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_KPOINTS_OPERATION(xvasp,"Xevenshift,Yevenshift,Zevenshift");   // KBIN::XVASP_KPOINTS_Kshift_Gamma_EVEN(_xvasp& xvasp);
      xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
    }
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("KSHIFT_GAMMA_ODD")) {            /*************** KPOINTS **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]KPOINTS=KSHIFT_GAMMA_ODD - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_KPOINTS_OPERATION(xvasp,"Xoddshift,Yoddshift,Zoddshift");   // KBIN::XVASP_KPOINTS_Kshift_Gamma_ODD(xvasp);
      xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
    }
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("KSCHEME_MONKHORST_PACK")) {      /*************** KPOINTS **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]KPOINTS=KSCHEME_MONKHORST_PACK - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_KPOINTS_OPERATION(xvasp,"Monkhorst-Pack");
      xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
    }
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("KSCHEME_GAMMA")) {               /*************** KPOINTS **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]KPOINTS=KSCHEME_GAMMA - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_KPOINTS_OPERATION(xvasp,"Gamma");
      xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
    }
    if(Krun && (vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE && vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("KSCHEME_AUTO"))) {              /*************** KPOINTS **************/
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]KPOINTS=KSCHEME_AUTO  Calculating structure lattice" << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.str.GetLatticeType();
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]KPOINTS=KSCHEME_AUTO  Found Lattice=" << xvasp.str.bravais_lattice_variation_type << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      xvasp.str.kpoints_kscheme="Monkhorst-Pack";
      //     if(xvasp.str.bravais_lattice_variation_type=="HEX" || xvasp.str.bravais_lattice_variation_type=="FCC") xvasp.str.kpoints_kscheme="Gamma"; // add RHL
      if(xvasp.str.bravais_lattice_variation_type=="HEX" || xvasp.str.bravais_lattice_variation_type=="FCC" || xvasp.str.bravais_lattice_variation_type=="RHL") xvasp.str.kpoints_kscheme="Gamma";  //CO+DX20200404 - 3-fold symmetries REQUIRE Gamma-centered
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]KPOINTS=KSCHEME_" << xvasp.str.kpoints_kscheme << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_KPOINTS_OPERATION(xvasp,xvasp.str.kpoints_kscheme);
      xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
    }
    if(vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("GAMMA")) {
      vflags.KBIN_VASP_WRITE_KPOINTS=TRUE;
      xvasp.str.kpoints_k1=1;xvasp.str.kpoints_k2=1;xvasp.str.kpoints_k3=1;
      xvasp.str.kpoints_kmax=max(xvasp.str.kpoints_k1,xvasp.str.kpoints_k2,xvasp.str.kpoints_k3);
      xvasp.str.kpoints_kppra=xvasp.str.kpoints_k1*xvasp.str.kpoints_k2*xvasp.str.kpoints_k3*xvasp.str.atoms.size();
      KBIN::XVASP_numbers2string(xvasp);
      xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]KPOINTS=GAMMA - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << " " << xvasp.str.kpoints_kmax << "] - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    if(vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("IBZKPT")) {
      aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]KPOINTS=IBZKPT - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    // ------------------------------------
    // end
    xvasp.aopts.flag("FLAG::XVASP_KPOINTS_generated",TRUE);
    return Krun;
  }; // KBIN::VASP_Modify_KPOINTS
}

namespace KBIN {
  bool VASP_Reread_KPOINTS(_xvasp& xvasp,ofstream &FileMESSAGE,_aflags &aflags) { // AFLOW_FUNCTION_IMPLEMENTATION
    ostringstream aus;
    bool Krun=TRUE;
    if(!aurostd::FileExist(xvasp.Directory+"/KPOINTS")) {
      aus << "EEEEE  KBIN::VASP_Reread_KPOINTS: KPOINTS not present in directory: " << xvasp.Directory << " - "  << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      Krun=FALSE;
      return Krun;
    }
    xvasp.KPOINTS_orig.str(std::string()); xvasp.KPOINTS_orig << xvasp.KPOINTS.str();
    xvasp.KPOINTS.str(std::string()); xvasp.KPOINTS << aurostd::file2string(xvasp.Directory+"/KPOINTS"); // DID REREAD
    xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
    return Krun;
  }
}

// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// POTCAR
namespace KBIN {
  bool VASP_Find_DATA_POTCAR(const string& species_pp,string &FilePotcar,string &DataPotcar,string &AUIDPotcar) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string pseudopotential="nothing";
    FilePotcar="";DataPotcar="";
    if(aurostd::substring2bool(species_pp,"pot_LDA") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POT_LDA))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POT_LDA;
    if(aurostd::substring2bool(species_pp,"pot_GGA") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POT_GGA))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POT_GGA;
    if(aurostd::substring2bool(species_pp,"pot_PBE") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POT_PBE))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POT_PBE;
    if(aurostd::substring2bool(species_pp,"potpaw_LDA") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA;
    if(aurostd::substring2bool(species_pp,"potpaw_GGA") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA;
    if(aurostd::substring2bool(species_pp,"potpaw_PBE") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE;
    if(aurostd::substring2bool(species_pp,"potpaw_LDA_KIN") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN;
    if(aurostd::substring2bool(species_pp,"potpaw_PBE_KIN") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN;

    if(LDEBUG) cerr << "[VASP_Find_DATA_POTCAR]: species_pp=" << species_pp << endl; 
    if(LDEBUG) cerr << "[VASP_Find_DATA_POTCAR]: pseudopotential=" << pseudopotential << endl;

    bool found=FALSE;

    string string2load_list="AFLOW_PSEUDOPOTENTIALS_LIST_TXT:"+species_pp+"/POTCAR";
    string string2load_data="AFLOW_PSEUDOPOTENTIALS_TXT:"+species_pp+"/POTCAR";
    if(aurostd::substring2bool(init::InitLoadString(string2load_list,FALSE),"POTCAR")) {
      found=TRUE;
      FilePotcar=init::InitLoadString(string2load_list,FALSE);
      DataPotcar=init::InitLoadString(string2load_data,FALSE);
      string file_tmp=aurostd::TmpFileCreate("DataPotcar");
      AUIDPotcar=aurostd::file2auid(file_tmp);
      if(aurostd::FileExist(file_tmp)) aurostd::RemoveFile(file_tmp);

      if(LDEBUG) cerr << "[VASP_Find_DATA_POTCAR]: species_pp=" << species_pp << endl; 
      if(LDEBUG) cerr << "[VASP_Find_DATA_POTCAR]: pseudopotential=" << pseudopotential << endl;
      if(LDEBUG) cerr << "[VASP_Find_DATA_POTCAR]: AUIDPotcar=" << AUIDPotcar << endl; 
    }
    return found;
  }
}

namespace KBIN {
  bool VASP_Find_FILE_POTCAR(const string& species_pp,string &FilePotcar,string &DataPotcar,string &AUIDPotcar) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string pseudopotential="nothing";
    FilePotcar="";DataPotcar="";
    if(aurostd::substring2bool(species_pp,"pot_LDA") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POT_LDA))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POT_LDA;
    if(aurostd::substring2bool(species_pp,"pot_GGA") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POT_GGA))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POT_GGA;
    if(aurostd::substring2bool(species_pp,"pot_PBE") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POT_PBE))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POT_PBE;
    if(aurostd::substring2bool(species_pp,"potpaw_LDA") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA;
    if(aurostd::substring2bool(species_pp,"potpaw_GGA") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA;
    if(aurostd::substring2bool(species_pp,"potpaw_PBE") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE;
    if(aurostd::substring2bool(species_pp,"potpaw_LDA_KIN") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN;
    if(aurostd::substring2bool(species_pp,"potpaw_PBE_KIN") || aurostd::substring2bool(species_pp,DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN))
      pseudopotential=DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN;

    if(LDEBUG) cerr << "[VASP_Find_FILE_POTCAR]: species_pp=" << species_pp << endl; 
    if(LDEBUG) cerr << "[VASP_Find_FILE_POTCAR]: pseudopotential=" << pseudopotential << endl;

    bool found=FALSE;
    FilePotcar="";
    for(uint j=0;j<vVASP_POTCAR_DIRECTORIES.size()&&!found;j++) {   // cycle through possible directories
      for(uint k=0;k<=9&&!found;k++) {                     // cycle through fixing current before/after
        string FilePotcark;
        if(k==0) FilePotcark=species_pp;      // k==0 dont touch
        if(k>0) FilePotcark=vVASP_POTCAR_DIRECTORIES.at(j)+"/"+species_pp;
        //if(k==1) FilePotcark=FilePotcark; // nothing done //CO20190329 - compiler doesn't like FilePotcark=FilePotcark
        if(k==2) aurostd::StringSubst(FilePotcark,pseudopotential,pseudopotential+"/"+DEFAULT_VASP_POTCAR_DATE);  // current POST
        if(k==3) aurostd::StringSubst(FilePotcark,pseudopotential,DEFAULT_VASP_POTCAR_DATE+"/"+pseudopotential);  // current PRE
        if(k==4) aurostd::StringSubst(FilePotcark,pseudopotential,pseudopotential+"/dateP1");   // dateP1 POST
        if(k==5) aurostd::StringSubst(FilePotcark,pseudopotential,"dateP1/"+pseudopotential);   // dateP1 PRE
        if(k==6) aurostd::StringSubst(FilePotcark,pseudopotential,pseudopotential+"/dateP2");   // dateP2 POST
        if(k==7) aurostd::StringSubst(FilePotcark,pseudopotential,"dateP2/"+pseudopotential);   // dateP2 PRE
        if(k==8) aurostd::StringSubst(FilePotcark,pseudopotential,pseudopotential+"/dateP3");   // dateP3 POST
        if(k==9) aurostd::StringSubst(FilePotcark,pseudopotential,"dateP3/"+pseudopotential);   // dateP3 PRE
        for(uint l=0;l<=1&&!found;l++) {                  // cycle through POTCAR as postfix
          FilePotcar=FilePotcark;                      // default
          if(l==0) FilePotcar=FilePotcark;             // current NO
          if(l==1) FilePotcar=FilePotcark+"/POTCAR";   // POTCAR after
          // clean up to avoid comments //
          // TESTING POTCARS
          if(!aurostd::substring2bool(FilePotcar,"POTCAR")) FilePotcar=FilePotcar+"/POTCAR"; // some fix, _AFLOWIN_ might have some personality issues
          // DEBUG=TRUE;
          FilePotcar=aurostd::CleanFileName(FilePotcar);
          if(aurostd::FileExist(FilePotcar) && !aurostd::FileEmpty(FilePotcar)) found=TRUE;
          if(LDEBUG) cout << "DDDDD  POTCAR  (" << species_pp << ",j=" << j << ",k=" << k << ",l=" << l << ")  [FilePotcar=" <<  FilePotcar << "] found=" << found << " aurostd::FileExist(FilePotcar)=" << aurostd::FileExist(FilePotcar) << endl;
        } // l-cycle through POTCAR as postfix
      } // k-cycle through fixing current before/after
    } // j-cycle through possible directories

    if(found) {
      stringstream aus;
      ifstream FileINPUT;
      FileINPUT.clear();
      FileINPUT.open(FilePotcar.c_str(),std::ios::in);
      char c;
      while (FileINPUT.get(c)) aus.put(c);
      FileINPUT.clear();FileINPUT.close();
      DataPotcar=aus.str();

      AUIDPotcar=aurostd::file2auid(FilePotcar);
      if(LDEBUG) cerr << "[VASP_Find_FILE_POTCAR]: species_pp=" << species_pp << endl; 
      if(LDEBUG) cerr << "[VASP_Find_FILE_POTCAR]: pseudopotential=" << pseudopotential << endl;
      if(LDEBUG) cerr << "[VASP_Find_FILE_POTCAR]: AUIDPotcar=" << AUIDPotcar << endl; 
    }

    return found;
  }
}

namespace KBIN {
  bool VASP_Produce_POTCAR(_xvasp& xvasp,const string& AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags,_vflags &vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy = XPID + "KBIN::VASP_Produce_POTCAR():";
    if(AflowIn.length()==0) {throw aurostd::xerror(_AFLOW_FILE_NAME_,"KBIN::VASP_Produce_POTCAR():","empty AflowIn",_FILE_CORRUPT_);} //CO20200624
    if(!kflags.AFLOW_MODE_VASP) {cerr << XPID << "KBIN::VASP_Produce_POTCAR: should kflags.AFLOW_MODE_VASP be set ??" << endl;}
    string::size_type sub_size1,sub_size2;
    string subS,subS1,subS2,subSDIR="$POTCARDIR";
    ostringstream aus;
    bool Krun=TRUE;
    xvasp.POTCAR.str(std::string());
    xvasp.POTCAR_AUID.clear();
    xvasp.aopts.flag("FLAG::XVASP_POTCAR_generated",FALSE);
    xvasp.aopts.flag("FLAG::XVASP_POTCAR_changed",FALSE);
    xvasp.POTCAR_POTENTIALS.str(std::string());

    // IMPLICIT or EXPLICIT or EXTERNAL for POTCAR
    Krun=(Krun && (vflags.KBIN_VASP_POTCAR_MODE.flag("IMPLICIT") || vflags.KBIN_VASP_POTCAR_MODE.flag("EXPLICIT") || vflags.KBIN_VASP_POTCAR_MODE.flag("EXTERNAL")));
    if(!Krun) {
      aurostd::StringstreamClean(aus);
      aus << "EEEEE  [VASP_POTCAR_MODE_IMPLICIT] or [VASP_POTCAR_MODE_EXPLICIT] or [VASP_POTCAR_MODE_EXPLICIT] must be specified " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      Krun=FALSE;
      return Krun;
    }
    // IMPLICIT **************************************************
    if(Krun && vflags.KBIN_VASP_POTCAR_MODE.flag("IMPLICIT") && vflags.KBIN_VASP_POTCAR_FILE.flag("KEYWORD")) { // [VASP_POTCAR_MODE_IMPLICIT] construction
      // Prepare POTCAR
      aus << "00000  MESSAGE POTCAR  generation IMPLICIT file from " << _AFLOWIN_ << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.isentry) {         /*************** POTCAR **************/
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]AUTO_PSEUDOPOTENTIALS - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        aus << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]AUTO_PSEUDOPOTENTIALS=\"" << vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.xscheme << "\" - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      // Prepare POTCAR
      ifstream FileINPUT;
      vector<string> FilePotcars;
      string FilePotcar,DataPotcar,AUIDPotcar;
      if(!vflags.KBIN_VASP_POTCAR_FILE.flag("SYSTEM_AUTO")) {
        subS2="[VASP_POTCAR_FILE]"; // STRING TO SEARCH
        subS1=AflowIn.substr(AflowIn.find(subS2));
        while (aurostd::substring2bool(subS1,subS2)) {
          subS=subS1;
          sub_size1=subS.find(subS2);
          sub_size2=(subS.substr(sub_size1)).find("\n");
          subS=subS.substr(sub_size1+subS2.length(),sub_size2-subS2.length());
          subS1=subS1.substr(sub_size1+1);
          if(aurostd::substring2bool(subS,subSDIR)) {
            string subSpre,subSpost;
            sub_size1=subS.find(subSDIR)+subSDIR.length();
            sub_size2=(subS.substr(sub_size1)).find("\n");
            subSpost=subS.substr(sub_size1,sub_size2);
            sub_size1=0;
            sub_size2=(subS.substr(sub_size1)).find(subSDIR);
            subSpre=subS.substr(sub_size1,sub_size2);
            subS.clear();
            subS=subSpre+vVASP_POTCAR_DIRECTORIES.at(0)+subSpost;
          }
          FilePotcars.push_back(subS);
        }
      }
      if(vflags.KBIN_VASP_POTCAR_FILE.flag("SYSTEM_AUTO")) {
        for(uint i=0;i<xvasp.str.species.size();i++) {
          if(xvasp.str.species.at(i)!="") {
            subS="";
            if(!vflags.KBIN_VASP_POTCAR_FILE.flag("PREFIX")) {
              subS+=vVASP_POTCAR_DIRECTORIES.at(0);
              subS+=DEFAULT_VASP_PSEUDOPOTENTIAL_TYPE+"/"+DEFAULT_VASP_POTCAR_DATE+"/";
            }
            if(vflags.KBIN_VASP_POTCAR_FILE.flag("PREFIX")) {
              subS=vflags.KBIN_VASP_POTCAR_FILE.getattachedscheme("PREFIX");  //aurostd::substring2string(AflowIn,"[VASP_POTCAR_FILE]PREFIX=",TRUE);  //CO20181226
              if(LDEBUG) cerr << "DEBUG subS=" << subS << endl;
              if(aurostd::substring2bool(subS,subSDIR)) {
                subS=vVASP_POTCAR_DIRECTORIES.at(0)+aurostd::substring2string(subS,subSDIR,TRUE);
              }
              if(LDEBUG) cerr << "DEBUG subS=" << subS << endl;
            }
            subS+="/"+xvasp.str.species.at(i);
            if(!vflags.KBIN_VASP_POTCAR_FILE.flag("SUFFIX"))
              subS+=DEFAULT_VASP_POTCAR_SUFFIX;
            if(vflags.KBIN_VASP_POTCAR_FILE.flag("SUFFIX"))
              subS+="/"+vflags.KBIN_VASP_POTCAR_FILE.getattachedscheme("SUFFIX"); //aurostd::substring2string(AflowIn,"[VASP_POTCAR_FILE]SUFFIX=",TRUE);  //CO20181226
            FilePotcars.push_back(subS);
          }
        }
      }
      // FIX PSEUDOPOTENTIALS
      //CO20181226 - originally this was applied regardless of AUTO_PSEUDOPOTENTIALS
      //this causes problems, as AUTO_PSEUDOPOTENTIALS strips PP information ([VASP_POTCAR_FILE]) in automatic aflow generation (see avasp)
      //therefore, Mn_pv becomes Mn (WRONG)
      //I believe the patch is to do this ONLY if !AUTO_PSEUDOPOTENTIALS
      //if this guard is removed, please patch fixEmptyAtomNames() to be force_fix=false (APL/aflow_apl_kphonons)
      //ME20190211 - This patch causes species_pp to be unoccupied if AUTO_PSEUDOPOTENTIALS (fixed below)
      //ME20190313 - This breaks when [VASP_POTCAR_FILE] is a path - need to break the path apart
      if(!vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.isentry) {  //CO20181226
        vector<string> tokens;  //ME20190313
        for(uint i=0;i<FilePotcars.size() && i<xvasp.str.species_pp.size();i++) {
          aurostd::string2tokens(FilePotcars.at(i), tokens, "/");  //ME20190313
          xvasp.str.species_pp.at(i) = tokens.back();  //ME20190313
          //xvasp.str.species_pp.at(i)=FilePotcars.at(i);
        }
      }
      // IF AUTO_PSEUDOPOTENTIALS
      if(LDEBUG) {cerr << soliloquy << " vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.isentry=" << vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.isentry << endl;}
      if(vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.isentry) {
        for(uint i=0;i<FilePotcars.size();i++) {
          string cleanname=KBIN::VASP_PseudoPotential_CleanName(FilePotcars.at(i));
          aurostd::StringSubst(cleanname,"/POTCAR","");
          vector<string> tokens;aurostd::string2tokens(cleanname,tokens,"/");
          if(tokens.size()>0) cleanname=tokens.at(tokens.size()-1);
          if(LDEBUG) cerr << "[VASP_Produce_POTCAR]: vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.xscheme=" << vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.xscheme << endl;
          if(vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.xscheme=="potpaw_PBE_KIN")
            FilePotcars.at(i)=DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN+"/"+AVASP_Get_PseudoPotential_PAW_PBE_KIN(cleanname);
          if(vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.xscheme=="potpaw_LDA_KIN")
            FilePotcars.at(i)=DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN+"/"+AVASP_Get_PseudoPotential_PAW_LDA_KIN(cleanname);

          if(vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.xscheme=="potpaw_PBE")
            FilePotcars.at(i)=DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE+"/"+AVASP_Get_PseudoPotential_PAW_PBE(cleanname);
          if(vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.xscheme=="potpaw_GGA")
            FilePotcars.at(i)=DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA+"/"+AVASP_Get_PseudoPotential_PAW_GGA(cleanname);
          if(vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.xscheme=="potpaw_LDA")
            FilePotcars.at(i)=DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA+"/"+AVASP_Get_PseudoPotential_PAW_LDA(cleanname);

          if(vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.xscheme=="pot_PBE")
            FilePotcars.at(i)=DEFAULT_VASP_POTCAR_DIR_POT_PBE+"/"+AVASP_Get_PseudoPotential_PBE(cleanname);
          if(vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.xscheme=="pot_GGA")
            FilePotcars.at(i)=DEFAULT_VASP_POTCAR_DIR_POT_GGA+"/"+AVASP_Get_PseudoPotential_GGA(cleanname);
          if(vflags.KBIN_VASP_FORCE_OPTION_AUTO_PSEUDOPOTENTIALS.xscheme=="pot_LDA")
            FilePotcars.at(i)=DEFAULT_VASP_POTCAR_DIR_POT_LDA+"/"+AVASP_Get_PseudoPotential_LDA(cleanname);

          //ME20190211 - Occupy species_pp
          aurostd::string2tokens(FilePotcars[i], tokens, "/");
          if (tokens.size() > 0) xvasp.str.species_pp[i] = tokens.back();

          if(LDEBUG) cerr << "[VASP_Produce_POTCAR]: cleanname=" << cleanname << endl;
          if(LDEBUG) cerr << "[VASP_Produce_POTCAR]: FilePotcars.at(" << i << ")=" << FilePotcars.at(i) << endl;
        }
      }

      // PRINT POTCARS
      for(uint i=0;i<FilePotcars.size();i++)
        aus << "00000  MESSAGE POTCAR  Loading potcar file = [" << FilePotcars.at(i) << "]" << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      // LOAD POTCARS (test some fixes because of different machines)

      string FilePotcarNaked;
      vector<string> tokens;

      bool AFLOW_PSEUDOPOTENTIALS=FALSE;
      if(aurostd::substring2bool(init::InitGlobalObject("AFLOW_PSEUDOPOTENTIALS"),"1")) {
        AFLOW_PSEUDOPOTENTIALS=TRUE;
        aus << "00000  MESSAGE POTCAR  AFLOW_PSEUDOPOTENTIALS=TRUE - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      }
      //	AFLOW_PSEUDOPOTENTIALS=FALSE;

      for(uint i=0;i<FilePotcars.size();i++) {                        // cycle through all POTCARS
        aurostd::StringSubst(FilePotcars.at(i),"//","/");
        // strip everything
        tokens.clear();
        FilePotcarNaked=FilePotcars.at(i);
        for(uint iext=0;iext<XHOST.vext.size();iext++) { 
          FilePotcarNaked=aurostd::RemoveSubStringFirst(FilePotcarNaked,"/POTCAR"+XHOST.vext.at(iext));
        }
        aurostd::string2tokens(FilePotcarNaked,tokens,"/");
        xvasp.POTCAR_POTENTIALS << tokens.at(tokens.size()-1);
        // now build it

        bool found_DATA=FALSE;
        bool found_FILE=FALSE;
        if(AFLOW_PSEUDOPOTENTIALS) {
          found_DATA=KBIN::VASP_Find_DATA_POTCAR(FilePotcars.at(i),FilePotcar,DataPotcar,AUIDPotcar);
          if(found_DATA) {
            aus << "00000  MESSAGE POTCAR  DATA: Found potcar FilePotcar=" << FilePotcar << " - DataPotcar.size()=" << DataPotcar.size() << " - AUIDPotcar=" << AUIDPotcar << " " << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            xvasp.POTCAR << DataPotcar << endl;
            xvasp.POTCAR_AUID.push_back(AUIDPotcar);
            // [OBSOLETE] aus << "00000  MESSAGE POTCAR  DATA: Found potcar xvasp.POTCAR.str().size()=" << xvasp.POTCAR.str().size() << " " << endl;
            // [OBSOLETE] aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	      
          }
        } else {
          found_FILE=KBIN::VASP_Find_FILE_POTCAR(FilePotcars.at(i),FilePotcar,DataPotcar,AUIDPotcar);
          if(found_FILE) {
            // FileINPUT.clear();FileINPUT.open(FilePotcar.c_str(),std::ios::in);
            // char c;while (FileINPUT.get(c)) xvasp.POTCAR.put(c);
            // FileINPUT.clear();FileINPUT.close();
            aus << "00000  MESSAGE POTCAR  FILE: Found potcar FilePotcar=" << FilePotcar << " - DataPotcar.size()=" << DataPotcar.size() << " - AUIDPotcar=" << AUIDPotcar << " " << endl;
            aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
            xvasp.POTCAR << DataPotcar;// << endl; file has already new line
            xvasp.POTCAR_AUID.push_back(AUIDPotcar);
            // [OBSOLETE] aus << "00000  MESSAGE POTCAR  FILE: Found potcar xvasp.POTCAR.str().size()=" << xvasp.POTCAR.str().size() << " " << endl;
            // [OBSOLETE] aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);	      
          } else {
            aus << "EEEEE  POTCAR [" << FilePotcars.at(i).c_str() << "] not found! (aflow_ivasp.cpp) " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
            aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
            Krun=FALSE;                                         // DO NOT RUN ANYMORE
            return Krun;
          }
        }
      }
    }
    // EXPLICIT **************************************************
    if(Krun && vflags.KBIN_VASP_POTCAR_MODE.flag("EXPLICIT")) {   // [VASP_POTCAR_MODE_EXPLICIT] construction
      // Prepare POTCAR
      aus << "00000  MESSAGE POTCAR  generation EXPLICIT " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      {
        ifstream FileINPUT;
        string FileNameINPUT=xvasp.Directory+"/"+_AFLOWIN_;
        string line;bool lflag=FALSE;
        FileINPUT.open(FileNameINPUT.c_str(),std::ios::in);
        while (getline(FileINPUT,line)) {
          if(lflag==TRUE) xvasp.POTCAR << line << endl;
          //	if(lflag==TRUE)  xvasp.POTCAR+=line+"\n";
          if(aurostd::substring2bool(line,"[VASP_POTCAR_MODE_EXPLICIT]")) lflag=TRUE;
        }
        FileINPUT.clear();FileINPUT.close();
      }
    }
    // EXTERNAL **************************************************
    if(Krun && vflags.KBIN_VASP_POTCAR_MODE.flag("EXTERNAL")) {  // [VASP_POTCAR_MODE_EXTERNAL] construction
      string file;
      aus << "00000  MESSAGE POTCAR  generation EXTERNAL file from " << _AFLOWIN_ << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      if(vflags.KBIN_VASP_POTCAR_FILE.flag("COMMAND") && vflags.KBIN_VASP_POTCAR_FILE.flag("FILE")) {
        aus << "EEEEE   [VASP_POTCAR_MODE]FILE=  and  [VASP_POTCAR_MODE]COMMAND=  can not be used together " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
        Krun=FALSE;
        return Krun;
      }
      if(!vflags.KBIN_VASP_POTCAR_FILE.flag("COMMAND") && (vflags.KBIN_VASP_POTCAR_FILE.flag("FILE") || !vflags.KBIN_VASP_POTCAR_FILE.flag("FILE"))) {
        if(vflags.KBIN_VASP_POTCAR_FILE.flag("FILE")) {
          file=vflags.KBIN_VASP_POTCAR_FILE.getattachedscheme("FILE"); //aurostd::substring2string(AflowIn,"[VASP_POTCAR_FILE]FILE=",TRUE); //CO20181226
          aus << "00000  MESSAGE POTCAR  generation from file=" << file << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        } else {
          file=DEFAULT_VASP_EXTERNAL_POTCAR;
          aus << "00000  MESSAGE POTCAR  generation from DEFAULT file=" << DEFAULT_VASP_EXTERNAL_POTCAR << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        }
        if(!aurostd::FileExist(file)) {
          aus << "EEEEE  ERROR POTCAR file=" << file << " does not exist! " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        if(aurostd::FileEmpty(file)) {
          aus << "EEEEE  ERROR POTCAR file=" << file << " is empty! " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        xvasp.POTCAR << aurostd::file2string(file);
      }
      if(vflags.KBIN_VASP_POTCAR_FILE.flag("COMMAND") && !vflags.KBIN_VASP_POTCAR_FILE.flag("FILE")) {
        file=vflags.KBIN_VASP_POTCAR_FILE.getattachedscheme("COMMAND"); //aurostd::substring2string(AflowIn,"[VASP_POTCAR_FILE]COMMAND=",FALSE);  //CO20181226
        aus << "00000  MESSAGE POTCAR  generation from command= '" << file << "' " << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
        file=file+" > ./_aflow_POTCAR."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";    // create temp //CO20200502 - threadID
        aurostd::execute(file);                           // create temp
        file="./_aflow_POTCAR."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp";            // file name //CO20200502 - threadID
        if(!aurostd::FileExist(file)) {  // could not write (directory protected)
          aus << "EEEEE  ERROR POTCAR file=" << file << " does not exist! " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        if(aurostd::FileEmpty(file)) {  // contains nothing good
          aus << "EEEEE  ERROR POTCAR file=" << file << " is empty! " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
          Krun=FALSE;
          return Krun;
        }
        xvasp.POTCAR << aurostd::file2string(file);       // load POTCAR
        aurostd::RemoveFile("./_aflow_POTCAR."+XHOST.ostrPID.str()+"."+XHOST.ostrTID.str()+".tmp");     // remove temp  //CO20200502 - threadID
      }
    }

    // get xvasp.POTCAR_ENMAX **************************************************
    xPOTCAR potcar;potcar.GetProperties(xvasp.POTCAR);
    xvasp.POTCAR_ENMAX=potcar.ENMAX;
    xvasp.POTCAR_ENMIN=potcar.ENMIN;
    aus << "00000  MESSAGE POTCAR  ENMAX max   = " << xvasp.POTCAR_ENMAX << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
    aus << "00000  MESSAGE POTCAR  ENMIN min   = " << xvasp.POTCAR_ENMIN << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    // get xvasp.POTCAR_PAW **************************************************
    xvasp.POTCAR_PAW=potcar.POTCAR_PAW;
    xvasp.POTCAR_TYPE=potcar.POTCAR_TYPE;
    aus << "00000  MESSAGE POTCAR  PAW = " << (xvasp.POTCAR_PAW ? "TRUE":"FALSE") << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    aus << "00000  MESSAGE POTCAR  POTCAR_TYPE=" <<  xvasp.POTCAR_TYPE << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    xvasp.POTCAR_orig << xvasp.POTCAR.str();
    xvasp.aopts.flag("FLAG::XVASP_POTCAR_generated",TRUE);

    // done now return
    return Krun;
  };  // KBIN::VASP_Produce_POTCAR
}



namespace KBIN {
  bool VASP_Modify_POTCAR(_xvasp& xvasp,ofstream &FileMESSAGE,_aflags &aflags,_vflags &vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    ostringstream aus;
    bool Krun=TRUE;
    // ------------------------------------
    if(Krun && vflags.KBIN_VASP_FORCE_OPTION_NOTUNE.isentry==FALSE) {
      if(0) {
        aus << "00000  MESSAGE-OPTION  XXXXX" << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      }
    }
    // end
    xvasp.aopts.flag("FLAG::XVASP_POTCAR_generated",TRUE);
    return Krun;
  }; // KBIN::VASP_Produce_POTCAR
}


namespace KBIN {
  bool VASP_Reread_POTCAR(_xvasp& xvasp,ofstream &FileMESSAGE,_aflags &aflags) { // AFLOW_FUNCTION_IMPLEMENTATION
    ostringstream aus;
    bool Krun=TRUE;
    if(!aurostd::FileExist(xvasp.Directory+"/POTCAR")) {
      aus << "EEEEE  KBIN::VASP_Reread_POTCAR: POTCAR not present in directory: " << xvasp.Directory << " - "  << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintErrorStream(FileMESSAGE,aus,XHOST.QUIET);
      Krun=FALSE;
      return Krun;
    }
    xvasp.POTCAR_orig.str(std::string()); xvasp.POTCAR_orig << xvasp.POTCAR.str();
    xvasp.POTCAR.str(std::string()); xvasp.POTCAR << aurostd::file2string(xvasp.Directory+"/POTCAR"); // DID REREAD
    xvasp.aopts.flag("FLAG::XVASP_POTCAR_changed",TRUE);
    return Krun;
  }
}


// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// INCAR MODIFICATIONS

// ***************************************************************************
// KBIN::XVASP_MPI_Autotune
// ***************************************************************************
namespace KBIN {
  void VASP_MPI_Autotune(_xvasp& xvasp,_aflags &aflags,bool VERBOSE) {
    string FileContent,strline;
    int imax;

    FileContent=xvasp.INCAR.str();
    xvasp.INCAR.str(std::string());
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i);
      if(aurostd::substring2bool(strline,"IALGO",TRUE)  || aurostd::substring2bool(strline,"#IALGO",TRUE) ||
          aurostd::substring2bool(strline,"LPLANE",TRUE) || aurostd::substring2bool(strline,"#LPLANE",TRUE) ||
          aurostd::substring2bool(strline,"NPAR",TRUE)   || aurostd::substring2bool(strline,"#NPAR",TRUE)   ||
          aurostd::substring2bool(strline,"LSCALU",TRUE) || aurostd::substring2bool(strline,"#LSCALU",TRUE) ||
          aurostd::substring2bool(strline,"NSIM",TRUE)   || aurostd::substring2bool(strline,"#NSIM",TRUE)) {
        if(VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::VASP_MPI_Autotune) " << endl;
      } else {
        if(!VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
        if(VERBOSE) xvasp.INCAR << strline << endl;
      }
    }
    // xvasp.INCAR << endl;
    if(VERBOSE) xvasp.INCAR << "# Performing MPI-AUTOTUNE [AFLOW] - begin" << endl;
    if(VERBOSE) xvasp.INCAR << "# MPI tuning for Linux cluster with CPUS=" << xvasp.NCPUS << endl;
    if(VERBOSE) xvasp.INCAR << "# look at the VASP manual " << endl;
    if(VERBOSE) xvasp.INCAR << "# http://www.teragridforum.org/mediawiki/images/7/7e/Performance_on_Ranger_of_the_medium_size_problem.pdf  " << endl;
    xvasp.INCAR << aurostd::PaddedPOST("IALGO=48",_incarpad_) << "# MPI algorithm " << endl;
    xvasp.INCAR << aurostd::PaddedPOST("LPLANE=.TRUE.",_incarpad_) << "# parallelization " << endl;
    // xvasp.INCAR << aurostd::PaddedPOST("NPAR="+aurostd::utype2string(xvasp.NCPUS),_incarpad_) << "# number of CPUs" << endl;
    // xvasp.INCAR << aurostd::PaddedPOST("NPAR="+aurostd::utype2string(sqrt(xvasp.NCPUS)),_incarpad_) << "# number of CPUs" << endl;
    int NPAR=0,NCORE=1;
    if(xvasp.NCPUS >   0 && xvasp.NCPUS <=  32) {NPAR=2;NCORE=xvasp.NCPUS/NPAR;}
    if(xvasp.NCPUS >  32 && xvasp.NCPUS <=  64) {NPAR=4;NCORE=xvasp.NCPUS/NPAR;}
    if(xvasp.NCPUS >  64 && xvasp.NCPUS <= 128) {NPAR=8;NCORE=xvasp.NCPUS/NPAR;}
    if(xvasp.NCPUS > 128 && xvasp.NCPUS <= 256) {NPAR=16;NCORE=xvasp.NCPUS/NPAR;}
    if(xvasp.NCPUS > 256 && xvasp.NCPUS <= 512) {NPAR=32;NCORE=xvasp.NCPUS/NPAR;}
    if(xvasp.NCPUS > 512 && xvasp.NCPUS <=1024) {NPAR=64;NCORE=xvasp.NCPUS/NPAR;}
    if(xvasp.NCPUS >1024) {NPAR=128;NCORE=xvasp.NCPUS/NPAR;}
    if(NPAR==0) {NPAR=4;NCORE=1;}
    if(xvasp.NCPUS==32)  {NPAR=4;NCORE=32;} // test for DX conrad and gordon
    if(xvasp.NCPUS==48)  {NPAR=4;NCORE=48;} // test for DX gaffney, koehr, and mustang
    if(xvasp.NCPUS==44)  {NPAR=4;NCORE=44;} // test for DX onyx

    // marylou fulton super computer center
    if(xvasp.NCPUS==4)  { // best 4 times
      //  NCPUS04 NPAR01 NCORE04 18576.588 20.6407
      //  NCPUS04 NPAR01 NCORE03 18289.541 20.3217
      //  NCPUS04 NPAR01 NCORE01 18223.492 20.2483
      //  NCPUS04 NPAR01 NCORE02 18176.841 20.1965
      NPAR=1;NCORE=4;  // trying to make NPAR*NCORE=NCPUS 
    }
    if(xvasp.NCPUS==6)  { // best 4 times
      //  NCPUS06 NPAR03 NCORE02 12758.774 21.2646
      //  NCPUS06 NPAR03 NCORE04 12724.680 21.2078
      //  NCPUS06 NPAR03 NCORE05 12678.754 21.1313
      //  NCPUS06 NPAR03 NCORE01 12669.103 21.1152
      NPAR=3;NCORE=2;  // trying to make NPAR*NCORE=NCPUS 
    }
    if(xvasp.NCPUS==8)  { // best 4 times   
      //  NCPUS08 NPAR02 NCORE06 10757.264 23.905
      //  NCPUS08 NPAR02 NCORE04 10731.053 23.8468
      //  NCPUS08 NPAR02 NCORE07 10728.917 23.842
      //  NCPUS08 NPAR02 NCORE01 10707.747 23.795
      NPAR=2;NCORE=4;  // trying to make NPAR*NCORE=NCPUS if possible
    }
    if(xvasp.NCPUS==12)  { // best times
      //  NCPUS12 NPAR03 NCORE07 7256.892 24.1896
      //  NCPUS12 NPAR03 NCORE11 7253.081 24.1769
      //  NCPUS12 NPAR03 NCORE04 7251.042 24.1701
      //  NCPUS12 NPAR03 NCORE09 7248.894 24.163
      //  NCPUS12 NPAR03 NCORE05 7245.983 24.1533
      //  NCPUS12 NPAR03 NCORE06 7242.303 24.141
      //  NCPUS12 NPAR03 NCORE02 7239.272 24.1309
      //  NCPUS12 NPAR03 NCORE12 7239.034 24.1301
      //  NCPUS12 NPAR03 NCORE03 7234.666 24.1156
      //  NCPUS12 NPAR03 NCORE10 7230.139 24.1005
      //  NCPUS12 NPAR03 NCORE08 7175.204 23.9173
      //  NCPUS12 NPAR03 NCORE01 7169.198 23.8973
      NPAR=3;NCORE=2; // trying to make NPAR*NCORE=NCPUS if possible
    }
    if(xvasp.NCPUS==14)  { // best times
      // NCPUS14 NPAR02 NCORE13 6647.652 25.852
      // NCPUS14 NPAR02 NCORE04 6625.906 25.7674
      // NCPUS14 NPAR02 NCORE01 6622.861 25.7556
      // NCPUS14 NPAR02 NCORE12 6606.536 25.6921
      // NCPUS14 NPAR02 NCORE14 6605.233 25.687
      // NCPUS14 NPAR02 NCORE06 6602.645 25.677
      // NCPUS14 NPAR02 NCORE02 6583.709 25.6033
      // NCPUS14 NPAR02 NCORE08 6581.525 25.5948
      // NCPUS14 NPAR02 NCORE05 6577.889 25.5807
      // NCPUS14 NPAR02 NCORE03 6573.488 25.5636
      // NCPUS14 NPAR02 NCORE07 6562.956 25.5226
      // NCPUS14 NPAR02 NCORE10 6536.529 25.4198
      // NCPUS14 NPAR02 NCORE11 6518.928 25.3514
      NPAR=2;NCORE=7; // trying to make NPAR*NCORE=NCPUS if possible
    }
    if(xvasp.NCPUS==16)  { // best times
      //  NCPUS16 NPAR04 NCORE04 6423.776 28.5501
      //  NCPUS16 NPAR04 NCORE02 6418.438 28.5264
      //  NCPUS16 NPAR04 NCORE13 6416.928 28.5197
      //  NCPUS16 NPAR04 NCORE05 6416.292 28.5169
      //  NCPUS16 NPAR04 NCORE14 6415.565 28.5136
      //  NCPUS16 NPAR04 NCORE10 6407.344 28.4771
      //  NCPUS16 NPAR04 NCORE01 6405.943 28.4709
      //  NCPUS16 NPAR04 NCORE16 6403.422 28.4597
      //  NCPUS16 NPAR04 NCORE06 6323.368 28.1039
      //  NCPUS16 NPAR04 NCORE08 6323.029 28.1024
      NPAR=4;NCORE=8; // trying to make NPAR*NCORE=NCPUS if possible
    }
    if(xvasp.NCPUS==24)  { // best times
      //  NCPUS24 NPAR03 NCORE02 5293.432 35.2895
      //  NCPUS24 NPAR03 NCORE05 5286.059 35.2404
      //  NCPUS24 NPAR03 NCORE12 5277.060 35.1804
      //  NCPUS24 NPAR03 NCORE09 5274.808 35.1654
      //  NCPUS24 NPAR03 NCORE17 5267.672 35.1178
      //  NCPUS24 NPAR03 NCORE15 5265.254 35.1017
      //  NCPUS24 NPAR03 NCORE01 5262.930 35.0862
      //  NCPUS24 NPAR03 NCORE23 5254.012 35.0267
      //  NCPUS24 NPAR03 NCORE06 5247.212 34.9814
      //  NCPUS24 NPAR03 NCORE24 5219.655 34.7977
      //  NCPUS24 NPAR03 NCORE07 5185.569 34.5705
      NPAR=3;NCORE=12;  // trying to make NPAR*NCORE=NCPUS if possible
    }

    if(xvasp.NCPUS==28)  { // best times  (delay %)
      //  NCPUS28 NPAR04 NCORE12 4817.701 37.471  1.97962
      //  NCPUS28 NPAR04 NCORE07 4809.654 37.4084 1.80929
      //  NCPUS28 NPAR04 NCORE13 4799.799 37.3318 1.60068
      //  NCPUS28 NPAR04 NCORE05 4799.391 37.3286 1.59204
      //  NCPUS28 NPAR04 NCORE20 4797.978 37.3176 1.56213
      //  NCPUS28 NPAR04 NCORE22 4794.220 37.2884 1.48259
      //  NCPUS28 NPAR04 NCORE01 4785.839 37.2232 1.30518
      //  NCPUS28 NPAR04 NCORE18 4769.633 37.0971 0.962135
      //  NCPUS28 NPAR04 NCORE23 4767.426 37.08   0.915418
      //  NCPUS28 NPAR04 NCORE06 4724.180 36.7436
      NPAR=4;NCORE=6;  // trying to make NPAR*NCORE=NCPUS if possible
    }

    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_EOS") && xvasp.NCPUS==32)  {
      xvasp.INCAR << aurostd::PaddedPOST("# Override for MACHINE::MPCDF_EOS and xvasp.NCPUS==32",_incarpad_) << endl; 
      NPAR=4;NCORE=32;
    }
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_DRACO") && xvasp.NCPUS==32)  {
      xvasp.INCAR << aurostd::PaddedPOST("# Override for MACHINE::MPCDF_DRACO and xvasp.NCPUS==32",_incarpad_) << endl; 
      NPAR=4;NCORE=32;
    }
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_COBRA") && xvasp.NCPUS==32)  {
      xvasp.INCAR << aurostd::PaddedPOST("# Override for MACHINE::MPCDF_COBRA and xvasp.NCPUS==32",_incarpad_) << endl; 
      NPAR=4;NCORE=32;
    }
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_COBRA") && xvasp.NCPUS==40)  {
      xvasp.INCAR << aurostd::PaddedPOST("# Override for MACHINE::MPCDF_COBRA and xvasp.NCPUS==40",_incarpad_) << endl; 
      NPAR=4;NCORE=40;
    }
    if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_HYDRA") && xvasp.NCPUS==32)  {
      xvasp.INCAR << aurostd::PaddedPOST("# Override for MACHINE::MPCDF_HYDRA and xvasp.NCPUS==32",_incarpad_) << endl; 
      NPAR=4;NCORE=32;
    }

    //  cerr << "xvasp.NCPUS=" << xvasp.NCPUS << endl;;
    // cerr << "NPAR=" << NPAR << endl;

    xvasp.INCAR << aurostd::PaddedPOST("NPAR="+aurostd::utype2string(NPAR),_incarpad_) << "# lookup table from NCPUS=" << xvasp.NCPUS << endl;    // 6 times faster than NPAR=cpus.
    xvasp.INCAR << aurostd::PaddedPOST("NCORE="+aurostd::utype2string(NCORE),_incarpad_) << "# lookup table from NCPUS=" << xvasp.NCPUS << endl;    
    //  xvasp.INCAR << aurostd::PaddedPOST("NPAR="+aurostd::utype2string(NPAR),_incarpad_) << "# seems to give good results with all MPI" << endl;    // 6 times faster than NPAR=cpus.
    //  xvasp.INCAR << aurostd::PaddedPOST("NPAR=4",_incarpad_) << "# seems to give good results with all MPI" << endl;    // 6 times faster than NPAR=cpus..
    xvasp.INCAR << aurostd::PaddedPOST("NSIM=4",_incarpad_) << "# tune size vector/matrix product" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("LSCALU=.FALSE.",_incarpad_) << "# neglect scalapack " << endl;
    if(VERBOSE) xvasp.INCAR << "# Performing MPI-AUTOTUNE [AFLOW] - end" << endl;
  }
}


// ***************************************************************************
// KBIN::XVASP_INCAR_System_Auto
// ***************************************************************************
namespace KBIN {
  void XVASP_INCAR_System_Auto(_xvasp& xvasp,bool VERBOSE) {        // AFLOW_FUNCTION_IMPLEMENTATION
    string FileContent,strline;
    FileContent=xvasp.INCAR.str();
    xvasp.INCAR.str(std::string());
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    int imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i);
      if(aurostd::substring2bool(strline,"SYSTEM",TRUE)    || aurostd::substring2bool(strline,"#SYSTEM",TRUE) ||
          aurostd::substring2bool(strline,"INFO",TRUE)      || aurostd::substring2bool(strline,"#INFO",TRUE) ||
          aurostd::substring2bool(strline,"PROTOTYPE",TRUE) || aurostd::substring2bool(strline,"#PROTOTYPE",TRUE)) {
        if(VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_System_Auto) " << endl;
      } else {
        if(!VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
        if(VERBOSE) xvasp.INCAR << strline << endl;
      }
    }
    FileContent=xvasp.INCAR.str();
    xvasp.INCAR.str(std::string());
    xvasp.INCAR << "#AFLOW INCAR automatically generated" << endl;
    xvasp.INCAR << "SYSTEM=" << xvasp.str.title << endl;
    xvasp.INCAR << "#PROTOTYPE=" << xvasp.str.prototype << endl;
    xvasp.INCAR << "#INFO=" << xvasp.str.info << endl;
    xvasp.INCAR << FileContent;
  }
}

// ***************************************************************************
// KBIN::XVASP_INCAR_Relax_ON
// ***************************************************************************
namespace KBIN {
  void XVASP_INCAR_Relax_ON(_xvasp& xvasp,_vflags& vflags,int number) {        // AFLOW_FUNCTION_IMPLEMENTATION
    string FileContent,strline;
    int imax,isif=3;
    if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS_CELL_VOLUME")==FALSE) {  // whatever is the number
      if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS")==TRUE  && vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_SHAPE")==FALSE && vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_VOLUME")==FALSE) isif=2; // 0,1,2 FIX
      if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS")==TRUE  && vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_SHAPE")==TRUE  && vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_VOLUME")==TRUE)  isif=3;
      if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS")==TRUE  && vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_SHAPE")==TRUE  && vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_VOLUME")==FALSE) isif=4;
      if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS")==FALSE && vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_SHAPE")==TRUE  && vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_VOLUME")==FALSE) isif=5;
      if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS")==FALSE && vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_SHAPE")==TRUE  && vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_VOLUME")==TRUE)  isif=6;
      if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS")==FALSE && vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_SHAPE")==FALSE && vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("CELL_VOLUME")==TRUE)  isif=7;
      if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("ALL")) isif=3;
      if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS_CELL_SHAPE")) isif=4;//AS20201208
    }
    if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS_CELL_VOLUME")) {  // whatever is the number
      if(aurostd::_isodd(number)) isif=7;  // VOLUME
      if(aurostd::_iseven(number)) isif=2; // IONS
    }

    FileContent=xvasp.INCAR.str();
    xvasp.INCAR.str(std::string());
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    imax=aurostd::GetNLinesString(FileContent);
    // RELAX_MODE=ENERGY mode
    if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.xscheme=="ENERGY") {
      for(int i=1;i<=imax;i++) {
        strline=aurostd::GetLineString(FileContent,i);
        if(aurostd::substring2bool(strline,"IBRION",TRUE) || aurostd::substring2bool(strline,"#IBRION",TRUE) ||
            aurostd::substring2bool(strline,"NSW",TRUE)    || aurostd::substring2bool(strline,"#NSW",TRUE)    ||
            aurostd::substring2bool(strline,"LORBIT",TRUE)    || aurostd::substring2bool(strline,"#LORBIT",TRUE)    || //CO20180130 get spinD
            aurostd::substring2bool(strline,"ISIF",TRUE)   || aurostd::substring2bool(strline,"#ISIF",TRUE)) {
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_Relax_ON) " << endl;
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }
      // xvasp.INCAR << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing RELAX=ENERGY [AFLOW] begin" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("IBRION=2",_incarpad_) << "# relax with Conjugate Gradient  " << endl;
      xvasp.INCAR << aurostd::PaddedPOST("NSW=51",_incarpad_) << "# relax for long              " << endl;
      xvasp.INCAR << aurostd::PaddedPOST("ISIF="+aurostd::utype2string(isif),_incarpad_) << "# relax appropriately         " << endl;
      xvasp.INCAR << aurostd::PaddedPOST("LORBIT=10",_incarpad_) << "# get spin decomposition      " << endl; //CO20180130 get spinD
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing RELAX=ENERGY [AFLOW] end " << endl;
    }
    // RELAX=MODE=FORCES mode
    if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_MODE.xscheme=="FORCES") {
      for(int i=1;i<=imax;i++) {
        strline=aurostd::GetLineString(FileContent,i);
        if(aurostd::substring2bool(strline,"IBRION",TRUE) || aurostd::substring2bool(strline,"#IBRION",TRUE) ||
            aurostd::substring2bool(strline,"NELMIN",TRUE)    || aurostd::substring2bool(strline,"#NELMIN",TRUE)    ||
            aurostd::substring2bool(strline,"ADDGRID",TRUE)    || aurostd::substring2bool(strline,"#ADDGRID",TRUE)    ||
            aurostd::substring2bool(strline,"EDIFFG",TRUE)    || aurostd::substring2bool(strline,"#EDIFFG",TRUE)    ||
            aurostd::substring2bool(strline,"LORBIT",TRUE)    || aurostd::substring2bool(strline,"#LORBIT",TRUE)    || //CO20180130 get spinD
            aurostd::substring2bool(strline,"NSW",TRUE)    || aurostd::substring2bool(strline,"#NSW",TRUE)    ||
            aurostd::substring2bool(strline,"ISIF",TRUE)   || aurostd::substring2bool(strline,"#ISIF",TRUE)) {
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_Relax_ON) " << endl;
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }
      // xvasp.INCAR << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing RELAX=FORCES [AFLOW] begin" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("NELMIN=4",_incarpad_) << "# The forces have to be well converged " << endl;
      xvasp.INCAR << aurostd::PaddedPOST("ADDGRID=.TRUE.",_incarpad_) << "# To support finer forces calculation " << endl;
      xvasp.INCAR << aurostd::PaddedPOST("EDIFFG="+aurostd::utype2string(DEFAULT_VASP_PREC_EDIFFG,12),_incarpad_) << "# The final structure has to have zero forces! " << endl;
      //    xvasp.INCAR << aurostd::PaddedPOST("EDIFFG=-1E-5",_incarpad_) << "# The final structure has to have zero forces! " << endl;
      xvasp.INCAR << aurostd::PaddedPOST("IBRION=1",_incarpad_) << "# More stable algorithm " << endl;
      xvasp.INCAR << aurostd::PaddedPOST("NSW=100",_incarpad_) << "# relax for very long " << endl;
      // xvasp.INCAR << "ISIF=3         # relax everything " << endl;
      xvasp.INCAR << aurostd::PaddedPOST("ISIF="+aurostd::utype2string(isif),_incarpad_)     << "# relax appropriately         " << endl;
      xvasp.INCAR << aurostd::PaddedPOST("LORBIT=10",_incarpad_) << "# get spin decomposition      " << endl; //CO20180130 get spinD
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing RELAX=FORCES [AFLOW] end " << endl;
    }
    // done now write if necessary
    if(vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS_CELL_VOLUME") && number>1) {  // whatever is the number
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",FALSE);
      aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"));
    }
  }
}

// ***************************************************************************
namespace KBIN {
  void XVASP_INCAR_Relax_ON(_xvasp& xvasp,bool VERBOSE) {        // AFLOW_FUNCTION_IMPLEMENTATION
    _vflags vflags;
    vflags.KBIN_VASP_INCAR_VERBOSE=VERBOSE;
    vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("ALL");
    vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("IONS");
    vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("CELL_SHAPE");
    vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.push("CELL_VOLUME");
    vflags.KBIN_VASP_FORCE_OPTION_RELAX_TYPE.flag("IONS_CELL_VOLUME",FALSE);
    KBIN::XVASP_INCAR_Relax_ON(xvasp,vflags,1);
  }
}

// ***************************************************************************
// KBIN::XVASP_INCAR_Static_ON
namespace KBIN {
  void XVASP_INCAR_Static_ON(_xvasp& xvasp,_vflags& vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    string FileContent,strline;
    FileContent=xvasp.INCAR.str();
    xvasp.INCAR.str(std::string());
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    int imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i); // ibrion = -1,5,6,7,8 dont need to be removed
      if(aurostd::substring2bool(strline,"IBRION=0",TRUE) || aurostd::substring2bool(strline,"#IBRION=0",TRUE) || // std md
          aurostd::substring2bool(strline,"IBRION=1",TRUE) || aurostd::substring2bool(strline,"#IBRION=1",TRUE) || // quasi-newton algo to relax ions
          aurostd::substring2bool(strline,"IBRION=2",TRUE) || aurostd::substring2bool(strline,"#IBRION=2",TRUE) || // conj-grad algo to relax
          aurostd::substring2bool(strline,"IBRION=3",TRUE) || aurostd::substring2bool(strline,"#IBRION=3",TRUE) || // damping-factor algo to relax
          aurostd::substring2bool(strline,"NSW",TRUE)    || aurostd::substring2bool(strline,"#NSW",TRUE)    ||
          aurostd::substring2bool(strline,"LORBIT",TRUE)    || aurostd::substring2bool(strline,"#LORBIT",TRUE)    || //CO20180130 get spinD
          (aurostd::substring2bool(strline,"LAECHG",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry) || (aurostd::substring2bool(strline,"#LAECHG",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry)||
          (aurostd::substring2bool(strline,"LELF",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry) || (aurostd::substring2bool(strline,"#LELF",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry)||
          aurostd::substring2bool(strline,"IALGO",TRUE)    || aurostd::substring2bool(strline,"#IALGO",TRUE) || //SC from ICSD
          (aurostd::substring2bool(strline,"ALGO",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved==FALSE) ||  //SC from ICSD
          (aurostd::substring2bool(strline,"#ALGO",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved==FALSE) || //SC from ICSD
          aurostd::substring2bool(strline,"ISIF",TRUE)   || aurostd::substring2bool(strline,"#ISIF",TRUE)) {
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_Static_ON) " << endl;
      } else {
        if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
      }
    }
    // xvasp.INCAR << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing STATIC [AFLOW]" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry && vflags.KBIN_VASP_FORCE_OPTION_BADER.option) xvasp.INCAR << aurostd::PaddedPOST("LAECHG=.TRUE.",_incarpad_) << "# Performing STATIC  (Bader ON)" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry && !vflags.KBIN_VASP_FORCE_OPTION_BADER.option) xvasp.INCAR << aurostd::PaddedPOST("LAECHG=.FALSE.",_incarpad_) << "# Performing STATIC  (Bader OFF)" << endl;
    // if(vflags.KBIN_VASP_FORCE_OPTION_BADER.option) xvasp.INCAR << aurostd::PaddedPOST("ADDGRID=.TRUE.",_incarpad_) << "# Performing STATIC  (Bader ON)" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry && vflags.KBIN_VASP_FORCE_OPTION_ELF.option) xvasp.INCAR << aurostd::PaddedPOST("LELF=.TRUE.",_incarpad_) << "# Performing STATIC  (Elf ON)" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry && !vflags.KBIN_VASP_FORCE_OPTION_ELF.option) xvasp.INCAR << aurostd::PaddedPOST("LELF=.FALSE.",_incarpad_) << "# Performing STATIC  (Elf OFF)" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved==FALSE) xvasp.INCAR << "ALGO=Normal      # Performing STATIC" << endl;
    xvasp.INCAR << "LORBIT=10        # Performing STATIC" << endl; //CO20180130 get spinD
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "#removing IBRION, NSW, ISIF" << endl;
  }
}


// ***************************************************************************
// KBIN::XVASP_INCAR_Relax_Static_ON
namespace KBIN {
  void XVASP_INCAR_Relax_Static_ON(_xvasp& xvasp,_vflags& vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    string FileContent,strline;
    FileContent=xvasp.INCAR.str();
    xvasp.INCAR.str(std::string());
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    int imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i);
      if(
          ((aurostd::substring2bool(strline,"PREC",TRUE) && !aurostd::substring2bool(strline,"SYMPREC",TRUE)) && vflags.KBIN_VASP_FORCE_OPTION_PREC.preserved==FALSE) || //CO20171003 - don't confuse PREC and SYMPREC
          (aurostd::substring2bool(strline,"#PREC",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_PREC.preserved==FALSE) || // Aleksey
          aurostd::substring2bool(strline,"EMIN",TRUE) || aurostd::substring2bool(strline,"#EMIN",TRUE) || // Aleksey
          aurostd::substring2bool(strline,"EMAX",TRUE) || aurostd::substring2bool(strline,"#EMAX",TRUE) || // Aleksey
          aurostd::substring2bool(strline,"NEDOS",TRUE) || aurostd::substring2bool(strline,"#NEDOS",TRUE) || // Aleksey
          aurostd::substring2bool(strline,"IALGO",TRUE)    || aurostd::substring2bool(strline,"#IALGO",TRUE) || //SC from ICSD
          (aurostd::substring2bool(strline,"ALGO",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved==FALSE) ||  //SC from ICSD
          (aurostd::substring2bool(strline,"#ALGO",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved==FALSE) || //SC from ICSD
          aurostd::substring2bool(strline,"ISMEAR",TRUE) || aurostd::substring2bool(strline,"#ISMEAR",TRUE) ||
          aurostd::substring2bool(strline,"SIGMA",TRUE) || aurostd::substring2bool(strline,"#SIGMA",TRUE) ||
          aurostd::substring2bool(strline,"ISIF",TRUE) || aurostd::substring2bool(strline,"#ISIF",TRUE) ||
          aurostd::substring2bool(strline,"IBRION",TRUE) || aurostd::substring2bool(strline,"#IBRION",TRUE) ||
          aurostd::substring2bool(strline,"NSW",TRUE) || aurostd::substring2bool(strline,"#NSW",TRUE) ||
          aurostd::substring2bool(strline,"NELM",TRUE) || aurostd::substring2bool(strline,"#NELM",TRUE) ||
          aurostd::substring2bool(strline,"NELMIN",TRUE) || aurostd::substring2bool(strline,"#NELMIN",TRUE) ||
          aurostd::substring2bool(strline,"LORBIT",TRUE) || aurostd::substring2bool(strline,"#LORBIT",TRUE) ||
          aurostd::substring2bool(strline,"RWIGS",TRUE) || aurostd::substring2bool(strline,"#RWIGS",TRUE) ||
          aurostd::substring2bool(strline,"LCHARG",TRUE) || aurostd::substring2bool(strline,"#LCHARG",TRUE) ||
          (aurostd::substring2bool(strline,"LAECHG",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry) || (aurostd::substring2bool(strline,"#LAECHG",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry)||
          aurostd::substring2bool(strline,"LWAVE",TRUE) || aurostd::substring2bool(strline,"#LWAVE",TRUE) ||
          (aurostd::substring2bool(strline,"LELF",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry) || (aurostd::substring2bool(strline,"#LELF",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry)||
          aurostd::substring2bool(strline,"LPARD",TRUE) || aurostd::substring2bool(strline,"#LPARD",TRUE) ||
          aurostd::substring2bool(strline,"IBAND",TRUE) || aurostd::substring2bool(strline,"#IBAND",TRUE) ||
          aurostd::substring2bool(strline,"LSEPB",TRUE) || aurostd::substring2bool(strline,"#LSEPB",TRUE)) {
            if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_Relax_Static_ON)" << endl;
          } else {
            if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
            if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
          }
    }
    // xvasp.INCAR << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing RELAX_STATIC [AFLOW]" << endl;
    // xvasp.INCAR << aurostd::PaddedPOST("ISMEAR=0",_incarpad_) << "# Performing RELAX_STATIC" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("ISMEAR=-5",_incarpad_) << "# Performing RELAX_STATIC (tetrahedron method with Blochl corrections : wahyu/stefano 2009/10/09)" << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << aurostd::PaddedPOST("",_incarpad_) << "# with -5 forces are not correct in metals, so it is good " << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << aurostd::PaddedPOST("",_incarpad_) << "# for semiconductors/insulator or for NO relaxations" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("SIGMA=0.05",_incarpad_) << "# Performing RELAX_STATIC (so the DOS will not spill too much from band edges)" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved==FALSE) xvasp.INCAR << aurostd::PaddedPOST("ALGO=Normal",_incarpad_) << "# Performing RELAX_STATIC (stefano from ICSD)" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("ISIF=2",_incarpad_) << "# Performing RELAX_STATIC (relax ions, but single step so no relax)" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("IBRION=2",_incarpad_) << "# Performing RELAX_STATIC (relax conj. grad, but single step so no relax)" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("NSW=0",_incarpad_) << "# Performing RELAX_STATIC (zero ionic steps, so no relax, just static)" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("NELM=120",_incarpad_) << "# Performing RELAX_STATIC" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("NELMIN=2",_incarpad_) << "# Performing RELAX_STATIC" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("LORBIT=10",_incarpad_) << "# Performing RELAX_STATIC" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("LCHARG=.TRUE.",_incarpad_) << "# Performing RELAX_STATIC" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry && vflags.KBIN_VASP_FORCE_OPTION_BADER.option) xvasp.INCAR << aurostd::PaddedPOST("LAECHG=.TRUE.",_incarpad_) << "# Performing RELAX_STATIC  (Bader ON)" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry && !vflags.KBIN_VASP_FORCE_OPTION_BADER.option) xvasp.INCAR << aurostd::PaddedPOST("LAECHG=.FALSE.",_incarpad_) << "# Performing RELAX_STATIC  (Bader OFF)" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry && vflags.KBIN_VASP_FORCE_OPTION_ELF.option) xvasp.INCAR << aurostd::PaddedPOST("LELF=.TRUE.",_incarpad_) << "# Performing RELAX_STATIC  (Elf ON)" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry && !vflags.KBIN_VASP_FORCE_OPTION_ELF.option) xvasp.INCAR << aurostd::PaddedPOST("LELF=.FALSE.",_incarpad_) << "# Performing RELAX_STATIC  (Elf OFF)" << endl;
    // if(vflags.KBIN_VASP_RUN_DIELECTRIC_STATIC || vflags.KBIN_VASP_RUN_DIELECTRIC_DYNAMIC) {
    //   xvasp.INCAR << aurostd::PaddedPOST("LWAVE=.TRUE.",_incarpad_) << "# Performing RELAX_STATIC (forced TRUE by DIELECTRIC)" << endl;
    // } else {
    //   xvasp.INCAR << aurostd::PaddedPOST("LWAVE=.FALSE.",_incarpad_) << "# Performing RELAX_STATIC" << endl;
    // }
    xvasp.INCAR << aurostd::PaddedPOST("LWAVE=.FALSE.",_incarpad_) << "# Performing RELAX_STATIC" << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << aurostd::PaddedPOST("#LPARD=.TRUE.",_incarpad_) << "# Performing RELAX_STATIC" << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << aurostd::PaddedPOST("#IBAND=8 9",_incarpad_) << "# Performing RELAX_STATIC" << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << aurostd::PaddedPOST("#LSEPB=.TRUE.",_incarpad_) << "# Performing RELAX_STATIC" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_PREC.preserved==FALSE) xvasp.INCAR << aurostd::PaddedPOST("PREC=ACCURATE",_incarpad_) << "# Performing RELAX_STATIC (aleksey)" << endl;
    // xvasp.INCAR << aurostd::PaddedPOST("EMIN= -14.0",_incarpad_) << "# Performing RELAX_STATIC (aleksey)" << endl;
    // xvasp.INCAR << aurostd::PaddedPOST("EMAX=  12.0",_incarpad_) << "# Performing RELAX_STATIC (aleksey)" << endl;
    // xvasp.INCAR << aurostd::PaddedPOST("EMIN= -35.0",_incarpad_) << "# Performing RELAX_STATIC (aleksey)" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("EMIN= -30.0",_incarpad_) << "# Performing RELAX_STATIC (aleksey) force search for EMIN" << endl; // was 35
    xvasp.INCAR << aurostd::PaddedPOST("EMAX=  45.0",_incarpad_) << "# Performing RELAX_STATIC (aleksey)" << endl;  // was 35
    xvasp.INCAR << aurostd::PaddedPOST("NEDOS= 5000",_incarpad_) << "# Performing RELAX_STATIC (aleksey)" << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "#adjusting ISMEAR,SIGMA,ISIF,IBRION,NSW,NELM,NELMIN,LORBIT,LCHARG,LAECHG,LWAVE,LELF,LPARD,IBAND,LSEPB" << endl;

    // check for LDA/GGA
    if(xvasp.POTCAR_TYPE=="LDA" || xvasp.POTCAR_TYPE=="GGA") {
      xvasp.INCAR << "#Fixing RWIGS/LORBIG for LDA/GGA" << endl;
      ofstream oss;KBIN::XVASP_INCAR_RWIGS_Static(xvasp,vflags,oss,ON);}
  }
}


// ***************************************************************************
// KBIN::XVASP_INCAR_Relax_Static_Bands_ON
namespace KBIN {
  void XVASP_INCAR_Relax_Static_Bands_ON(_xvasp& xvasp,_vflags& vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    string FileContent,strline;
    FileContent=xvasp.INCAR.str();
    xvasp.INCAR.str(std::string());
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    int imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i);
      if(
          ((aurostd::substring2bool(strline,"PREC",TRUE) && !aurostd::substring2bool(strline,"SYMPREC",TRUE)) && vflags.KBIN_VASP_FORCE_OPTION_PREC.preserved==FALSE) || //CO20171003 - don't confuse PREC and SYMPREC
          (aurostd::substring2bool(strline,"#PREC",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_PREC.preserved==FALSE) || // Aleksey
          aurostd::substring2bool(strline,"ISMEAR",TRUE) || aurostd::substring2bool(strline,"#ISMEAR",TRUE) ||
          aurostd::substring2bool(strline,"SIGMA",TRUE) || aurostd::substring2bool(strline,"#SIGMA",TRUE) ||
          aurostd::substring2bool(strline,"IALGO",TRUE)    || aurostd::substring2bool(strline,"#IALGO",TRUE) || //SC from ICSD
          (aurostd::substring2bool(strline,"ALGO",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved==FALSE) ||  //SC from ICSD
          (aurostd::substring2bool(strline,"#ALGO",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved==FALSE) || //SC from ICSD
          aurostd::substring2bool(strline,"ISIF",TRUE) || aurostd::substring2bool(strline,"#ISIF",TRUE) ||
          aurostd::substring2bool(strline,"IBRION",TRUE) || aurostd::substring2bool(strline,"#IBRION",TRUE) ||
          aurostd::substring2bool(strline,"NSW",TRUE) || aurostd::substring2bool(strline,"#NSW",TRUE) ||
          aurostd::substring2bool(strline,"NELM",TRUE) || aurostd::substring2bool(strline,"#NELM",TRUE) ||
          aurostd::substring2bool(strline,"NELMIN",TRUE) || aurostd::substring2bool(strline,"#NELMIN",TRUE) ||
          aurostd::substring2bool(strline,"EMIN",TRUE) || aurostd::substring2bool(strline,"#EMIN",TRUE) || // Aleksey
          aurostd::substring2bool(strline,"EMAX",TRUE) || aurostd::substring2bool(strline,"#EMAX",TRUE) || // Aleksey
          aurostd::substring2bool(strline,"NEDOS",TRUE) || aurostd::substring2bool(strline,"#NEDOS",TRUE) || // Aleksey
          aurostd::substring2bool(strline,"LORBIT",TRUE) || aurostd::substring2bool(strline,"#LORBIT (Wahyu REMOVED, CORRECT)",TRUE) ||
          aurostd::substring2bool(strline,"RWIGS",TRUE) || aurostd::substring2bool(strline,"#RWIGS",TRUE) ||
          aurostd::substring2bool(strline,"LELF",TRUE) || aurostd::substring2bool(strline,"#LELF (Wahyu REMOVED, CORRECT)",TRUE) ||
          aurostd::substring2bool(strline,"ICHARG",TRUE) || aurostd::substring2bool(strline,"#ICHARG",TRUE) ||
          aurostd::substring2bool(strline,"LCHARG",TRUE) || aurostd::substring2bool(strline,"#LCHARG",TRUE) ||
          (aurostd::substring2bool(strline,"LAECHG",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry) || (aurostd::substring2bool(strline,"#LAECHG",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry)||
          (aurostd::substring2bool(strline,"LELF",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry) || (aurostd::substring2bool(strline,"#LELF",TRUE) && vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry)||
          aurostd::substring2bool(strline,"LWAVE",TRUE) || aurostd::substring2bool(strline,"#LWAVE",TRUE)) {
            if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_Relax_Static_Bands_ON)" << endl;
          } else {
            if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
            if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
          }
    }
    // xvasp.INCAR << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing RELAX_STATIC_BANDS [AFLOW]" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_PREC.preserved==FALSE) xvasp.INCAR << aurostd::PaddedPOST("PREC=ACCURATE",_incarpad_) << "# Performing RELAX_STATIC_BANDS (aleksey)" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("ISMEAR=0",_incarpad_) << "# Performing RELAX_STATIC_BANDS (put back Gauss, use 1 if metals)" << endl;
    // xvasp.INCAR << aurostd::PaddedPOST("ISMEAR=-5",_incarpad_) << "# Performing RELAX_STATIC (tetrahedron method with Blochl corrections : WSETYAWAN+SC20091009)" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("SIGMA=0.05",_incarpad_) << "# Performing RELAX_STATIC_BANDS (so the DOS will not spill too much from band edges)" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_ALGO.preserved==FALSE) xvasp.INCAR << aurostd::PaddedPOST("ALGO=Normal",_incarpad_) << "# Performing RELAX_STATIC_BANDS (stefano from ICSD)" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("ISIF=2",_incarpad_) << "# Performing RELAX_STATIC_BANDS" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("IBRION=2",_incarpad_) << "# Performing RELAX_STATIC_BANDS" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("NSW=0 ",_incarpad_) << "# Performing RELAX_STATIC_BANDS" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("NELM=120",_incarpad_) << "# Performing RELAX_STATIC_BANDS" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("NELMIN=2",_incarpad_) << "# Performing RELAX_STATIC_BANDS" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("LORBIT=10",_incarpad_) << "# Performing RELAX_STATIC_BANDS" << endl;  //CO20180130 get spinD
    xvasp.INCAR << aurostd::PaddedPOST("ICHARG=11",_incarpad_) << "# Performing RELAX_STATIC_BANDS" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("LCHARG=.FALSE.",_incarpad_) << "# Performing RELAX_STATIC_BANDS" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry && vflags.KBIN_VASP_FORCE_OPTION_BADER.option) xvasp.INCAR << aurostd::PaddedPOST("LAECHG=.TRUE.",_incarpad_) << "# Performing RELAX_STATIC_BANDS (Bader ON)" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry && !vflags.KBIN_VASP_FORCE_OPTION_BADER.option) xvasp.INCAR << aurostd::PaddedPOST("LAECHG=.FALSE.",_incarpad_) << "# Performing RELAX_STATIC_BANDS (Bader OFF)" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry && vflags.KBIN_VASP_FORCE_OPTION_ELF.option) xvasp.INCAR << aurostd::PaddedPOST("LELF=.TRUE.",_incarpad_) << "# Performing RELAX_STATIC_BANDS (Elf ON)" << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry && !vflags.KBIN_VASP_FORCE_OPTION_ELF.option) xvasp.INCAR << aurostd::PaddedPOST("LELF=.FALSE.",_incarpad_) << "# Performing RELAX_STATIC_BANDS (Elf OFF)" << endl;

    // if(vflags.KBIN_VASP_RUN_DIELECTRIC_STATIC || vflags.KBIN_VASP_RUN_DIELECTRIC_DYNAMIC) {
    //   xvasp.INCAR << aurostd::PaddedPOST("LWAVE=.TRUE.",_incarpad_) << "# Performing RELAX_STATIC_BANDS (forced TRUE by DIELECTRIC)" << endl;
    // } else {
    //   xvasp.INCAR << aurostd::PaddedPOST("LWAVE=.FALSE.",_incarpad_) << "# Performing RELAX_STATIC_BANDS " << endl;
    // }
    xvasp.INCAR << aurostd::PaddedPOST("LWAVE=.FALSE.",_incarpad_) << "# Performing RELAX_STATIC_BANDS " << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "#adjusting ISMEAR,SIGMA,ISIF,IBRION,NSW,NELM,NELMIN,ICHARG,LCHARG,LAECHG,LWAVE" << endl;
    // cerr << "FIX KBIN::XVASP_INCAR_Relax_Static_Bands_Bands_ON(_xvasp& xvasp,bool vflags.KBIN_VASP_INCAR_VERBOSE)" << endl;

    // check for LDA/GGA
    // check for LDA/GGA
    if(xvasp.POTCAR_TYPE=="LDA" || xvasp.POTCAR_TYPE=="GGA") {
      xvasp.INCAR << "#Fixing RWIGS/LORBIG for LDA/GGA" << endl;
      ofstream oss;KBIN::XVASP_INCAR_RWIGS_Static(xvasp,vflags,oss,ON);}
  }
}


// ***************************************************************************
// KBIN::XVASP_INCAR_RWIGS_Static
namespace KBIN {
  void XVASP_INCAR_RWIGS_Static(_xvasp& xvasp,_vflags& vflags,ofstream &FileMESSAGE,bool OPERATION) {        // AFLOW_FUNCTION_IMPLEMENTATION
    // cerr << "DEBUG RWIGS" << endl;
    ostringstream aus;
    string FileContent,strline;
    int imax;
    if(OPERATION==ON)  aus << "00000  MESSAGE FORCE RWIGS_STATIC  ON : [VASP_FORCE_OPTION]RWIGS_STATIC " << " - " << Message(_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
    if(OPERATION==OFF) aus << "00000  MESSAGE FORCE RWIGS_STATIC  OFF: [VASP_FORCE_OPTION]RWIGS_STATIC " << " - " << Message(_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);

    FileContent=xvasp.INCAR.str();
    xvasp.INCAR.str(std::string());
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i);
      if(aurostd::substring2bool(strline,"RWIGS",TRUE) || aurostd::substring2bool(strline,"#RWIGS",TRUE) ||
          aurostd::substring2bool(strline,"LORBIT",TRUE)   || aurostd::substring2bool(strline,"#LORBIT",TRUE)) {
        if(vflags.KBIN_VASP_INCAR_VERBOSE && OPERATION==ON ) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_RWIGS_Static_ON) " << endl;
        if(vflags.KBIN_VASP_INCAR_VERBOSE && OPERATION==OFF) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_RWIGS_Static_OFF) " << endl;
      } else {
        if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
      }
    }
    // xvasp.INCAR << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE && OPERATION==ON ) xvasp.INCAR << "# Performing RWIGS_STATIC=ON [AFLOW]" << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE && OPERATION==OFF) xvasp.INCAR << "# Performing RWIGS_STATIC=OFF [AFLOW]" << endl;	
    if(OPERATION==ON) {
      xvasp.INCAR << "RWIGS=";
      vector<string> tokens,tokens2;
      aurostd::string2tokens(xvasp.POTCAR.str(),tokens,"\n");
      for(uint i=0;i<tokens.size();i++)
        if(aurostd::substring2bool(tokens.at(i),"RWIGS")) {
          aurostd::string2tokens(tokens.at(i),tokens2);
          xvasp.INCAR << tokens2.at(5) << " ";
        }
      xvasp.INCAR << "# Performing Density of states" << endl;
      if(xvasp.POTCAR_TYPE!="LDA" && xvasp.POTCAR_TYPE!="GGA")
        xvasp.INCAR << aurostd::PaddedPOST("LORBIT=0",_incarpad_) << "# Performing Density of states (PAW)" << endl;
      if(xvasp.POTCAR_TYPE=="LDA" || xvasp.POTCAR_TYPE=="GGA")
        xvasp.INCAR << aurostd::PaddedPOST("LORBIT=1",_incarpad_) << "# Performing Density of states (LDA/GGA)" << endl;
    }
    if(OPERATION==OFF) {;} // dummy nothng to do
    if(vflags.KBIN_VASP_INCAR_VERBOSE && OPERATION==ON ) xvasp.INCAR << "# adjusting RWIGS/LORBIT " << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE && OPERATION==OFF) xvasp.INCAR << "# removing RWIGS/LORBIT " << endl;

    if(vflags.KBIN_VASP_FORCE_OPTION_SPIN.option) {;} // dummy
  }
}


// double round(double x,double eps) {
// double out=x;
// out=out/eps;
// out=round(out)*eps;
// return out;
// }

// ***************************************************************************
// KBIN::XVASP_INCAR_Precision
namespace KBIN {
  void XVASP_INCAR_Precision(_xvasp& xvasp,_vflags& vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    string FileContent,strline;
    int imax;
    //  cerr << xvasp.INCAR.str() << endl;
    FileContent=xvasp.INCAR.str();
    xvasp.INCAR.str(std::string());
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i);
      if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="LOW" || vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="MEDIUM" || vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="NORMAL") {
        if((aurostd::substring2bool(strline,"PREC",TRUE) && !aurostd::substring2bool(strline,"SYMPREC",TRUE)) || aurostd::substring2bool(strline,"#PREC",TRUE) || //CO20171003 - don't confuse PREC and SYMPREC
            aurostd::substring2bool(strline,"ENMAX",TRUE) || aurostd::substring2bool(strline,"#ENMAX",TRUE)) {
          if(vflags.KBIN_VASP_INCAR_VERBOSE) {
            if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="LOW") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR PREC=LOW) " << endl;
            if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="MEDIUM") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR PREC=MEDIUM) " << endl;
            if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="NORMAL") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR PREC=NORMAL) " << endl;
            if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="HIGH") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR PREC=HIGH) " << endl;
            if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="ACCURATE") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR PREC=ACCURATE) " << endl;
          }
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }
      if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="ACCURATE" || vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="HIGH") {
        if(((aurostd::substring2bool(strline,"PREC",TRUE) && !aurostd::substring2bool(strline,"SYMPREC",TRUE)) || aurostd::substring2bool(strline,"#PREC",TRUE) ||  //CO20171003 - don't confuse PREC and SYMPREC
              aurostd::substring2bool(strline,"ENMAX",TRUE) || aurostd::substring2bool(strline,"#ENMAX",TRUE) ||
              aurostd::substring2bool(strline,"LREAL",TRUE) || aurostd::substring2bool(strline,"#LREAL",TRUE) ||
              aurostd::substring2bool(strline,"EDIFF",TRUE) || aurostd::substring2bool(strline,"#EDIFF",TRUE) ||
              aurostd::substring2bool(strline,"ALGO",TRUE) || aurostd::substring2bool(strline,"#ALGO",TRUE)) &&
            (!aurostd::substring2bool(strline,"IALGO",TRUE)  && !aurostd::substring2bool(strline,"#IALGO",TRUE)) &&
            (!aurostd::substring2bool(strline,"EDIFFG",TRUE)  && !aurostd::substring2bool(strline,"#EDIFFG",TRUE))) {
          if(vflags.KBIN_VASP_INCAR_VERBOSE) {
            if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="HIGH") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR PREC=HIGH) " << endl;
            if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="ACCURATE") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR PREC=ACCURATE) " << endl;
          }
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }
      //BEGIN JJPR
      if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="PHONONS") {
        if(((aurostd::substring2bool(strline,"PREC",TRUE) && !aurostd::substring2bool(strline,"SYMPREC",TRUE)) || aurostd::substring2bool(strline,"#PREC",TRUE) ||  //CO20171003 - don't confuse PREC and SYMPREC
              aurostd::substring2bool(strline,"ENMAX",TRUE) || aurostd::substring2bool(strline,"#ENMAX",TRUE) ||
              aurostd::substring2bool(strline,"LREAL",TRUE) || aurostd::substring2bool(strline,"#LREAL",TRUE) ||
              aurostd::substring2bool(strline,"EDIFF",TRUE) || aurostd::substring2bool(strline,"#EDIFF",TRUE) ||
              aurostd::substring2bool(strline,"ALGO",TRUE) || aurostd::substring2bool(strline,"#ALGO",TRUE)) &&
            (!aurostd::substring2bool(strline,"IALGO",TRUE)  && !aurostd::substring2bool(strline,"#IALGO",TRUE)) &&
            (!aurostd::substring2bool(strline,"EDIFFG",TRUE)  && !aurostd::substring2bool(strline,"#EDIFFG",TRUE))) {
          if(vflags.KBIN_VASP_INCAR_VERBOSE) {
            if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="PHONONS") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR PREC=PHONONS) " << endl;
          }
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }
      //END JJPR 
    }
    // xvasp.INCAR << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing Precision PREC=" << vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme << " [AFLOW] begin " << endl;

    if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="LOW") {
      xvasp.INCAR << aurostd::PaddedPOST("PREC=Low",_incarpad_) << "# low and fast" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("ENMAX="+aurostd::utype2string(xvasp.POTCAR_ENMAX*DEFAULT_VASP_PREC_ENMAX_LOW,_IVASP_DOUBLE2STRING_PRECISION_),_incarpad_) << "# " << DEFAULT_VASP_PREC_ENMAX_LOW << "*ENMAX (" << xvasp.POTCAR_ENMAX << ") of pseudopotentials " << endl;
    };
    if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="MEDIUM") {
      xvasp.INCAR << aurostd::PaddedPOST("PREC=Medium",_incarpad_) << "# medium resonable" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("ENMAX="+aurostd::utype2string(xvasp.POTCAR_ENMAX*DEFAULT_VASP_PREC_ENMAX_MEDIUM,_IVASP_DOUBLE2STRING_PRECISION_),_incarpad_) << "# " << DEFAULT_VASP_PREC_ENMAX_MEDIUM << "*ENMAX (" << xvasp.POTCAR_ENMAX << ") of pseudopotentials " << endl;
    };
    if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="NORMAL") {
      xvasp.INCAR << aurostd::PaddedPOST("PREC=Normal",_incarpad_) << "# avoid wrap around errors" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("ENMAX="+aurostd::utype2string(xvasp.POTCAR_ENMAX*DEFAULT_VASP_PREC_ENMAX_NORMAL,_IVASP_DOUBLE2STRING_PRECISION_),_incarpad_) << "# " << DEFAULT_VASP_PREC_ENMAX_NORMAL << "*ENMAX (" << xvasp.POTCAR_ENMAX << ") of pseudopotentials " << endl; 
    };
    if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="HIGH") {
      xvasp.INCAR << aurostd::PaddedPOST("PREC=High",_incarpad_) << "# avoid wrap around errors" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("ENMAX="+aurostd::utype2string(xvasp.POTCAR_ENMAX*DEFAULT_VASP_PREC_ENMAX_HIGH,_IVASP_DOUBLE2STRING_PRECISION_),_incarpad_) << "# " << DEFAULT_VASP_PREC_ENMAX_HIGH << "*ENMAX (" << xvasp.POTCAR_ENMAX << ") of pseudopotentials " << endl;
    };
    if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="ACCURATE") // || vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="HIGH")
    { //CO20200106 - patching for auto-indenting
      xvasp.INCAR << aurostd::PaddedPOST("PREC=Accurate",_incarpad_) << "# avoid wrap around errors" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("ENMAX="+aurostd::utype2string(xvasp.POTCAR_ENMAX*DEFAULT_VASP_PREC_ENMAX_ACCURATE,_IVASP_DOUBLE2STRING_PRECISION_),_incarpad_) << "# " << DEFAULT_VASP_PREC_ENMAX_ACCURATE << "*ENMAX (" << xvasp.POTCAR_ENMAX << ") of pseudopotentials " << endl;
      xvasp.INCAR << aurostd::PaddedPOST("LREAL=.FALSE.",_incarpad_) << "# reciprocal space projection technique " << endl;
      xvasp.INCAR << aurostd::PaddedPOST("EDIFF=1E-6",_incarpad_) << "# high accuracy required          " << endl;
      xvasp.INCAR << aurostd::PaddedPOST("ALGO=Fast",_incarpad_) << "# fast determination of ground state " << endl;
    };
    if(vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme=="PHONONS") { 
      xvasp.INCAR << aurostd::PaddedPOST("PREC=Accurate",_incarpad_) << "# avoid wrap around errors" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("ENMAX="+aurostd::utype2string(xvasp.POTCAR_ENMAX*DEFAULT_VASP_PREC_ENMAX_ACCURATE,_IVASP_DOUBLE2STRING_PRECISION_),_incarpad_) << "# " << DEFAULT_VASP_PREC_ENMAX_ACCURATE << "*ENMAX (" << xvasp.POTCAR_ENMAX << ") of pseudopotentials " << endl;
      xvasp.INCAR << aurostd::PaddedPOST("LREAL=.FALSE.",_incarpad_) << "# reciprocal space projection technique " << endl;
      xvasp.INCAR << aurostd::PaddedPOST("EDIFF=1E-8",_incarpad_) << "# high accuracy required          " << endl;
      xvasp.INCAR << aurostd::PaddedPOST("EDIFFG=-0.001",_incarpad_) << "# high accuracy required          " << endl;
      xvasp.INCAR << aurostd::PaddedPOST("NELMIN = 4",_incarpad_) << "# The forces have to be well converged" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("NELM = 120",_incarpad_) << "# Many electronic steps" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("ADDGRID = .TRUE.",_incarpad_) << "# For finer forces" << endl;
    };
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing Precision PREC=" << vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme << " [AFLOW] end " << endl;
    //  cerr << xvasp.INCAR.str() << endl;
  }
}

// ***************************************************************************
// KBIN::XVASP_INCAR_Metagga
namespace KBIN {
  void XVASP_INCAR_Metagga(_xvasp& xvasp,_vflags& vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION   TPSS | RTPSS | M06L | MBJL | SCAN | MS0 | MS1 | MS2 | NONE
    string FileContent,strline;
    int imax;
    // cerr << xvasp.INCAR.str() << endl;
    FileContent=xvasp.INCAR.str();
    xvasp.INCAR.str(std::string());
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i);
      if(aurostd::substring2bool(strline,"METAGGA",TRUE) || aurostd::substring2bool(strline,"#METAGGA",TRUE)) {
        if(vflags.KBIN_VASP_INCAR_VERBOSE) {
          if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="TPSS") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR METAGGA=TPSS) " << endl;
          if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="RTPSS") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR METAGGA=RTPSS) " << endl;
          if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="M06L") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR METAGGA=M06L) " << endl;
          if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="MBJL") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR METAGGA=MBJL) " << endl;
          if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="SCAN") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR METAGGA=SCAN) " << endl;
          if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="MS0") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR METAGGA=MS0) " << endl;
          if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="MS1") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR METAGGA=MS1) " << endl;
          if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="MS2") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR METAGGA=MS2) " << endl;
          if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="NONE") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR METAGGA=NONE) " << endl;
        }
      } else {
        if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
      }
    }
    // xvasp.INCAR << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme!="NONE") {
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing Metagga METAGGA=" << vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme << " [AFLOW] begin " << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="TPSS") xvasp.INCAR << aurostd::PaddedPOST("METAGGA=TPSS",_incarpad_) << "# METAGGA = TPSS  " << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="RTPSS") xvasp.INCAR << aurostd::PaddedPOST("METAGGA=RTPSS",_incarpad_) << "# METAGGA = RTPSS  " << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="M06L") xvasp.INCAR << aurostd::PaddedPOST("METAGGA=M06L",_incarpad_) << "# METAGGA = M06L  " << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="MBJL") xvasp.INCAR << aurostd::PaddedPOST("METAGGA=MBJL",_incarpad_) << "# METAGGA = MBJL  " << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="SCAN") xvasp.INCAR << aurostd::PaddedPOST("METAGGA=SCAN",_incarpad_) << "# METAGGA = SCAN  " << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="MS0") xvasp.INCAR << aurostd::PaddedPOST("METAGGA=MS0",_incarpad_) << "# METAGGA = MS0  " << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="MS1") xvasp.INCAR << aurostd::PaddedPOST("METAGGA=MS1",_incarpad_) << "# METAGGA = MS1  " << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="MS2") xvasp.INCAR << aurostd::PaddedPOST("METAGGA=MS2",_incarpad_) << "# METAGGA = MS2  " << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme=="NONE") xvasp.INCAR << aurostd::PaddedPOST("METAGGA=NONE",_incarpad_) << "# METAGGA = NONE  " << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing Metagga METAGGA=" << vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme << " [AFLOW] end " << endl;
    }
  }
}

// ***************************************************************************
// KBIN::XVASP_INCAR_Ivdw
namespace KBIN {
  void XVASP_INCAR_Ivdw(_xvasp& xvasp,_vflags& vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION   number_for_VASP_see_manual_for_IVDW | 0
    string FileContent,strline;
    int imax;
    // cerr << xvasp.INCAR.str() << endl;
    FileContent=xvasp.INCAR.str();
    xvasp.INCAR.str(std::string());
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i);
      if(aurostd::substring2bool(strline,"IVDW",TRUE) || aurostd::substring2bool(strline,"#IVDW",TRUE)) {
        if(vflags.KBIN_VASP_INCAR_VERBOSE) {
          xvasp.INCAR << "# " << strline << " # AFLOW REMOVED " << endl;
        }
      } else {
        if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
      }
    }
    // xvasp.INCAR << endl;
    if(vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme!="0") {
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing Ivdw IVDW=" << vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme << " [AFLOW] begin " << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme=="0" || vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme=="OFF" ||
          vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme=="NONE" || vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme=="") {
        xvasp.INCAR << aurostd::PaddedPOST("#IVDW=0",_incarpad_) << "# IVDW = 0  " << endl;
      } else {	
        xvasp.INCAR << aurostd::PaddedPOST("IVDW="+vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme,_incarpad_) << "# IVDW = " << vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme << "  " << endl;
      }
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing Ivdw IVDW=" << vflags.KBIN_VASP_FORCE_OPTION_IVDW.xscheme << " [AFLOW] end " << endl;
    }
  }
}

// ***************************************************************************
// KBIN::XVASP_INCAR_ABMIX
namespace KBIN {
  void XVASP_INCAR_ABMIX(_xvasp& xvasp,_vflags& vflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    if(!vflags.KBIN_VASP_FORCE_OPTION_ABMIX.isentry) return;
    string FileContent,strline;
    int imax;
    // cerr << xvasp.INCAR.str() << endl;
    FileContent=xvasp.INCAR.str();
    xvasp.INCAR.str(std::string());
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i);
      if(aurostd::substring2bool(strline,"AMIX",TRUE) || aurostd::substring2bool(strline,"#AMIX",TRUE) ||
          aurostd::substring2bool(strline,"BMIX",TRUE) || aurostd::substring2bool(strline,"#BMIX",TRUE)) {
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_ABMIX)" << endl;
      } else {
        if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
      }
    }

    // xvasp.INCAR << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing KBIN::XVASP_INCAR_ABMIX [AFLOW]" << endl;
    uint isUS=FALSE,isPAW=FALSE;
    if(vflags.KBIN_VASP_FORCE_OPTION_ABMIX.xscheme.at(0)=='A') { // AUTO
      // needs to check PPs
      //   cerr << "xvasp.POTCAR_TYPE=" << xvasp.POTCAR_TYPE << endl;
      if(xvasp.POTCAR_TYPE=="LDA") isUS=TRUE;
      if(xvasp.POTCAR_TYPE=="GGA") isUS=TRUE;
      if(xvasp.POTCAR_TYPE=="PAW_LDA") isPAW=TRUE;
      if(xvasp.POTCAR_TYPE=="PAW_GGA") isPAW=TRUE;
      if(xvasp.POTCAR_TYPE=="PAW_PBE") isPAW=TRUE;
      if(xvasp.POTCAR_TYPE=="PAW_LDA_KIN") isPAW=TRUE;
      if(xvasp.POTCAR_TYPE=="PAW_PBE_KIN") isPAW=TRUE;
      if(xvasp.POTCAR_TYPE=="PAW_RPBE") isPAW=TRUE;
      if(!isUS && !isPAW) {
        isPAW=TRUE; // some default
        xvasp.INCAR << "# ABMIX=AUTO => POTCAR type not found, switching to PAW (aflow_ivasp.cpp)." << endl;   
      }      
    }
    if(vflags.KBIN_VASP_FORCE_OPTION_ABMIX.xscheme=="US" || isUS) { // US
      if(isUS) xvasp.INCAR << "# ABMIX=AUTO => Found US pseudopotentials" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("AMIX=0.5",_incarpad_) << "# Performing KBIN::XVASP_INCAR_ABMIX (US)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("BMIX=1.0",_incarpad_) << "# Performing KBIN::XVASP_INCAR_ABMIX (US)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("AMIX_MAG=2.0",_incarpad_) << "# Performing KBIN::XVASP_INCAR_ABMIX (US)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("BMIX_MAG=1.0",_incarpad_) << "# Performing KBIN::XVASP_INCAR_ABMIX (US)" << endl;
      // needs to check PPs
    }
    if(vflags.KBIN_VASP_FORCE_OPTION_ABMIX.xscheme=="PAW" || isPAW) { // PAW
      if(isPAW) xvasp.INCAR << "# ABMIX=AUTO => Found PAW pseudopotentials" << endl;   
      xvasp.INCAR << aurostd::PaddedPOST("AMIX=0.2",_incarpad_) << "# Performing KBIN::XVASP_INCAR_ABMIX (PAW)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("BMIX=0.00001",_incarpad_) << "# Performing KBIN::XVASP_INCAR_ABMIX (PAW)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("AMIX_MAG=0.8",_incarpad_) << "# Performing KBIN::XVASP_INCAR_ABMIX (PAW)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("BMIX_MAG=0.00001",_incarpad_) << "# Performing KBIN::XVASP_INCAR_ABMIX (PAW)" << endl;
      // needs to check PPs
    }
    // check for LDA/GGA
    if(xvasp.POTCAR_TYPE=="LDA" || xvasp.POTCAR_TYPE=="GGA") {
      //    xvasp.INCAR << "#Fixing RWIGS/LORBIG for LDA/GGA" << endl;
      ofstream oss;KBIN::XVASP_INCAR_RWIGS_Static(xvasp,vflags,oss,ON);
    }
  }
}


// ***************************************************************************
// KBIN::XVASP_INCAR_GetNBANDS
namespace KBIN {
  int XVASP_INCAR_GetNBANDS(_xvasp& xvasp,bool ispin) {
    vector<double> vZVAL;
    double ZVAL=GetZVAL(xvasp.POTCAR,vZVAL);
    // // cout << XPID << "00000  MESSAGE POTCAR ZVAL max   = " << ZVAL << endl;
    int nbands;
    nbands=GetNBANDS((int) ZVAL,(int) xvasp.str.atoms.size(),5,ispin);
    nbands=GetNBANDS((int) ZVAL,(int) xvasp.str.atoms.size(),5,TRUE);  // SAFETY
    // nbands=nbands+20+nbands/5; // MORE SAFETY
    // cout << XPID << "00000  MESSAGE POTCAR NBANDS = " << nbands << endl;
    if(!XHOST.QUIET) cout << XPID << "00000  MESSAGE-OPTION  [VASP_FORCE_OPTION]NBANDS  = " << nbands << endl;
    return nbands;
  }
}

// ***************************************************************************
// KBIN::XVASP_INCAR_PSTRESS
namespace KBIN {
  bool XVASP_INCAR_PREPARE_GENERIC(string command,_xvasp& xvasp,_vflags& vflags,string svalue,int ivalue,double dvalue,bool OPTION)
    //  bool XVASP_INCAR_PREPARE_GENERIC(string command,_xvasp& xvasp,_kflags& kflags,_vflags& vflags,string svalue,int ivalue,double dvalue,bool OPTION)
  { //CO20200106 - patching for auto-indenting
    bool DONE=FALSE;
    bool LDEBUG=FALSE;
    string FileContent,strline;
    FileContent=xvasp.INCAR.str();
    xvasp.INCAR.str(std::string());
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    int imax=aurostd::GetNLinesString(FileContent);

    // ***************************************************************************
    // GENERIC GENERIC GENERIC GENERIC GENERIC GENERIC
    if(command=="GENERIC") {
      DONE=TRUE;
    }

    // ***************************************************************************
    // ALGO ALGO ALGO ALGO ALGO ALGO
    if(command=="ALGO") {
      for(int i=1;i<=imax;i++) {
        strline=aurostd::GetLineString(FileContent,i);
        if(aurostd::substring2bool(strline,"ALGO",TRUE) || aurostd::substring2bool(strline,"#ALGO",TRUE) ||
            aurostd::substring2bool(strline,"IALGO",TRUE) || aurostd::substring2bool(strline,"#IALGO",TRUE)) {
          if(vflags.KBIN_VASP_INCAR_VERBOSE) {
            if(vflags.KBIN_VASP_FORCE_OPTION_ALGO.xscheme=="NORMAL") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR ALGO=NORMAL) " << endl;
            if(vflags.KBIN_VASP_FORCE_OPTION_ALGO.xscheme=="VERYFAST") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR ALGO=VERYFAST) " << endl;
            if(vflags.KBIN_VASP_FORCE_OPTION_ALGO.xscheme=="FAST") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR ALGO=FAST) " << endl;
            if(vflags.KBIN_VASP_FORCE_OPTION_ALGO.xscheme=="ALL") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR ALGO=ALL) " << endl;
            if(vflags.KBIN_VASP_FORCE_OPTION_ALGO.xscheme=="DAMPED") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR ALGO=DAMPED) " << endl;
          }
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }
      // xvasp.INCAR << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing Algo ALGO=" << vflags.KBIN_VASP_FORCE_OPTION_ALGO.xscheme << " [AFLOW] begin " << endl;

      if(vflags.KBIN_VASP_FORCE_OPTION_ALGO.xscheme=="NORMAL") xvasp.INCAR << aurostd::PaddedPOST("ALGO=Normal",_incarpad_) << "# ALGO = Normal: will select, IALGO=38 (blocked Davidson block iteration scheme)" << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_ALGO.xscheme=="VERYFAST") xvasp.INCAR << aurostd::PaddedPOST("ALGO=Veryfast",_incarpad_) << "# ALGO = Very_Fast: will select IALGO=48 (RMM-DIIS)" << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_ALGO.xscheme=="FAST") xvasp.INCAR << aurostd::PaddedPOST("ALGO=Fast",_incarpad_) << "# ALGO = Fast: IALGO=38 is used for the initial phase, and then VASP switches to IALGO=48" << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_ALGO.xscheme=="ALL") xvasp.INCAR << aurostd::PaddedPOST("ALGO=All",_incarpad_) << "# ALGO = All: The all band simultaneous update of wavefunctions can be selected using (IALGO=58)" << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_ALGO.xscheme=="DAMPED") xvasp.INCAR << aurostd::PaddedPOST("ALGO=Damped",_incarpad_) << "# ALGO = Damped: A damped velocity friction algorithm is selected (IALGO=53)" << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing Algo ALGO=" << vflags.KBIN_VASP_FORCE_OPTION_ALGO.xscheme << " [AFLOW] end " << endl;
      DONE=TRUE;
    }

    // ***************************************************************************
    // ENMAX_MULTIPLY ENMAX_MULTIPLY ENMAX_MULTIPLY ENMAX_MULTIPLY ENMAX_MULTIPLY ENMAX_MULTIPLY
    if(command=="ENMAX_MULTIPLY") {
      for(int i=1;i<=imax;i++) {
        strline=aurostd::GetLineString(FileContent,i);
        if(aurostd::substring2bool(strline,"ENMAX",TRUE) || aurostd::substring2bool(strline,"#ENMAX",TRUE)) {
          if(vflags.KBIN_VASP_INCAR_VERBOSE) {
            xvasp.INCAR << "# " << strline << " # AFLOW REMOVED " << endl;
          }
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }
      // xvasp.INCAR << endl;     
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing ENMAX_MULTIPLY=" << vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme << " [AFLOW] begin " << endl;
      xvasp.INCAR << aurostd::PaddedPOST("ENMAX="+aurostd::utype2string(xvasp.POTCAR_ENMAX*vflags.KBIN_VASP_FORCE_OPTION_ENMAX_MULTIPLY_EQUAL.content_double,_IVASP_DOUBLE2STRING_PRECISION_),_incarpad_) << "# " << vflags.KBIN_VASP_FORCE_OPTION_ENMAX_MULTIPLY_EQUAL.content_double << "*ENMAX (" << xvasp.POTCAR_ENMAX << ") of pseudopotentials " << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing ENMAX_MULTIPLY=" << vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme << " [AFLOW] end " << endl;
      //  cerr << xvasp.INCAR.str() << endl;
      DONE=TRUE;
    }

    // ***************************************************************************
    // IMIX IMIX IMIX IMIX IMIX IMIX
    if(command=="IMIX") {
      for(int i=1;i<=imax;i++) {
        strline=aurostd::GetLineString(FileContent,i);
        if(aurostd::substring2bool(strline,"IMIX",TRUE) || aurostd::substring2bool(strline,"#IMIX",TRUE)) {
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC(IMIX)) " << endl;
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }
      // xvasp.INCAR << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing IMIX=X [AFLOW] begin" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("IMIX="+aurostd::utype2string(ivalue),_incarpad_) << "# IMIX=X" << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing IMIX=X [AFLOW] end " << endl;
      DONE=TRUE;
    }

    // ***************************************************************************
    // IALGO IALGO IALGO IALGO IALGO IALGO
    if(command=="IALGO") {
      // // aus_exec << "cat INCAR | grep -v IALGO > aflow.tmp && mv aflow.tmp INCAR" << endl;
      // aus_exec << "cat INCAR | grep -v IALGO | grep -v IMIX > aflow.tmp && mv aflow.tmp INCAR" << endl;
      // aus_exec << "echo \"IALGO=48\"  >> INCAR " << endl;
      for(int i=1;i<=imax;i++) {
        strline=aurostd::GetLineString(FileContent,i);
        if(aurostd::substring2bool(strline,"IALGO",TRUE) || aurostd::substring2bool(strline,"#IALGO",TRUE) ||
            aurostd::substring2bool(strline,"ALGO",TRUE) || aurostd::substring2bool(strline,"#ALGO",TRUE) ||
            aurostd::substring2bool(strline,"IMIX",TRUE) || aurostd::substring2bool(strline,"#IMIX",TRUE)) {
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC(IALGO)) " << endl;
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }
      // xvasp.INCAR << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing IALGO=X [AFLOW] begin" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("IALGO="+aurostd::utype2string(ivalue),_incarpad_) << "# IALGOX1" << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing IALGO=X [AFLOW] end " << endl;
      DONE=TRUE;
    }

    // ***************************************************************************
    // TYPE TYPE TYPE TYPE TYPE TYPE
    if(command=="TYPE") {
      for(int i=1;i<=imax;i++) {
        strline=aurostd::GetLineString(FileContent,i);
        if(aurostd::substring2bool(strline,"ISMEAR",TRUE) || aurostd::substring2bool(strline,"#ISMEAR",TRUE) ||
            aurostd::substring2bool(strline,"SIGMA",TRUE) || aurostd::substring2bool(strline,"#SIGMA",TRUE)) {
          if(vflags.KBIN_VASP_INCAR_VERBOSE) {
            if(svalue=="DEFAULT") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC(TYPE) TYPE=" << svalue << ") " << endl;
            if(svalue=="METAL") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC(TYPE) TYPE=" << svalue << ") " << endl;
            if(svalue=="SEMICONDUCTOR") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC(TYPE) TYPE=" << svalue << ") " << endl;
            if(svalue=="INSULATOR") xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC(TYPE) TYPE=" << svalue << ") " << endl;
          }
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC(TYPE) TYPE=" << svalue << ") " << endl;
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }
      // xvasp.INCAR << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing TYPE=" << svalue << " [AFLOW] begin " << endl;
      if(svalue=="DEFAULT") {
        xvasp.INCAR << aurostd::PaddedPOST("ISMEAR=1",_incarpad_) << "# for default (as metal)" << endl;
        xvasp.INCAR << aurostd::PaddedPOST("SIGMA=0.1",_incarpad_) << "# for default (as metal)" << endl;
        //    xvasp.INCAR << aurostd::PaddedPOST("ISMEAR=2",_incarpad_) << "# default" << endl;
        //   xvasp.INCAR << aurostd::PaddedPOST("SIGMA=0.2",_incarpad_) << "# default" << endl;
      }
      if(svalue=="METAL") {
        xvasp.INCAR << aurostd::PaddedPOST("ISMEAR=1",_incarpad_) << "# for metal" << endl;
        xvasp.INCAR << aurostd::PaddedPOST("SIGMA=0.1",_incarpad_) << "# for metal" << endl;
      }
      if(svalue=="SEMICONDUCTOR" || svalue=="INSULATOR") {
        xvasp.INCAR << aurostd::PaddedPOST("ISMEAR=0",_incarpad_) << "# for insulators/semiconductors" << endl;
        xvasp.INCAR << aurostd::PaddedPOST("SIGMA=0.05",_incarpad_) << "# for insulators/semiconductors" << endl;
      }
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing TYPE=" << svalue << " [AFLOW] end " << endl;
      DONE=TRUE;
    }

    // ***************************************************************************
    // ISMEAR ISMEAR ISMEAR ISMEAR ISMEAR ISMEAR
    if(command=="ISMEAR") { //CO20181129
      for(int i=1;i<=imax;i++) {
        strline=aurostd::GetLineString(FileContent,i);
        if(aurostd::substring2bool(strline,"ISMEAR",TRUE) || aurostd::substring2bool(strline,"#ISMEAR",TRUE) || 
            aurostd::substring2bool(strline,"SIGMA",TRUE) || aurostd::substring2bool(strline,"#SIGMA",TRUE)) {
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC(ISMEAR) ISMEAR=" << aurostd::utype2string(vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_EQUAL.content_int) << ") " << endl;
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }
      // xvasp.INCAR << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing ISMEAR=" << aurostd::utype2string(vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_EQUAL.content_int) << " [AFLOW] begin " << endl;
      xvasp.INCAR << aurostd::PaddedPOST("ISMEAR="+aurostd::utype2string(vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_EQUAL.content_int),_incarpad_) << "# KBIN::XVASP_INCAR_PREPARE_GENERIC(ISMEAR)" << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing ISMEAR=" << aurostd::utype2string(vflags.KBIN_VASP_FORCE_OPTION_ISMEAR_EQUAL.content_int) << " [AFLOW] end " << endl;
      DONE=TRUE;
    }
    // ***************************************************************************

    // ***************************************************************************
    // SIGMA SIGMA SIGMA SIGMA SIGMA SIGMA
    if(command=="SIGMA") { //CO20181129
      for(int i=1;i<=imax;i++) {
        strline=aurostd::GetLineString(FileContent,i);
        if(aurostd::substring2bool(strline,"SIGMA",TRUE) || aurostd::substring2bool(strline,"#SIGMA",TRUE)) {
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC(SIGMA) SIGMA=" << aurostd::utype2string(vflags.KBIN_VASP_FORCE_OPTION_SIGMA_EQUAL.content_double,_IVASP_DOUBLE2STRING_PRECISION_) << ") " << endl;
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }
      // xvasp.INCAR << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing SIGMA=" << aurostd::utype2string(vflags.KBIN_VASP_FORCE_OPTION_SIGMA_EQUAL.content_double,_IVASP_DOUBLE2STRING_PRECISION_) << " [AFLOW] begin " << endl;
      xvasp.INCAR << aurostd::PaddedPOST("SIGMA="+aurostd::utype2string(vflags.KBIN_VASP_FORCE_OPTION_SIGMA_EQUAL.content_double,_IVASP_DOUBLE2STRING_PRECISION_),_incarpad_) << "# KBIN::XVASP_INCAR_PREPARE_GENERIC(SIGMA)" << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing SIGMA=" << aurostd::utype2string(vflags.KBIN_VASP_FORCE_OPTION_SIGMA_EQUAL.content_double,_IVASP_DOUBLE2STRING_PRECISION_) << " [AFLOW] end " << endl;
      DONE=TRUE;
    }
    // ***************************************************************************

    // PAW_CORRECTIONS PAW_CORRECTIONS PAW_CORRECTIONS PAW_CORRECTIONS PAW_CORRECTIONS PAW_CORRECTIONS
    if(command=="PAW_CORRECTIONS") {
      for(int i=1;i<=imax;i++) {
        strline=aurostd::GetLineString(FileContent,i);
        if(aurostd::substring2bool(strline,"LASPH",TRUE) || aurostd::substring2bool(strline,"#LASPH",TRUE)) {
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC(PAW_CORRECTIONS)) " << endl;
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }
      // xvasp.INCAR << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE)  xvasp.INCAR << "# Performing PAW=ON [AFLOW] begin" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("LASPH=.TRUE.",_incarpad_) << "# LASPH  http://cms.mpi.univie.ac.at/vasp/vasp/LASPH_tag.html#incar-lasph" << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE)  xvasp.INCAR << "# Performing PAW=ON [AFLOW] end " << endl;
      DONE=TRUE;
    }

    // ***************************************************************************
    // NBANDS NBANDS NBANDS NBANDS NBANDS NBANDS
    if(command=="NBANDS") {
      for(int i=1;i<=imax;i++) {
        strline=aurostd::GetLineString(FileContent,i);
        if(aurostd::substring2bool(strline,"NBANDS",TRUE) || aurostd::substring2bool(strline,"#NBANDS",TRUE)) {
          if(vflags.KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry && vflags.KBIN_VASP_INCAR_VERBOSE)
            xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC(NBANDS) <= KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry) " << endl;
          if(vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.isentry && vflags.KBIN_VASP_INCAR_VERBOSE)
            xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC(NBANDS) <= KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.isentry) " << endl;
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }
      // xvasp.INCAR << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry) {
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing NBANDS [AFLOW] begin (KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry)" << endl;
        bool ispin=TRUE;
        ispin=vflags.KBIN_VASP_FORCE_OPTION_SPIN.option;
        // xvasp.INCAR << "NBANDS=" << KBIN::XVASP_INCAR_GetNBANDS(xvasp,ispin) << " \t#NBANDS" << endl;
        xvasp.INCAR << aurostd::PaddedPOST("NBANDS="+aurostd::utype2string(KBIN::XVASP_INCAR_GetNBANDS(xvasp,ispin)),_incarpad_) << "# KBIN::XVASP_INCAR_GetNBANDS" << endl;
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing NBANDS [AFLOW] end (KBIN_VASP_FORCE_OPTION_NBANDS_AUTO_isentry)" << endl;
      }
      if(vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.isentry) {
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing NBANDS [AFLOW] begin (KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.isentry)" << endl;
        xvasp.INCAR << aurostd::PaddedPOST("NBANDS="+aurostd::utype2string(vflags.KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.content_int),_incarpad_) << "# KBIN::XVASP_INCAR_GetNBANDS" << endl;
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing NBANDS [AFLOW] end (KBIN_VASP_FORCE_OPTION_NBANDS_EQUAL.isentry)" << endl;
      }
      DONE=TRUE;
    }

    // ***************************************************************************
    // PSTRESS PSTRESS PSTRESS PSTRESS PSTRESS PSTRESS
    if(command=="PSTRESS") {
      for(int i=1;i<=imax;i++) {
        strline=aurostd::GetLineString(FileContent,i);
        if(aurostd::substring2bool(strline,"PSTRESS",TRUE) || aurostd::substring2bool(strline,"#PSTRESS",TRUE)) {
          if(vflags.KBIN_VASP_FORCE_OPTION_PSTRESS_EQUAL.isentry && vflags.KBIN_VASP_INCAR_VERBOSE)
            xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC(PSTRESS) <= KBIN_VASP_FORCE_OPTION_PSTRESS_EQUAL.isentry) " << endl;
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }
      // xvasp.INCAR << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_PSTRESS_EQUAL.isentry) {
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing PSTRESS [AFLOW] begin (KBIN_VASP_FORCE_OPTION_PSTRESS_EQUAL.isentry)" << endl;
        xvasp.INCAR << aurostd::PaddedPOST("PSTRESS="+aurostd::utype2string(vflags.KBIN_VASP_FORCE_OPTION_PSTRESS_EQUAL.content_double,_IVASP_DOUBLE2STRING_PRECISION_),_incarpad_) << "# KBIN::XVASP_INCAR_PREPARE_GENERIC(PSTRESS)" << endl;
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing PSTRESS [AFLOW] end (KBIN_VASP_FORCE_OPTION_PSTRESS_EQUAL.isentry)" << endl;
      }
      DONE=TRUE;
    }
    // ***************************************************************************
    // EDIFFG EDIFFG EDIFFG EDIFFG EDIFFG EDIFFG
    if(command=="EDIFFG") {
      for(int i=1;i<=imax;i++) {
        strline=aurostd::GetLineString(FileContent,i);
        if(aurostd::substring2bool(strline,"EDIFFG",TRUE) || aurostd::substring2bool(strline,"#EDIFFG",TRUE)) {
          if(vflags.KBIN_VASP_FORCE_OPTION_EDIFFG_EQUAL.isentry && vflags.KBIN_VASP_INCAR_VERBOSE)
            xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC(EDIFFG) <= KBIN_VASP_FORCE_OPTION_EDIFFG_EQUAL.isentry) " << endl;
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }
      // xvasp.INCAR << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_EDIFFG_EQUAL.isentry) {
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing EDIFFG [AFLOW] begin (KBIN_VASP_FORCE_OPTION_EDIFFG_EQUAL.isentry)" << endl;
        xvasp.INCAR << aurostd::PaddedPOST("EDIFFG="+aurostd::utype2string(vflags.KBIN_VASP_FORCE_OPTION_EDIFFG_EQUAL.content_double,_IVASP_DOUBLE2STRING_PRECISION_),_incarpad_) << "# KBIN::XVASP_INCAR_PREPARE_GENERIC(EDIFFG)" << endl;
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing EDIFFG [AFLOW] end (KBIN_VASP_FORCE_OPTION_EDIFFG_EQUAL.isentry)" << endl;
      }
      DONE=TRUE;
    }

    // ***************************************************************************
    // POTIM POTIM POTIM POTIM POTIM POTIM
    if(command=="POTIM") {
      for(int i=1;i<=imax;i++) {
        strline=aurostd::GetLineString(FileContent,i);
        if(aurostd::substring2bool(strline,"POTIM",TRUE) || aurostd::substring2bool(strline,"#POTIM",TRUE)) {
          if(vflags.KBIN_VASP_FORCE_OPTION_POTIM_EQUAL.isentry && vflags.KBIN_VASP_INCAR_VERBOSE)
            xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC(POTIM) <= KBIN_VASP_FORCE_OPTION_POTIM_EQUAL.isentry) " << endl;
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }
      // xvasp.INCAR << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_POTIM_EQUAL.isentry) {
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing POTIM [AFLOW] begin (KBIN_VASP_FORCE_OPTION_POTIM_EQUAL.isentry)" << endl;
        xvasp.INCAR << aurostd::PaddedPOST("POTIM="+aurostd::utype2string(vflags.KBIN_VASP_FORCE_OPTION_POTIM_EQUAL.content_double,_IVASP_DOUBLE2STRING_PRECISION_),_incarpad_) << "# KBIN::XVASP_INCAR_PREPARE_GENERIC(POTIM)" << endl;
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing POTIM [AFLOW] end (KBIN_VASP_FORCE_OPTION_POTIM_EQUAL.isentry)" << endl;
      }
      DONE=TRUE;
    }

    // ***************************************************************************
    // SPIN SPIN SPIN SPIN SPIN SPIN
    if(command=="SPIN") {
      bool MAGMOM_ALREADY_SPECIFIED=FALSE;    //CO
      for(int i=1;i<=imax;i++) {
        strline=aurostd::GetLineString(FileContent,i);
        if(aurostd::substring2bool(strline,"ISPIND",TRUE) || aurostd::substring2bool(strline,"#ISPIND",TRUE) ||
            aurostd::substring2bool(strline,"ISPIN",TRUE) || aurostd::substring2bool(strline,"#ISPIN",TRUE) ||
            ((vflags.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.isentry && vflags.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.option)    //CO
             && (aurostd::substring2bool(strline,"MAGMOM",TRUE) || aurostd::substring2bool(strline,"#MAGMOM",TRUE)))) { //CO
          if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==ON)  xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC(SPIN+ON) " << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==OFF) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC(SPIN+OFF) " << endl;
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
          //if(aurostd::substring2bool(strline,"MAGMOM")) {cerr << strline << endl;}        //CO
          if(aurostd::substring2bool(strline,"MAGMOM") && !aurostd::substring2bool(strline,"#MAGMOM")) {MAGMOM_ALREADY_SPECIFIED=TRUE;}       //CO
        }
      }
      // xvasp.INCAR << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==ON)  xvasp.INCAR << "# Performing SPIN=ON [AFLOW] begin" << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==OFF) xvasp.INCAR << "# Performing SPIN=OFF [AFLOW] begin" << endl;
      if(OPTION==ON)  xvasp.INCAR << aurostd::PaddedPOST("ISPIND=2",_incarpad_) << "# SPIN=ON" << endl;
      if(OPTION==ON)  xvasp.INCAR << aurostd::PaddedPOST("ISPIN=2",_incarpad_) << "# SPIN=ON" << endl;
      if(!MAGMOM_ALREADY_SPECIFIED && OPTION==ON && vflags.KBIN_VASP_FORCE_OPTION_LSCOUPLING.option==FALSE) {// && vflags.KBIN_VASP_FORCE_OPTION_LSCOUPLING.option_OFF==FALSE)    //CO
        if(vflags.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.isentry && vflags.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.option)
          xvasp.INCAR << aurostd::PaddedPOST("MAGMOM= "+aurostd::utype2string(xvasp.str.atoms.size())+"*5",_incarpad_) << "# " << xvasp.str.atoms.size() << " atoms " << endl;
        if(!vflags.KBIN_VASP_FORCE_OPTION_AUTO_MAGMOM.isentry)
          xvasp.INCAR << aurostd::PaddedPOST("MAGMOM= "+aurostd::utype2string(xvasp.str.atoms.size())+"*1",_incarpad_) << "# SPIN=ON - " << xvasp.str.atoms.size() << " atoms (default) " << endl;
      }
      if(OPTION==OFF) xvasp.INCAR << aurostd::PaddedPOST("ISPIN=1",_incarpad_) << "# SPIN=OFF" << endl;
      // with spin, makes the mag mom +5....
      if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==ON)  xvasp.INCAR << "# Performing SPIN=ON [AFLOW] end " << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==OFF) xvasp.INCAR << "# Performing SPIN=OFF [AFLOW] end " << endl;

      DONE=TRUE;
    }

    // ***************************************************************************
    // LS_COUPLING LS_COUPLING LS_COUPLING LS_COUPLING LS_COUPLING LS_COUPLING
    if(command=="LS_COUPLING") {
      // strip magmom from incar and make the mabmom from the order parameter
      for(int i=1;i<=imax;i++) {
        strline=aurostd::GetLineString(FileContent,i);
        if(aurostd::substring2bool(strline,"LSORBIT",TRUE) || aurostd::substring2bool(strline,"#LSORBIT",TRUE) ||
            aurostd::substring2bool(strline,"LNONCOLLINEAR",TRUE) || aurostd::substring2bool(strline,"#LNONCOLLINEAR",TRUE) ||
            aurostd::substring2bool(strline,"MAGMOM",TRUE) || aurostd::substring2bool(strline,"#MAGMOM",TRUE)) {
          if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==ON)  xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC(LS_COUPLING+ON)) " << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==OFF) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC(LS_COUPLING+OFF)) " << endl;
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }
      // xvasp.INCAR << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==ON) xvasp.INCAR << "# Performing LSCOUPLING=ON [AFLOW] begin" << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==OFF) xvasp.INCAR << "# Performing LSCOUPLING=OFF [AFLOW] begin" << endl;
      if(OPTION==ON ) xvasp.INCAR << aurostd::PaddedPOST("LSORBIT=.TRUE.",_incarpad_) << "# LSORBIT=ON" << endl;
      if(OPTION==ON ) xvasp.INCAR << aurostd::PaddedPOST("LNONCOLLINEAR=.TRUE.",_incarpad_) << "# LNONCOLLINEAR=ON" << endl;
      if(OPTION==ON ) {
        xvasp.INCAR << "MAGMOM= ";
        for(uint i=0;i<xvasp.str.atoms.size();i++) xvasp.INCAR << " 0 0 5";
        xvasp.INCAR << " \t# " << xvasp.str.atoms.size() << " 3*atoms " << endl;
      }
      if(OPTION==OFF) xvasp.INCAR << aurostd::PaddedPOST("LSORBIT=.FALSE.",_incarpad_) << "# LSORBIT=OFF" << endl;
      if(OPTION==OFF) xvasp.INCAR << aurostd::PaddedPOST("LNONCOLLINEAR=.FALSE.",_incarpad_) << "# LNONCOLLINEAR=OFF" << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==ON)  xvasp.INCAR << "# Performing LSCOUPLING=ON [AFLOW] end " << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==OFF) xvasp.INCAR << "# Performing LSCOUPLING=OFF [AFLOW] end " << endl;

      if(vflags.KBIN_VASP_FORCE_OPTION_LSCOUPLING.option==FALSE) {;} // dummy

      DONE=TRUE;
    }

    // ***************************************************************************
    // AUTO_MAGMOM AUTO_MAGMOM AUTO_MAGMOM AUTO_MAGMOM AUTO_MAGMOM AUTO_MAGMOM
    if(command=="AUTO_MAGMOM") {
      for(int i=1;i<=imax;i++) {
        strline=aurostd::GetLineString(FileContent,i);
        if(OPTION && (aurostd::substring2bool(strline,"MAGMOM",TRUE) || aurostd::substring2bool(strline,"#MAGMOM",TRUE))) { //CO
          if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==ON) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC(AUTO_MAGMOM+ON)) " << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==OFF) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC(AUTO_MAGMOM+OFF)) " << endl;
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }
      // xvasp.INCAR << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==ON) xvasp.INCAR << "# Performing AUTO_MAGMOM=ON [AFLOW] begin" << endl;
      // if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==OFF) xvasp.INCAR << "# Performing AUTO_MAGMOM=OFF [AFLOW] begin" << endl;
      if(OPTION==ON) {
        xvasp.INCAR << "MAGMOM= ";
        for(uint i=0;i<xvasp.str.atoms.size();i++) {
          if(xvasp.str.atoms.at(i).order_parameter_atom==FALSE) xvasp.INCAR << " 0";
          if(xvasp.str.atoms.at(i).order_parameter_atom==TRUE) xvasp.INCAR << " " << xvasp.str.atoms.at(i).order_parameter_value;
        }
        xvasp.INCAR << "   # " << xvasp.str.atoms.size() << " atoms " << endl;
      }
      if(OPTION==OFF) {;} // nothing to be done
      if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==ON) xvasp.INCAR << "# Performing AUTO_MAGMOM=ON [AFLOW] end " << endl;
      // if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==OFF) xvasp.INCAR << "# Performing AUTO_MAGMOM=OFF [AFLOW] end " << endl;

      DONE=TRUE;
    }

    // ***************************************************************************
    // NSW NSW NSW NSW NSW NSW
    if(command=="NSW") {
      for(int i=1;i<=imax;i++) {
        strline=aurostd::GetLineString(FileContent,i);
        if(aurostd::substring2bool(strline,"NSW",TRUE) || aurostd::substring2bool(strline,"#NSW",TRUE)) {
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC(NSW)) " << endl;
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }    // xvasp.INCAR << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing NSW [AFLOW] begin" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("NSW="+aurostd::utype2string(vflags.KBIN_VASP_FORCE_OPTION_NSW_EQUAL_VALUE),_incarpad_) << "# nsw " << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing NSW [AFLOW] end " << endl;

      DONE=TRUE;
    }


    // ***************************************************************************
    // SYM SYM SYM SYM SYM SYM
    if(command=="SYM") {
      for(int i=1;i<=imax;i++) {
        strline=aurostd::GetLineString(FileContent,i);
        if(aurostd::substring2bool(strline,"ISYM",TRUE) || aurostd::substring2bool(strline,"#ISYM",TRUE)) {
          if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==ON)  xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC(SYM+ON)) " << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==OFF) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC(SYM+OFF)) " << endl;
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }   // xvasp.INCAR << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==ON)  xvasp.INCAR << "# Performing SYM=ON [AFLOW] begin" << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==OFF) xvasp.INCAR << "# Performing SYM=OFF [AFLOW] begin" << endl;
      if(OPTION==ON)  xvasp.INCAR << aurostd::PaddedPOST("ISYM=2",_incarpad_) << "# SYMMETRY=ON" << endl; //CO20181113
      if(OPTION==OFF) xvasp.INCAR << aurostd::PaddedPOST("ISYM=0",_incarpad_) << "# SYMMETRY=OFF" << endl; //CO20181113
      if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==ON)  xvasp.INCAR << "# Performing SYM=ON [AFLOW] end " << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==OFF) xvasp.INCAR << "# Performing SYM=OFF [AFLOW] end " << endl;

      DONE=TRUE;
    }

    // ***************************************************************************
    // WAVECAR WAVECAR WAVECAR WAVECAR WAVECAR WAVECAR
    if(command=="WAVECAR") {
      for(int i=1;i<=imax;i++) {
        strline=aurostd::GetLineString(FileContent,i);
        if(aurostd::substring2bool(strline,"LWAVE",TRUE) || aurostd::substring2bool(strline,"#LWAVE",TRUE)) {
          if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==ON ) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC(WAVECAR+ON)) " << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==OFF) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC(WAVECAR+OFF)) " << endl;
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }    // xvasp.INCAR << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==ON ) xvasp.INCAR << "# Performing WAVECAR=ON [AFLOW] begin" << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==OFF) xvasp.INCAR << "# Performing WAVECAR=OFF [AFLOW] begin" << endl;
      if(OPTION==ON ) xvasp.INCAR << aurostd::PaddedPOST("LWAVE=.TRUE.",_incarpad_) << "# WAVECAR=ON" << endl;
      if(OPTION==OFF) {
        // if(vflags.KBIN_VASP_RUN_DIELECTRIC_STATIC || vflags.KBIN_VASP_RUN_DIELECTRIC_DYNAMIC) {
        //   xvasp.INCAR << aurostd::PaddedPOST("LWAVE=.TRUE.",_incarpad_) << "# WAVECAR=ON (forced TRUE by DIELECTRIC)" << endl;
        // } else {
        //   xvasp.INCAR << aurostd::PaddedPOST("LWAVE=.FALSE.",_incarpad_) << "# WAVECAR=OFF" << endl;
        // }
        xvasp.INCAR << aurostd::PaddedPOST("LWAVE=.FALSE.",_incarpad_) << "# WAVECAR=OFF" << endl;
      }
      if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==ON ) xvasp.INCAR << "# Performing WAVECAR=ON [AFLOW] end " << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==OFF) xvasp.INCAR << "# Performing WAVECAR=OFF [AFLOW] end " << endl;

      DONE=TRUE;
    }

    // ***************************************************************************
    // CHGCAR CHGCAR CHGCAR CHGCAR CHGCAR CHGCAR
    if(command=="CHGCAR") {
      for(int i=1;i<=imax;i++) {
        strline=aurostd::GetLineString(FileContent,i);
        if(aurostd::substring2bool(strline,"LCHARG",TRUE) || aurostd::substring2bool(strline,"#LCHARG",TRUE)) {
          if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==ON ) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC(CHGCAR+ON)) " << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==OFF) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC(CHGCAR+OFF)) " << endl;
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }   // xvasp.INCAR << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==ON ) xvasp.INCAR << "# Performing CHGCAR=ON [AFLOW] begin" << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==OFF) xvasp.INCAR << "# Performing CHGCAR=OFF [AFLOW] begin" << endl;
      if(OPTION==ON ) xvasp.INCAR << aurostd::PaddedPOST("LCHARG=.TRUE.",_incarpad_) << "# CHGCAR=ON" << endl;
      if(OPTION==OFF) xvasp.INCAR << aurostd::PaddedPOST("LCHARG=.FALSE.",_incarpad_) << "# CHGCAR=OFF" << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==ON ) xvasp.INCAR << "# Performing CHGCAR=ON [AFLOW] end " << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE && OPTION==OFF) xvasp.INCAR << "# Performing CHGCAR=OFF [AFLOW] end " << endl;

      DONE=TRUE;
    }

    //ME20191028
    // ICHARG ICHARG ICHARG ICHARG ICHARG ICHARG
    if (command == "ICHARG") {
      const string& chgcar = vflags.KBIN_VASP_FORCE_OPTION_CHGCAR_FILE.content_string;
      for (int i = 1; i <= imax; i++) {
        strline = aurostd::GetLineString(FileContent, i);
        if (aurostd::substring2bool(strline, "ICHARG", TRUE) || aurostd::substring2bool(strline, "#ICHARG", TRUE)) {
          if (vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_PREPARE_GENERIC)" << std::endl;
        } else {
          if (!vflags.KBIN_VASP_INCAR_VERBOSE && !strline.empty()) xvasp.INCAR << strline << std::endl;
          if (vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << std::endl;
        }
      }
      // Use negative values to just remove ICHARG
      if (ivalue >= 0) {
        if (vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing CHGCAR_FILE=" << chgcar << " [AFLOW] begin" << std::endl;
        xvasp.INCAR << "ICHARG=" << ivalue << std::endl;
        if (vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing CHGCAR_FILE=" << chgcar << " [AFLOW] end" << std::endl;
      }
      DONE = true;
    }

    // ***************************************************************************
    // GENERIC GENERIC GENERIC GENERIC GENERIC GENERIC
    if(command=="GENERIC") {
      DONE=TRUE;
    }

    // ***************************************************************************
    // GENERIC GENERIC GENERIC GENERIC GENERIC GENERIC
    if(command=="GENERIC") {
      DONE=TRUE;
    }

    // ***************************************************************************

    if(LDEBUG) cout << XPID << "KBIN::XVASP_INCAR_PREPARE_GENERIC: RETURNING CORRECTLY" << endl;
    if(LDEBUG) cout << XPID << "KBIN::XVASP_INCAR_PREPARE_GENERIC: command=" << command << endl;
    if(LDEBUG) cout << XPID << "KBIN::XVASP_INCAR_PREPARE_GENERIC: svalue=" << svalue << endl;
    if(LDEBUG) cout << XPID << "KBIN::XVASP_INCAR_PREPARE_GENERIC: ivalue=" << ivalue << endl;
    if(LDEBUG) cout << XPID << "KBIN::XVASP_INCAR_PREPARE_GENERIC: dvalue=" << dvalue << endl;
    if(LDEBUG) cout << XPID << "KBIN::XVASP_INCAR_PREPARE_GENERIC: OPTION=" << OPTION << endl;

    if(DONE) return TRUE;

    cerr << XPID << "KBIN::XVASP_INCAR_PREPARE_GENERIC: ERROR" << endl;
    cerr << XPID << "KBIN::XVASP_INCAR_PREPARE_GENERIC: command=" << command << " not found " << endl;
    cerr << XPID << "KBIN::XVASP_INCAR_PREPARE_GENERIC: svalue=" << svalue << endl;
    cerr << XPID << "KBIN::XVASP_INCAR_PREPARE_GENERIC: ivalue=" << ivalue << endl;
    cerr << XPID << "KBIN::XVASP_INCAR_PREPARE_GENERIC: dvalue=" << dvalue << endl;
    cerr << XPID << "KBIN::XVASP_INCAR_PREPARE_GENERIC: OPTION=" << OPTION << endl;

    return FALSE; //previously an exit
  }
} 

// ***************************************************************************
// KBIN::XVASP_INCAR_ADJUST_ICHARG
namespace KBIN {
  //ME20191028
  // When a CHGCAR file is specified in the aflow.in file, it is really only
  // useful in the first relxation calculations. For all other calculations,
  // it should not read from that CHGCAR file, so set ICHARG = 2 (default, but do not write) and comment
  // out the original CHGCAR file. Exception: if a CHGCAR file is output
  // after the relaxation, assume that the user wants to reuse it.
  void XVASP_INCAR_ADJUST_ICHARG(_xvasp& xvasp, _vflags& vflags, _aflags& aflags, int step, ofstream& FileMESSAGE) {
    if ((step == 1) && vflags.KBIN_VASP_FORCE_OPTION_CHGCAR_FILE.isentry) {
      // Do not set ICHARG when a CHGCAR file is output
      if (!vflags.KBIN_VASP_FORCE_OPTION_CHGCAR.option) {
        ostringstream aus;
        aus << "00000  MESSAGE ICHARG: Removing ICHARG - " << Message(aflags, _AFLOW_MESSAGE_DEFAULTS_) << std::endl;
        aurostd::PrintMessageStream(FileMESSAGE, aus, XHOST.QUIET);
        KBIN::XVASP_INCAR_PREPARE_GENERIC("ICHARG", xvasp, vflags, "", -1, 0.0, false); //remove any ICHARG from INCAR, leaving default ==2
        xvasp.aopts.flag("FLAG::XVASP_INCAR_changed", true);
      }
      if (xvasp.aopts.flag("FLAG::XVASP_INCAR_changed")) {
        xvasp.aopts.flag("FLAG::XVASP_INCAR_generated", true);
        xvasp.INCAR_orig.str(std::string());
        xvasp.INCAR_orig << xvasp.INCAR.str();
        if (step < xvasp.NRELAX) {  // Do not write when at the last step or there will be an extra INCAR
          aurostd::stringstream2file(xvasp.INCAR, string(xvasp.Directory+"/INCAR"));
        }
      }

      // Comment out CHGCAR file in aflow.in
      if (vflags.KBIN_VASP_FORCE_OPTION_CHGCAR_FILE.isentry) {
        stringstream aflowin_fixed;
        string filename = aurostd::CleanFileName(xvasp.Directory + "/" + _AFLOWIN_);
        string filecontent = aurostd::file2string(aurostd::CleanFileName(filename));
        int nlines = aurostd::GetNLinesString(filecontent);
        string line = "";
        for (int l = 1; l <= nlines; l++) {
          line = aurostd::GetLineString(filecontent, l);
          if (aurostd::substring2bool("[VASP_FORCE_OPTION]CHGCAR_FILE=", line)) {
            string line_fixed = aurostd::RemoveWhiteSpacesFromTheFront(line);
            if (line_fixed[0] != '#') line = "#" + line_fixed;
          }
          aflowin_fixed << line << std::endl;
        }
        aurostd::stringstream2file(aflowin_fixed, filename);
      }
    }
  }
}  // namespace KBIN

// ***************************************************************************
// KBIN::XVASP_INCAR_SPIN_REMOVE_RELAX
namespace KBIN {
  void XVASP_INCAR_SPIN_REMOVE_RELAX(_xvasp& xvasp,_aflags &aflags,_vflags& vflags,int step,ofstream &FileMESSAGE) {        // AFLOW_FUNCTION_IMPLEMENTATION
    ostringstream aus;
    bool fix_aflowlin=FALSE;
    _kflags kflags;
    if(step==1 && vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1 && abs(xvasp.str.qm_mag_atom)<DEFAULT_VASP_SPIN_REMOVE_CUTOFF) {
      aus << "00000  MESSAGE FORCE SPIN_OFF: [VASP_FORCE_OPTION]SPIN_REMOVE_RELAX_1  SPIN= " << xvasp.str.qm_mag_atom << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("SPIN",xvasp,vflags,"",0,0.0,OFF);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
      fix_aflowlin=TRUE;
    }
    if(step==2 && vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2 && abs(xvasp.str.qm_mag_atom)<DEFAULT_VASP_SPIN_REMOVE_CUTOFF) {
      aus << "00000  MESSAGE FORCE SPIN_OFF: [VASP_FORCE_OPTION]SPIN_REMOVE_RELAX_2  SPIN= " << xvasp.str.qm_mag_atom << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      KBIN::XVASP_INCAR_PREPARE_GENERIC("SPIN",xvasp,vflags,"",0,0.0,OFF);
      xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
      fix_aflowlin=TRUE;
    }
    if(xvasp.aopts.flag("FLAG::XVASP_INCAR_changed")) {
      xvasp.aopts.flag("FLAG::XVASP_INCAR_generated",TRUE);
      xvasp.INCAR_orig.str(std::string()); xvasp.INCAR_orig << xvasp.INCAR.str();
      if (step < xvasp.NRELAX) {  //ME20200107 - do not write when at the last step or there will be an extra INCAR
        aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"));
      }
      // xvasp.INCAR << aurostd::file2string(xvasp.Directory+"/INCAR"); // DID REREAD
    }
    if(fix_aflowlin) {
      // fix aflowlin
      stringstream aus_exec;
      aus_exec << "cd " << xvasp.Directory << endl;
      aus_exec << "cat " << _AFLOWIN_ << " | sed \"s/\\[VASP_FORCE_OPTION\\]SPIN=/#\\[VASP_FORCE_OPTION\\]SPIN=/g\" | sed \"s/##\\[/#\\[/g\" > aflow.tmp && mv aflow.tmp " << _AFLOWIN_ << "" << endl;
      aus_exec << "echo \"[VASP_FORCE_OPTION]SPIN=OFF" << "      // Self Correction\"" << " >> " << _AFLOWIN_ << " " << endl;
      aurostd::execute(aus_exec);
      //ME20181117 - Also change vflags
      vflags.KBIN_VASP_FORCE_OPTION_SPIN.option = false;
      vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_1 = false;
      vflags.KBIN_VASP_FORCE_OPTION_SPIN_REMOVE_RELAX_2 = false;
    }
  }
}


// ***************************************************************************
// KBIN::XVASP_KPOINTS_IBZKPT_UPDATE
namespace KBIN {
  void XVASP_KPOINTS_IBZKPT_UPDATE(_xvasp& xvasp,_aflags &aflags,_vflags& vflags,int step,ofstream &FileMESSAGE) {        // AFLOW_FUNCTION_IMPLEMENTATION
    ostringstream aus;
    // aus << "00000  MESSAGE IBZKPT.relax1 - step=" << step << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
    // aus << "00000  MESSAGE IBZKPT.relax1 - vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("IBZKPT")=" << vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("IBZKPT") << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
    // aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    if(step>=1 && vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.flag("IBZKPT")) {
      aus << "00000  MESSAGE FORCE KPOINTS from IBZKPT.relax1 - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
      xvasp.KPOINTS.str(std::string());
      xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
      aurostd::file2stringstream(string(xvasp.Directory+"/IBZKPT.relax1"),xvasp.KPOINTS);
      xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
    }
    if(xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed")) {
      xvasp.aopts.flag("FLAG::XVASP_KPOINTS_generated",TRUE);
      xvasp.KPOINTS_orig.str(std::string()); xvasp.KPOINTS_orig << xvasp.KPOINTS.str();
      if (step < xvasp.NRELAX) {  //ME20200107 - do not write when at the last step or there will be an extra INCAR
        aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS"));
      }
      // xvasp.KPOINTS << aurostd::file2string(xvasp.Directory+"/KPOINTS"); // DID REREAD
    }
  }
}


// ***************************************************************************
// KBIN::XVASP_INCAR_LDAU_OFF
namespace KBIN {
  void XVASP_INCAR_LDAU_OFF(_xvasp& xvasp,bool VERBOSE) {        // AFLOW_FUNCTION_IMPLEMENTATION
    string FileContent,strline;
    FileContent=xvasp.INCAR.str();
    xvasp.INCAR.str(std::string());
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    int imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i);
      if(!aurostd::substring2bool(strline,"SYSTEM", TRUE) &&
          (aurostd::substring2bool(strline,"LDAU",TRUE) || aurostd::substring2bool(strline,"#LDAU",TRUE) ||
           aurostd::substring2bool(strline,"LMAXMIX",TRUE) || aurostd::substring2bool(strline,"#LMAXMIX",TRUE))) {
        if(VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_LDAU_OFF) " << endl;
      } else {
        if(!VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
        if(VERBOSE) xvasp.INCAR << strline << endl;
      }
    } // xvasp.INCAR << endl;
    if(VERBOSE) xvasp.INCAR << "# Performing LDAU=OFF [AFLOW] begin" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("LDAU=.FALSE.",_incarpad_) << "# LDAU=OFF" << endl;
    if(VERBOSE) xvasp.INCAR << "# Performing LDAU=OFF [AFLOW] end " << endl;
  }
}


// ***************************************************************************
// KBIN::XVASP_INCAR_LDAU_ON
namespace KBIN {
  void XVASP_INCAR_LDAU_ON(_xvasp& xvasp,_vflags& vflags,uint type) {        // AFLOW_FUNCTION_IMPLEMENTATION
    string soliloquy=XPID+"KBIN::XVASP_INCAR_LDAU_ON():";
    stringstream message;
    /// is LDAU necessary ?
    if(type!=1 && type!=2) return; // NO LDAU
    // YES
    string FileContent,strline;
    FileContent=xvasp.INCAR.str();
    xvasp.INCAR.str(std::string());
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    int imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i);
      if(!aurostd::substring2bool(strline,"SYSTEM", TRUE) &&
          !aurostd::substring2bool(strline,"LDAUL",TRUE) && !aurostd::substring2bool(strline,"#LDAUL",TRUE) &&
          !aurostd::substring2bool(strline,"LDAUU",TRUE) && !aurostd::substring2bool(strline,"#LDAUU",TRUE) &&
          !aurostd::substring2bool(strline,"LDAUJ",TRUE) && !aurostd::substring2bool(strline,"#LDAUJ",TRUE) &&
          (aurostd::substring2bool(strline,"LDAU",TRUE) || aurostd::substring2bool(strline,"#LDAU",TRUE) ||
           aurostd::substring2bool(strline,"LMAXMIX",TRUE) || aurostd::substring2bool(strline,"#LMAXMIX",TRUE))) {
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_LDAU" << type << "_ON) " << endl;
      } else {
        if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
        if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
      }
    } // xvasp.INCAR << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing LDAU" << type << "=ON [AFLOW] begin" << endl;
    int LMAXMIX=6;
    xvasp.INCAR << aurostd::PaddedPOST("LDAU=.TRUE.",_incarpad_) << "# AFLOW LSDA+U" << endl;
    // LDAU generation

    vector<string> vLDAUspecies;vector<uint> vLDAUtype;vector<int> vLDAUL;vector<double> vLDAUU,vLDAUJ;
    bool LDAU=FALSE;vLDAUtype.clear();vLDAUL.clear();vLDAUU.clear();vLDAUJ.clear();
    uint nspecies=0;
    vector<string> tokens_group,tokens;

    if(vflags.KBIN_VASP_LDAU_AFLOW_AUTO_flag==TRUE) { // get them from the names
      // cerr << soliloquy << " vflags.KBIN_VASP_LDAU_AFLOW_AUTO_flag==TRUE => get them from the names " << endl;
      if(vflags.KBIN_VASP_LDAU_SPECIES.length()>1) {
        aurostd::string2tokens(vflags.KBIN_VASP_LDAU_SPECIES,tokens," ");
        nspecies=tokens.size();
        if(xvasp.str.species.size()!=nspecies) {
          throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"xvasp.str.species.size() != nspecies=tokens.size()",_INDEX_MISMATCH_); //CO20200624
        }
        for(uint i=0;i<xvasp.str.species.size();i++)
          AVASP_Get_LDAU_Parameters(aurostd::CleanStringASCII(KBIN::VASP_PseudoPotential_CleanName(tokens.at(i))),LDAU,vLDAUspecies,vLDAUtype,vLDAUL,vLDAUU,vLDAUJ); // parameters for LDAU2
      }
    } else { // get them from PARAMETERS
      // cerr << soliloquy << " vflags.KBIN_VASP_LDAU_AFLOW_AUTO_flag==FALSE => get them from the parameters:= " << vflags.KBIN_VASP_LDAU_PARAMETERS << endl;
      LDAU=TRUE;
      aurostd::string2tokens(vflags.KBIN_VASP_LDAU_PARAMETERS,tokens_group,";");
      if(tokens_group.size()<3) {
        throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"must specify 4 entries in "+vflags.KBIN_VASP_LDAU_PARAMETERS,_INPUT_NUMBER_); //CO20200624
      }
      aurostd::string2tokens(tokens_group.at(0),tokens,",");
      nspecies=tokens.size();
      for(uint j=0;j<nspecies;j++) { vLDAUspecies.push_back("");vLDAUtype.push_back(2);vLDAUL.push_back(-1);vLDAUU.push_back(0.0);vLDAUJ.push_back(0.0);} // make space
      for(uint i=0;i<tokens_group.size();i++) {
        // cerr << i << " " << tokens_group.at(i) << endl;
        aurostd::string2tokens(tokens_group.at(i),tokens,",");
        if(nspecies!=tokens.size()) {
          throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"must specify "+aurostd::utype2string(nspecies)+" entries in "+tokens_group.at(i),_INPUT_NUMBER_); //CO20200624
        }
        //	for(uint j=0;j<nspecies;j++) cerr << tokens.at(j) << " "; cerr << endl;
        if(i==0) for(uint j=0;j<nspecies;j++) vLDAUspecies.at(j)=tokens.at(j);
        if(i==1) for(uint j=0;j<nspecies;j++) vLDAUL.at(j)=aurostd::string2utype<int>(tokens.at(j));
        if(i==2) for(uint j=0;j<nspecies;j++) vLDAUU.at(j)=aurostd::string2utype<double>(tokens.at(j));
        if(i==3) for(uint j=0;j<nspecies;j++) vLDAUJ.at(j)=aurostd::string2utype<double>(tokens.at(j));
      }
    }
    // get type
    if(xvasp.str.species.size()!=nspecies) {
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"xvasp.str.species.size() != nspecies=tokens.size()",_INPUT_NUMBER_); //CO20200624
    }

    bool type2=FALSE,type1=FALSE;
    for(uint i=0;i<xvasp.str.species.size();i++) {
      if(vLDAUtype.at(i)==0) {;} // do nothing
      if(vLDAUtype.at(i)==1) type1=TRUE;
      if(vLDAUtype.at(i)==2) type2=TRUE;
    }
    if(type1==FALSE && type2==FALSE)             // no LDAU
      cerr << soliloquy << " no need for LDAU, all species do not need or are not in the table : " << xvasp.str.title << endl;
    if(type1==TRUE && type2==FALSE && type!=1)   // LDAU type 1
      cerr << soliloquy << " your type=" << type << " is incompatible with type1 of table : " << xvasp.str.title << endl;
    if(type1==FALSE && type2==TRUE && type!=2)   // LDAU type 2
      cerr << soliloquy << " your type=" << type << " is incompatible with type2 of table : " << xvasp.str.title << endl;
    if(type1==TRUE && type2==TRUE) {
      message << "Your type=" << type << " is incompatible with type1 and type2 of table : " << xvasp.str.title << endl;
      message << "  You can NOT HAVE type1 and type2 simultaneously from table.... My suggestion is that you remove LDAU_SPECIES from " << _AFLOWIN_ << "" << endl;
      message << "  and specify LDAUL, LDAUU, LDAUJ in the INCAR part of aflow,in manually" << endl;
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_INPUT_ILLEGAL_); //CO20200624
    }

    LMAXMIX=0;
    for(uint i=0;i<vLDAUL.size();i++) if(2*vLDAUL.at(i)>LMAXMIX) LMAXMIX=2*vLDAUL.at(i);        // get the MAX MIX for the orbitals
    // this is to max all the orbitals... if necessary..
    if(0) { //WSETYAWAN DISAGREES // do I need to max the orbitals ?
      for(uint i=0;i<vLDAUL.size();i++) if(vLDAUL.at(i)>vLDAUL.at(0)) vLDAUL.at(0)=vLDAUL.at(i);  // MAX the orbitals
      for(uint i=0;i<vLDAUL.size();i++) vLDAUL.at(i)=vLDAUL.at(0);                                // MAX the orbitals
    }
    // FIX NEGATIVES
    for(uint i=0;i<vLDAUL.size();i++) if(vLDAUL.at(i)<0) vLDAUL.at(i)=0;

    if(1) {
      stringstream straus;
      straus.clear();straus.str(std::string());
      straus << "#LDAU_SPECIES=";for(uint i=0;i<nspecies;i++) straus <<  vLDAUspecies.at(i) << " ";
      xvasp.INCAR << aurostd::PaddedPOST(straus.str(),_incarpad_) << "# LDAU species " << endl;
      if(vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.content_int==0) {
        straus.clear();straus.str(std::string());
        straus << "LDAUL=";for(uint i=0;i<nspecies;i++) straus <<  vLDAUL.at(i) << " ";
        xvasp.INCAR << aurostd::PaddedPOST(straus.str(),_incarpad_) << "# l-quantum number for which the on site interaction is added (Default 2) automatic LDAUL table" << endl;
        straus.clear();straus.str(std::string());
        straus << "LDAUU=";for(uint i=0;i<nspecies;i++) straus <<  vLDAUU.at(i) << " ";
        xvasp.INCAR << aurostd::PaddedPOST(straus.str(),_incarpad_) << "# UEFF parameter. Automatic LDAUU table" << endl;
        straus.clear();straus.str(std::string());
        straus << "LDAUJ=";for(uint i=0;i<nspecies;i++) straus <<  vLDAUJ.at(i) << " ";
        xvasp.INCAR << aurostd::PaddedPOST(straus.str(),_incarpad_) << "# J parameter (if used). Automatic LDAUJ table" << endl;
      } else {
        xvasp.INCAR << "#ADIABATIC_LDAU table = {nsteps,nspecies,{{LUJ}_nspecies}_steps}" << endl;
        straus.clear();straus.str(std::string());
        straus << "#ADIABATIC_LDAU="<< vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.content_int << "," << nspecies << ",";
        for(int j=1;j<=vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.content_int;j++) {
          for(uint i=0;i<nspecies;i++) {
            double cff=0.0;
            if(vLDAUL.at(i)==0) cff=0.0;                                                        // s orbital no U
            if(vLDAUL.at(i)==1) cff=0.0;                                                        // p orbital no U
            if(vLDAUL.at(i)==2) cff=(double) ((j+1.5)/(double) (vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.content_int));  // d orbital U
            if(vLDAUL.at(i)==3) cff=(double) ((j+1.5)/(double) (vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.content_int));  // f orbital U
            if(cff<0.0) cff=0.0;
            if(cff>1.0) cff=1.0;
            straus << vLDAUL.at(i) << "," << vLDAUU.at(i)*cff << "," <<  vLDAUJ.at(i)*cff;
            if(j<vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.content_int || i<nspecies-1) straus << ",";
          }
        }
        xvasp.INCAR << straus.str() << " " << endl; //  << "# ADIABATIC LDAU table = nsteps,nspecies,{{LUJ}_nspecies}_steps" << endl;
      }
    }

    if(0) {
      xvasp.INCAR << "#LDAU_SPECIES=";for(uint i=0;i<nspecies;i++) xvasp.INCAR << vLDAUspecies.at(i) << " ";xvasp.INCAR << " # LDAU species " << endl;
      xvasp.INCAR << "LDAUL=";for(uint i=0;i<nspecies;i++) xvasp.INCAR << vLDAUL.at(i) << " ";xvasp.INCAR << " # l-quantum number for which the on site interaction is added (Default 2) automatic LDAUL table" << endl;
      xvasp.INCAR << "LDAUU=";for(uint i=0;i<nspecies;i++) xvasp.INCAR << vLDAUU.at(i) << " ";xvasp.INCAR << " # UEFF parameter. Automatic LDAUU table" << endl;
      xvasp.INCAR << "LDAUJ=";for(uint i=0;i<nspecies;i++) xvasp.INCAR << vLDAUJ.at(i) << " ";xvasp.INCAR << " # J parameter (if used). Automatic LDAUJ table" << endl;
    }

    xvasp.INCAR << aurostd::PaddedPOST("LDAUTYPE="+aurostd::utype2string(type),_incarpad_) << "# Type of LDA+U. " << endl;
    xvasp.INCAR << aurostd::PaddedPOST("LMAXMIX="+aurostd::utype2string(LMAXMIX),_incarpad_) << "# Controls up to which l-quantum number the onsite PAW charge densities are passed through the charge density mixer. " << endl;
    // xvasp.INCAR << aurostd::PaddedPOST("LDAUPRINT=1",_incarpad_) << "# Controls verbosity of the L(S)DA+U module. (Default 0) # AFLOW LSDA+U" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("LDAUPRINT=0",_incarpad_) << "# Controls verbosity of the L(S)DA+U module. (Default 0) # AFLOW LSDA+U" << endl;

    // cerr << xvasp.INCAR.str() << endl;
    if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing LDAU" << type << "=ON [AFLOW] end " << endl;
  }
}


// ***************************************************************************
// KBIN::XVASP_INCAR_LDAU_ADIABATIC
namespace KBIN {
  void XVASP_INCAR_LDAU_ADIABATIC(_xvasp& xvasp,int step) {        // AFLOW_FUNCTION_IMPLEMENTATION
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy=XPID+"KBIN::XVASP_INCAR_LDAU_ADIABATIC():";
    // reload incar
    xvasp.INCAR_orig.clear(); xvasp.INCAR_orig.str(xvasp.INCAR.str());
    xvasp.INCAR.clear(); xvasp.INCAR.str(aurostd::file2string(xvasp.Directory+"/INCAR"));
    //  xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    // operate
    string FileContent=xvasp.INCAR.str(),strline;
    int imax=aurostd::GetNLinesString(FileContent);
    stringstream straus;
    straus.clear();straus.str(std::string());

    for(int iline=1;iline<=imax;iline++) {
      strline=aurostd::GetLineString(FileContent,iline);
      if(!aurostd::substring2bool(strline,"LDAUL=",TRUE) && !aurostd::substring2bool(strline,"LDAUU=",TRUE) && !aurostd::substring2bool(strline,"LDAUJ=",TRUE)) {
        xvasp.INCAR << strline << endl;
        if(strline.find("#ADIABATIC_LDAU=")!=string::npos) {
          // FOUND write LDAUL LDAUU LDAUJ
          if(LDEBUG) cerr << "[KBIN::XVASP_INCAR_LDAU_ADIABATIC] LDAU_ADIABATIC" << endl;
          vector<string> tokens1,tokens;
          aurostd::string2tokens(strline,tokens,"=");
          aurostd::string2tokens(tokens.at(1),tokens1,"#");
          aurostd::string2tokens(tokens1.at(0),tokens,",");
          uint t=0; // start of tokens. j
          uint nsteps=aurostd::string2utype<int>(tokens.at(t++));
          uint nspecies=aurostd::string2utype<int>(tokens.at(t++));
          if(LDEBUG) cerr << "[KBIN::XVASP_INCAR_LDAU_ADIABATIC] LDAU_STEP=" << step << endl;
          if(LDEBUG) cerr << "[KBIN::XVASP_INCAR_LDAU_ADIABATIC] LDAU_ADIABATIC_NSTEPS=" << nsteps << endl;
          if(LDEBUG) cerr << "[KBIN::XVASP_INCAR_LDAU_ADIABATIC] LDAU_ADIABATIC_NSPECIES=" << nspecies << endl;
          if(tokens.size() != (uint) nsteps*nspecies*3+2) {
            throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"wrong number of ADIABATIC_LDAU entries ("+aurostd::utype2string(tokens.size())+")",_INPUT_NUMBER_); //CO20200624
          }
          vector<int> vLDAUL;vector<double> vLDAUU,vLDAUJ;
          for(uint k=1;k<=nsteps;k++) {
            for(uint j=0;j<nspecies;j++) {
              if((int) k==step) {
                vLDAUL.push_back(aurostd::string2utype<int>(tokens.at(t++)));
                vLDAUU.push_back(aurostd::string2utype<double>(tokens.at(t++)));
                vLDAUJ.push_back(aurostd::string2utype<double>(tokens.at(t++)));
              } else {
                t+=3;
              }
            }
          }
          straus.clear();straus.str(std::string());
          straus << "LDAUL=";for(uint j=0;j<nspecies;j++) straus <<  vLDAUL.at(j) << " ";
          xvasp.INCAR << aurostd::PaddedPOST(straus.str(),_incarpad_) << "# l-quantum number for which the on site interaction is added (Default 2) automatic LDAUL table" << endl;
          if(LDEBUG) cerr << "[KBIN::XVASP_INCAR_LDAU_ADIABATIC] " << aurostd::PaddedPOST(straus.str(),_incarpad_) << endl;
          straus.clear();straus.str(std::string());
          straus << "LDAUU=";for(uint j=0;j<nspecies;j++) straus <<  vLDAUU.at(j) << " ";
          xvasp.INCAR << aurostd::PaddedPOST(straus.str(),_incarpad_) << "# UEFF parameter. Automatic LDAUU table" << endl;
          if(LDEBUG) cerr << "[KBIN::XVASP_INCAR_LDAU_ADIABATIC] " << aurostd::PaddedPOST(straus.str(),_incarpad_) << endl;
          straus.clear();straus.str(std::string());
          straus << "LDAUJ=";for(uint j=0;j<nspecies;j++) straus <<  vLDAUJ.at(j) << " ";
          xvasp.INCAR << aurostd::PaddedPOST(straus.str(),_incarpad_) << "# J parameter (if used). Automatic LDAUJ table" << endl;
          if(LDEBUG) cerr << "[KBIN::XVASP_INCAR_LDAU_ADIABATIC] " << aurostd::PaddedPOST(straus.str(),_incarpad_) << endl;
        }
      }
    }
    // xvasp.INCAR << endl;
    // rewrite incar
    aurostd::RemoveFile(string(xvasp.Directory+"/INCAR"));
    aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"));
  }
}


// ***************************************************************************
// KBIN::XVASP_INCAR_LDAU_CUTOFF
namespace KBIN {
  void XVASP_INCAR_LDAU_CUTOFF(_xvasp& xvasp,bool VERBOSE) {        // AFLOW_FUNCTION_IMPLEMENTATION
    string FileContent,strline;
    FileContent=xvasp.INCAR.str();
    xvasp.INCAR.str(std::string());
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    int imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i);
      if(!aurostd::substring2bool(strline,"SYSTEM", TRUE) &&
          (aurostd::substring2bool(strline,"LDAU",TRUE) || aurostd::substring2bool(strline,"#LDAU",TRUE))) {
        if(VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_LDAU_CUTOFF) " << endl;
      } else {
        if(!VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
        if(VERBOSE) xvasp.INCAR << strline << endl;
      }
    }
    // xvasp.INCAR << endl;
    if(VERBOSE) xvasp.INCAR << "# Performing LDAU=CUTOFF [AFLOW] begin" << endl;
    xvasp.INCAR << aurostd::PaddedPOST("LDAU=.FALSE.",_incarpad_) << "# LDAU=CUTOFF" << endl;
    if(VERBOSE) xvasp.INCAR << "# Performing LDAU=CUTOFF [AFLOW] end " << endl;
  }
}

// ***************************************************************************
// KBIN::XVASP_INCAR_EFIELD_PEAD
namespace KBIN {
  void XVASP_INCAR_EFIELD_PEAD(_xvasp& xvasp,double dvalue,bool VERBOSE) {        // AFLOW_FUNCTION_IMPLEMENTATION
    string FileContent,strline;
    FileContent=xvasp.INCAR.str();
    xvasp.INCAR.str(std::string());
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    int imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i);
      if(aurostd::substring2bool(strline,"EFIELD",TRUE) || aurostd::substring2bool(strline,"#EFIELD",TRUE)) {
        if(VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_EFIELD_PEAD) " << endl;
      } else {
        if(!VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
        if(VERBOSE) xvasp.INCAR << strline << endl;
      }
    }
    xvector<double> E(3);
    E(1)=DEFAULT_EFIELD_PEAD;
    E(2)=DEFAULT_EFIELD_PEAD;
    E(3)=DEFAULT_EFIELD_PEAD;
    vector<string> tokens;
    // from vasp.out
    if(0) for(uint i=1;i<=3;i++) {
      aurostd::string2tokens(aurostd::execute2string("cat "+xvasp.Directory+"/vasp.out | grep E_g/N_"+aurostd::utype2string(i)),tokens,">");
      if(tokens.size()>0) {
        aurostd::string2tokens(string(tokens.at(1)),tokens,"=");
        if(tokens.size()>0) E(i)=aurostd::string2utype<double>(tokens.at(1));
        if(E(i)<=1.0E-6) E(i)=DEFAULT_EFIELD_PEAD;
      }
      //    cerr << "E(" << i << ")=" << E(i) << endl;
    }
    // from OUTCAR
    aurostd::string2tokens(aurostd::execute2string("cat "+xvasp.Directory+"/OUTCAR | grep EFIELD_PEAD"),tokens,"=");
    if(tokens.size()>0) {
      //   cerr << tokens.at(1) << endl;
      if(tokens.size()>0) aurostd::string2tokens(string(tokens.at(1)),tokens);
      if(tokens.size()>2) 
        for(uint i=1;i<=3;i++) {
          E(i)=aurostd::string2utype<double>(tokens.at(i-1));
          if(E(i)<=1.0E-6) E(i)=DEFAULT_EFIELD_PEAD;
          //	cerr << "E(" << i << ")=" << E(i) << endl;
        }
    }

    // if(i<tokens.size()-1 && aurostd::substring2bool(tokens.at(i),"NBANDS")) NBANDS_OUTCAR=aurostd::string2utype<int>(tokens.at(i+1));

    // xvasp.INCAR << endl;
    E(1)*=dvalue;E(2)*=dvalue;E(3)*=dvalue;
    // for(uint i=1;i<=3;i++) cerr << "E(" << i << ")=" << E(i) << endl;

    if(VERBOSE) xvasp.INCAR << "# Performing EFIELD_PEAD [AFLOW] begin dvalue=" << dvalue << endl;
    xvasp.INCAR << aurostd::PaddedPOST("EFIELD_PEAD="+aurostd::utype2string(E(1),6)+" "+aurostd::utype2string(E(2),6)+" "+aurostd::utype2string(E(3),6)+" ",_incarpad_) << "# EFIELD_PEAD" << endl;
    if(VERBOSE) xvasp.INCAR << "# Performing EFIELD_PEAD [AFLOW] end " << endl;
  }
}


// ***************************************************************************
// KBIN::XVASP_INCAR_KPOINTS_Dielectric  DIELECTRIC
namespace KBIN {
  void XVASP_INCAR_KPOINTS_Dielectric_SET(_xvasp& xvasp,_kflags &kflags,_vflags& vflags,string svalue) {        // AFLOW_FUNCTION_IMPLEMENTATION
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    // STATIC
    if(svalue=="STATIC" || svalue=="static") {
      //.  Retain the following static run entries and their values: ALGO, LREAL, NSIM, ISYM, IBRION, NSW, NELM, NELMIN, ENMAX, ISPIN, ISMEAR, SIGMA, and everything LDA+U related.
      // d.  Eliminate PSTRESS, EMIN, EMAX, LORBIT, ISIF, NEDOS.
      // NELM = 0

      string FileContent,strline;
      FileContent=xvasp.INCAR.str();
      xvasp.INCAR.str(std::string());xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
      int imax=aurostd::GetNLinesString(FileContent);
      for(int i=1;i<=imax;i++) {
        strline=aurostd::GetLineString(FileContent,i);
        if(aurostd::substring2bool(strline,"NEDOS",TRUE) || aurostd::substring2bool(strline,"#NEDOS",TRUE) // Camilo
            || aurostd::substring2bool(strline,"ISIF",TRUE) || aurostd::substring2bool(strline,"#ISIF",TRUE) // Camilo
            || aurostd::substring2bool(strline,"EMIN",TRUE) || aurostd::substring2bool(strline,"#EMIN",TRUE) // Camilo
            || aurostd::substring2bool(strline,"EMAX",TRUE) || aurostd::substring2bool(strline,"#EMAX",TRUE) // Camilo
            // ||  aurostd::substring2bool(strline,"NELM",TRUE) || aurostd::substring2bool(strline,"#NELM",TRUE) // Camilo
            || aurostd::substring2bool(strline,"ICHARG",TRUE) || aurostd::substring2bool(strline,"#ICHARG",TRUE) // Camilo
            || aurostd::substring2bool(strline,"LCHARG",TRUE) || aurostd::substring2bool(strline,"#LCHARG",TRUE) // Camilo
            || aurostd::substring2bool(strline,"LEPSILON",TRUE) || aurostd::substring2bool(strline,"#LEPSILON",TRUE) // Camilo
            || aurostd::substring2bool(strline,"LCALCEPS",TRUE) || aurostd::substring2bool(strline,"#LCALCEPS",TRUE) // Camilo
            || aurostd::substring2bool(strline,"LRPA",TRUE) || aurostd::substring2bool(strline,"#LRPA",TRUE) // Camilo
            || aurostd::substring2bool(strline,"LORBIT",TRUE) || aurostd::substring2bool(strline,"#LORBIT",TRUE)  // Camilo
            || aurostd::substring2bool(strline,"LWAVE",TRUE) || aurostd::substring2bool(strline,"#LWAVE",TRUE)  // Camilo
            || aurostd::substring2bool(strline,"NPAR",TRUE) || aurostd::substring2bool(strline,"#NPAR",TRUE)   // NPAR UNCOMPATIBLE WITH Dielectric
            // hope...	 || aurostd::substring2bool(strline,"NBANDS",TRUE) || aurostd::substring2bool(strline,"#NBANDS",TRUE) || aurostd::substring2bool(strline,"#fixed nbands",TRUE)   // NPAR UNCOMPATIBLE WITH Dielectric
          ) { // Camilo
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_KPOINTS_Dielectric_SET)" << endl;
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }
      // xvasp.INCAR << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing DIELECTRIC_STATIC [AFLOW] begin" << endl;
      // xvasp.INCAR << aurostd::PaddedPOST("NELM=0",_incarpad_) << "# Performing DIELECTRIC_STATIC (Camilo)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("LEPSILON=.TRUE.",_incarpad_) << "# Performing DIELECTRIC_STATIC (Camilo)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("LRPA=.TRUE.",_incarpad_) << "# Performing DIELECTRIC_STATIC (Camilo)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("LWAVE=.TRUE.",_incarpad_) << "# Performing DIELECTRIC_STATIC (Camilo)" << endl;
      // FIX ISMEAR=0 and SIGMA=0.01 ?
      if(0) { // NO NPAR
        if(kflags.KBIN_MPI)
          //      xvasp.INCAR << aurostd::PaddedPOST("NPAR="+aurostd::utype2string(kflags.KBIN_MPI_NCPUS),_incarpad_) << "# Performing DIELECTRIC_STATIC (Camilo)" << endl;
          xvasp.INCAR << aurostd::PaddedPOST("NPAR=1",_incarpad_) << "# Performing DIELECTRIC_STATIC (Camilo)" << endl;
      }
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing DIELECTRIC_STATIC [AFLOW] end" << endl;

      // now to the KPOINTS
      FileContent=xvasp.KPOINTS.str();
      imax=aurostd::GetNLinesString(FileContent);
      xvasp.KPOINTS.str(std::string());xvasp.aopts.flag("FLAG::XVASP_KPOINTS_changed",TRUE);
      xvasp.KPOINTS << "KPOINTS File, automatically generated by KPPRA_DELTA - aflow with DK=" << DIELECTRIC_DK << endl;
      xvasp.KPOINTS << aurostd::GetLineString(FileContent,2) << endl;
      KPPRA_DELTA(xvasp.str,DIELECTRIC_DK);  // COUT if you want...
      if(aurostd::FileExist(xvasp.Directory+string("/KPOINTS.bands"))) {  // with some LUCK I`ve already nailed the FCC/HEX
        xvasp.str.kpoints_kscheme=aurostd::GetLineString(FileContent,3);
        string KPOINTS_bands=aurostd::file2string(xvasp.Directory+"/KPOINTS.bands"); // cerr << KPOINTS_bands << endl;
        bool isFCC=aurostd::substring2bool(KPOINTS_bands,"FCC","fcc"); // cerr << XPID << "KBIN::XVASP_INCAR_KPOINTS_Dielectric_SET: isFCC=" << isFCC << endl;
        bool isHEX=aurostd::substring2bool(KPOINTS_bands,"HEX","hex"); // cerr << XPID << "KBIN::XVASP_INCAR_KPOINTS_Dielectric_SET: isHEX=" << isHEX << endl;
        bool isRHL=aurostd::substring2bool(KPOINTS_bands,"RHL","rhl"); // cerr << XPID << "KBIN::XVASP_INCAR_KPOINTS_Dielectric_SET: isRHL=" << isRHL << endl;
        if(abs(xvasp.str.kpoints_s1)+abs(xvasp.str.kpoints_s2)+abs(xvasp.str.kpoints_s3)>0.1) {isFCC=FALSE;isHEX=FALSE;isRHL=FALSE;} // no shift if already shifted
        //  if(isFCC || isHEX) {cerr << XPID << "KBIN::XVASP_INCAR_KPOINTS_Dielectric_SET switching to GAMMA" << endl;}
        if(isFCC || isHEX || isRHL) {
          xvasp.str.kpoints_kscheme="Gamma";
          if(_iseven(xvasp.str.kpoints_k1)) {xvasp.str.kpoints_k1++;xvasp.str.kpoints_s1=0.0;if(LDEBUG) cout << "Xeven k1=" << xvasp.str.kpoints_k1 << endl;}
          if(_iseven(xvasp.str.kpoints_k2)) {xvasp.str.kpoints_k2++;xvasp.str.kpoints_s2=0.0;if(LDEBUG) cout << "Yeven k2=" << xvasp.str.kpoints_k2 << endl;}
          if(_iseven(xvasp.str.kpoints_k3)) {xvasp.str.kpoints_k3++;xvasp.str.kpoints_s3=0.0;if(LDEBUG) cout << "Zeven k3=" << xvasp.str.kpoints_k3 << endl;}
        }  // else dont touch it
      }
      xvasp.KPOINTS << xvasp.str.kpoints_kscheme << endl;
      xvasp.KPOINTS << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << " " << endl;
      xvasp.KPOINTS << xvasp.str.kpoints_s1 << " " << xvasp.str.kpoints_s2 << " " << xvasp.str.kpoints_s3 << " " << endl;
      //    xvasp.KPOINTS << aurostd::GetLineString(FileContent,5) << endl;
    }
    // DYNAMIC
    if(svalue=="DYNAMIC" || svalue=="dynamic") {
      cerr << XPID << "KBIN::XVASP_INCAR_KPOINTS_Dielectric_SET DYNAMIC" << endl;
      // a.  Reuse the STEP 01 WAVECAR
      // b.  ALGO=EXACT NELM=1 LOPTICS=.TRUE.   CSHIFT=0.15  OMEGAMAX=25   NEDOS=12500
      // i.  Remove LEPSILON and LRPA
      string FileContent,strline;
      FileContent=xvasp.INCAR.str();
      xvasp.INCAR.str(std::string());xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
      int imax=aurostd::GetNLinesString(FileContent);
      for(int i=1;i<=imax;i++) {
        strline=aurostd::GetLineString(FileContent,i);
        if(aurostd::substring2bool(strline,"ALGO",TRUE) || aurostd::substring2bool(strline,"#ALGO",TRUE) || // Camilo
            aurostd::substring2bool(strline,"NELM",TRUE) || aurostd::substring2bool(strline,"#NELM",TRUE) || // Camilo
            aurostd::substring2bool(strline,"LOPTICS",TRUE) || aurostd::substring2bool(strline,"#LOPTICS",TRUE) || // Camilo
            aurostd::substring2bool(strline,"CSHIFT",TRUE) || aurostd::substring2bool(strline,"#CSHIFT",TRUE) || // Camilo
            aurostd::substring2bool(strline,"OMEGAMAX",TRUE) || aurostd::substring2bool(strline,"#OMEGAMAX",TRUE) || // Camilo
            aurostd::substring2bool(strline,"NEDOS",TRUE) || aurostd::substring2bool(strline,"#NEDOS",TRUE) || // Camilo
            aurostd::substring2bool(strline,"NBANDS",TRUE) || aurostd::substring2bool(strline,"#NBANDS",TRUE) || // Camilo
            aurostd::substring2bool(strline,"LEPSILON",TRUE) || aurostd::substring2bool(strline,"#LEPSILON",TRUE) || // Camilo
            aurostd::substring2bool(strline,"LCALCEPS",TRUE) || aurostd::substring2bool(strline,"#LCALCEPS",TRUE) || // Camilo
            aurostd::substring2bool(strline,"LRPA",TRUE) || aurostd::substring2bool(strline,"#LRPA",TRUE) || // Camilo
            aurostd::substring2bool(strline,"LORBIT",TRUE) || aurostd::substring2bool(strline,"#LORBIT",TRUE)) { // Camilo
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_KPOINTS_Dielectric_SET)" << endl;
        } else {
          if(!vflags.KBIN_VASP_INCAR_VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
          if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << strline << endl;
        }
      }
      // get NBANDS from OUTCAR.dielectric_static
      int NBANDS_OUTCAR=0;
      xOUTCAR OUTCAR_NBANDS(xvasp.Directory+"/OUTCAR.dielectric_static");
      NBANDS_OUTCAR=OUTCAR_NBANDS.NBANDS;

      // xvasp.INCAR << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing DIELECTRIC_DYNAMIC [AFLOW] begin" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("ALGO=EXACT",_incarpad_) << "# Performing DIELECTRIC_DYNAMIC (Camilo)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("NELM=1",_incarpad_) << "# Performing DIELECTRIC_DYNAMIC (Camilo)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("LOPTICS=.TRUE.",_incarpad_) << "# Performing DIELECTRIC_DYNAMIC (Camilo)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("CSHIFT=0.15",_incarpad_) << "# Performing DIELECTRIC_DYNAMIC (Camilo)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("OMEGAMAX=25",_incarpad_) << "# Performing DIELECTRIC_DYNAMIC (Camilo)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("NEDOS=12500",_incarpad_) << "# Performing DIELECTRIC_DYNAMIC (Camilo)" << endl;
      xvasp.INCAR << aurostd::PaddedPOST("NBANDS="+aurostd::utype2string((double) 4*NBANDS_OUTCAR),_incarpad_) << "# Performing DIELECTRIC_DYNAMIC NBANDS_OUTCAR="  << NBANDS_OUTCAR << " (Camilo)" << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) xvasp.INCAR << "# Performing DIELECTRIC_DYNAMIC [AFLOW] end" << endl;
    }
  }
}


// ***************************************************************************
// KBIN::XVASP_INCAR_REMOVE_ENTRY
namespace KBIN {
  void XVASP_INCAR_REMOVE_ENTRY(_xvasp& xvasp,string ENTRY,string COMMENT,bool VERBOSE) {        // AFLOW_FUNCTION_IMPLEMENTATION
    string FileContent,strline;
    FileContent=xvasp.INCAR.str();
    xvasp.INCAR.str(std::string());
    xvasp.aopts.flag("FLAG::XVASP_INCAR_changed",TRUE);
    int imax=aurostd::GetNLinesString(FileContent);
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(FileContent,i);
      if(aurostd::substring2bool(strline,ENTRY,TRUE) || aurostd::substring2bool(strline,"#"+ENTRY,TRUE)) {
        if(VERBOSE) xvasp.INCAR << "# " << strline << " # AFLOW REMOVED (KBIN::XVASP_INCAR_REMOVE_ENTRY) " << COMMENT << endl;
      } else {
        if(!VERBOSE && strline.length()) xvasp.INCAR << strline << endl;
        if(VERBOSE) xvasp.INCAR << strline << endl;
      }
    }
  }
}


// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// KPOINTS MODIFICATIONS
namespace KBIN {
  bool XVASP_KPOINTS_KPOINTS(_xvasp &xvasp,ofstream &FileMESSAGE,bool VERBOSE) {
    ostringstream aus;
    xvasp.KPOINTS.str(std::string()); xvasp.KPOINTS.clear();
    xvasp.KPOINTS_orig.str(std::string());xvasp.KPOINTS_orig.clear();
    xvasp.KPOINTS << "KPOINTS File, automatically generated by KPPRA - aflow with Kpoints=" << xvasp.str.kpoints_kppra << endl;
    xvasp.KPOINTS << xvasp.str.kpoints_mode << endl;  // MODE AUTO 0
    xvasp.KPOINTS << xvasp.str.kpoints_kscheme << endl;
    xvasp.KPOINTS << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << endl;
    xvasp.KPOINTS.precision(3);
    xvasp.KPOINTS << xvasp.str.kpoints_s1 << " " << xvasp.str.kpoints_s2 << " " << xvasp.str.kpoints_s3 << endl;

    aurostd::StringstreamClean(aus);
    if(VERBOSE) {
      aus << "00000  MESSAGE KPOINTS K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " ";
      aus << " Kshift=[" << xvasp.str.kpoints_s1 << " " << xvasp.str.kpoints_s2 << " " << xvasp.str.kpoints_s3 << "]" << " ";
      aus << " Kpoints=" << xvasp.str.kpoints_kppra  << " with " << xvasp.str.kpoints_kscheme << "   " << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);			
    }
    // KPOINTS done
    xvasp.KPOINTS_orig << xvasp.KPOINTS.str();
    return TRUE;
  }
}

namespace KBIN {
  bool XVASP_KPOINTS_KPOINTS(_xvasp &xvasp) {
    xvasp.KPOINTS.str(std::string()); xvasp.KPOINTS.clear();
    xvasp.KPOINTS_orig.str(std::string());xvasp.KPOINTS_orig.clear();
    xvasp.KPOINTS << "KPOINTS File, automatically generated by KBIN::XVASP_KPOINTS_KPOINTS " << endl;
    xvasp.KPOINTS << xvasp.str.kpoints_mode << endl;
    xvasp.KPOINTS << xvasp.str.kpoints_kscheme << endl;
    xvasp.KPOINTS << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << endl;
    xvasp.KPOINTS.precision(3);
    xvasp.KPOINTS << xvasp.str.kpoints_s1 << " " << xvasp.str.kpoints_s2 << " " << xvasp.str.kpoints_s3 << endl;
    xvasp.KPOINTS_orig << xvasp.KPOINTS.str();
    return TRUE;
  }
}

namespace KBIN {
  bool XVASP_KPOINTS_OPERATION(_xvasp& xvasp,string operation) {
    bool LDEBUG=(FALSE || XHOST.DEBUG); // TRUE;
    if(LDEBUG) cout << XPID << "KBIN::XVASP_KPOINTS_OPERATION  operation=" << operation << endl;
    if(LDEBUG) cout << XPID << "KBIN::XVASP_KPOINTS_OPERATION  xvasp.str.kpoints_k*=(" << xvasp.str.kpoints_k1 << "," << xvasp.str.kpoints_k2 << "," << xvasp.str.kpoints_k3 << ")" << endl;
    if(LDEBUG) cout << XPID << "KBIN::XVASP_KPOINTS_OPERATION  xvasp.str.kpoints_s*=(" << xvasp.str.kpoints_s1 << "," << xvasp.str.kpoints_s2 << "," << xvasp.str.kpoints_s3 << ")" << endl;
    if(LDEBUG) cout << XPID << "KBIN::XVASP_KPOINTS_OPERATION  xvasp.str.kpoints_kscheme=*=(" << xvasp.str.kpoints_kscheme << ")" << endl;
    uint i=0;
    bool change_made=false; //CO20210201
    if(operation.find("X--")) {i++;xvasp.str.kpoints_k1--;change_made=true;if(LDEBUG) cout << "X-- k1=" << xvasp.str.kpoints_k1 << endl;}
    if(operation.find("Y--")) {i++;xvasp.str.kpoints_k2--;change_made=true;if(LDEBUG) cout << "Y-- k2=" << xvasp.str.kpoints_k2 << endl;}
    if(operation.find("Z--")) {i++;xvasp.str.kpoints_k3--;change_made=true;if(LDEBUG) cout << "Z-- k3=" << xvasp.str.kpoints_k3 << endl;}
    if(operation.find("X++")) {i++;xvasp.str.kpoints_k1++;change_made=true;if(LDEBUG) cout << "X++ k1=" << xvasp.str.kpoints_k1 << endl;}
    if(operation.find("Y++")) {i++;xvasp.str.kpoints_k2++;change_made=true;if(LDEBUG) cout << "Y++ k2=" << xvasp.str.kpoints_k2 << endl;}
    if(operation.find("Z++")) {i++;xvasp.str.kpoints_k3++;change_made=true;if(LDEBUG) cout << "Z++ k3=" << xvasp.str.kpoints_k3 << endl;}
    change_made=false;  //reset for Xeven/Xodd/etc.
    if(operation.find("Xeven")) {i++;if(_isodd(xvasp.str.kpoints_k1)) {xvasp.str.kpoints_k1++;xvasp.str.kpoints_s1=0.0;change_made=true;if(LDEBUG) cout << "Xeven k1=" << xvasp.str.kpoints_k1 << endl;}}
    if(operation.find("Yeven")) {i++;if(_isodd(xvasp.str.kpoints_k2)) {xvasp.str.kpoints_k2++;xvasp.str.kpoints_s2=0.0;change_made=true;if(LDEBUG) cout << "Yeven k2=" << xvasp.str.kpoints_k2 << endl;}}
    if(operation.find("Zeven")) {i++;if(_isodd(xvasp.str.kpoints_k3)) {xvasp.str.kpoints_k3++;xvasp.str.kpoints_s3=0.0;change_made=true;if(LDEBUG) cout << "Zeven k3=" << xvasp.str.kpoints_k3 << endl;}}
    if(!change_made&&(operation.find("Xeven")||operation.find("Yeven")||operation.find("Zeven"))){  //CO20210201 - was already even before, increment by 2
      if(operation.find("Xeven")) {i++;xvasp.str.kpoints_k1+=2;xvasp.str.kpoints_s1=0.0;change_made=true;if(LDEBUG) cout << "Xeven k1=" << xvasp.str.kpoints_k1 << endl;}
      if(operation.find("Yeven")) {i++;xvasp.str.kpoints_k2+=2;xvasp.str.kpoints_s2=0.0;change_made=true;if(LDEBUG) cout << "Yeven k2=" << xvasp.str.kpoints_k2 << endl;}
      if(operation.find("Zeven")) {i++;xvasp.str.kpoints_k3+=2;xvasp.str.kpoints_s3=0.0;change_made=true;if(LDEBUG) cout << "Zeven k3=" << xvasp.str.kpoints_k3 << endl;}
    }
    change_made=false;  //reset for Xeven/Xodd/etc.
    if(operation.find("Xodd")) {i++;if(_iseven(xvasp.str.kpoints_k1)) {xvasp.str.kpoints_k1++;xvasp.str.kpoints_s1=0.0;change_made=true;if(LDEBUG) cout << "Xeven k1=" << xvasp.str.kpoints_k1 << endl;}}
    if(operation.find("Yodd")) {i++;if(_iseven(xvasp.str.kpoints_k2)) {xvasp.str.kpoints_k2++;xvasp.str.kpoints_s2=0.0;change_made=true;if(LDEBUG) cout << "Yeven k2=" << xvasp.str.kpoints_k2 << endl;}}
    if(operation.find("Zodd")) {i++;if(_iseven(xvasp.str.kpoints_k3)) {xvasp.str.kpoints_k3++;xvasp.str.kpoints_s3=0.0;change_made=true;if(LDEBUG) cout << "Zeven k3=" << xvasp.str.kpoints_k3 << endl;}}
    if(!change_made&&(operation.find("Xodd")||operation.find("Yodd")||operation.find("Zodd"))){  //CO20210201 - was already odd before, increment by 2
      if(operation.find("Xodd")) {i++;xvasp.str.kpoints_k1+=2;xvasp.str.kpoints_s1=0.0;change_made=true;if(LDEBUG) cout << "Xodd k1=" << xvasp.str.kpoints_k1 << endl;}
      if(operation.find("Yodd")) {i++;xvasp.str.kpoints_k2+=2;xvasp.str.kpoints_s2=0.0;change_made=true;if(LDEBUG) cout << "Yodd k2=" << xvasp.str.kpoints_k2 << endl;}
      if(operation.find("Zodd")) {i++;xvasp.str.kpoints_k3+=2;xvasp.str.kpoints_s3=0.0;change_made=true;if(LDEBUG) cout << "Zodd k3=" << xvasp.str.kpoints_k3 << endl;}
    }
    if(operation.find("Xevenshift")) {i++;if(_iseven(xvasp.str.kpoints_k1)) xvasp.str.kpoints_s1=0.5; else xvasp.str.kpoints_s1=0.0;change_made=true;if(LDEBUG) cout << "Xevenshift s1=" << xvasp.str.kpoints_s1 << endl;}
    if(operation.find("Yevenshift")) {i++;if(_iseven(xvasp.str.kpoints_k2)) xvasp.str.kpoints_s2=0.5; else xvasp.str.kpoints_s2=0.0;change_made=true;if(LDEBUG) cout << "Yevenshift s2=" << xvasp.str.kpoints_s2 << endl;}
    if(operation.find("Zevenshift")) {i++;if(_iseven(xvasp.str.kpoints_k3)) xvasp.str.kpoints_s3=0.5; else xvasp.str.kpoints_s3=0.0;change_made=true;if(LDEBUG) cout << "Zevenshift s3=" << xvasp.str.kpoints_s3 << endl;}
    if(operation.find("Xoddshift")) {i++;if(_isodd(xvasp.str.kpoints_k1)) xvasp.str.kpoints_s1=0.5; else xvasp.str.kpoints_s1=0.0;change_made=true;if(LDEBUG) cout << "Xoddshift s1=" << xvasp.str.kpoints_s1 << endl;}
    if(operation.find("Yoddshift")) {i++;if(_isodd(xvasp.str.kpoints_k2)) xvasp.str.kpoints_s2=0.5; else xvasp.str.kpoints_s2=0.0;change_made=true;if(LDEBUG) cout << "Yoddshift s2=" << xvasp.str.kpoints_s2 << endl;}
    if(operation.find("Zoddshift")) {i++;if(_isodd(xvasp.str.kpoints_k3)) xvasp.str.kpoints_s3=0.5; else xvasp.str.kpoints_s3=0.0;change_made=true;if(LDEBUG) cout << "Zoddshift s3=" << xvasp.str.kpoints_s3 << endl;}
    if(operation.at(0)=='M' || operation.at(0)=='m') {i++;xvasp.str.kpoints_kscheme="Monkhorst-Pack";change_made=true;if(LDEBUG) cout << "Monkhorst-Pack" << endl;} // "Monkhorst-Pack";
    if(operation.at(0)=='G' || operation.at(0)=='g') {i++;xvasp.str.kpoints_kscheme="Gamma";change_made=true;if(LDEBUG) cout << "Gamma" << endl;} // "Gamma";
    if(operation.find("Monkhorst-Pack")) {i++;xvasp.str.kpoints_kscheme="Monkhorst-Pack";change_made=true;if(LDEBUG) cout << "Monkhorst-Pack" << endl;} // "Monkhorst-Pack";
    if(operation.find("Gamma")) {i++;xvasp.str.kpoints_kscheme="Gamma";change_made=true;if(LDEBUG) cout << "Gamma" << endl;} // "Gamma";
    if(xvasp.str.kpoints_k1<1) xvasp.str.kpoints_k1=1;
    if(xvasp.str.kpoints_k2<1) xvasp.str.kpoints_k2=1;
    if(xvasp.str.kpoints_k3<1) xvasp.str.kpoints_k3=1;
    if(LDEBUG) cout << XPID << "KBIN::XVASP_KPOINTS_OPERATION  i=" << i << endl;
    if(LDEBUG) cout << XPID << "KBIN::XVASP_KPOINTS_OPERATION  xvasp.str.kpoints_k*=(" << xvasp.str.kpoints_k1 << "," << xvasp.str.kpoints_k2 << "," << xvasp.str.kpoints_k3 << ")" << endl;
    if(LDEBUG) cout << XPID << "KBIN::XVASP_KPOINTS_OPERATION  xvasp.str.kpoints_s*=(" << xvasp.str.kpoints_s1 << "," << xvasp.str.kpoints_s2 << "," << xvasp.str.kpoints_s3 << ")" << endl;
    if(LDEBUG) cout << XPID << "KBIN::XVASP_KPOINTS_OPERATION  xvasp.str.kpoints_kscheme=*=(" << xvasp.str.kpoints_kscheme << ")" << endl;
    xvasp.str.kpoints_kmax=max(xvasp.str.kpoints_k1,xvasp.str.kpoints_k2,xvasp.str.kpoints_k3);
    xvasp.str.kpoints_kppra=xvasp.str.kpoints_k1*xvasp.str.kpoints_k2*xvasp.str.kpoints_k3*xvasp.str.atoms.size();
    KBIN::XVASP_KPOINTS_KPOINTS(xvasp);
    return TRUE;
  }
}

// bool XVASP_KPOINTS_Kscheme(_xvasp& xvasp,string kscheme) {
//   xvasp.str.kpoints_kscheme=kscheme;
//   if(kscheme.at(0)=='M' || kscheme.at(0)=='m') xvasp.str.kpoints_kscheme="Monkhorst-Pack";// "Monkhorst-Pack";
//   if(kscheme.at(0)=='G' || kscheme.at(0)=='g') xvasp.str.kpoints_kscheme="Gamma";// "Gamma";
//   KBIN::XVASP_KPOINTS_KPOINTS(xvasp);
//   cerr << XPID << "KBIN::XVASP_KPOINTS_Kscheme" << endl << xvasp.KPOINTS.str() << endl;
//   return TRUE;
// }


namespace KBIN {
  bool XVASP_KPOINTS_Fix_KPPRA(_xvasp &xvasp,int NK,ofstream &FileMESSAGE,bool VERBOSE) {
    ostringstream aus;
    string stringKPPRA;
    stringKPPRA=KPPRA(xvasp.str,NK);
    if(VERBOSE) FileMESSAGE << stringKPPRA;
    // VERBOSE on the SCREEN
    if(VERBOSE) {
      aus.precision(5);
      aus << "00000  MESSAGE KPOINTS KPPRA routine ["
        <<  xvasp.str.kpoints_k1 << "," << xvasp.str.kpoints_k2 << "," << xvasp.str.kpoints_k3 << "]=" << xvasp.str.kpoints_k1*xvasp.str.kpoints_k2*xvasp.str.kpoints_k3 << "=["
        <<  modulus(xvasp.str.klattice(1)/((double) xvasp.str.kpoints_k1)) << "," << modulus(xvasp.str.klattice(2)/((double) xvasp.str.kpoints_k2)) << ","
        <<  modulus(xvasp.str.klattice(3)/((double) xvasp.str.kpoints_k3)) << "]   " << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    }
    return TRUE;
  }
}

namespace KBIN {
  bool XVASP_KPOINTS_Fix_KSHIFT(_xvasp &xvasp,_xvasp &rxvasp,bool KAUTOSHIFT,bool VERBOSE) {
    stringstream aus;
    if(KAUTOSHIFT) {
      if(_isodd(xvasp.str.kpoints_k1+rxvasp.str.kpoints_k1)) { xvasp.str.kpoints_s1=0.5; } else { xvasp.str.kpoints_s1=0.0; } 
      if(_isodd(xvasp.str.kpoints_k2+rxvasp.str.kpoints_k2)) { xvasp.str.kpoints_s2=0.5; } else { xvasp.str.kpoints_s2=0.0; } 
      if(_isodd(xvasp.str.kpoints_k3+rxvasp.str.kpoints_k3)) { xvasp.str.kpoints_s3=0.5; } else { xvasp.str.kpoints_s3=0.0; } 
    }
    if(VERBOSE) {
      aus << "00000  MESSAGE KPOINTS KSHIFT=[" << xvasp.str.kpoints_s1 << " " << xvasp.str.kpoints_s2 << " " << xvasp.str.kpoints_s3 << "]" << " " << endl;
      // aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);			
      cout << aus.str();
    }
    return TRUE;
  }
}

namespace KBIN {
  bool XVASP_KPOINTS_Fix_KPOINTS(_xvasp &xvasp,int NK,ofstream &FileMESSAGE,bool VERBOSE) {
    stringstream aus;
    // create KPOINTS stringstream from KPPRA
    xvasp.KPOINTS.str(std::string()); xvasp.KPOINTS.clear();
    xvasp.KPOINTS_orig.str(std::string());xvasp.KPOINTS_orig.clear();
    xvasp.KPOINTS << "KPOINTS File, automatically generated by KPPRA - aflow with Kpoints=" << xvasp.str.kpoints_kppra  << "  [KPPRA=" << NK << "]" << endl;
    xvasp.KPOINTS << xvasp.str.kpoints_mode << endl;  // MODE AUTO 0
    xvasp.KPOINTS << xvasp.str.kpoints_kscheme << endl;
    xvasp.KPOINTS << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << endl;
    xvasp.KPOINTS.precision(3);
    xvasp.KPOINTS << xvasp.str.kpoints_s1 << " " << xvasp.str.kpoints_s2 << " " << xvasp.str.kpoints_s3 << endl;

    aurostd::StringstreamClean(aus);
    if(VERBOSE) {
      aus << "00000  MESSAGE KPOINTS K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " " << endl;
      aus << "00000  MESSAGE KPOINTS Kpoints=" << xvasp.str.kpoints_kppra  << "  [KPPRA=" << NK  << "]    with " << xvasp.str.kpoints_kscheme << "   " << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);			
    }
    // KPOINTS done
    xvasp.KPOINTS_orig << xvasp.KPOINTS.str();
    return TRUE;
  }
}

namespace KBIN {
  void XVASP_string2numbers(_xvasp& xvasp) {
    string strline;
    int imax=aurostd::GetNLinesString(xvasp.KPOINTS.str());
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(xvasp.KPOINTS,i);
      if(i==3) xvasp.str.kpoints_kscheme=strline;
      if(i==4) stringstream(strline) >> xvasp.str.kpoints_k1 >> xvasp.str.kpoints_k2 >> xvasp.str.kpoints_k3;
    }
  }
}
namespace KBIN {
  void XVASP_numbers2string(_xvasp& xvasp) {
    string strline,stringKPOINTS=xvasp.KPOINTS.str();
    int imax=aurostd::GetNLinesString(stringKPOINTS);
    xvasp.KPOINTS.str(std::string());
    for(int i=1;i<=imax;i++) {
      strline=aurostd::GetLineString(stringKPOINTS,i);
      if(i!=3 && i!=4 && i!=imax) xvasp.KPOINTS << strline << endl;
      if(i==3) xvasp.KPOINTS << xvasp.str.kpoints_kscheme << endl;
      if(i==4) xvasp.KPOINTS << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << endl;
      if(i==imax && strline.length()!=0) xvasp.KPOINTS << strline << endl;
    }
  }
}
// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// FIX MODIFICATIONS

namespace KBIN {
  void XVASP_Afix_Clean(_xvasp& xvasp,string preserve_name) {
    aurostd::execute("cat "+xvasp.Directory+"/vasp.out >> "+xvasp.Directory+"/"+preserve_name);
    deque<string> vrem;aurostd::string2tokens("aflow.tmp,CHG,CONTCAR,DOSCAR,EIGENVAL,IBZKPT,OUTCAR,OSZICAR,PCDAT,WAVECAR,aflow.qsub*,XDATCAR,vasprun.xml,vasp.out,core*",vrem,",");
    for (uint i=0;i<vrem.size();i++)
      aurostd::RemoveFile(xvasp.Directory+"/"+vrem.at(i));
    if(!xvasp.aopts.flag("FLAG::CHGCAR_PRESERVED")) aurostd::RemoveFile(xvasp.Directory+"/CHGCAR"); 
  }
}
namespace KBIN {
  void XVASP_Afix_ROTMAT(_xvasp& xvasp,int mode,bool verbose,_aflags &aflags,ofstream &FileMESSAGE) {
    // cerr << XPID << "KBIN::XVASP_Afix_ROTMAT mode=" << mode << endl;
    // mode=1 put gamma and even
    // mode=2 put gamma and odd
    // mode=3 put max
    // mode=4 put odd
    // mode=5 put even
    if(xvasp.aopts.flag("FLAG::KPOINTS_PRESERVED")) return; // don`t touch kpoints
    // AFIX *************************************
    int kmax=xvasp.str.kpoints_kmax;
    // int k1=xvasp.str.kpoints_k1,k2=xvasp.str.kpoints_k2,k3=xvasp.str.kpoints_k3,kppra=xvasp.str.kpoints_kppra;
    // double s1=xvasp.str.kpoints_s1,s2=xvasp.str.kpoints_s2,s3=xvasp.str.kpoints_s3;
    bool LQUIET=!verbose;
    stringstream aus_exec;
    ostringstream aus;
    bool DEBUG_KBIN_XVASP_Afix_ROTMAT=FALSE;
    // mode *************************************
    if(mode<=0) {	
      KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run Reciprocal Rotation Matrix problems (XXXX) ");
      aus << "00000  FIX RECIPROCAL - KLATTICE SYMMETRY MISMATCH FIX (XXXX) - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,LQUIET);
    }
    // mode *************************************
    if(mode==1) {
      KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run Reciprocal Rotation Matrix problems (KPOINTS_KSHIFT_EVEN) ");
      aus << "00000  FIX RECIPROCAL - KLATTICE SYMMETRY MISMATCH FIX (KPOINTS_KSHIFT_EVEN) - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,LQUIET);
      KBIN::XVASP_KPOINTS_OPERATION(xvasp,"Xevenshift,Yevenshift,Zevenshift");   // KBIN::XVASP_KPOINTS_Kshift_Gamma_EVEN(_xvasp& xvasp);
      aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS"));
      if(DEBUG_KBIN_XVASP_Afix_ROTMAT) aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS_KBIN_XVASP_KPOINTS_Kshift_Gamma_EVEN"));
    }
    // mode *************************************
    if(mode==2) {
      KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run Reciprocal Rotation Matrix problems (KPOINTS_KSHIFT_ODD) ");
      aus << "00000  FIX RECIPROCAL - KLATTICE SYMMETRY MISMATCH FIX (KPOINTS_KSHIFT_ODD) - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,LQUIET);
      KBIN::XVASP_KPOINTS_OPERATION(xvasp,"Xoddshift,Yoddshift,Zoddshift");   // KBIN::XVASP_KPOINTS_Kshift_Gamma_ODD(xvasp);
      aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS"));
      if(DEBUG_KBIN_XVASP_Afix_ROTMAT) aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS_KBIN_XVASP_KPOINTS_Kshift_Gamma_ODD"));
    }
    // mode *************************************
    if(mode==3) {
      KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run Reciprocal Rotation Matrix problems (KPOINTS_MAX) ");
      aus << "00000  FIX RECIPROCAL - KLATTICE SYMMETRY MISMATCH FIX (KPOINTS_MAX) - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,LQUIET);
      xvasp.str.kpoints_k1=kmax;xvasp.str.kpoints_k2=kmax;xvasp.str.kpoints_k3=kmax;
      xvasp.str.kpoints_kppra=kmax*kmax*kmax*xvasp.str.atoms.size();
      xvasp.str.kpoints_s1=0.0;xvasp.str.kpoints_s2=0.0;xvasp.str.kpoints_s3=0.0;
      // if(_isodd(xvasp.str.kpoints_k1))  {xvasp.str.kpoints_s1=0.0;xvasp.str.kpoints_s2=0.0;xvasp.str.kpoints_s3=0.0;}
      // if(_iseven(xvasp.str.kpoints_k1)) {xvasp.str.kpoints_s1=0.5;xvasp.str.kpoints_s2=0.5;xvasp.str.kpoints_s3=0.5;}
      KBIN::XVASP_KPOINTS_KPOINTS(xvasp);
      aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS"));
      if(DEBUG_KBIN_XVASP_Afix_ROTMAT) aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS_KBIN_XVASP_KPOINTS_KPOINTS_MAX"));
    }
    // mode *************************************
    if(mode==4) {
      KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run Reciprocal Rotation Matrix problems (KPOINTS_ODD) ");
      aus << "00000  FIX RECIPROCAL - KLATTICE SYMMETRY MISMATCH FIX (KPOINTD_ODD) - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,LQUIET);
      KBIN::XVASP_KPOINTS_OPERATION(xvasp,"Xodd,Yodd,Zodd");  // KBIN::XVASP_KPOINTS_ODD(xvasp);
      aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS"));
      if(DEBUG_KBIN_XVASP_Afix_ROTMAT) aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS_KBIN_XVASP_KPOINTS_ODD"));
    }
    // mode *************************************
    if(mode==5) {
      KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run Reciprocal Rotation Matrix problems (KPOINTS_EVEN) ");
      aus << "00000  FIX RECIPROCAL - KLATTICE SYMMETRY MISMATCH FIX (KPOINTD_EVEN) - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,LQUIET);
      KBIN::XVASP_KPOINTS_OPERATION(xvasp,"Xeven,Yeven,Zeven"); //    KBIN::XVASP_KPOINTS_EVEN(xvasp);
      aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS"));
      if(DEBUG_KBIN_XVASP_Afix_ROTMAT) aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS_KBIN_XVASP_KPOINTS_EVEN"));
    }
    // mode *************************************
    if(mode==6) { // SYM is the desperate case
      KBIN::VASP_Error(xvasp,"WWWWW  ERROR KBIN::VASP_Run Reciprocal Rotation Matrix problems (ISYM=0) ");
      aus << "00000  FIX RECIPROCAL - KLATTICE SYMMETRY MISMATCH FIX (ISYM=0) - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,LQUIET);
      _kflags kflags;_vflags vflags; // temporary
      KBIN::XVASP_INCAR_PREPARE_GENERIC("SYM",xvasp,vflags,"",0,0.0,OFF);
      // [OBSOLETE]      KBIN::XVASP_INCAR_SYM(xvasp,OFF,TRUE);
      aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"));
      KBIN::XVASP_KPOINTS_OPERATION(xvasp,"Xodd,Yodd,Zodd");  // KBIN::XVASP_KPOINTS_ODD(xvasp);
      aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS"));
      if(DEBUG_KBIN_XVASP_Afix_ROTMAT) aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS_KBIN_XVASP_KPOINTS_ODD"));
    }
    // reload to restart ---------------------------------

    xvasp.KPOINTS_orig.str(std::string()); xvasp.KPOINTS_orig << xvasp.KPOINTS.str();
    xvasp.KPOINTS.str(std::string()); xvasp.KPOINTS << aurostd::file2string(xvasp.Directory+"/KPOINTS");
    xvasp.INCAR_orig.str(std::string()); xvasp.INCAR_orig << xvasp.INCAR.str();
    xvasp.INCAR.str(std::string()); xvasp.INCAR << aurostd::file2string(xvasp.Directory+"/INCAR");
    // clean to restart ----------------------------------    
    aurostd::execute("cat "+xvasp.Directory+"/vasp.out >> "+xvasp.Directory+"/aflow.error.rotmat");
    deque<string> vrem;aurostd::string2tokens("aflow.tmp,CHG,CONTCAR,DOSCAR,EIGENVAL,IBZKPT,OUTCAR,OSZICAR,PCDAT,WAVECAR,aflow.qsub*,XDATCAR,vasprun.xml,vasp.out,core*",vrem,",");
    for (uint i=0;i<vrem.size();i++)
      aurostd::RemoveFile(xvasp.Directory+"/"+vrem.at(i));
    if(!xvasp.aopts.flag("FLAG::CHGCAR_PRESERVED")) aurostd::RemoveFile(xvasp.Directory+"/CHGCAR");
  }
}

namespace KBIN {
  void XVASP_Afix_NBANDS(_xvasp& xvasp,int& nbands,bool VERBOSE) {
    string soliloquy=XPID+"XVASP_Afix_NBANDS():";
    stringstream aus_exec;
    // [OBSOLETE] char buffer[BUFFER_MAXLEN];
    ifstream FileAUS;
    string FileNameAUS;
    // [OBSOLETE] string word;

    // get NBANDS from OUTCAR
    aurostd::RemoveFile(xvasp.Directory+"/aflow.tmp");
    nbands=0;
    xOUTCAR OUTCAR_NBANDS(xvasp.Directory+"/OUTCAR");
    nbands=OUTCAR_NBANDS.NBANDS;

    if(nbands==0) {                                            // GET BANDS FROM INCAR
      if(VERBOSE) cerr << XPID << "KBIN::XVASP_Afix_NBANDS: GET NBANDS FROM INCAR" << endl;
      //aus_exec << "cat INCAR | grep -v \"\\#\" | grep NBANDS | sed \"s/=//g\" | sed \"s/NBANDS//g\" > aflow.tmp " << endl;
      aurostd::execute("cat "+xvasp.Directory+"/INCAR | grep -v \"\\#\" | grep NBANDS | sed \"s/=//g\" | sed \"s/NBANDS//g\" > "+xvasp.Directory+"/aflow.tmp");
      FileNameAUS=xvasp.Directory+"/aflow.tmp";
      FileAUS.open(FileNameAUS.c_str(),std::ios::in);
      if(!FileAUS) {
        throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"cannot get NBANDS from INCAR",_FILE_CORRUPT_); //CO20200624
      }
      FileAUS >> nbands;
      FileAUS.clear();FileAUS.close();
      aurostd::RemoveFile(xvasp.Directory+"/aflow.tmp");
    }
    //  cerr << XPID << "KBIN::XVASP_Afix_NBANDS: nbands=" << nbands << endl;
    if(nbands<5) {
      nbands=KBIN::XVASP_INCAR_GetNBANDS(xvasp,TRUE);
    } else {
      nbands=nbands+20+nbands/5; // why did I remove it ?
    }
    //  cerr << XPID << "KBIN::XVASP_Afix_NBANDS: nbands=" << nbands << endl;

    // nbands++;
    //    cout << "DEBUG KBIN::XVASP_Afix_NBANDS nbands=" << nbands << " (20150602) directory=" << xvasp.Directory << endl;
    //    cerr << "DEBUG KBIN::XVASP_Afix_NBANDS nbands=" << nbands << " (20150602) directory=" << xvasp.Directory << endl;
    aus_exec << "cd " << xvasp.Directory << endl;
    // aus_exec << "cat INCAR | grep -v NBANDS > aflow.tmp && mv aflow.tmp INCAR" << endl;
    aus_exec << "cat INCAR | sed \"s/NBANDS/#NBANDS/g\" > aflow.tmp && mv aflow.tmp INCAR" << endl;
    // aus_exec << "echo \"# Performing KBIN::XVASP_Afix_NBANDS (new) [AFLOW] begin\" >> INCAR " << endl;
    aus_exec << "echo \"#fixed nbands=\"" << nbands << " >> INCAR " << endl;
    aus_exec << "echo \"NBANDS=\"" << nbands << " >> INCAR " << endl;
    // aus_exec << "echo \"# Performing KBIN::XVASP_Afix_NBANDS (new) [AFLOW] end\" >> INCAR " << endl;
    aurostd::execute(aus_exec);
    // fix aflowlin
    aus_exec << "cd " << xvasp.Directory << endl;
    aus_exec << "cat " << _AFLOWIN_ << " | sed \"s/\\[VASP_FORCE_OPTION\\]NBANDS/#\\[VASP_FORCE_OPTION\\]NBANDS/g\" | sed \"s/##\\[/#\\[/g\" > aflow.tmp && mv aflow.tmp " << _AFLOWIN_ << "" << endl;
    aus_exec << "echo \"[VASP_FORCE_OPTION]NBANDS=" << nbands << "      // Self Correction\"" << " >> " << _AFLOWIN_ << " " << endl;
    aurostd::execute(aus_exec);
    // reload to restart ---------------------------------
    xvasp.INCAR_orig.str(std::string()); xvasp.INCAR_orig << xvasp.INCAR.str();
    xvasp.INCAR.str(std::string()); xvasp.INCAR << aurostd::file2string(xvasp.Directory+"/INCAR");
    // clean to restart ----------------------------------
    KBIN::XVASP_Afix_Clean(xvasp,"aflow.error.nbands");
  }

  void XVASP_Afix_POTIM(_xvasp& xvasp,double& potim,bool VERBOSE) {
    string soliloquy=XPID+"KBIN::XVASP_Afix_POTIM():";
    stringstream aus_exec;
    ifstream FileAUS;
    string FileNameAUS;

    // get POTIM from OUTCAR
    aurostd::RemoveFile(xvasp.Directory+"/aflow.tmp");
    potim=0;
    xOUTCAR OUTCAR_POTIM(xvasp.Directory+"/OUTCAR");
    potim=OUTCAR_POTIM.POTIM;
    if(potim==0) {                                            // GET BANDS FROM INCAR
      if(VERBOSE) cerr << XPID << "KBIN::XVASP_Afix_POTIM: GET POTIM FROM INCAR" << endl;
      aurostd::execute("cat "+xvasp.Directory+"/INCAR | grep -v \"\\#\" | grep POTIM | sed \"s/=//g\" | sed \"s/POTIM//g\" > "+xvasp.Directory+"/aflow.tmp");
      FileNameAUS=xvasp.Directory+"/aflow.tmp";
      FileAUS.open(FileNameAUS.c_str(),std::ios::in);
      if(!FileAUS) {
        throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"cannot get POTIM from INCAR",_FILE_CORRUPT_); //CO20200624
      }
      FileAUS >> potim;
      FileAUS.clear();FileAUS.close();
      aurostd::RemoveFile(xvasp.Directory+"/aflow.tmp");
    }
    if(VERBOSE) cerr << XPID << "KBIN::XVASP_Afix_POTIM: potim=" << potim << endl;
    potim=potim/2.0;
    if(potim<0.01) potim=0.01;
    if(VERBOSE) cerr << XPID << "KBIN::XVASP_Afix_POTIM: potim=" << potim << endl;

    aus_exec << "cd " << xvasp.Directory << endl;
    aus_exec << "cat INCAR | sed \"s/POTIM/#POTIM/g\" > aflow.tmp && mv aflow.tmp INCAR" << endl;
    aus_exec << "echo \"#fixed potim=\"" << potim << " >> INCAR " << endl;
    aus_exec << "echo \"POTIM=\"" << potim << " >> INCAR " << endl;
    aurostd::execute(aus_exec);
    // fix aflowlin
    aus_exec << "cd " << xvasp.Directory << endl;
    aus_exec << "cat " << _AFLOWIN_ << " | sed \"s/\\[VASP_FORCE_OPTION\\]POTIM/#\\[VASP_FORCE_OPTION\\]POTIM/g\" | sed \"s/##\\[/#\\[/g\" > aflow.tmp && mv aflow.tmp " << _AFLOWIN_ << "" << endl;
    aus_exec << "echo \"[VASP_FORCE_OPTION]POTIM=" << potim << "      // Self Correction\"" << " >> " << _AFLOWIN_ << " " << endl;
    aurostd::execute(aus_exec);
    // reload to restart ---------------------------------
    xvasp.INCAR_orig.str(std::string()); xvasp.INCAR_orig << xvasp.INCAR.str();
    xvasp.INCAR.str(std::string()); xvasp.INCAR << aurostd::file2string(xvasp.Directory+"/INCAR");
  }
}

namespace KBIN {
  double XVASP_Afix_GENERIC(string mode,_xvasp& xvasp,_kflags& kflags,_vflags& vflags,double param_double,int param_int) {
    string soliloquy=XPID+"KBIN::XVASP_Afix_GENERIC():";
    string file_error,file_stream;
    stringstream aus_exec;
    double out=0.0;
    bool rewrite_incar=FALSE,rewrite_poscar=FALSE,rewrite_kpoints=FALSE;
    bool reload_incar=FALSE,reload_poscar=FALSE,reload_kpoints=FALSE;
    //  cerr << XPID << "KBIN::XVASP_Afix_GENERIC=" << mode << endl;

    vflags.KBIN_VASP_INCAR_VERBOSE=TRUE; 
    aus_exec << "cd " << xvasp.Directory << endl; //CO
    aus_exec << "echo \"# Performing XVASP_Afix_GENERIC (" << mode << ") [AFLOW] end\" >> INCAR " << endl;
    aurostd::execute(aus_exec);

    if(mode=="SYMPREC" || mode=="SGRCON" || mode=="INVGRP") {
      file_error="aflow.error.symprec";
      reload_incar=TRUE;
      aus_exec << "cd " << xvasp.Directory << endl;
      aus_exec << "cat INCAR | sed \"s/SYMPREC/#SYMPREC/g\" > aflow.tmp && mv aflow.tmp INCAR" << endl; // remove SYMPREC
      if(vflags.KBIN_VASP_INCAR_VERBOSE) aus_exec << "echo \"# Performing KBIN::XVASP_Afix (" << mode << ") [AFLOW] begin\" >> INCAR " << endl;
      // aus_exec << "echo \"#fixed symprec \" >> INCAR " << endl;
      // aus_exec << "echo \"SYMPREC=1e-6\" >> INCAR " << endl;
      aus_exec << "echo \"SYMPREC=1e-7    # (fix=" << mode << "\" >> INCAR " << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) aus_exec << "echo \"# Performing KBIN::XVASP_Afix (" << mode << ") [AFLOW] end\" >> INCAR " << endl;
      aurostd::execute(aus_exec);
    }
    if(mode=="READ_KPOINTS_RD_SYM") {
      file_error="aflow.error.read_kpoints_rd_sym";
      reload_incar=TRUE;
      aus_exec << "cd " << xvasp.Directory << endl;
      aus_exec << "cat INCAR | sed \"s/ISYM/#ISYM/g\" > aflow.tmp && mv aflow.tmp INCAR" << endl; // remove SYMPREC
      if(vflags.KBIN_VASP_INCAR_VERBOSE) aus_exec << "echo \"# Performing KBIN::XVASP_Afix (" << mode << ") [AFLOW] begin\" >> INCAR " << endl;
      aus_exec << "echo \"ISYM=0    # (fix=" << mode << "\" >> INCAR " << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) aus_exec << "echo \"# Performing KBIN::XVASP_Afix (" << mode << ") [AFLOW] end\" >> INCAR " << endl;
      aurostd::execute(aus_exec);
    }
    if(mode=="IBZKPT") {
      file_error="aflow.error.ibzkpt";
      reload_kpoints=TRUE;
      if(xvasp.aopts.flag("FLAG::KPOINTS_PRESERVED")) return 0.0; // don`t touch kpoints
      if(vflags.KBIN_VASP_INCAR_VERBOSE) aus_exec << "KBIN::XVASP_Afix_GENERIC MAKING (" << mode << ") BEGIN" << endl;
      KBIN::XVASP_KPOINTS_OPERATION(xvasp,"Xodd,Yodd,Zodd");  // KBIN::XVASP_KPOINTS_ODD(xvasp);  // this should put the origin in GAMMA ??
      if(vflags.KBIN_VASP_INCAR_VERBOSE) aus_exec << "KBIN::XVASP_Afix_GENERIC MAKING (" << mode << ") END" << endl;
      xvasp.str.kpoints_s1=0.0;xvasp.str.kpoints_s2=0.0;xvasp.str.kpoints_s3=0.0;
      rewrite_kpoints=TRUE;
    }
    if(mode=="GAMMA_SHIFT") {
      file_error="aflow.error.gamma_shift";
      reload_kpoints=TRUE;
      if(xvasp.aopts.flag("FLAG::KPOINTS_PRESERVED")) return 0.0; // don`t touch kpoints
      if(vflags.KBIN_VASP_INCAR_VERBOSE) aus_exec << "KBIN::XVASP_Afix_GENERIC MAKING (" << mode << ") BEGIN" << endl;
      KBIN::XVASP_KPOINTS_OPERATION(xvasp,"Gamma");
      if(vflags.KBIN_VASP_INCAR_VERBOSE) aus_exec << "KBIN::XVASP_Afix_GENERIC MAKING (" << mode << ") END" << endl;
      //    xvasp.str.kpoints_s1=0.0;xvasp.str.kpoints_s2=0.0;xvasp.str.kpoints_s3=0.0;
      rewrite_kpoints=TRUE;
    }
    if(mode=="MPICH11") {
      file_error="aflow.error.mpich11";
      kflags.KBIN_MPI_OPTIONS=string("ulimit -s unlimited ");//+string(" && ")+kflags.KBIN_MPI_OPTIONS;
    }
    if(mode=="MPICH139") {
      file_error="aflow.error.mpich139";
      reload_kpoints=TRUE;
      if(xvasp.aopts.flag("FLAG::KPOINTS_PRESERVED")) return 0.0; // don`t touch kpoints
      if(vflags.KBIN_VASP_INCAR_VERBOSE) aus_exec << "KBIN::XVASP_Afix_GENERIC MAKING (" << mode << ") BEGIN" << endl;
      KBIN::XVASP_KPOINTS_OPERATION(xvasp,"X--,Y--,Z--");  // reduce KPOINTS withoug bugging the origin
      KBIN::XVASP_KPOINTS_OPERATION(xvasp,"X--,Y--,Z--");  // reduce KPOINTS withoug bugging the origin
      if(vflags.KBIN_VASP_INCAR_VERBOSE) aus_exec << "KBIN::XVASP_Afix_GENERIC MAKING (" << mode << ") END" << endl;
      // xvasp.str.kpoints_s1=0.0;xvasp.str.kpoints_s2=0.0;xvasp.str.kpoints_s3=0.0;
      rewrite_kpoints=TRUE;
      kflags.KBIN_MPI_OPTIONS=string("ulimit -s unlimited ");//+string("\n")+kflags.KBIN_MPI_OPTIONS;
    }
    if(mode=="NKXYZ_IKPTD") {
      file_error="aflow.error.nkxyz_ikptd";
      reload_kpoints=TRUE;
      // if(xvasp.aopts.flag("FLAG::KPOINTS_PRESERVED")) return 0.0; // don`t touch kpoints MUST FIX ANYWAY
      KBIN::XVASP_KPOINTS_OPERATION(xvasp,"X--,Y--,Z--");  // this should put the origin in GAMMA ??
      rewrite_kpoints=TRUE;
    }
    if(mode=="EDDRMM") {
      file_error="aflow.error.eddrmm";
      reload_kpoints=TRUE;
      if(xvasp.aopts.flag("FLAG::KPOINTS_PRESERVED")) return 0.0; // don`t touch kpoints
      KBIN::XVASP_KPOINTS_OPERATION(xvasp,"Xodd,Yodd,Zodd");  // KBIN::XVASP_KPOINTS_ODD(xvasp);  // this should put the origin in GAMMA ??
      xvasp.str.kpoints_s1=0.0;xvasp.str.kpoints_s2=0.0;xvasp.str.kpoints_s3=0.0; // you go on gamma ONLY if you put shift =0 !!
      rewrite_kpoints=TRUE;
    }
    if(mode=="LREAL") {
      file_error="aflow.error.lreal";
      reload_incar=TRUE;
      aus_exec << "cd " << xvasp.Directory << endl;
      aus_exec << "cat INCAR | sed \"s/LREAL/#LREAL/g\" > aflow.tmp && mv aflow.tmp INCAR" << endl; // remove LREAL
      if(vflags.KBIN_VASP_INCAR_VERBOSE) aus_exec << "echo \"# Performing KBIN::XVASP_Afix (" << mode << ") [AFLOW] begin\" >> INCAR " << endl;
      aus_exec << "echo \"LREAL=.TRUE.    # (fix=" << mode << "\" >> INCAR " << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) aus_exec << "echo \"# Performing KBIN::XVASP_Afix (" << mode << ") [AFLOW] end\" >> INCAR " << endl;
      aurostd::execute(aus_exec);
    }
    if(mode=="EXCCOR") {
      file_error="aflow.error.exccor";
      reload_poscar=TRUE;
      if(xvasp.aopts.flag("FLAG::POSCAR_PRESERVED")) return 0.0; // don`t touch poscar
      aurostd::stringstream2file(xvasp.POSCAR,string(xvasp.Directory+"/POSCAR.orig"));
      xvasp.str.scale=xvasp.str.scale*1.0627;   // cubic root of 1.2
      xvasp.str.scale=xvasp.str.scale*1.0627;   // cubic root of 1.2
      xvasp.str.scale=xvasp.str.scale*1.0627;   // cubic root of 1.2
      KBIN::VASP_Produce_POSCAR(xvasp);
      rewrite_poscar=TRUE;
    }
    if(mode=="BRMIX") {
      file_error="aflow.error.brmix";
      reload_incar=TRUE;
      vflags.KBIN_VASP_FORCE_OPTION_ALGO.xscheme="VERYFAST";
      KBIN::XVASP_INCAR_PREPARE_GENERIC("ALGO",xvasp,vflags,"",0,0.0,FALSE);
      // [OBSOLETE] KBIN::XVASP_INCAR_ALGO(xvasp,vflags);
      rewrite_incar=TRUE;
    }
    if(mode=="DAV") {
      file_error="aflow.error.dav";
      reload_incar=TRUE;
      KBIN::XVASP_INCAR_PREPARE_GENERIC("IALGO",xvasp,vflags,"",48,0.0,FALSE);
      rewrite_incar=TRUE;
    }
    if(mode=="EDDDAV") {
      file_error="aflow.error.edddav";
      reload_incar=TRUE;
      KBIN::XVASP_INCAR_PREPARE_GENERIC("IALGO",xvasp,vflags,"",48,0.0,FALSE);   // same as DAV
      rewrite_incar=TRUE;
    }
    if(mode=="EFIELD_PEAD") {
      file_error="aflow.error.efield_pead";
      reload_incar=TRUE;
      KBIN::XVASP_INCAR_EFIELD_PEAD(xvasp,0.20,vflags.KBIN_VASP_INCAR_VERBOSE);   // same as DAV
      rewrite_incar=TRUE;
    }
    if(mode=="ZPOTRF") {
      file_error="aflow.error.zpotrf";
      reload_poscar=TRUE;
      aus_exec << "cd " << xvasp.Directory << endl;
      aus_exec << "cat " << _AFLOWIN_ << " | grep -v \"CONVERT_UNIT_CELL\" > aflow.tmp" << endl;
      aus_exec << "echo \"[VASP_FORCE_OPTION]CONVERT_UNIT_CELL=SCONV   // KBIN::XVASP_Afix_ZPOTRF \" >> aflow.tmp && mv aflow.tmp " << _AFLOWIN_ << " " << endl;
      aurostd::execute(aus_exec);
    }
    if(mode=="ZPOTRF_POTIM") {
      file_error="aflow.error.zpotrf_potim";
      reload_incar=TRUE;
      double potim;
      KBIN::XVASP_Afix_POTIM(xvasp,potim,vflags.KBIN_VASP_INCAR_VERBOSE);
      rewrite_incar=TRUE;
    }
    if(mode=="NATOMS") {
      file_error="aflow.error.natoms";
      reload_poscar=TRUE;
      if(xvasp.aopts.flag("FLAG::POSCAR_PRESERVED")) return 0.0; // don`t touch poscar
      aurostd::stringstream2file(xvasp.POSCAR,string(xvasp.Directory+"/POSCAR.orig"));
      xvasp.str.scale=xvasp.str.scale*pow(2.0,1.0/3.0);   // cubic root of 1.2
      KBIN::VASP_Produce_POSCAR(xvasp);
      rewrite_poscar=TRUE;
    }
    if(mode=="PSMAXN") {
      double enmax=param_double;
      file_error="aflow.error.psmaxn";
      reload_incar=TRUE;
      aurostd::file2string(xvasp.Directory+"/INCAR",file_stream);
      enmax=0.0;
      if(aurostd::substring2bool(file_stream,"ENMAX=",TRUE)) {  // GET ENMAX FROM OUTCAR
        enmax=aurostd::substring2utype<double>(file_stream,"ENMAX=",TRUE);
      } else {enmax=500.0;}
      // enmax=enmax*0.99; // reduce 1%.... if enough
      enmax=enmax*0.97; // reduce 3%.... if enough
      // cerr << "enmax=" << enmax << endl;
      aus_exec << "cd " << xvasp.Directory << endl;
      aus_exec << "cat INCAR | sed \"s/ENMAX/#ENMAX/g\" > aflow.tmp && mv aflow.tmp INCAR" << endl; // remove ENMAX
      aus_exec << "cat INCAR | sed \"s/LREAL/#LREAL/g\" > aflow.tmp && mv aflow.tmp INCAR" << endl; // remove LREAL
      if(vflags.KBIN_VASP_INCAR_VERBOSE) aus_exec << "echo \"# Performing KBIN::XVASP_Afix (" << mode << ") [AFLOW] begin\" >> INCAR " << endl;
      aus_exec << "echo \"#fixed enmax=\"" << enmax << " >> INCAR " << endl;
      aus_exec << "echo \"ENMAX=\"" << enmax << " >> INCAR " << endl;
      aus_exec << "echo \"LREAL=.TRUE.\" >> INCAR " << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) aus_exec << "echo \"# Performing KBIN::XVASP_Afix (" << mode << ") [AFLOW] end\" >> INCAR " << endl;
      aurostd::execute(aus_exec);
      out=enmax;
    }
    if(mode=="NPAR") {
      file_error="aflow.error.npar";
      reload_incar=TRUE;
      // fix aflowlin
      aus_exec << "cd " << xvasp.Directory << endl;
      aus_exec << "cat INCAR | sed \"s/NPAR/#NPAR/g\" > aflow.tmp && mv aflow.tmp INCAR" << endl; // remove NPAR
      if(vflags.KBIN_VASP_INCAR_VERBOSE) aus_exec << "echo \"# Performing KBIN::XVASP_Afix (" << mode << ") [AFLOW] begin\" >> INCAR " << endl;
      aus_exec << "echo \"NPAR=1 \" >> INCAR " << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) aus_exec << "echo \"# Performing KBIN::XVASP_Afix (" << mode << ") [AFLOW] end\" >> INCAR " << endl;
      aurostd::execute(aus_exec);
    }
    if(mode=="NPARC") {
      file_error="aflow.error.nparc";
      reload_incar=TRUE;
      // fix aflowlin
      aus_exec << "cd " << xvasp.Directory << endl;
      aus_exec << "cat INCAR | sed \"s/NPAR/#NPAR/g\" > aflow.tmp && mv aflow.tmp INCAR" << endl; // remove NPAR
      if(vflags.KBIN_VASP_INCAR_VERBOSE) aus_exec << "echo \"# Performing KBIN::XVASP_Afix (" << mode << ") [AFLOW] begin\" >> INCAR " << endl;
      // remove NPAR
      //    if(param_int>1) {aus_exec << "echo \"NPAR=" << param_int << " \" >> INCAR " << endl;}
      // else { aus_exec << "echo \"NPAR=4 \" >> INCAR " << endl;}
      aus_exec << "echo \"NPAR=4 \" >> INCAR " << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) aus_exec << "echo \"# Performing KBIN::XVASP_Afix (" << mode << ") [AFLOW] end\" >> INCAR " << endl;
      aurostd::execute(aus_exec);
    }
    if(mode=="NPARN") {
      file_error="aflow.error.nparn";
      reload_incar=TRUE;
      // fix aflowlin
      aus_exec << "cd " << xvasp.Directory << endl;
      aus_exec << "cat INCAR | sed \"s/NPAR/#NPAR/g\" > aflow.tmp && mv aflow.tmp INCAR" << endl; // remove NPAR
      if(vflags.KBIN_VASP_INCAR_VERBOSE) aus_exec << "echo \"# Performing KBIN::XVASP_Afix (" << mode << ") [AFLOW] begin\" >> INCAR " << endl;
      // remove NPAR
      //    if(param_int>1) {aus_exec << "echo \"NPAR=" << param_int << " \" >> INCAR " << endl;}
      // else { aus_exec << "echo \"NPAR=4 \" >> INCAR " << endl;}
      aus_exec << "echo \"NPAR=4 \" >> INCAR " << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) aus_exec << "echo \"# Performing KBIN::XVASP_Afix (" << mode << ") [AFLOW] end\" >> INCAR " << endl;
      aurostd::execute(aus_exec);
    }
    if(mode=="NPAR_REMOVE") {
      file_error="aflow.error.npar_remove";
      reload_incar=TRUE;
      // fix aflowlin
      aus_exec << "cd " << xvasp.Directory << endl;
      aus_exec << "cat INCAR | sed \"s/NPAR/#REMOVED NPAR/g\" > aflow.tmp && mv aflow.tmp INCAR" << endl; // remove NPAR
      if(vflags.KBIN_VASP_INCAR_VERBOSE) aus_exec << "echo \"# Performing KBIN::XVASP_Afix (" << mode << ") [AFLOW] begin\" >> INCAR " << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) aus_exec << "echo \"# Performing KBIN::XVASP_Afix (" << mode << ") [AFLOW] end\" >> INCAR " << endl;
      aurostd::execute(aus_exec);
    }
    if(mode=="CSLOSHING") {
      file_error="aflow.error.csloshing";
      reload_incar=TRUE;
      vflags.KBIN_VASP_FORCE_OPTION_ALGO.xscheme="VERYFAST";
      vflags.KBIN_VASP_FORCE_OPTION_ALGO.xscheme="NORMAL";
      KBIN::XVASP_INCAR_PREPARE_GENERIC("ALGO",xvasp,vflags,"",0,0.0,FALSE);
      // [OBSOLETE] KBIN::XVASP_INCAR_ALGO(xvasp,vflags);
      aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"));
      // fix aflowlin
      aus_exec << "cd " << xvasp.Directory << endl;
      aus_exec << "cat " << _AFLOWIN_ << " | sed \"s/\\[VASP_FORCE_OPTION\\]ALGO/#\\[VASP_FORCE_OPTION\\]ALGO/g\" | sed \"s/##\\[/#\\[/g\" > aflow.tmp && mv aflow.tmp " << _AFLOWIN_ << "" << endl;
      aus_exec << "echo \"[VASP_FORCE_OPTION]ALGO=NORMAL      // Self Correction\"" << " >> " << _AFLOWIN_ << " " << endl;
      aurostd::execute(aus_exec);
    }
    if(mode=="DENTET") {
      file_error="aflow.error.dentet";
      reload_incar=TRUE;
      // fix aflowlin
      aus_exec << "cd " << xvasp.Directory << endl;
      aus_exec << "cat INCAR | sed \"s/ISMEAR/#ISMEAR/g\" > aflow.tmp && mv aflow.tmp INCAR" << endl; // remove ISMEAR
      if(vflags.KBIN_VASP_INCAR_VERBOSE) aus_exec << "echo \"# Performing KBIN::XVASP_Afix (" << mode << ") [AFLOW] begin\" >> INCAR " << endl;
      aus_exec << "echo \"ISMEAR=2                                        # Performing RELAX_STATIC (Methfessel-Paxton order 2)\" >> INCAR " << endl;
      if(vflags.KBIN_VASP_INCAR_VERBOSE) aus_exec << "echo \"# Performing KBIN::XVASP_Afix (" << mode << ") [AFLOW] end\" >> INCAR " << endl;
      aurostd::execute(aus_exec);
    }
    // rewrite to restart ---------------------------------
    if(rewrite_incar) {aurostd::stringstream2file(xvasp.INCAR,string(xvasp.Directory+"/INCAR"));}
    if(rewrite_poscar) {aurostd::stringstream2file(xvasp.POSCAR,string(xvasp.Directory+"/POSCAR"));}
    if(rewrite_kpoints) {aurostd::stringstream2file(xvasp.KPOINTS,string(xvasp.Directory+"/KPOINTS"));}
    // reload to restart ---------------------------------
    if(reload_incar) {
      xvasp.INCAR_orig.str(std::string()); xvasp.INCAR_orig << xvasp.INCAR.str();
      xvasp.INCAR.str(std::string()); xvasp.INCAR << aurostd::file2string(xvasp.Directory+"/INCAR");
    }
    if(reload_poscar) {
      xvasp.POSCAR_orig.str(std::string()); xvasp.POSCAR_orig << xvasp.POSCAR.str();
      xvasp.POSCAR.str(std::string()); xvasp.POSCAR << aurostd::file2string(xvasp.Directory+"/POSCAR");
    }
    if(reload_kpoints) {
      xvasp.KPOINTS_orig.str(std::string()); xvasp.KPOINTS_orig << xvasp.KPOINTS.str();
      xvasp.KPOINTS.str(std::string()); xvasp.KPOINTS << aurostd::file2string(xvasp.Directory+"/KPOINTS");
    }

    // clean to restart ----------------------------------
    if(file_error!="") KBIN::XVASP_Afix_Clean(xvasp,file_error);
    if(file_error.empty()) {
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"cannot log error",_FILE_CORRUPT_); //CO20200624
    }
    // return
    if(param_int==0) {;} // dummy load
    if(mode=="PSMAXN") return out;
    return 0.0;
  }
}
// ---------------------------------------------------------------------------------------------------------------------------------------------------------

// ***************************************************************************//
// KBIN::GetMostRelaxedStructure
// ***************************************************************************//
namespace KBIN {
  xstructure GetMostRelaxedStructure(string directory) {
    string soliloquy="KBIN::GetMostRelaxedStructure():";  //CO20200404
    string POSCARfile;

    if(XHOST.vext.size()!=XHOST.vcat.size()) {
      throw aurostd::xerror(_AFLOW_FILE_NAME_, soliloquy, "XHOST.vext.size()!=XHOST.vcat.size(), aborting", _RUNTIME_ERROR_); //CO20200404
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    //READ POSCAR.orig
    string file_poscar_tmp=aurostd::TmpFileCreate("POSCAR.tmp");

    //GET its lattice vectors and reciprocal lattice vectors
    bool found_POSCAR=FALSE;
    if(!found_POSCAR && aurostd::FileExist(directory+"/POSCAR.bands")) {
      found_POSCAR=TRUE;
      POSCARfile=directory+"/POSCAR.bands";
    }
    if(!found_POSCAR && aurostd::FileExist(directory+"/POSCAR.static")) {
      found_POSCAR=TRUE;
      POSCARfile=directory+"/POSCAR.static";
    }
    if(!found_POSCAR && aurostd::FileExist(directory+"/POSCAR.relax2")) {
      found_POSCAR=TRUE;
      POSCARfile=directory+"/POSCAR.relax2";
    }
    if(!found_POSCAR && aurostd::FileExist(directory+"/POSCAR.relax1")) {
      found_POSCAR=TRUE;
      POSCARfile=directory+"/POSCAR.relax1";
    }
    if(!found_POSCAR && aurostd::FileExist(directory+"/POSCAR")) { 
      found_POSCAR=TRUE;
      POSCARfile=directory+"/POSCAR";
    }
    for(uint iext=0;iext<XHOST.vext.size();iext++) {
      if(!found_POSCAR && aurostd::FileExist(directory+"/POSCAR.bands"+XHOST.vext.at(iext))) {
        found_POSCAR=TRUE;
        aurostd::execute(XHOST.vcat.at(iext)+" "+directory+"/POSCAR.bands"+XHOST.vext.at(iext)+" > "+file_poscar_tmp);
        POSCARfile=file_poscar_tmp;
      }
    }
    for(uint iext=0;iext<XHOST.vext.size();iext++) {
      if(!found_POSCAR && aurostd::FileExist(directory+"/POSCAR.static"+XHOST.vext.at(iext))) {
        found_POSCAR=TRUE;
        aurostd::execute(XHOST.vcat.at(iext)+" "+directory+"/POSCAR.static"+XHOST.vext.at(iext)+" > "+file_poscar_tmp);
        POSCARfile=file_poscar_tmp;
      }
    }
    for(uint iext=0;iext<XHOST.vext.size();iext++) {
      if(!found_POSCAR && aurostd::FileExist(directory+"/POSCAR.relax1"+XHOST.vext.at(iext))) {
        found_POSCAR=TRUE;
        aurostd::execute(XHOST.vcat.at(iext)+" "+directory+"/POSCAR.relax1"+XHOST.vext.at(iext)+" > "+file_poscar_tmp);
        POSCARfile=file_poscar_tmp;
      }
    }
    for(uint iext=0;iext<XHOST.vext.size();iext++) {
      if(!found_POSCAR && aurostd::FileExist(directory+"/POSCAR.relax2"+XHOST.vext.at(iext))) {
        found_POSCAR=TRUE;
        aurostd::execute(XHOST.vcat.at(iext)+" "+directory+"/POSCAR.relax2"+XHOST.vext.at(iext)+" > "+file_poscar_tmp);
        POSCARfile=file_poscar_tmp;
      }
    }
    if(!found_POSCAR)  {
      throw aurostd::xerror(_AFLOW_FILE_NAME_, soliloquy, "No POSCAR[.bands|.static|.relax2|.relax1][.EXT] found in the directory, aborting", _INPUT_MISSING_); //CO20200404
    }

    xstructure xstr_name(POSCARfile, IOVASP_POSCAR); //import xstructure from filename
    if(aurostd::FileExist(file_poscar_tmp)) aurostd::RemoveFile(file_poscar_tmp);
    return xstr_name;
  }
}

// ***************************************************************************//
// KBIN::ExtractAtomicSpecies
// ***************************************************************************//
namespace KBIN {
  vector<string> ExtractAtomicSpecies(const string& directory,ostream& oss) {ofstream FileMESSAGE;return ExtractAtomicSpecies(directory,FileMESSAGE,oss);} //CO20200404 - added ofstream
  vector<string> ExtractAtomicSpecies(const string& directory,ofstream& FileMESSAGE,ostream& oss) { //CO20200404 - added ofstream
    string soliloquy="KBIN::ExtractAtomicSpecies():"; //CO20200404
    stringstream message;
    string OUTCARfile; //, POSCARfile;
    string file_outcar_tmp=aurostd::TmpFileCreate("OUTCAR.tmp");

    if(XHOST.vext.size()!=XHOST.vcat.size()) {
      throw aurostd::xerror(_AFLOW_FILE_NAME_, soliloquy, "XHOST.vext.size()!=XHOST.vcat.size(), aborting", _RUNTIME_ERROR_); //CO20200404
    }

    xstructure xstr_name=GetMostRelaxedStructure(directory); //CO20180626

    //WARNING!!!!
    if(aurostd::FileExist(directory+"/POSCAR.orig")) {
      xstructure xstr_poscar_orig(directory+"/POSCAR.orig", IOVASP_POSCAR);
      if(xstr_poscar_orig.atoms.size()!=xstr_name.atoms.size()) {
        //[CO20200404 - OBSOLETE]cerr << "!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!" << endl;
        //[CO20200404 - OBSOLETE]cerr << "I tell you--'POSCAR.orig' has different atoms number from 'POSCAR.bands', though I can still produce PEDOS plots for you!" << endl;
        pflow::logger(_AFLOW_FILE_NAME_, soliloquy, "POSCAR.orig has different atoms number from POSCAR.bands", directory, FileMESSAGE, oss, _LOGGER_WARNING_);
      }
    }

    //Extract the atom label (name) from POSCAR.relax1
    vector<string> AtomicName;
    for(unsigned int i=0; i<xstr_name.species.size(); i++) {
      if(xstr_name.species.at(i).compare("")!=0) {
        for (int j=0; j<xstr_name.num_each_type.at(i); j++) {
          xstr_name.species.at(i)=KBIN::VASP_PseudoPotential_CleanName(xstr_name.species.at(i));
          AtomicName.push_back(xstr_name.species.at(i));
        }
      }
      else { // cout << "no name given taking from OUTCAR" << endl;
        // ********************************************************************************  
        //READ  OUTCAR
        bool found_OUTCAR=FALSE;
        if(aurostd::FileExist(directory+"./OUTCARfile.tmp")) aurostd::RemoveFile(directory+"./OUTCARfile.tmp");
        if(!found_OUTCAR && aurostd::FileExist(directory+"/OUTCAR")) {
          found_OUTCAR=TRUE;
          OUTCARfile=directory+"/OUTCAR";
        }
        if(!found_OUTCAR && aurostd::FileExist(directory+"/OUTCAR.static")) {
          found_OUTCAR=TRUE;
          OUTCARfile=directory+"/OUTCAR.static";
        }
        for(uint iext=0;iext<XHOST.vext.size();iext++) {
          if(!found_OUTCAR && aurostd::FileExist(directory+"/OUTCAR.static"+XHOST.vext.at(iext))) {
            found_OUTCAR=TRUE;	  
            aurostd::execute(XHOST.vcat.at(iext)+" "+directory+"/OUTCAR.static"+XHOST.vext.at(iext)+" > "+file_outcar_tmp);
            OUTCARfile=file_outcar_tmp;
          }
        }
        if(!found_OUTCAR) {
          throw aurostd::xerror(_AFLOW_FILE_NAME_, soliloquy, "No OUTCAR[.static][.EXT] found in the directory, aborting", _INPUT_MISSING_);  //CO20200404
        }

        vector<string> vposcars, vpp;
        vector<string> AtomSpeciesName;
        string str_grep=aurostd::execute2string("grep TITEL "+OUTCARfile);
        aurostd::string2tokens(str_grep, vposcars, "\n");

        for(unsigned int k=0;k<vposcars.size();k++) {
          aurostd::string2tokens(vposcars.at(k),vpp," ");
          AtomSpeciesName.push_back(vpp.at(3));
          AtomSpeciesName.at(k)=KBIN::VASP_PseudoPotential_CleanName(AtomSpeciesName.at(k));
        }
        // ********************************************************************************  
        xstr_name.species.at(i)=AtomSpeciesName.at(i);
        for (int j=0; j<xstr_name.num_each_type.at(i); j++) {
          AtomicName.push_back(xstr_name.species.at(i));
        }
      }
    }

    //Double Check
    for (unsigned int i=0; i<AtomicName.size(); i++) {
      AtomicName.at(i)=KBIN::VASP_PseudoPotential_CleanName(AtomicName.at(i));
    }

    if(aurostd::FileExist(file_outcar_tmp)) aurostd::RemoveFile(file_outcar_tmp);

    /****************************************************************************************************
      for (unsigned int i=0; i<AtomicName.size(); i++) {
      cout <<AtomicName.at(i) << endl;
      }
      */
    return AtomicName;
  }
} // namespace KBIN

// ***************************************************************************//
// KBIN::ExtractEfermiOUTCAR
// ***************************************************************************//
namespace KBIN {
  double ExtractEfermiOUTCAR(string directory) {
    string soliloquy="KBIN::ExtractEfermiOUTCAR():";  //CO20200404
    double Efermi=0;
    string OUTCARfile, stmp, line;
    stringstream  strline;

    if(XHOST.vext.size()!=XHOST.vcat.size()) {
      throw aurostd::xerror(_AFLOW_FILE_NAME_, soliloquy, "XHOST.vext.size()!=XHOST.vcat.size(), aborting", _RUNTIME_ERROR_); //CO20200404
    }

    bool found_OUTCAR=FALSE;

    //READ OUTCAR, get the number of IONS, Fermi Level	
    string file_outcar_tmp=aurostd::TmpFileCreate("OUTCAR.tmp");
    if(aurostd::FileExist(directory+"./OUTCARfile.tmp")) aurostd::RemoveFile(directory+"./OUTCARfile.tmp");
    if(!found_OUTCAR && aurostd::FileExist(directory+"/OUTCAR")) {
      found_OUTCAR=TRUE;
      OUTCARfile=directory+"/OUTCAR";
    }
    if(!found_OUTCAR && aurostd::FileExist(directory+"/OUTCAR.static")) {
      found_OUTCAR=TRUE;
      OUTCARfile=directory+"/OUTCAR.static";
    }
    for(uint iext=0;iext<XHOST.vext.size();iext++) {
      if(!found_OUTCAR && aurostd::FileExist(directory+"/OUTCAR.static"+XHOST.vext.at(iext))) {
        found_OUTCAR=TRUE;
        aurostd::execute(XHOST.vcat.at(iext)+" "+directory+"/OUTCAR.static"+XHOST.vext.at(iext)+" > "+file_outcar_tmp);
        OUTCARfile=file_outcar_tmp;
      }
    }
    if(!found_OUTCAR) {
      throw aurostd::xerror(_AFLOW_FILE_NAME_, soliloquy, "No OUTCAR[.static][.EXT] found in the directory, aborting", _INPUT_MISSING_);  //CO20200404
    }

    ////GET Number of IONS and Fermi
    string anchor_word_Efermi="E-fermi";
    vector<string> vlines;
    aurostd::file2vectorstring(OUTCARfile, vlines);
    for(unsigned int i=0; i<vlines.size(); i++) {
      if(vlines.at(i).find(anchor_word_Efermi) !=string::npos) {
        strline.clear();
        strline.str(vlines.at(i));
        strline >> stmp >> stmp >> Efermi;
      }
    }
    if(aurostd::FileExist(file_outcar_tmp)) aurostd::RemoveFile(file_outcar_tmp);
    //cout << Efermi << "fermi" << endl;
    return Efermi;
  }
} // namespace KBIN

// ***************************************************************************
// KBIN::ExtractSystemName(string directory)
// ***************************************************************************
namespace KBIN {
  //ME20200217 - the old ExtractSystemName function (now ExtractSystemNameVASP)
  // breaks for long system names since VASP cuts them off. The system in the
  // aflow.in file should be the canonical system name since all VASP output
  // ultimately depends on it.
  string ExtractSystemName(const string& directory) {
    string system_name = ExtractSystemNameFromAFLOWIN(directory);
    if (system_name.empty()) {
      string function = "KBIN::ExtractSystemName()";
      string message = "Could not extract system.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_CORRUPT_);
    }
    aurostd::StringSubst(system_name,":.","."); //CO20200731 - patching bug in system name generation
    aurostd::StringSubst(system_name,"@ARUN.A",":ARUN.A");  //CO20200731 - patching for AEL/AGL ARUNs - @ in filenames are NOT good
    return system_name;
  }

  string ExtractSystemNameFromAFLOWIN(const string& _directory) {
    //CO20200731 - adjusting for old aflow.ins that only contain 'SYSTEM=' in INCAR section, or
    //[AFLOW] SYSTEM = Mo_pvRu_pv.A30B0
    string directory=_directory;  //CO20200624
    if (directory.empty()) directory = ".";
    string system_name = "";
    string aflowin_path = directory + "/" + _AFLOWIN_;
    if(!aurostd::FileExist(aflowin_path)) {
      string function = "KBIN::ExtractSystemNameFromAFLOWIN():";
      string message = "Could not find file " + aflowin_path + ".";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_NOT_FOUND_);
    }
    // ME20200525 - Need to take potential white space between [AFLOW],
    // SYSTEM, and = into account.
    string AflowIn=aurostd::RemoveComments(aurostd::file2string(aflowin_path)); // for safety //CO20180502
    vector<string> vAflowIn;aurostd::string2vectorstring(AflowIn,vAflowIn);
    uint nlines = vAflowIn.size();
    for (uint iline = 0; iline < nlines && system_name.empty(); iline++) {
      if (vAflowIn[iline].find("SYSTEM") != string::npos) {
        const string& line = vAflowIn[iline];
        string::size_type t = line.find_first_of("=");
        system_name = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(line.substr(t + 1));
        if(!system_name.empty()){return system_name;}
      }
    }

    return system_name;
  }

  string ExtractSystemNameFromVASP(const string& _directory) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy=XPID+"KBIN::ExtractSystemNameFromVASP():";
    string directory="", SystemName="", stmp="", DOSCARfile="";
    stringstream strline;

    directory=_directory; //Get directory
    if(directory=="") directory="."; // here

    if(LDEBUG) {cerr << XPID << "KBIN::ExtractSystemName: directory=" << directory << endl;}

    if(XHOST.vext.size()!=XHOST.vcat.size()) {
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"XHOST.vext.size()!=XHOST.vcat.size(), aborting",_RUNTIME_ERROR_); //CO20200404
    }

    bool found=FALSE;

    //READ DOSCAR.static; Firstly, we try to get SystemName from DOSCAR
    if(LDEBUG) {cerr << XPID << "KBIN::ExtractSystemName: [1]" << endl;}
    string doscarfile_tmp=aurostd::TmpFileCreate("DOSCAR.tmp");
    if(aurostd::FileExist(directory+"/DOSCARfile.tmp")) aurostd::RemoveFile(directory+"/DOSCARfile.tmp");
    if(!found && aurostd::FileExist(directory+"/DOSCAR")) {
      found=TRUE;
      DOSCARfile=directory+"/DOSCAR";
    }
    if(!found && aurostd::FileExist(directory+"/DOSCAR.static")) {
      found=TRUE;
      DOSCARfile=directory+"/DOSCAR.static";
    }
    for(uint iext=0;iext<XHOST.vext.size();iext++) {
      if(!found && aurostd::FileExist(directory+"/DOSCAR.static"+XHOST.vext.at(iext))) {
        found=TRUE;
        aurostd::execute(XHOST.vcat.at(iext)+" \""+directory+"/DOSCAR.static"+XHOST.vext.at(iext)+"\""+" > "+doscarfile_tmp);
        DOSCARfile=doscarfile_tmp;found=TRUE;
      } 
    }    
    if(found) {
      //OPEN DOSCAR file  
      //
      vector<string> vlines;
      aurostd::file2vectorstring(DOSCARfile, vlines);
      strline.clear();
      strline.str(vlines.at(4));
      strline >> SystemName;
      //cout <<SystemName << " Is this a systemname?" << endl;
      if(aurostd::FileExist(doscarfile_tmp)) aurostd::RemoveFile(doscarfile_tmp);

      //Get the working directory
      //[CO20191112 - OBSOLETE]const int PATH_LENGTH_MAX=1024;
      //[CO20191112 - OBSOLETE]char work_dir[PATH_LENGTH_MAX];
      //[CO20191112 - OBSOLETE]getcwd(work_dir, PATH_LENGTH_MAX);  
      string work_dir=aurostd::getPWD();  //CO20191112

      //Get the data directory
      chdir(directory.c_str());  //Jump into the data directory
      //[CO20191112 - OBSOLETE]char data_dir[PATH_LENGTH_MAX];
      //[CO20191112 - OBSOLETE]getcwd(data_dir, PATH_LENGTH_MAX);  
      string data_dir=aurostd::getPWD();  //CO20191112
      chdir(work_dir.c_str()); //Jump back into the work directory  //CO20191112

      //Check whether it is MAGNETIC DATABASE
      bool FLAG_MAGNETIC=FALSE;
      string abs_dir(data_dir);  //convert the char type of data_dir into string type
      if(abs_dir.find("LIB3")!=string::npos) FLAG_MAGNETIC=TRUE;

      if(FLAG_MAGNETIC) {
        return SystemName;
      } else {
        vector<string> data;
        aurostd::string2tokens(data_dir.c_str(), data, "/");  //CO20191112
        int datasize=data.size();
        bool FLAG_DIRECTORY_ICSD=FALSE;
        for (uint i=0; i<data.size();i++) {
          if((data.at(i).find("ICSD")!=string::npos)) {
            FLAG_DIRECTORY_ICSD =TRUE;
          }
        }

        //Check whether the last entry is "BANDS"
        string Last_Entry=data.back();
        bool FLAG_LAST_ENTRY_BANDS=FALSE;
        if(Last_Entry.compare("BANDS")==0) FLAG_LAST_ENTRY_BANDS=TRUE;

        //Check whether the SystemName is "unknown", to fix the bug of aflow
        bool FLAG_SYSTEMNAME_UNKNOWN=FALSE;
        if(SystemName.compare("unknown")==0) FLAG_SYSTEMNAME_UNKNOWN=TRUE;
        //if ICSD is not found in system name
        if((SystemName.find("ICSD")==string::npos)) {
          if(datasize >=2) {
            if(FLAG_LAST_ENTRY_BANDS) {
              if(FLAG_SYSTEMNAME_UNKNOWN && FLAG_DIRECTORY_ICSD) {
                SystemName=data.at(datasize-2);
              } else SystemName=data.at(datasize-3)+"."+data.at(datasize-2);
            } else {
              if(FLAG_SYSTEMNAME_UNKNOWN && FLAG_DIRECTORY_ICSD) {
                SystemName=Last_Entry;
              } else SystemName=data.at(datasize-2)+"."+data.back();
            }
          } else {
            SystemName=data.back();
          }
        }
        return SystemName;
      }
    }

    // TRY OUTCAR.relax2
    for(uint iext=0;iext<XHOST.vext.size();iext++) {
      if(!found && aurostd::FileExist(directory+"/OUTCAR.relax2"+XHOST.vext.at(iext))) {
        found=TRUE;
        xOUTCAR outcar(directory+"/OUTCAR.relax2"+XHOST.vext.at(iext));
        SystemName=outcar.SYSTEM;
        return SystemName;
      }
    }
    // TRY OUTCAR.relax1
    for(uint iext=0;iext<XHOST.vext.size();iext++) {
      if(!found && aurostd::FileExist(directory+"/OUTCAR.relax1"+XHOST.vext.at(iext))) {
        found=TRUE;
        xOUTCAR outcar(directory+"/OUTCAR.relax1"+XHOST.vext.at(iext));
        SystemName=outcar.SYSTEM;
        return SystemName;
      }
    }
    // NOTHING FOUND !!!
    if(!found) {
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"No DOSCAR[.static][.EXT], OUTCAR[.relax1|.relax2]][.EXT] found in the directory, aborting",_FILE_CORRUPT_); //CO20200624
    }
    return "";
  }
} // namespace KBIN

#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2021           *
// *                                                                         *
// ***************************************************************************

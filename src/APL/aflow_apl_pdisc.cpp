

#include "aflow_apl.h"

// Some parts are written within the C++0x support in GCC, especially the std::thread,
// which is implemented in gcc 4.4 and higher.... For multithreads with std::thread see:
// http://www.justsoftwaresolutions.co.uk/threading/multithreading-in-c++0x-part-1-starting-threads.html
#if GCC_VERSION >= 40400  // added two zeros
#define AFLOW_APL_MULTITHREADS_ENABLE 1
#include <thread>
#else
#warning "The multithread parts of APL will be not included, since they need gcc 4.4 and higher (C++0x support)."
#endif

using namespace std;

namespace apl {

// ///////////////////////////////////////////////////////////////////////////

PhononDispersionCalculator::PhononDispersionCalculator(IPhononCalculator& pc, Logger& l) : _pc(pc), _logger(l) {
  _system = _pc.getSystemName();  // ME190614
}

// ///////////////////////////////////////////////////////////////////////////

PhononDispersionCalculator::~PhononDispersionCalculator() {
  clear();
}

// ///////////////////////////////////////////////////////////////////////////

void PhononDispersionCalculator::clear() {
  _qpoints.clear();
  _freqs.clear();
  _temperature = 0.0;  // ME190614
}

//////////////////////////////////////////////////////////////////////////////

void PhononDispersionCalculator::initPathCoords(  //CO 180406
    const string& USER_DC_INITCOORDS,
    const string& USER_DC_INITLABELS,
    int USER_DC_NPOINTS, 
    bool CARTESIAN_COORDS) {
  if(USER_DC_INITCOORDS.empty() || USER_DC_INITLABELS.empty()){throw APLRuntimeError("apl::PhononDispersionCalculator::initPathCoords; Inputs are empty.");}
  _pb.defineCustomPoints(USER_DC_INITCOORDS,USER_DC_INITLABELS,_pc.getSupercell(),CARTESIAN_COORDS);
  _pb.setDensity(USER_DC_NPOINTS);
  //_qpoints = _pb.getPath(); // Get points // OBSOLETE ME190429 - this function should just define points; there is no path to set or get
}

void PhononDispersionCalculator::initPathLattice(const string& USER_DC_INITLATTICE,int USER_DC_NPOINTS){
  string lattice = USER_DC_INITLATTICE;
  if (lattice.empty()) {
    xstructure a(_pc.getInputCellStructure());
    //CO - START
    if (a.bravais_lattice_variation_type == "") {
      if (a.spacegroup == "") {
        if(a.space_group_ITC<1 || a.space_group_ITC>230){a.space_group_ITC = a.SpaceGroup_ITC();} //if (a.space_group_ITC == 0) { //CO 180214 - if not set then it could be 32767
        a.spacegroup = GetSpaceGroupName(a.space_group_ITC) + " #" + aurostd::utype2string(a.space_group_ITC);  //will break here if spacegroup is bad
      }
      // Use PLATON to get space group number if user did not get it...
      //a.platon2sg(_PLATON_P_EQUAL_DEFAULT,    //corey
      //	      _PLATON_P_EXACT_DEFAULT,
      //	      _PLATON_P_ANG_DEFAULT,
      //	      _PLATON_P_D1_DEFAULT,
      //	      _PLATON_P_D2_DEFAULT,
      //            _PLATON_P_D3_DEFAULT);
      //if( a.spacegroup.empty() )
      //throw apl::APLRuntimeError("apl::PhononDispersionCalculator::initPath(); The PLATON call to get spacegroup number failed. You have to specify it by DCINITSG in "+_AFLOWIN_);

      //vector<string> tokens;
      //aurostd::string2tokens(a.spacegroup,tokens,"#");
      //int spacegroupNumber = aurostd::string2utype<int>(tokens[1]);
      //tokens.clear();

      //lattice = LATTICE_Lattice_Variation_SpaceGroup(spacegroupNumber,_pc.getInputCellStructure());
      //lattice = LATTICE::SpaceGroup2LatticeVariation(spacegroupNumber,_pc.getInputCellStructure());
      lattice = LATTICE::SpaceGroup2LatticeVariation(a.space_group_ITC, a);
    } else {lattice = a.bravais_lattice_variation_type;}
    _logger << "The phonon dispersion curves will be generated for lattice variation " << lattice << "." << apl::endl;
  }
  //CO - END

  // cerr << "LATTICE=" << lattice << std::endl;
  // Suck point definition from the electronic structure part of AFLOW...
  _pb.takeAflowElectronicPath(lattice,_pc.getSupercell());             //CO 180406
                              //_pc.getInputCellStructure(),        //CO 180406
                              //_pc.getSuperCellStructure());       //CO 180406
  
  _pb.setDensity(USER_DC_NPOINTS);
  _qpoints = _pb.getPath(); // Get points
}

//////////////////////////////////////////////////////////////////////////////

void PhononDispersionCalculator::setPath(const string& USER_DC_OWNPATH) {
  // Get user's path...
  if (!USER_DC_OWNPATH.empty()) {
    if (USER_DC_OWNPATH.find('|') != string::npos) {
      // ME190614 - START
      // This breaks "mixed" paths such as G-X-W-L|K-U (interprets as G-X|W-L|K-U
      // _qpoints = _pb.getPath(apl::PathBuilder::COUPLE_POINT_MODE, USER_DC_OWNPATH);
      vector<string> tokens, tokens_pt;
      aurostd::string2tokens(USER_DC_OWNPATH, tokens, "-");
      string path;
      if (tokens[0].find('|') != string::npos) {
        string function = "PhononDispersionCalculator::setPath()";
        string message = "Cannot have | in the first path coordinate";
        throw aurostd::xerror(function, message, _INPUT_ILLEGAL_);
      } else {
        path = tokens[0];
      }
      for (uint i = 1; i < tokens.size(); i++) {
        if ((tokens[i].find('|') != string::npos) || (i == tokens.size() - 1)) {
          path += '-' + tokens[i];
        } else {
          path += '-' + tokens[i] + '|' + tokens[i];
        }
      }
      _qpoints = _pb.getPath(apl::PathBuilder::COUPLE_POINT_MODE, path);
      // ME190614 - END
    } else {
      _qpoints = _pb.getPath(apl::PathBuilder::SINGLE_POINT_MODE, USER_DC_OWNPATH);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void PhononDispersionCalculator::calculateInOneThread(int startIndex, int endIndex) {
  //cout << "Thread: from " << startIndex << " to " <<  endIndex << std::endl;
  for (int iqp = startIndex; iqp < endIndex; iqp++) {
    _logger.updateProgressBar(iqp, _qpoints.size());
    _freqs[iqp] = _pc.getFrequency(_qpoints[iqp], _frequencyFormat);
    //std::this_thread::yield();
  }
}

//////////////////////////////////////////////////////////////////////////////

void PhononDispersionCalculator::calc(const IPCFreqFlags frequencyFormat) {
  // Save
  _frequencyFormat = frequencyFormat;

  // Maybe there was some error and the list of q-points is empty, hence bye-bye...
  if (_qpoints.empty())
    throw apl::APLRuntimeError("There are no points for calculation.");

// Compute frequencies for each q-point

#ifdef AFLOW_APL_MULTITHREADS_ENABLE

  // Get the number of CPUS
  int ncpus; //= sysconf(_SC_NPROCESSORS_ONLN);  // AFLOW_MachineNCPUs;  //CO 180214
  _pc.get_NCPUS(ncpus);  //CO 180214
  if (ncpus < 1) ncpus = 1;
//  int qpointsPerCPU = _qpoints.size() / ncpus;  OBSOLETE ME180801

  // Show info
  if (ncpus == 1)
    _logger.initProgressBar("Calculating frequencies for PDIS");
  else
    _logger.initProgressBar("Calculating frequencies for PDIS (" + stringify(ncpus) + " threads)");

  // Prepare storage
  _freqs.clear();
  xvector<double> zero(_pc.getNumberOfBranches());
  for (uint i = 0; i < _qpoints.size(); i++)
    _freqs.push_back(zero);

  // Distribute the calculation
  int startIndex, endIndex;
  std::vector<std::thread*> threads;
  vector<vector<int> > thread_dist = getThreadDistribution((int) _qpoints.size(), ncpus);
  for (int icpu = 0; icpu < ncpus; icpu++) {
    startIndex = thread_dist[icpu][0];
    endIndex = thread_dist[icpu][1];
    threads.push_back(new std::thread(&PhononDispersionCalculator::calculateInOneThread, this, startIndex, endIndex));
  }

/* OBSOLETE ME 180801
  for (int icpu = 0; icpu < ncpus; icpu++) {
    startIndex = icpu * qpointsPerCPU;
    endIndex = startIndex + qpointsPerCPU;
    if (((uint)endIndex > _qpoints.size()) ||
        ((icpu == ncpus - 1) && ((uint)endIndex < _qpoints.size())))
      endIndex = _qpoints.size();
    threads.push_back(new std::thread(&PhononDispersionCalculator::calculateInOneThread, this, startIndex, endIndex));
  }
*/

  // Wait to finish all threads here!
  for (uint i = 0; i < threads.size(); i++) {
    threads[i]->join();
    delete threads[i];
  }
  threads.clear();

  // Done
  _logger.finishProgressBar();

#else

  _logger.initProgressBar("Calculating frequencies for PDIS");
  for (uint iqp = 0; iqp < _qpoints.size(); iqp++) {
    _logger.updateProgressBar(iqp, _qpoints.size());
    _freqs.push_back(_pc.getFrequency(_qpoints[iqp], _frequencyFormat));
  }
  _logger.finishProgressBar();

#endif
}

//////////////////////////////////////////////////////////////////////////////

void PhononDispersionCalculator::writePDIS() {
  string filename = DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_PDIS_FILE; //ME181226
  _logger << "Writing dispersion curves into file " << filename << "." << apl::endl; //ME181226

  //CO - START
  //ofstream outfile("PDIS",ios_base::out);
  stringstream outfile;
  //if( !outfile.is_open() ) {
  //  throw apl::APLRuntimeError("Cannot open output PDIS file.");
  //}
  //CO - END

  // Write header
  outfile << "# Phonon dispersion curves calculated by Aflow" << std::endl;
  outfile << "#" << std::endl;
  outfile << "# <system>    \"" << _system << "\"" << std::endl;  // ME190614 - use system name, not structure title
  outfile << "#" << std::endl;
  outfile << "# <units>     " << _frequencyFormat << std::endl;
  outfile << "# <nbranches> " << _pc.getNumberOfBranches() << std::endl;
  outfile << "# <npoints>   " << _freqs.size() << std::endl;
  outfile << "# <nsubpathp> " << _pb.getDensity() + 1 << std::endl;
  outfile << "#" << std::endl;

  // Write table of label points
  outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
  outfile << setprecision(8);
  double x = 0.0;
  double wholePathLength = _pb.getPathLength();
  map<double, string> labelMap;
  for (uint i = 1; i < _pb.getPointSize(); i++) {
    outfile << "# <label>     " << x << " "
            << setw(5) << _pb.getPointLabel(i)
            << std::endl;
    labelMap.insert(std::pair<double, string>(x, _pb.getPointLabel(i)));
    x += _pb.getPathLength(i) / wholePathLength;
  }
  outfile << "# <label>     " << 1.0 << " " << setw(5) << _pb.getPointLabel(_pb.getPointSize()) << std::endl;
  labelMap.insert(std::pair<double, string>(1.0, _pb.getPointLabel(_pb.getPointSize())));
  outfile << "#" << std::endl;

    //writing high-symmetry qpoints [PINKU] //PN180705
    stringstream ouths;  //PN180705
    ouths << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);  //PN180705
    ouths << setprecision(8);  //PN180705

  // Write table of exact _qpoints + use label map to identify its labels
  x = 0.0;
  int subpath = 0;
  double xstep = 0.0;
  int p = 0;
  vector<double> exactPointPositions;
  for (uint i = 0; i < _qpoints.size(); i++) {
    // Check it
    if (isExactQPoint(_qpoints[i], _pc.getSuperCellStructure().lattice)) {
      // Is it new exact points
      uint j = 0;
      for (; j < exactPointPositions.size(); j++)
        if (exactPointPositions[j] == x) break;

      // If yes, add it....
      if (j == exactPointPositions.size()) {
        exactPointPositions.push_back(x);
        string name = "-";
        std::map<double, string>::iterator iter = labelMap.begin();
        for (; iter != labelMap.end(); iter++)
          if (fabs(iter->first - x) < _AFLOW_APL_EPS_) break;
        if (iter != labelMap.end())
          name = iter->second;
        outfile << "# <exact>     " << x << " "
                << setw(5) << name
                << std::endl;
        ouths << "# <exact>     " << x << " "  //PN180705
                << setw(5) << name  //PN180705
                << setw(15) << _qpoints[i][1]<< setw(15) << _qpoints[i][2]<< setw(15) << _qpoints[i][3]  //PN180705
                << '\n';  //PN180705
      }
    }

    // Step of x will change in each subpath
    if (i % (_pb.getDensity() + 1) == 0) {
      if (i + 1 != _freqs.size())
        xstep = _pb.getPathLength(++subpath) / wholePathLength / (_pb.getDensity());
    }
    x += xstep;
    if (p != 0 && (p % _pb.getDensity()) == 0) {
      x -= xstep;
      p = -1;
    }
    p++;
  }
  outfile << "#" << std::endl;
  labelMap.clear();
  exactPointPositions.clear();

  // Write frequencies
  //[OBSOLETE PN180705]path_segment.clear();  //[PINKU]
  //[OBSOLETE PN180705]path.clear();          //[PINKU]
  outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
  outfile << setprecision(8);
  x = 0.0;
  subpath = 0;
  xstep = 0.0;
  p = 0;
  for (uint i = 0; i < _freqs.size(); i++) {
      ouths<<setw(15)<<_qpoints[i][1]<<setw(15) //PN180705
           <<_qpoints[i][2]<<setw(15)<<_qpoints[i][3]<<setw(15)<<p<<setw(15)<<x<<"\n"; //PN180705

    outfile << setw(4) << p << " ";
    //[OBSOLETE PN180705]path_segment.push_back(p);  //[PINKU]
    outfile << setw(15) << x << " ";
    //[OBSOLETE PN180705]path.push_back(x);  //[PINKU]
    for (uint j = 1; j <= _pc.getNumberOfBranches(); j++)
      outfile << setw(15) << _freqs[i](j) << " ";
    outfile << std::endl;

    // Step of x will change in each subpath
    if (i % (_pb.getDensity() + 1) == 0) {
      if (i + 1 != _freqs.size())
        xstep = _pb.getPathLength(++subpath) / wholePathLength / (_pb.getDensity());
    }
    x += xstep;
    if (p != 0 && (p % _pb.getDensity()) == 0) {
      x -= xstep;
      p = -1;
    }
    p++;
  }

  //CO - START
  aurostd::stringstream2file(outfile, filename); //ME181226
  if (!aurostd::FileExist(filename)) { //ME181226
    string function = "PhononDispersionCalculator::writePDIS()";
    string message = "Cannot open output file " + filename + "."; //ME181226
    throw aurostd::xerror(function, message, _FILE_ERROR_);
//    throw apl::APLRuntimeError("Cannot open output PDIS file.");
  }
  //
  //outfile.clear();
  //outfile.close();
  //CO - END

  //PINKU //PN180705
  string hskptsfile = DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_HSKPTS_FILE; //ME181226
  aurostd::stringstream2file(ouths, hskptsfile); //ME181226
  if (!aurostd::FileExist(hskptsfile)) { //ME181226
    string function = "PhononDispersionCalculator::writePDIS()";
    string message = "Cannot open output file " + hskptsfile + "."; //ME181226
    throw aurostd::xerror(function, message, _FILE_ERROR_);
//    throw apl::APLRuntimeError("Cannot open output aflow.apl_hskpoints.out file.");
  }
  //PINKU
}

//////////////////////////////////////////////////////////////////////////////

bool PhononDispersionCalculator::isExactQPoint(const xvector<double>& qpoint,
                                               const xmatrix<double>& lattice) {
  xcomplex<double> iONE(0.0, 1.0);

  bool isExact = false;
  for (_AFLOW_APL_REGISTER_ int i = 1; i <= 1; i++) {
    for (_AFLOW_APL_REGISTER_ int j = 1; j <= 1; j++) {
      for (_AFLOW_APL_REGISTER_ int k = 1; k <= 1; k++) {
        xvector<double> L = (((double)i) * lattice(1) +
                             ((double)j) * lattice(2) +
                             ((double)k) * lattice(3));
        xcomplex<double> p = exp(iONE * scalar_product(qpoint, L));
        if ((fabs(p.imag()) < _AFLOW_APL_EPS_) &&
            (fabs(p.real() - 1.0) < _AFLOW_APL_EPS_)) {
          isExact = true;
          break;
        }
      }
      if (isExact) break;
    }
    if (isExact) break;
  }

  return (isExact);
}

// ///////////////////////////////////////////////////////////////////////////

// ME190614 - START
// Write the eigenvalues into a VASP EIGENVAL-formatted file
void PhononDispersionCalculator::writePHEIGENVAL() {
  string filename = DEFAULT_APL_PHEIGENVAL_FILE;
  _logger << "Writing phonon eigenvalues file " << filename << "." << apl::endl;
  stringstream outfile;
  // Header
  int nions = (int) _pc.getInputCellStructure().atoms.size();
  outfile << std::setw(4) << nions
          << std::setw(4) << nions
          << std::setw(4) << 0
          << std::setw(4) << 1 << std::endl;
  double vol, a, b, c, potim;
  vol = GetVolume(_pc.getInputCellStructure())/nions;
  a = _pc.getInputCellStructure().a * 1E-10;
  b = _pc.getInputCellStructure().b * 1E-10;
  c = _pc.getInputCellStructure().c * 1E-10;
  potim = 0.5E-15;
  outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
  outfile << std::fixed;
  outfile << std::setprecision(7) << std::scientific;
  outfile << std::setw(15) << vol
          << std::setw(15) << a
          << std::setw(15) << b
          << std::setw(15) << c
          << std::setw(15) << potim << std::endl;
  
  outfile << std::setw(23) << std::setprecision(15) << _temperature << std::endl;
  outfile << "  PHON" << std::endl;
  outfile << " " << _pc.getSystemName() << std::endl;

  uint nbranch = _pc.getNumberOfBranches();
  uint nqpts = _freqs.size();
  outfile << std::dec << std::setw(4) << 0 << "  " << nqpts << std::setw(4) << nbranch << std::endl;
  
  double weight = 1.0/_freqs.size();
  double factorTHz2Raw = _pc.getFrequencyConversionFactor(apl::THZ, apl::RAW);
  double factorRaw2meV = _pc.getFrequencyConversionFactor(apl::RAW, apl::MEV);
  double conv = factorTHz2Raw * factorRaw2meV/1000;
  apl::PathBuilder::StoreEnumType store = _pb.getStore();
  xvector<double> qpt(3);
  xmatrix<double> c2f = inverse(trasp(ReciprocalLattice(_pc.getInputCellStructure().lattice)));

  for (uint q = 0; q < nqpts; q++) {
    if (store == apl::PathBuilder::CARTESIAN_LATTICE) qpt = c2f * _qpoints[q];
    else qpt = _qpoints[q];
    outfile << " " << std::endl;  // Space MUST be there or the xEIGENVAL reader breaks!
    for (int i = 1; i < 4; i++) {
      outfile << std::scientific << std::setprecision(8) << "  " << qpt[i];
    }
    outfile << std::scientific << std::setprecision(8) << "  " << weight << std::endl;
    for (uint br = 1; br <= nbranch; br++) {
      outfile << std::dec << std::setprecision(0) << std::setw(4) << br
              << std::setprecision(8) << std::fixed << std::setw(15) << (conv * _freqs[q][br]) << std::endl;
    }
  }
  
  // Write file
  aurostd::stringstream2file(outfile, filename);
  if (!aurostd::FileExist(filename)) {
    string function = "PhononDispersionCalculator::writePHEIGENVAL()";
    string message = "Cannot open output file " + filename + ".";
    throw aurostd::xerror(function, message, _FILE_ERROR_);
  }
  // Also write PHKPOINTS and PHPOSCAR file
  writePHPOSCAR();  // Structure required for plotting
  writePHKPOINTS();
}

// ///////////////////////////////////////////////////////////////////////////

// Write the structure into a VASP POSCAR-formatted file
void PhononDispersionCalculator::writePHPOSCAR() {
  string filename = DEFAULT_APL_PHPOSCAR_FILE;
  stringstream poscar;
  poscar << _pc.getInputCellStructure();
  aurostd::stringstream2file(poscar, filename);
  if (!aurostd::FileExist(filename)) {
    string function = "PhononDispersionCalculator::writePHPOSCAR()";
    string message = "Cannot open output file " + filename + ".";
    throw aurostd::xerror(function, message, _FILE_ERROR_);
  }
}

// ///////////////////////////////////////////////////////////////////////////

// Write the k-point path into a VASP KPOINTS-formatted file
void PhononDispersionCalculator::writePHKPOINTS() {
  _pb.writePHKPOINTS(_pc.getSupercell());
}

// ///////////////////////////////////////////////////////////////////////////
// ME190614 - END

}  // namespace apl

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *                                                                         *
// ***************************************************************************
// Written by Andriy Smolyanyuk (andriy.smolyanyuk@duke.edu)
//
//    This file is part of Aflow.
//
//    Aflow is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    Aflow is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with Aflow.  If not, see <https://www.gnu.org/licenses/>.

#include "aflow_apl.h"

#if __cplusplus >= 201103L
template <typename T>
using auto_ptr = std::unique_ptr<T>;
#else
using std::auto_ptr;
#endif

// Field width modifiers: used to format text during the output via std::setw function
#define SW 5  // width of columns with blank space separator
#define TW 15 // width of columns containing label/number
#define DCOEFF 1e-2 // this coefficient is used in numerical differentiation
// (should not be too small). Usage dT = DCOEFF*T

#define DEBUG_QHA false // toggles QHA-related debug output
#define DEBUG_QHA_GP_FIT false // toggles debug output related to the fit functionality
// in the calcGrueneisen function

//=============================================================================
//              Definitions of the NonlinearFit class members

namespace apl
{
  NonlinearFit::NonlinearFit()
  {
    free();
  }

  NonlinearFit::NonlinearFit(const NonlinearFit &nlf)
  {
    if (this==&nlf) return;
    free(); copy(nlf);
  }

  NonlinearFit::NonlinearFit(xvector<double> &x_in, xvector<double> &y_in,
      xvector<double> &guess_in,
      double foo(const double x, const xvector<double> &p, xvector<double> &dydp),
      double tol_in, double tau_in, int max_iter_in)
  {
    free();
    Npoints = x_in.rows;
    Nparams = guess_in.rows;
    tol = tol_in;
    tau = tau_in;
    max_iter = max_iter_in;
    x = x_in;
    y = y_in;
    residuals = xvector<double> (Npoints);
    guess = guess_in;
    p = xvector<double> (Nparams);
    dydp = xvector<double> (Nparams);
    A = xmatrix<double> (Nparams, Nparams);
    J = xmatrix<double> (Npoints, Nparams);
    f = foo;
  }

  NonlinearFit::~NonlinearFit() { free(); }

  const NonlinearFit& NonlinearFit::operator=(const NonlinearFit &nlf)
  {
    copy(nlf);
    return *this;
  }

  void NonlinearFit::clear() { free(); }
  void NonlinearFit::free()
  {
    Npoints = 0;
    Nparams = 0;
    tol = 0.0;
    tau = 0.0;
    max_iter = 0;
    x.clear();
    y.clear();
    residuals.clear();
    guess.clear();
    p.clear();
    dydp.clear();
    A.clear();
    J.clear();
    f = NULL;
  }

  void NonlinearFit::copy(const NonlinearFit &nlf)
  {
    if (this==&nlf) return;

    Npoints = nlf.Npoints;
    Nparams = nlf.Nparams;
    tol = nlf.tol;
    tau = nlf.tau;
    max_iter = nlf.max_iter;
    x = nlf.x;
    y = nlf.y;
    residuals = nlf.residuals;
    guess = nlf.guess;
    p = nlf.p;
    dydp = nlf.dydp;
    A = nlf.A;
    J = nlf.J;
    f = nlf.f;
  }

  /// Calculates the residual sum of squares of a model function w.r.t the given data
  ///
  double NonlinearFit::calculateResidualSquareSum(const xvector<double> &params)
  {
    double chi_sqr = 0.0;

    calculateResiduals(params); // calculate residuals and save them to residuals
    for (int i=1; i<= Npoints; i++) chi_sqr += pow(residuals[i], 2);

    return chi_sqr;
  }

  /// Calculates residuals of a given model function and stores the results
  ///
  void NonlinearFit::calculateResiduals(const xvector<double> &params)
  {
    for (int i=1; i<=Npoints; i++) residuals[i] = y[i] - f(x[i], params, dydp);
  }

  /// Calculates the Jacobian of a given model function for the given "guess" parameters
  ///
  void NonlinearFit::Jacobian(const xvector<double> &guess)
  {
    for (int i=1; i<=Npoints; i++){
      f(x[i], guess, dydp); 
      for (int j=1; j<=Nparams; j++){
        J[i][j] = -dydp[j]; // derivative of a given function w.r.t fit parameters
      }
    }
  }

  /// Nonlinear fit using the Levenberg-Marquardt algorithm.
  /// The implementation here is based on the ideas from Numerical Recipes and
  /// K.Madsen et al. Methods For Non-linear Least Squares Problems
  /// http://www2.imm.dtu.dk/pubdb/views/edoc_download.php/3215/pdf/imm3215.pdf
  ///
  bool NonlinearFit::fitLevenbergMarquardt()
  {
    string function = "NonlinearFit::fit(): ";
    bool LDEBUG = (FALSE || DEBUG_QHA || XHOST.DEBUG);
    if (LDEBUG) cerr << function << "begin"  << std::endl;

    static const double step_scaling_factor = 10.0;

    xvector<double> pnew(Nparams);
    xmatrix<double> G(Nparams, 1), M(Nparams, Nparams);

    int iter = 0;
    double chi_sqr = 0.0, new_chi_sqr =0.0; // old and new residual square sums

    p = guess; Jacobian(p); calculateResiduals(p); // p are the fitted parameters

    xmatrix<double> A = trasp(J)*J;
    xvector<double> g = -trasp(J)*residuals;

    // determine initial step
    double lambda = tau*maxDiagonalElement(A);
    if (LDEBUG) cerr << function << "lambda = " << lambda << std::endl;

    while (iter<max_iter){
      iter++;
      if (LDEBUG) cerr << function << "iteration: " << iter << std::endl;

      // M = A + lambda*diag(A)
      M = A; for (int i=1; i<=Nparams; i++) M[i][i] += lambda*A[i][i];

      if (LDEBUG) cerr << function << " M:" << std::endl << M << std::endl;

      // transformation: xvector(N) => xmatrix(N,1) to use GaussJordan function
      for (int i=1; i<=Nparams; i++) G[i][1] = g[i];
      aurostd::GaussJordan(M,G);

      pnew = p + M*g;

      chi_sqr = calculateResidualSquareSum(p);
      new_chi_sqr = calculateResidualSquareSum(pnew);

      if (abs(new_chi_sqr - chi_sqr) < tol){
        p = pnew;
        break;
      }

      // update only if step leads to a smaller residual sum of squares
      if (new_chi_sqr < chi_sqr){
        p = pnew; Jacobian(p); calculateResiduals(p);

        A = trasp(J)*J;
        g = -trasp(J)*residuals;

        lambda /= step_scaling_factor;
      }
      else{
        lambda *= step_scaling_factor;
      }
      if (LDEBUG) cerr << function << "pnew = " << p << std::endl;
    }

    if (LDEBUG) cerr << function << "end"  << std::endl;
    return iter < max_iter;
  }
}

//================================================================================
//                    EOS related

/// Murnaghan equation of state
/// https://en.wikipedia.org/wiki/Murnaghan_equation_of_state
/// Ref: Proceedings of the National Academy of Sciences of the United States of 
/// America, 30 (9): 244–247
/// 
double Murnaghan(const double x, const xvector<double> &p, xvector<double> &dydp)
{
  double Eeq = p[1];
  double Veq = p[2];
  double B   = p[3];
  double Bp  = p[4];
  double V = x;

  dydp[1] = 1.0;
  dydp[2] = (B*(Veq-V*pow(Veq/V,Bp)))/(Veq-Bp*Veq);
  dydp[3] = -(Veq/(Bp-1.0))+(V*(1+pow(Veq/V,Bp)/(Bp-1.0)))/Bp;
  dydp[4] = (B*(-V+2.0*Bp*V-pow(Bp,2)*V+pow(Bp,2)*Veq+V*pow(Veq/V,Bp)*(1.0-2.0*Bp
          +(Bp-1.0)*Bp*log(Veq/V))))/(pow(Bp-1.0,2)*pow(Bp,2));

  return Eeq-(B*Veq)/(Bp-1.0)+(B*V*(1+pow(Veq/V,Bp)/(Bp-1)))/Bp;
}

/// Birch-Murnaghan equation of state
/// https://en.wikipedia.org/wiki/Birch%E2%80%93Murnaghan_equation_of_state
/// Ref: Physical Review. 71 (11): 809–824
/// 
double BirchMurnaghan(const double x, const xvector<double> &p, xvector<double> &dydp)
{
  double Eeq = p[1];
  double Veq = p[2];
  double B   = p[3];
  double Bp  = p[4];
  double V = x;

  double kappa = pow(Veq/V, 2.0/3.0);

  dydp[1] = 1.0;
  dydp[2] = 3.0*B*kappa*(kappa-1)*(3.0*(Bp-6.0)*V*pow(Veq/V,1.0/3.0)+
      Veq*(62.0-36.0*kappa+3*Bp*(3*kappa-4)))/(16.0*Veq);
  dydp[3] = 9.0/16.0*Veq*pow(kappa-1,2.0)*(6.0-4.0*kappa+Bp*(kappa-1));
  dydp[4] = 9.0/16.0*B*Veq*pow(kappa-1,3.0);

  return Eeq + 9.0/16.0*B*Veq*(pow(kappa-1.0,3.0)*Bp + pow(kappa-1.0,2.0)*(6.0-4.0*kappa));
}


/// Checks if there is a minimum within a given data set 
/// (at least one internal point should be lower than the edge points)
bool isMinimumWithinBounds(const xvector<double> &y){
  for (int i=y.lrows+1; i<y.urows; i++){
    if (y[y.lrows] > y[i] && y[y.urows] > y[i]) return true;
  }

  return false;
}

/// The polynomial model that is used for equation of state fitting.
///
/// E = a + b*V^(-2/3) + c*V^(-4/3) + d*V^(-6/3) + e*V^(-8/3)
/// Ref: (eq (5.2) page 106)  https://doi.org/10.1017/CBO9781139018265
///
double EOSpoly(double V, const xvector<double> &p){
  if (p.rows != 5) return 0;
  return p[1] + p[2]*pow(V,-2/3.0) + p[3]*pow(V,-4/3.0)
    + p[4]*pow(V,-6/3.0) + p[5]*pow(V,-8/3.0);
}

/// First derivative w.r.t volume of the polynomial model
double dEOSpoly(double x, const xvector<double> &p){
  if (p.rows != 5) return 0;
  return -2/3.0*p[2]*pow(x,-5/3.0) - 4/3.0*p[3]*pow(x, -7/3.0)
    -6/3.0*p[4]*pow(x,-9/3.0) - 8/3.0*p[5]*pow(x,-11/3.0);
}

/// Second derivative w.r.t volume of the polynomial model
double d2EOSpoly(double x, const xvector<double> &p){
  if (p.rows != 5) return 0;
  return  10/9.0*p[2]*pow(x, -8/3.0) + 28/9.0*p[3]*pow(x,-10/3.0)
    +     6*p[4]*pow(x,-12/3.0) + 88/9.0*p[5]*pow(x,-14/3.0);
}

/// Third derivative w.r.t volume of the polynomial model
double d3EOSpoly(double x, const xvector<double> &p){
  if (p.rows != 5) return 0;
  return -(80.0/27.0*p[2]*pow(x, -11.0/3.0) + 280.0/27.0*p[3]*pow(x, -13.0/3.0) +
      24.0*p[4]*pow(x, -15.0/3.0) + 88*14/27.0*p[5]*pow(x, -17.0/3.0));
}

/// Calculates the equilibrium volume for the polynomial model.
///
/// Veq is determined as a solution of dE/dV=0 using the bisection method.
/// The function returns a negative number if an error has occurred.
///
/// It is assumed that only one minimum exists between Vmin and Vmax,
/// since it is the only physically meaningful situation.
///
/// The main idea: if a continuous function f(x) has the root f(x0)=0 in the interval [a,b], 
/// then its sign in the interval [a,x0) is opposite to its sign in the interval (x0,b].
/// The root searching algorithm is iterative: one bisects the interval [a,b] into two 
/// subintervals and for the next iteration picks a subinterval where f(x) has 
/// opposite sign on its ends.
/// This process continues until the length of the interval is less than some negligibly
/// small number and one picks the middle of the interval as x0.
/// In this situation, f(x0)~0 and this condition is used as a criterion to stop the
/// iteration procedure.
///
double EOSpolyGetEqVolume(const xvector<double> &p, double Vmin, double Vmax, 
    double tol=_mm_epsilon) {
  bool LDEBUG = (FALSE || DEBUG_QHA || XHOST.DEBUG);
  string function = "EOSpolyGetEqVolume(): ";
  if (LDEBUG) cerr << function << "begin" << std::endl;

  double left_end = Vmin; double right_end = Vmax;
  double middle = 0.5*(left_end + right_end);

  double f_at_left_end  = dEOSpoly(left_end, p);
  double f_at_right_end = dEOSpoly(right_end, p);
  double f_at_middle    = dEOSpoly(middle, p);

  // no root within a given interval
  if (sign(f_at_left_end) == sign(f_at_right_end)) return -1;

  if (LDEBUG){
    cerr << function << "left_end= "  << left_end  << "middle= ";
    cerr << middle  << "right_end= "  << right_end  << std::endl;
    cerr << function << "f_left= " << f_at_left_end;
    cerr << "f_middle= " << f_at_middle << "f_right= ";
    cerr << f_at_right_end << std::endl;
  }

  // Iterate until the convergence criterion is reached:
  // f(middle) is sufficiently close to zero.
  // Meanwhile, do a sanity check that the function has opposite signs at the interval ends.
  while ((sign(f_at_left_end) != sign(f_at_right_end)) && (abs(f_at_middle)>tol)){
    if (sign(f_at_left_end) == sign(f_at_middle)){
      std::swap(left_end, middle);
      std::swap(f_at_left_end, f_at_middle);
    }

    if (sign(f_at_right_end) == sign(f_at_middle)){
      std::swap(right_end, middle);
      std::swap(f_at_right_end, f_at_middle);
    }

    middle = 0.5*(left_end + right_end);
    f_at_middle = dEOSpoly(middle, p);

    if (LDEBUG){
      cerr << function << "left_end= "  << left_end  << "middle= ";
      cerr << middle  << "right_end= "  << right_end  << std::endl;
      cerr << function << "f_left= " << f_at_left_end;
      cerr << "f_middle= " << f_at_middle << "f_right= ";
      cerr << f_at_right_end << std::endl;
    }
  }

  // double-check that the convergence criterion was reached
  if (abs(f_at_middle) > tol) return -1;

  if (LDEBUG) cerr << function << "end" << std::endl;
  return middle;
}

/// Perform a fit to the polynomial EOS model for a given set of volumes and energies
/// using the linear least squares method.
///
xvector<double> fitEOSpoly(const xvector<double> &Volumes, xvector<double> &E){
  xmatrix<double> M(Volumes.rows,5);
  xvector<double> sigma(Volumes.rows);
  for (int i=sigma.lrows;i<=sigma.urows;i++) sigma[i]=1.0;

  for (int i=Volumes.lrows; i<=Volumes.urows; i++) {
    for (int j=0; j<=4; j++) M[i][j+1] = pow(Volumes[i],-2*j/3.0);
  }

  aurostd::cematrix M_ce(M);
  M_ce.LeastSquare(E, sigma);
  return M_ce.GetFitVector();
}

/// Calculates the bulk modulus for the polynomial EOS model.
double BulkModulus(double V0, const xvector<double> &fit_parameters){
  return eV2GPa*V0*d2EOSpoly(V0,fit_parameters);
}

/// Calculates the pressure derivative of bulk modulus for the polynomial EOS model.
double Bprime(double V0, const xvector<double> &fit_parameters){
  return -(V0*d3EOSpoly(V0, fit_parameters)/d2EOSpoly(V0, fit_parameters) + 1);
}

//=============================================================================
//                         Definitions of the QHA class members
namespace apl
{
  QHAN::QHAN(ostream& oss) : xStream(oss) { free(); }
  QHAN::QHAN(const QHAN &qha) : xStream(*qha.getOFStream(),*qha.getOSS()) {
    free(); copy(qha);
  }
  QHAN::~QHAN() { xStream::free(); free(); }
  const QHAN& QHAN::operator=(const QHAN &qha){
    copy(qha);
    return *this;
  }

  void QHAN::clear()
  {
    free();
  }

  /// Initializes all values to "zero" and attempts to clear all containers.
  ///
  void QHAN::free()
  {
    system_title = "";
    supercellopts.clear();
    isEOS = false; isGP_FD = false;
    ignore_imaginary = false;
    runQHA   = false; runQHA3P = false; runSCQHA = false;
    isInitialized = false;
    includeElectronicContribution = false;
    Ntemperatures = 0;
    N_GPvolumes = 3;
    N_EOSvolumes = 0;
    Nbranches = 0;
    NatomsOrigCell = 0;
    origStructure.clear();
    Temperatures.clear();
    GPvolumes.clear();
    EOSvolumes.clear();
    coefGPVolumes.clear();
    coefEOSVolumes.clear();
    DOS_Ef.clear();
    Efermi_V.clear();
    E0_V.clear();
    static_eigvals.clear();
    static_ibzkpts.clear();
    energies_V.clear();
    edos_V.clear();
    frequencies_V.clear();
    pdos_V.clear();
    qpWeights.clear();
    qPoints.clear();
    gp_fit_matrix.clear();
    omegaV_mesh.clear();
    omegaV_mesh_EOS.clear();
    gp_ph_dispersions.clear();
    eos_vib_thermal_properties.clear();
    subdirectories_apl_gp.clear();
    subdirectories_apl_eos.clear();
    subdirectories_static.clear();
    arun_runnames_static.clear();
    xinput.clear();
    currentDirectory = ".";
  }

  void QHAN::copy(const QHAN &qha)
  {
    if (this==&qha) return;

    apl_options       = qha.apl_options;
    system_title      = qha.system_title;
    supercellopts     = qha.supercellopts;
    isEOS             = qha.isEOS;
    isGP_FD           = qha.isGP_FD;
    ignore_imaginary  = qha.ignore_imaginary;
    runQHA            = qha.runQHA;
    runQHA3P          = qha.runQHA3P;
    runSCQHA          = qha.runSCQHA;
    isInitialized     = qha.isInitialized;
    includeElectronicContribution = qha.includeElectronicContribution;
    Ntemperatures     = qha.Ntemperatures;
    N_GPvolumes       = qha.N_GPvolumes;
    N_EOSvolumes      = qha.N_EOSvolumes;
    Nbranches         = qha.Nbranches;
    NatomsOrigCell    = qha.NatomsOrigCell;
    origStructure     = qha.origStructure;
    Temperatures      = qha.Temperatures;
    GPvolumes         = qha.GPvolumes;
    EOSvolumes        = qha.EOSvolumes;
    coefGPVolumes     = qha.coefGPVolumes;
    coefEOSVolumes    = qha.coefEOSVolumes;
    DOS_Ef            = qha.DOS_Ef;
    Efermi_V          = qha.Efermi_V;
    E0_V              = qha.E0_V;
    static_eigvals    = qha.static_eigvals;
    static_ibzkpts    = qha.static_ibzkpts;
    energies_V        = qha.energies_V;
    edos_V            = qha.edos_V;
    frequencies_V     = qha.frequencies_V;
    pdos_V            = qha.pdos_V;
    qpWeights         = qha.qpWeights;
    qPoints           = qha.qPoints;
    gp_fit_matrix     = qha.gp_fit_matrix;
    omegaV_mesh       = qha.omegaV_mesh;
    omegaV_mesh_EOS   = qha.omegaV_mesh_EOS;
    gp_ph_dispersions = qha.gp_ph_dispersions;
    eos_vib_thermal_properties = qha.eos_vib_thermal_properties;
    subdirectories_apl_eos  = qha.subdirectories_apl_eos;
    subdirectories_apl_gp   = qha.subdirectories_apl_gp;
    subdirectories_static   = qha.subdirectories_static;
    arun_runnames_static    = qha.arun_runnames_static;
    xinput            = qha.xinput;
    currentDirectory  = qha.currentDirectory;

    xStream::copy(qha);
  }

  ///////////////////////////////////////////////////////////////////////////////

  QHAN::QHAN(string &tpt, _xinput &xinput, _kflags &kflags, xoption &supercellopts,
      ofstream &messageFile, ostream &oss)
  {
    initialize(tpt, xinput, kflags, supercellopts, messageFile, oss);
  }

  /// Initializes the QHA class with all the necessary data.
  ///
  void QHAN::initialize(string &tpt, _xinput &xinput, _kflags &kflags,
      xoption &supercellopts, ofstream &messageFile, ostream &oss)
  {
    static const int REQUIRED_MIN_NUM_OF_DATA_POINTS_FOR_EOS_FIT = 5;
    static const int precision_format = 3;

    isInitialized = false;

    xStream::initialize(messageFile, oss);

    string function = "QHAN::initialize():";
    string msg  = "Initializing QHA.";
    pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *p_FileMESSAGE, *p_oss,
        _LOGGER_MESSAGE_);

    free();

    this->xinput = xinput;
    this->supercellopts = supercellopts;

    currentDirectory = xinput.xvasp.Directory; // remember the current directory

    // initialization of the data that is derived from the structure
    origStructure = xinput.xvasp.str;
    origStructure.ReScale(1.0);
    NatomsOrigCell = origStructure.atoms.size();
    Nbranches = NatomsOrigCell * 3;
    double Volume = origStructure.GetVolume();

    // parse QHA-related aflow.in options
    string dirname = "";
    double gp_distortion = 0.0;
    vector<double> eosrange(3);
    vector<string> tokens;
    vector<xoption>::iterator option;
    for (option  = kflags.KBIN_MODULE_OPTIONS.qhaflags.begin();
        option != kflags.KBIN_MODULE_OPTIONS.qhaflags.end(); ++option){
      if (option->keyword=="EOS") isEOS = option->option;
      if (option->keyword=="INCLUDE_ELE") includeElectronicContribution = option->option;
      if (option->keyword=="GP_FINITE_DIFF") isGP_FD = option->option;
      if (option->keyword=="IGNORE_IMAGINARY") ignore_imaginary = option->option;

      if (option->keyword=="EOS_DISTORTION_RANGE"){
        aurostd::string2tokens(option->content_string, tokens, string(" :"));

        if (tokens.size() != 3) {
          string message = "Wrong setting in the ";
          message += "[AFLOW_QHA]EOS_DISTORTION_RANGE.";
          message += "Specify as EOS_DISTORTION_RANGE=";
          message += AFLOWRC_DEFAULT_QHA_EOS_DISTORTION_RANGE;
          throw aurostd::xerror(_AFLOW_FILE_NAME_, QHA_ARUN_MODE, message,
              _INPUT_NUMBER_);
        }

        for (uint j=0; j<3; j++){
          eosrange[j] = aurostd::string2utype<double>(tokens[j]);
        }
      } // eos_distortion_range

      if (option->keyword=="GP_DISTORTION"){
        gp_distortion = option->content_double/100.0;
      }

      // note: QHA, QHA3P and SCQHA could run "simultaneously"
      if (option->keyword=="MODE"){
        tokens.clear();
        aurostd::string2tokens(option->content_string, tokens, ",");

        if ((tokens.size()<1) || (tokens.size() > 3)){
          string message = "Wrong setting in ";
          message += "[AFLOW_QHA]MODE. ";
          message += "Specify as MODE=QHA,QHA3P,SCQHA";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, QHA_ARUN_MODE, message,
              _INPUT_NUMBER_);
        }

        for (uint i=0; i<tokens.size(); i++){
          if (tokens[i].length()==3){
            if (tokens[i].find("QHA")!=std::string::npos) runQHA = true;
          }
          else if (tokens[i].length()==5){
            if (tokens[i].find("QHA3P")!=std::string::npos) runQHA3P = true;
            if (tokens[i].find("SCQHA")!=std::string::npos) runSCQHA = true;
          }
        }
      }
    }

    // determine the names for the directories for the calculation of the Grueneisen parameter
    // (calculated using finite differences method)
    vector<double> gprange(3);
    gprange[0] = 1.0-gp_distortion; gprange[1] = 1.0; gprange[2] = 1.0+gp_distortion;
    N_GPvolumes = 3;
    for (int j=0; j<N_GPvolumes; j++){
      subdirectories_apl_gp.push_back(ARUN_DIRECTORY_PREFIX + QHA_ARUN_MODE + 
          "_PHONON_" + aurostd::utype2string<double>(gprange[j],precision_format));

      coefGPVolumes.push_back(gprange[j]);
      GPvolumes.push_back(gprange[j]*Volume/NatomsOrigCell);
    }


    // determine the names for the directories used for the EOS calculation
    if (isEOS || runQHA3P || runSCQHA){
      eosrange[0] = 1.0 + eosrange[0]/100.0;
      eosrange[1] = 1.0 + eosrange[1]/100.0;
      eosrange[2] = eosrange[2]/100.0;

      N_EOSvolumes = floor((eosrange[1]-eosrange[0])/eosrange[2])+1;
      if (N_EOSvolumes < REQUIRED_MIN_NUM_OF_DATA_POINTS_FOR_EOS_FIT){
        isEOS = false;

        stringstream msg;
        msg << "QHA EOS calculation requires at least ";
        msg << aurostd::utype2string<int>(REQUIRED_MIN_NUM_OF_DATA_POINTS_FOR_EOS_FIT);
        msg << " APL calculations." << std::endl;
        msg << "The current choice of volume range and increment produces ";
        msg << aurostd::utype2string<int>(N_EOSvolumes);
        msg << " APL calculations." << std::endl;
        msg << "QHA EOS calculation will be skipped!";
        pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *p_FileMESSAGE,
            *p_oss, _LOGGER_ERROR_);
      }

      // get a set of volumes that would be used for the QHA-EOS calculation
      for (double i=eosrange[0]; i<=eosrange[1]; i+=eosrange[2]){
        subdirectories_apl_eos.push_back(ARUN_DIRECTORY_PREFIX + QHA_ARUN_MODE + 
            "_PHONON_" + aurostd::utype2string(i, precision_format));

        arun_runnames_static.push_back("STATIC_" + 
            aurostd::utype2string(i, precision_format));
        dirname = ARUN_DIRECTORY_PREFIX + QHA_ARUN_MODE + '_' +
          arun_runnames_static.back();
        subdirectories_static.push_back(dirname);

        coefEOSVolumes.push_back(i);
        EOSvolumes.push_back(i*Volume/NatomsOrigCell);
      }
    }

    tokens.clear();
    aurostd::string2tokens(tpt, tokens, string (" :"));
    if (tokens.size() != 3){
      stringstream msg;
      msg << "Wrong setting in ";
      msg << "[AFLOW_APL]TPT.";
      msg << "Specify as TPT="+AFLOWRC_DEFAULT_APL_TPT+"." << std::endl;
      msg << "See README_AFLOW_APL.TXT for the details.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, QHA_ARUN_MODE, msg,
          _INPUT_NUMBER_);
    }
    double tp_start = aurostd::string2utype<double>(tokens[0]);
    double tp_end   = aurostd::string2utype<double>(tokens[1]);
    double tp_step  = aurostd::string2utype<double>(tokens[2]);


    // define a set of temperatures for thermodynamic calculations
    Ntemperatures = floor((tp_end - tp_start)/tp_step) + 1;

    for (int T=tp_start; T<=tp_end; T+=tp_step) Temperatures.push_back(T);

    // this matrix will be used in the fit of the frequency-volume dependency
    // w(V) = a + b*V + c*V^2 + d*V^3
    gp_fit_matrix = xmatrix<double> (N_EOSvolumes, 4);
    for (int Vid=1; Vid<=N_EOSvolumes; Vid++){
      gp_fit_matrix[Vid][1] = 1;
      gp_fit_matrix[Vid][2] = EOSvolumes[Vid-1];
      gp_fit_matrix[Vid][3] = pow(EOSvolumes[Vid-1],2);
      gp_fit_matrix[Vid][4] = pow(EOSvolumes[Vid-1],3);
    }

    isInitialized = true;
  }

  /// Performs a QHA calculation.
  /// 
  /// For a regular QHA calculation, there are two possible choices:
  /// 1) calculate the Grueneisen parameter using the finite difference method.
  ///    This calculation requires 3 phonon calculations (for volumes V, V-dV and V+dV);
  /// 2) calculate temperature-dependent parameters (such as equilibrium volume,
  /// free energy, bulk modulus, thermal expansion, isochoric and isobaric specific heat,
  /// average Grueneisen parameter) employing a set of phonon calculations and
  /// making a fit to some model equation of state.
  /// It requires N_EOSvolumes static and N_EOSvolumes phonon calculations
  /// (N_EOSvolumes is determined by the EOS_DISTORTION_RANGE option).
  /// 
  /// QHA3P and SCQHA calculations require 3 phonon calculations and N_EOSvolumes static
  /// calculations.
  /// 
  void QHAN::run(_xflags &xflags, _aflags &aflags, _kflags &kflags, string &aflowin)
  {
    bool LDEBUG = (FALSE || DEBUG_QHA || XHOST.DEBUG);

    string function = "QHAN::run():";
    string msg = "Performing a QHA calculation.";
    pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *p_FileMESSAGE, *p_oss,
        _LOGGER_MESSAGE_);

    if (!isInitialized){
      msg = "QHA was not initialized properly and the QHA calculation will be aborted.";
      pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *p_FileMESSAGE,
          *p_oss, _LOGGER_ERROR_);
      return;
    }

    try{
      bool eos_data_available = false;
      bool gp_data_available = false;

      // QHA3P and SCQHA require a set of static EOS calculations.
      // But for a QHA calculation, the EOS flag is used to toggle these types of 
      // calculations.
      if (isEOS || runQHA3P || runSCQHA){
        if (checkStaticCalculations() == N_EOSvolumes){
          readStaticCalculationsData();
        }
        else{
          createSubdirectoriesStaticRun(xflags, aflags, kflags);
        }
      }

      // In a QHA calculation, the EOS flag performs APL calculations for a set of volumes.
      // This flag is used when one is interested in T-dependent properties.
      if (isEOS && runQHA){
        eos_data_available = runAPLcalculations(subdirectories_apl_eos,
            coefEOSVolumes, xflags, aflags, kflags, aflowin, false);

        if (eos_data_available){
          if (includeElectronicContribution) DOSatEf();
          if (LDEBUG) writeFrequencies();
          writeFVT();

          writeThermalProperties(EOS_POLYNOMIAL, QHA_CALC);
          writeThermalProperties(EOS_MURNAGHAN, QHA_CALC);
          writeThermalProperties(EOS_BIRCH_MURNAGHAN, QHA_CALC);
        }
      }

      // QHA3P and SCQHA require that mode-dependent Grueneisen parameters are
      // calculated via a finite difference numerical derivative.
      // For a QHA calculation, use the GP_FD flag to turn on this type of calculation.
      if (isGP_FD || runQHA3P || runSCQHA){
        gp_data_available = runAPLcalculations(subdirectories_apl_gp, coefGPVolumes,
            xflags, aflags, kflags, aflowin);

        if (gp_data_available){
          if (runQHA || runQHA3P){
            double V0 = origStructure.GetVolume()/NatomsOrigCell;
            writeGPpath(V0);
            writeGPmeshFD();
            writeAverageGPfiniteDifferences();

            if (runQHA3P){
              writeThermalProperties(EOS_POLYNOMIAL, QHA3P_CALC);
              writeThermalProperties(EOS_MURNAGHAN, QHA3P_CALC);
              writeThermalProperties(EOS_BIRCH_MURNAGHAN, QHA3P_CALC);
            }
          }

          if (runSCQHA){
            RunSCQHA(EOS_POLYNOMIAL, true);
            RunSCQHA(EOS_POLYNOMIAL, false);
          }
        }
      }
    }
    catch(aurostd::xerror &e){
      pflow::logger(e.whereFileName(),e.whereFunction(), e.what(), currentDirectory,
          *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
      return;
    }
  }

  /// Creates subdirectories with aflow.in for a set of static DFT calculations.
  ///
  void QHAN::createSubdirectoriesStaticRun(const _xflags &xflags, const _aflags &aflags,
      const _kflags &kflags)
  {
    // use static_bands calculations to get a reasonable electronic DOS
    xinput.xvasp.AVASP_flag_RUN_STATIC_BANDS       = true;
    xinput.xvasp.AVASP_flag_RUN_STATIC             = false;
    xinput.xvasp.AVASP_flag_RUN_RELAX_STATIC_BANDS = false;
    xinput.xvasp.AVASP_flag_RUN_RELAX_STATIC       = false;
    xinput.xvasp.AVASP_flag_RUN_RELAX              = false;
    xinput.xvasp.AVASP_flag_GENERATE               = false;
    xinput.xvasp.aopts.pop_attached("AFLOWIN_FLAG::MODULE");

    stringstream aflow;

    for (uint i=0; i<subdirectories_static.size(); i++){
      xinput.xvasp.str = origStructure;
      xinput.xvasp.str.InflateVolume(coefEOSVolumes[i]);

      xinput.setDirectory(currentDirectory);
      xinput.xvasp.AVASP_arun_mode = QHA_ARUN_MODE;
      xinput.xvasp.AVASP_arun_runname = arun_runnames_static[i];

      if (xinput.xvasp.AVASP_path_BANDS.empty()) xinput.xvasp.AVASP_path_BANDS = AFLOWRC_DEFAULT_BANDS_LATTICE;
      if (!xinput.xvasp.AVASP_value_BANDS_GRID)
        xinput.xvasp.AVASP_value_BANDS_GRID=DEFAULT_BANDS_GRID;
      AVASP_populateXVASP(aflags,kflags,xflags.vflags, xinput.xvasp);
      AVASP_MakeSingleAFLOWIN(xinput.xvasp, aflow, true);
    }
  }

  /// Checks if all required static calculations exist and returns the number of
  /// finished calculations.
  /// 
  /// If not all calculations are finished, it returns 0.
  /// 
  int QHAN::checkStaticCalculations()
  {
    string function = "QHAN::checkStaticCalculations():", msg = "";

    if (!isInitialized){
      msg = "QHA was not initialized properly and the QHA calculation will be aborted.";
      pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *p_FileMESSAGE,
          *p_oss, _LOGGER_ERROR_);
      return 0;
    }

    int count = 0;
    string dosfile = "", outcarfile = "";
    for (uint i=0; i<subdirectories_static.size(); i++){
      dosfile    = subdirectories_static[i]+'/'+"DOSCAR.static";
      outcarfile = subdirectories_static[i]+'/'+"OUTCAR.static";
      if ((aurostd::EFileExist(dosfile) || aurostd::FileExist(dosfile)) &&
          (aurostd::EFileExist(outcarfile) || aurostd::FileExist(outcarfile))){
        count++;
      }
      else{
        msg = "QHA is not able to proceed: the ";
        msg += subdirectories_static[i] + " directory is missing.";
        pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *p_FileMESSAGE,
            *p_oss, _LOGGER_WARNING_);
        return 0;
      }
    }

    return count;
  }

  /// Creates a set of EOS APL subdirectories with corresponding aflow.in or gathers
  /// and processes data from finished APL calculations.
  /// 
  bool QHAN::runAPLcalculations(const vector<string> &subdirectories,
      const vector<double> &coefVolumes, _xflags &xflags,
      _aflags &aflags, _kflags &kflags, string &aflowin, bool gp)
  {
    string function = "QHAN::runAPLcalculations():";
    string msg = "Reading phonon DOS and dispersion relations.";
    pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *p_FileMESSAGE, *p_oss,
        _LOGGER_MESSAGE_);

    if (!isInitialized){
      msg = "QHA was not initialized properly and the QHA calculation will be aborted.";
      pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *p_FileMESSAGE,
          *p_oss, _LOGGER_ERROR_);
      return false;
    }

    int Nqpoints = 0;
    bool apl_data_calculated = true;

    for (uint i=0; i<subdirectories.size(); i++){
      xinput.xvasp.str = origStructure;
      xinput.xvasp.str.InflateVolume(coefVolumes[i]);

      apl::PhononCalculator phcalc(*p_FileMESSAGE, *p_oss);
      phcalc.initialize_supercell(xinput.getXStr());
      phcalc.getSupercell().build(supercellopts);
      phcalc.setDirectory(subdirectories[i]);
      phcalc.setNCPUs(kflags);

      auto_ptr<apl::ForceConstantCalculator> fccalc;

      if (apl_options.getattachedscheme("ENGINE") == string("DM")){
        apl::DirectMethodPC* dmPC = new apl::DirectMethodPC(phcalc.getSupercell(),
            *p_FileMESSAGE, *p_oss);
        fccalc.reset(dmPC);

        // set options for direct method phonon calculations
        dmPC->setGeneratePlusMinus(apl_options.flag("AUTO_DIST"),apl_options.flag("DPM"));
        dmPC->setGenerateOnlyXYZ(apl_options.flag("XYZONLY"));
        dmPC->setDistortionSYMMETRIZE(apl_options.flag("SYMMETRIZE"));
        dmPC->setDistortionINEQUIVONLY(apl_options.flag("INEQUIVONLY"));
        dmPC->setDistortionMagnitude(aurostd::string2utype<double>(
              apl_options.getattachedscheme("DIST_MAGNITUDE")));
        dmPC->setCalculateZeroStateForces(apl_options.flag("ZEROSTATE"));
      }
      else{
        fccalc.reset(new apl::LinearResponsePC(phcalc.getSupercell(), *p_FileMESSAGE,
              *p_oss));
      }
      fccalc->setPolarMaterial(apl_options.flag("POLAR"));
      fccalc->setDirectory(subdirectories[i]);

      // set the name of the subdirectory
      xinput.setDirectory(subdirectories[i]);

      // if it is the first run, create APL subdirectories with corresponding aflow.in
      // files and skip to the next volume calculation
      if (fccalc->runVASPCalculations(xinput, aflags, kflags, xflags, aflowin)){
        apl_data_calculated = false;
        continue;
      }

      // check if APL calculation is ready and calculate the force constants
      if (!fccalc->run()){
        apl_data_calculated = false;
        continue;
      }
      phcalc.setHarmonicForceConstants(*fccalc);

      // calculate all phonon-related data: DOS, frequencies along the q-mesh and
      // phonon dispersions
      vector<xvector<double> > dummy_dos_projections; // do not need for QHA
      vector<int> dos_mesh(3);

      vector<string> tokens;

      aurostd::string2tokens(apl_options.getattachedscheme("DOS_MESH"), tokens, 
          string(" xX"));
      for (uint j=0; j<tokens.size(); j++){
        dos_mesh[j] = aurostd::string2utype<int>(tokens[j]);
      }

      phcalc.initialize_qmesh(dos_mesh);

      apl::DOSCalculator dosc(phcalc, apl_options.getattachedscheme("DOS_METHOD"),
          dummy_dos_projections);
      dosc.calc(aurostd::string2utype<double>(apl_options.getattachedscheme("DOS_NPOINTS")),
          aurostd::string2utype<double>(apl_options.getattachedscheme("DOS_SMEAR")));
      dosc.writePHDOSCAR(subdirectories[i]);

      if (dosc.hasNegativeFrequencies()){
        stringstream msg;
        msg << "Phonon dispersions of the APL calculation in the " << subdirectories[i];
        msg << " directory contain imaginary frequencies." << std::endl;
        if (ignore_imaginary){
          msg << "Imaginary part of the phonon dispersions and phonon DOS will be ignored.";
          msg << " Check if the results are still meaningful!" << std::endl;
          pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *p_FileMESSAGE,
              *p_oss, _LOGGER_WARNING_);
        }
        else{
          msg << "QHA calculation will be stopped after checking all available APL calculations." << std::endl;
          msg << "Please check the phonon DOS and phonon dispersions." << std::endl;
          msg << "Workaround: adjust the EOS_DISTORTION_RANGE option to exclude";
          msg << "the problematic calculations or ignore this error by setting "; 
          msg << "IGNORE_IMAGINARY=ON." << std::endl;

          pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *p_FileMESSAGE,
              *p_oss, _LOGGER_ERROR_);
          apl_data_calculated = false;
        }
      }

      if (phcalc.getSupercell().projectToPrimitive()){
        // if projection to primitive was successful update origStructure
        double Volume = origStructure.GetVolume();
        origStructure = phcalc.getSupercell().getInputStructure();
        origStructure.InflateVolume(Volume/origStructure.GetVolume());
      }
      else{
        msg = "Could not map the AFLOW standard primitive cell to the supercell. ";
        msg += "Phonon dispersions will be calculated using the original structure instead.";
        pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *p_FileMESSAGE,
            *p_oss, _LOGGER_ERROR_);
      }
      Nbranches = phcalc.getNumberOfBranches();

      string USER_DC_INITLATTICE="";
      int USER_DC_NPOINTS = aurostd::string2utype<int>(
          apl_options.getattachedscheme("BAND_NPOINTS"));
      apl::PhononDispersionCalculator pdisc(phcalc);
      pdisc.initPathLattice(USER_DC_INITLATTICE, USER_DC_NPOINTS);
      pdisc.calc(apl::THZ | apl::ALLOW_NEGATIVE);
      pdisc.writePHEIGENVAL(subdirectories[i]);

      // save all the data that is necessary for QHA calculations
      if (i==0){// qmesh data is the same for all volumes: need to store only once
        qpWeights = phcalc.getQMesh().getWeights();
        qPoints   = phcalc.getQMesh().getIrredQPointsFPOS();
        Nqpoints  = qPoints.size();
      }
      else{
        if (qpWeights.size() != phcalc.getQMesh().getWeights().size()){
          msg = "Inconsistent size of q-points weights for calculations at different volumes.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, QHA_ARUN_MODE, msg, _INDEX_MISMATCH_);
        }

        if (qPoints.size() != phcalc.getQMesh().getIrredQPointsFPOS().size()){
          msg = "Inconsistent number of irreducible q-points for calculations at different volumes.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, QHA_ARUN_MODE, msg, _INDEX_MISMATCH_);
        }
      }

      // we have two different sets of data: one for the finite difference Grueneisen
      // parameter calculation and one for the EOS APL calculation
      if (gp){
        gp_ph_dispersions.push_back(pdisc.createEIGENVAL());

        // allocate memory at the first run (entire chunk of memory is required because
        // data will be reordered)
        if (i==0){
          omegaV_mesh = vector<vector<vector<double> > > (Nqpoints,
              vector<vector<double> >(Nbranches, vector<double>(N_GPvolumes)));
        }

        vector<xvector<double> > freqs = dosc.getFreqs();

        for (int q=0; q<Nqpoints; q++){
          for (int branch=0; branch<Nbranches; branch++){
            // clean imaginary frequencies in case if user wants to ignore them
            omegaV_mesh[q][branch][i] = freqs[q][branch+1] < 0 ? 0 : freqs[q][branch+1];
          }
        }
      }
      else {
        eos_vib_thermal_properties.push_back(ThermalPropertiesCalculator(dosc,
              *p_FileMESSAGE));

        // allocate memory at the first run (an entire chunk of memory is required because
        // data will be reordered)
        if (i==0){
          omegaV_mesh_EOS = vector<vector<vector<double> > > (Nqpoints,
              vector<vector<double> >(Nbranches, vector<double>(N_EOSvolumes)));
        }

        vector<xvector<double> > freqs = dosc.getFreqs();

        for (int q=0; q<Nqpoints; q++){
          for (int branch=0; branch<Nbranches; branch++){
            // clean imaginary frequencies in case if user wants to ignore them
            omegaV_mesh_EOS[q][branch][i] = freqs[q][branch+1] < 0 ? 0 : freqs[q][branch+1];
          }
        }
      }
    }

    return apl_data_calculated;
  }

  /// Reads data from a set of static DFT calculations.
  void QHAN::readStaticCalculationsData()
  {
    string function = "QHAN::readStaticCalculationsData():";
    string msg = "";
    pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *p_FileMESSAGE, *p_oss,
        _LOGGER_MESSAGE_);

    if (!isInitialized){
      msg = "QHA was not initialized properly and the QHA calculation will be aborted.";
      pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *p_FileMESSAGE,
          *p_oss, _LOGGER_ERROR_);
      return;
    }

    xDOSCAR doscar;
    xOUTCAR outcar;
    string outcarfile = "", dosfile = "";
    for (uint i=0; i<subdirectories_static.size(); i++){
      msg = "Reading data from the static DFT calculation in the ";
      msg += subdirectories_static[i] + " directory.";
      pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *p_FileMESSAGE, 
          *p_oss, _LOGGER_MESSAGE_);

      outcarfile = subdirectories_static[i]+'/'+"OUTCAR.static";
      outcar.GetPropertiesFile(outcarfile);

      dosfile = subdirectories_static[i]+'/'+"DOSCAR.static";
      doscar.GetPropertiesFile(dosfile);

      vector<double> edos(doscar.number_energies);
      vector<double> energies = aurostd::deque2vector(doscar.venergy);

      for (uint j=0; j<doscar.number_energies; j++){
        // sum spin contributions if we deal with magnetic calculation
        for (uint s=0; s<doscar.vDOS.at(0).at(0).size(); s++){
          // if Methfessel-Paxton method was used some of the DOS values might
          // be negative => fix them to be zero
          if (doscar.vDOS[0][0][s][j]>=0){
            edos[j] += doscar.vDOS[0][0][s][j]/doscar.number_atoms;
          }
        }
      }

      edos_V.push_back(edos);
      energies_V.push_back(energies);
      Efermi_V.push_back(doscar.Efermi);
      E0_V.push_back(outcar.energy_cell/outcar.natoms);

      static_eigvals.push_back(xEIGENVAL(subdirectories_static[i]+'/'+"EIGENVAL.static"));
      static_ibzkpts.push_back(xIBZKPT(subdirectories_static[i]+'/'+"IBZKPT.static"));
    }
  }

  /// Calculates the Grueneisen parameter of an individual vibrational mode for a
  /// given volume.
  /// 
  /// gamma_i = -V/w*dw/dV,
  /// where w is a frequency of a given mode at a given volume
  /// 
  /// Volume dependence of the frequency is approximated by the following polynomial:
  /// w = a + b*V  + c*V**2 + d*V**3
  /// 
  /// @param V is a volume at which the Grueneisen parameter is calculated.
  /// @param xomega is an array of the frequency-volume dependence for a given
  /// phonon branch.
  /// 
  double QHAN::calcGrueneisen(double V, xvector<double> &xomega, double &w)
  {
    string function = "calcGrueneisen():";
    // set all weight in fit to 1
    xvector<double> s(xomega.rows); for (int i=s.lrows;i<=s.urows;i++) s[i]=1;
    aurostd::cematrix lsfit(gp_fit_matrix);

    double a=0.0, b=0.0, c=0.0, d=0.0; // fit parameters
    double gamma = 0; // mode-dependent Grueneisen parameter
    w = 0; // frequency of a given branch that corresponds to the input volume V

    /** Workaround: for w=0, gamma_i is assumed to be zero, since
      / the actual value depends on the path used to approach w->0. */
    bool freqs_are_nonzero = true;
    for (int i=xomega.lrows; i<=xomega.urows; i++){
      if (!(xomega[i] > AUROSTD_ROUNDOFF_TOL)){
        freqs_are_nonzero = false;
        break;
      }
    }

    // calculate gamma_i only if all frequencies in a given branch are nonzero
    if (freqs_are_nonzero){
      lsfit.LeastSquare(xomega,s);
      a = lsfit.AVec()[0];
      b = lsfit.AVec()[1];
      c = lsfit.AVec()[2];
      d = lsfit.AVec()[3];

#if DEBUG_QHA_GP_FIT
      double err = 0;
      xvector<double> tmp(4);
      tmp[1] = a; tmp[2] = b; tmp[3] = c; tmp[4] = d;

      // check the fit error
      for (int Vid=0; Vid<N_EOSvolumes; Vid++){
        err = abs(a+b*EOSvolumes[Vid]+c*pow(EOSvolumes[Vid],2)+d*pow(EOSvolumes[Vid],3)
            -xomega[Vid+1])/xomega[Vid+1];
        err *= 100.0;
        if (err>=10.0){
          string msg="Relative error of the log(w)=f(V) fit (used to ";
          msg+="determine mode-decomposed Grueneisen parameter) is larger than ";
          msg+="10\% for V="+aurostd::utype2string<double>(EOSvolumes[Vid]);

          pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory,
              *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
          cerr << function << " original frequencies:" << std::endl;
          cerr << function << xomega << std::endl;
          cerr << function << " fit to frequencies:" << std::endl;
          cerr << function << gp_fit_matrix * tmp << std::endl;
          break;
        }
      }
#endif

      w = a + b*V + c*pow(V,2) + d*pow(V,3);
      gamma = -V/w * (b + 2*c*V + 3*d*pow(V,2));
    }

    return gamma;
  }

  /// Calculates the Grueneisen parameter of an individual vibrational mode using
  /// the central finite difference method.
  /// 
  /// gamma = -V0/w0*dw/dV |V->V0 ~ -V0/w0*(w(V0+dV)-w(V0-dV))/(2*dV)
  /// 
  double QHAN::calcGrueneisenFD(const xvector<double> &xomega)
  {
    string function = "calcGrueneisenFD(): ", msg = "";
    if (xomega.rows != 3){
      msg = "Wrong size of the xomega array passed to calcGrueneisenFD function.";
      msg += "Expected size: 3. Actual size: ";
      msg += aurostd::utype2string<int>(xomega.rows);
      throw aurostd::xerror(_AFLOW_FILE_NAME_, QHA_ARUN_MODE, msg, _INDEX_BOUNDS_);
    }

    if (xomega[2] > _ZERO_TOL_){
      return -GPvolumes[1]/xomega[2]*(xomega[3]-xomega[1])/
        (GPvolumes[2]-GPvolumes[0]);
    }
    else{
      return 0.0;
    }
  }


  /// Calculates the free energy (without electronic contribution) as a function of 
  /// volume and temperature.
  /// The volume dependency is obtained via a fit to the model equation of state.
  /// @param qha_method defines what kind of QHA calculation is performed.
  ///
  double QHAN::FreeEnergyFit(double T, double V, EOSmethod eos_method,
      QHAmethod qha_method)
  {
    string function = "FreeEnergyFit():", msg = "";
    if (T<_ZERO_TOL_) return 0;

    xvector<double> E(N_EOSvolumes);
    switch(qha_method){
      case (QHA3P_CALC):
        for (int i=E.lrows; i<=E.urows; i++){
          E[i] = FreeEnergyTaylorExpansion(T, i-1);
        }
        break;
      case(QHA_CALC):
        for (int i=E.lrows; i<=E.urows; i++){
          E[i] = FreeEnergy(T, i-1);
        }
        break;
      default:
        msg = "Nonexistent QHA method was passed to " + function;
        throw aurostd::xerror(_AFLOW_FILE_NAME_, QHA_ARUN_MODE, msg, _INPUT_UNKNOWN_);
        break;
    }
    if (includeElectronicContribution) E += electronicFreeEnergySommerfeld(T);
    xvector<double> p = fitToEOSmodel(E, eos_method);
    return evalEOSmodel(V, p, eos_method);
  }

  /// Calculates the internal energy as a function of volume and temperature.
  ///
  double QHAN::InternalEnergyFit(double T, double V)
  {
    static xvector<double> U(N_EOSvolumes);
    static xvector<double> xvolumes = aurostd::vector2xvector(EOSvolumes);

    if (T<0) return 0;

    for (int i=0; i<N_EOSvolumes; i++){
      U[i+1] = eos_vib_thermal_properties[i].getInternalEnergy(T, apl::eV);
      U[i+1] /= NatomsOrigCell;
    }

    return EOSpoly(V, fitEOSpoly(xvolumes, U));
  }

  /// Calculates the free energy (without electronic contribution) for a calculation at
  /// a specific volume (given by its id in the EOSvolumes array).
  ///
  double QHAN::FreeEnergy(double T, int id)
  {
    return eos_vib_thermal_properties[id].getVibrationalFreeEnergy(T, apl::eV)
      /NatomsOrigCell + E0_V[id];
  }

  /// Fits the (free) energy-volume dependency to one of the following equation of state
  /// models:
  /// polynomial: (eq (5.2) page 106)  https://doi.org/10.1017/CBO9781139018265
  /// Murnaghan: Proceedings of the National Academy of Sciences of the United States of 
  /// America, 30 (9): 244–247
  /// Birch-Murnaghan: Physical Review. 71 (11): 809–824
  ///
  xvector<double> QHAN::fitToEOSmodel(xvector<double> &E, EOSmethod method)
  {
    string function = "EOSfit::fit():", msg = "";
    xvector<double> V = aurostd::vector2xvector(EOSvolumes);
    xvector<double> fit_params;
    xvector<double> guess(4);
    switch(method){
      case(EOS_POLYNOMIAL):
        fit_params = fitEOSpoly(V, E);
        EOS_volume_at_equilibrium = EOSpolyGetEqVolume(fit_params, min(V), max(V));
        EOS_energy_at_equilibrium = EOSpoly(EOS_volume_at_equilibrium, fit_params);
        EOS_bulk_modulus_at_equilibrium = BulkModulus(EOS_volume_at_equilibrium,fit_params);
        EOS_Bprime_at_equilibrium  = Bprime(EOS_volume_at_equilibrium, fit_params);
        break;
      case(EOS_BIRCH_MURNAGHAN):
        {
          guess[1] = min(E);
          guess[2] = (max(V)+min(V))/2;
          guess[3] = V[1]*(E[3]-2*E[2]+E[1])/pow(V[2]-V[1],2); // B from central differences
          guess[4] = 3.5; // a reasonable initial value for most materials

          NonlinearFit bmfit(V,E,guess,BirchMurnaghan);
          bmfit.fitLevenbergMarquardt();

          fit_params = bmfit.p;
          EOS_energy_at_equilibrium = bmfit.p[1];
          EOS_volume_at_equilibrium = bmfit.p[2];
          EOS_bulk_modulus_at_equilibrium   = bmfit.p[3]*eV2GPa;
          EOS_Bprime_at_equilibrium  = bmfit.p[4];
        }
        break;
      case(EOS_MURNAGHAN):
        {
          guess[1] = min(E);
          guess[2] = (max(V)+min(V))/2;
          guess[3] = V[1]*(E[3]-2*E[2]+E[1])/pow(V[2]-V[1],2); // B from central differences 
          guess[4] = 3.5; // a reasonable initial value for most materials

          NonlinearFit nlfit(V,E,guess,Murnaghan);
          nlfit.fitLevenbergMarquardt();

          fit_params = nlfit.p;
          EOS_energy_at_equilibrium = nlfit.p[1];
          EOS_volume_at_equilibrium = nlfit.p[2];
          EOS_bulk_modulus_at_equilibrium   = nlfit.p[3]*eV2GPa;
          EOS_Bprime_at_equilibrium  = nlfit.p[4];
        }
        break;
      default:
        msg = "Nonexistent EOS method was passed to " + function;
        throw aurostd::xerror(_AFLOW_FILE_NAME_, QHA_ARUN_MODE, msg, _INPUT_UNKNOWN_);
        break;
    }

    return fit_params;
  }

  /// Returns the (free) energy at a given volume for a chosen EOS model.
  double QHAN::evalEOSmodel(double V, const xvector<double> &p, EOSmethod method)
  {
    string function = "EOSfit::eval():", msg = "";
    double energy = 0;

    static xvector<double> dydp(4);
    switch(method){
      case(EOS_POLYNOMIAL):
        energy = EOSpoly(V, p);
        break;
      case(EOS_BIRCH_MURNAGHAN):
        energy = BirchMurnaghan(V, p, dydp);
        break;
      case(EOS_MURNAGHAN):
        energy = Murnaghan(V, p, dydp);
        break;
      default:
        msg = "Nonexistent EOS method was passed to " + function;
        throw aurostd::xerror(_AFLOW_FILE_NAME_, QHA_ARUN_MODE, msg, _INPUT_UNKNOWN_);
        break;
    }

    return energy;
  }

  /// Calculates the electronic free energy using integration over the electronic DOS.
  /// Currently this function does not work (to be fixed in the future)
  ///
  /// @param id selects the volume at which corresponding data was obtained 
  ///
  double QHAN::electronicFreeEnergy(double T, int id)
  {
    double EelecT = 0; double SelecT = 0; double FelecT = 0; double Eelec0K = 0;
    double f = 0.0, f0 = 0.0;

    double dE = (max(energies_V[id])-min(energies_V[id]))/(energies_V[id].size()-1);
    for (uint i=0; i<edos_V[id].size(); i++){
      f0 = aurostd::FermiDirac(energies_V[id][i], Efermi_V[id], 0);
      f  = aurostd::FermiDirac(energies_V[id][i], Efermi_V[id], T);

      Eelec0K += energies_V[id][i] * edos_V[id][i] * f0;
      EelecT  += energies_V[id][i] * edos_V[id][i] * f;
      // limit of values for f->0 or f->1 is 0
      if (f>0 && f<1) SelecT -= (f*log(f) + (1-f)*log(1-f)) * edos_V[id][i];
    }

    EelecT  *= dE;
    Eelec0K *= dE;
    SelecT  *= KBOLTZEV*dE;

    FelecT = (EelecT-Eelec0K) - T*SelecT;

    return FelecT;
  }

  /// Calculates the equilibrium volume at a given temperature.
  /// @param qha_method defines what kind of QHA calculation is performed.
  /// 
  double QHAN::getEqVolumeT(double T, EOSmethod eos_method, QHAmethod qha_method)
  {
    string function = "getEqVolumeT():", msg = "";
    xvector<double> E(N_EOSvolumes);
    switch(qha_method){
      case (QHA3P_CALC):
        for (int i=E.lrows; i<=E.urows; i++){
          E[i] = FreeEnergyTaylorExpansion(T, i-1);
        }
        break;
      case(QHA_CALC):
        for (int i=E.lrows; i<=E.urows; i++){
          E[i] = FreeEnergy(T, i-1);
        }
        break;
      default:
        msg = "Nonexistent QHA method was passed to " + function;
        throw aurostd::xerror(_AFLOW_FILE_NAME_, QHA_ARUN_MODE, msg, _INPUT_UNKNOWN_);
        break;
    }

    if (includeElectronicContribution) E += electronicFreeEnergySommerfeld(T);
    fitToEOSmodel(E, eos_method);

    return EOS_volume_at_equilibrium;
  }


  /// Calculates the average Grueneisen parameter (GP) and the isochoric specific 
  /// heat (CV) at a given temperature using a weighted sum over the k-points mesh.
  /// This function is used when the Grueneisen parameter is calculated using the finite
  /// difference method.
  ///
  void QHAN::calcCVandGP(double T, double &CV, double &GP)
  {
    CV = 0; GP = 0;
    if (T<_ZERO_TOL_) return;

    xvector<double> xomega;

    uint NIrQpoints = omegaV_mesh.size();
    int NQpoints = 0; //  total number of q-points

    double Cvi  = 0;   // mode-dependent specific heat at V=const
    double w    = 0;   // frequency for a given volume
    double expx = 0;   // temperature-dependent exponential factor
    double beta = 1.0/KBOLTZEV/T;

    for (uint q=0; q<NIrQpoints; q++){
      for (uint branch=0; branch<omegaV_mesh[q].size(); branch++){
        xomega = aurostd::vector2xvector(omegaV_mesh[q][branch]);

        w = xomega[2]*THz2Hz*PLANCKSCONSTANTEV_h; // [THz] -> [eV]
        if (w>0){
          expx = exp(beta*w);

          Cvi = pow(w,2)*expx/pow(expx-1.0,2) * qpWeights[q];

          GP += calcGrueneisenFD(xomega) * Cvi;
          CV += Cvi;
        }
      }
      NQpoints += qpWeights[q];
    }

    GP /= CV;
    CV /= NQpoints; CV /= Nbranches;
    CV *= 3*pow(beta,2);
  }

  /// Calculates the average Grueneisen parameter (GP) and the isochoric specific 
  /// heat (CV) at a given volume and a given temperature.
  /// This function is used when the Grueneisen parameter is calculated using a fit to the
  /// w(V) relation, that is obtained from a set of EOS calculations.
  ///
  void QHAN::calcCVandGPfit(double T, double V, double &CV, double &GP)
  {
    CV = 0; GP = 0;
    if (T<_ZERO_TOL_) return;

    xvector<double> xomega;

    uint NIrQpoints = omegaV_mesh_EOS.size();
    int NQpoints = 0; // the total number of q-points

    double Cvi  = 0;  // mode-dependent specific heat at V=const
    double w    = 0;  // frequency for a given volume
    double expx = 0;  // temperature-dependent exponential factor
    double gamma = 0;
    double beta = 1.0/KBOLTZEV/T;

    for (uint q=0; q<NIrQpoints; q++){
      for (int branch=0; branch<Nbranches; branch++){
        xomega = aurostd::vector2xvector(omegaV_mesh_EOS[q][branch]);

        gamma = calcGrueneisen(V, xomega, w);
        w *= THz2Hz*PLANCKSCONSTANTEV_h; // [THz] -> [eV]
        if (w > AUROSTD_ROUNDOFF_TOL){
          expx = exp(w*beta);

          Cvi = pow(w,2)*expx/pow(expx-1.0,2) * qpWeights[q];

          GP += gamma * Cvi;
          CV += Cvi;
        }
      }
      NQpoints += qpWeights[q];
    }

    GP /= CV;
    CV /= NQpoints; CV /= NatomsOrigCell;
    CV *= pow(beta,2);
  }

  /// Calculates the volumetric thermal expansion coefficient (beta).
  /// 
  /// Central finite differences are used to calculate beta = 1/V dV/dT.
  /// Values of V(T+dT) and V(T-dT) are taken from the EOS fit to the corresponding
  /// free energies.
  /// 
  /// For T=0 K, the function returns 0. Negative temperatures are ignored.
  /// 
  /// @param T is the temperature at which the calculation is done.
  /// @param eos_method defines which model is used for the EOS fit
  /// @return volumetric thermal expansion coefficient.
  /// 
  double QHAN::ThermalExpansion(double T, EOSmethod eos_method, QHAmethod qha_method)
  {
    if (!(T>0)) return 0;

    double dT = DCOEFF*T;
    return 0.5*(getEqVolumeT(T+dT,eos_method,qha_method)
        -getEqVolumeT(T-dT,eos_method,qha_method))/
      dT/getEqVolumeT(T,eos_method,qha_method);
  }

  /// Calculates the isochoric specific heat as a temperature derivative of the free
  /// energy using the central finite differences method.
  /// 
  /// @param eos_method defines which model is used for the EOS fit
  /// 
  double QHAN::IsochoricSpecificHeat(double T, double V, EOSmethod eos_method,
      QHAmethod qha_method)
  {
    double dT = DCOEFF*T;
    double CV = 0;
    if (T>0){
      CV = -(FreeEnergyFit(T+dT, V, eos_method, qha_method)
          -2*FreeEnergyFit(T, V, eos_method, qha_method)
          +FreeEnergyFit(T-dT ,V, eos_method, qha_method));
      CV *= T/pow(dT,2);
    }

    return CV;
  }

  /// Calculates the entropy as a temperature derivative of the free energy using
  /// the central finite differences method.
  /// 
  /// @param eos_method defines which model is used for EOS fit
  /// 
  double QHAN::Entropy(double T, double V, EOSmethod eos_method, QHAmethod qha_method)
  {
    double dT = DCOEFF*T;
    return -0.5*(FreeEnergyFit(T+dT,V,eos_method,qha_method)
        -FreeEnergyFit(T-dT,V,eos_method,qha_method))/dT;
  }

  /// Calculates the DOS value at the Fermi level using the linear tetrahedron method. 
  /// For the details check: https://doi.org/10.1103/PhysRevB.49.16223
  ///
  xvector<double> QHAN::DOSatEf()
  {
    xvector<double> energy_tetrahedron(4);//energy eigenvalues at the corners of the tetrahedron
    double weighted_volume = 0; // weighted volume of tetrahedron
    double DEf = 0.0; // DOS at Fermi level
    double Ef = 0.0;  // Fermi energy
    // Eij = energy_tetrahedron[i] - energy_tetrahedron[j]
    double E21 = 0.0, E31 = 0.0, E41 = 0.0, E32 = 0.0, E42 = 0.0, E43 = 0.0; 

    DOS_Ef =  xvector<double>(N_EOSvolumes);

    for (int id=0; id<N_EOSvolumes; id++){
      DEf = 0;
      for (uint band=0; band<static_eigvals[id].number_bands; band++){
        for (uint i=0; i<static_ibzkpts[id].ntetrahedra; i++){
          weighted_volume = static_ibzkpts[id].wtetrahedra * static_ibzkpts[id].vtetrahedra[i][1];

          for (uint s=0; s<static_eigvals[id].spin+1; s++){
            for (int j=1; j<=4; j++)
              energy_tetrahedron[j] = static_eigvals[id].venergy[static_ibzkpts[id].vtetrahedra[i][j+1]-1][band][s];

            energy_tetrahedron = aurostd::sort(energy_tetrahedron);
            Ef = Efermi_V[id];

            if (Ef < energy_tetrahedron[1] && energy_tetrahedron[4] < Ef) DEf += 0;

            if (energy_tetrahedron[1] < Ef && Ef <= energy_tetrahedron[2]){
              E21 = energy_tetrahedron[2] - energy_tetrahedron[1];
              E31 = energy_tetrahedron[3] - energy_tetrahedron[1];
              E41 = energy_tetrahedron[4] - energy_tetrahedron[1];

              DEf += weighted_volume*3*pow(Ef-energy_tetrahedron[1],2)/(E21*E31*E41);
            }
            else if (energy_tetrahedron[2] < Efermi_V[id] && Efermi_V[id] <= energy_tetrahedron[3]){
              E21 = energy_tetrahedron[2] - energy_tetrahedron[1];
              E31 = energy_tetrahedron[3] - energy_tetrahedron[1];
              E41 = energy_tetrahedron[4] - energy_tetrahedron[1];

              E32 = energy_tetrahedron[3] - energy_tetrahedron[2];
              E42 = energy_tetrahedron[4] - energy_tetrahedron[2];

              DEf += weighted_volume/(E31*E41)*(3*E21+6*(Ef-energy_tetrahedron[2])-3*(E31+E42)*pow(Ef-energy_tetrahedron[2],2)/(E32*E42));
            }
            else if (energy_tetrahedron[3] < Ef && Ef <= energy_tetrahedron[4]){
              E41 = energy_tetrahedron[4] - energy_tetrahedron[1];
              E42 = energy_tetrahedron[4] - energy_tetrahedron[2];
              E43 = energy_tetrahedron[4] - energy_tetrahedron[3];

              DEf += 3*weighted_volume*pow(energy_tetrahedron[4]-Ef,2)/(E41*E42*E43);
            }
          }
        }
      }
      DOS_Ef[id+1] = (2-static_eigvals[id].spin)*DEf; // factor 2 if non-magnetic and 1 otherwise
    }
    return DOS_Ef;
  }

  /// Calculates the electronic free energy using the Sommerfeld expansion.
  xvector<double> QHAN::electronicFreeEnergySommerfeld(double T)
  {
    xvector<double> F_Som = DOS_Ef;
    for (int i=F_Som.lrows; i<=F_Som.urows; i++){
      F_Som[i] *= -pow(M_PI*KBOLTZEV*T,2)/6.0;
    }
    return F_Som;
  }

  // QHA3P-related functions

  /// Returns the frequency extrapolated by a Taylor expansion.
  ///
  double QHAN::extrapolateFrequency(double V, const xvector<double> &xomega)
  {
    double dwdV = (xomega[1]-xomega[3])/(GPvolumes[0]-GPvolumes[2]);
    double d2wdV2 = (xomega[1]+xomega[3]-2.0*xomega[2])/
      pow(0.5*(GPvolumes[0]-GPvolumes[2]),2);

    double dV = (V-GPvolumes[1]);
    return xomega[2] + dwdV*dV + 0.5*d2wdV2*pow(dV,2);
  }

  /// Returns the Grueneisen parameter extrapolated by a Taylor expansion.
  ///
  double QHAN::extrapolateGamma(double V, const xvector<double> &xomega)
  {
    double gamma = 0;
    double V0 = GPvolumes[1];
    double w  = extrapolateFrequency(V, xomega);

    double dwdV = 0, d2wdV2 = 0;

    if (w>1e-6){
      dwdV = (xomega[1]-xomega[3])/(GPvolumes[0]-GPvolumes[2]);
      d2wdV2 = (xomega[1]+xomega[3]-2.0*xomega[2])/
        pow(0.5*(GPvolumes[0]-GPvolumes[2]),2);

      gamma = -V/w*(dwdV + d2wdV2 * (V-V0));
    }

    return gamma;
  }

  /// Calculates the free energy (total energy + vibrational contribution) obtained 
  /// using a Taylor expansion of the frequencies.
  /// This function is used in the QHA3P and SCQHA methods.
  /// 
  double QHAN::FreeEnergyTaylorExpansion(double T, int Vid)
  {
    if (T < _ZERO_TOL_) return E0_V[Vid];

    double w = 0.0; // extrapolated frequency at V_id volume
    double F = 0.0; // Free energy
    double fi = 0.0;
    double beta = 1.0/KBOLTZEV/T;
    int NQpoints = 0;
    xvector<double> xomega;
    for (uint q=0; q<omegaV_mesh.size(); q++){
      for (int branch=0; branch<Nbranches; branch++){
        xomega = aurostd::vector2xvector(omegaV_mesh[q][branch]);

        w = extrapolateFrequency(EOSvolumes[Vid], xomega) * THz2Hz *
          PLANCKSCONSTANTEV_h;
        fi = 0.5*w;

        if (w> AUROSTD_ROUNDOFF_TOL && T>_ZERO_TOL_) fi += KBOLTZEV*T*log(1-exp(-w*beta));

        fi *= qpWeights[q];
        F += fi;
      }
      NQpoints += qpWeights[q];
    }

    F /= NQpoints;
    F /= NatomsOrigCell;

    return F + E0_V[Vid];
  }

  /// Calculates the vibrational internal energy obtained using a Taylor expansion of
  /// the frequencies.
  /// This function is used in the QHA3P and SCQHA methods.
  double QHAN::InternalEnergyTaylorExpansion(double T, double V)
  {
    if (T<_ZERO_TOL_) return 0.0;

    double U = 0.0, ui = 0.0,  w = 0.0; 
    double beta = 1.0/KBOLTZEV/T;
    int NQpoints = 0;
    xvector<double> xomega;
    for (uint q=0; q<omegaV_mesh.size(); q++){
      for (int branch=0; branch<Nbranches; branch++){
        xomega = aurostd::vector2xvector(omegaV_mesh[q][branch]);

        w = extrapolateFrequency(V, xomega) * THz2Hz * PLANCKSCONSTANTEV_h;
        ui = 0.5*w;

        if (w>_mm_epsilon && T>_mm_epsilon) ui += w/(exp(w*beta)-1.0);

        ui *= qpWeights[q];
        U += ui;
      }
      NQpoints += qpWeights[q];
    }

    U /= NQpoints;
    U /= NatomsOrigCell;

    return U;
  }

  // Definition of the functions used by SCQHA.
  // Implementation is based on http://dx.doi.org/10.1103/PhysRevMaterials.3.073801
  // and https://doi.org/10.1016/j.commatsci.2016.04.012

  /// Calculates the phononic pressure multiplied by volume for a given temperature and
  /// a given volume.
  /// Check http://dx.doi.org/10.1103/PhysRevMaterials.3.073801
  /// for more details
  /// 
  double QHAN::VPgamma(double T, double V)
  {
    double VPgamma = 0.0, ui = 0.0,  w = 0.0; 
    double beta = 1.0/KBOLTZEV/T;
    int NQpoints = 0;
    xvector<double> xomega;
    for (uint q=0; q<omegaV_mesh.size(); q++){
      for (int branch=0; branch<Nbranches; branch++){
        xomega = aurostd::vector2xvector(omegaV_mesh[q][branch]);

        w = extrapolateFrequency(V, xomega) * THz2Hz * PLANCKSCONSTANTEV_h;
        ui = 0.5*w;

        if (w>_mm_epsilon && T>_mm_epsilon) ui += w/(exp(w*beta)-1.0);

        ui *= qpWeights[q];
        VPgamma += ui * extrapolateGamma(V, xomega);
      }
      NQpoints += qpWeights[q];
    }

    VPgamma /= NQpoints;
    VPgamma /= NatomsOrigCell;

    return VPgamma;
  }

  /// Performs SCQHA calculations.
  /// There are two implementations:
  /// 1) perform a self-consistent loop for the initial nonzero temperature and extrapolate
  /// the volume at the next temperature step using V(T+dT) ~ (1+beta dT)*V.
  /// Expect it to be inaccurate at high temperatures.
  /// 2) perform a self-consistent loop for each temperature.
  ///
  /// 1) corresponds to the original implementation as described by the following papers:
  /// http://dx.doi.org/10.1103/PhysRevMaterials.3.073801
  /// and https://doi.org/10.1016/j.commatsci.2016.04.012
  void QHAN::RunSCQHA(EOSmethod method, bool all_iterations_self_consistent)
  {
    const int max_scqha_iteration = 10000;
    const double dV = 1e-3;
    const double Vtol = 1e-5;

    string function = "QHAN::RunSCQHA():", msg = "";
    pflow::logger(QHA_ARUN_MODE, function, "", currentDirectory, *p_FileMESSAGE, *p_oss,
        _LOGGER_MESSAGE_);

    // get the equilibrium volume from the fit to the EOS model fitted to the set of 
    // energies obtained from static DFT calculations
    xvector<double> E = aurostd::vector2xvector<double>(E0_V);
    xvector<double> fit_params = fitToEOSmodel(E, method);
    double Pe = 0.0, VPg = 0.0; // electronic pressure and volume multiplied by phononic pressure
    // to avoid division by zero in the self-consistent loop,
    // the initial volume is taken to be 10% bigger
    double V = 1.1*EOS_volume_at_equilibrium;
    double Vnew = 0.0;

    // self-consistent loop for 0 K
    double V0K = 1.1*EOS_volume_at_equilibrium;
    for (int i=0; i<max_scqha_iteration; i++){
      Pe   = dEOSpoly(V, fit_params);
      VPg  = VPgamma(0.0,V);
      Vnew = VPg/Pe;
      if (abs(V0K - Vnew)/V > Vtol) V0K += (Vnew - V0K) * dV; else break;
    }

    // the name of the output file depends on the EOS fit method
    stringstream file;
    string filename = DEFAULT_QHA_FILE_PREFIX + "scqha.";
    string sc = all_iterations_self_consistent ? "sc." : "";
    switch(method){
      case(EOS_POLYNOMIAL):
        filename += sc+"polynomial.";
        break;
      case(EOS_BIRCH_MURNAGHAN):
        filename += sc+"birch-murnaghan.";
        break;
      case(EOS_MURNAGHAN):
        filename += sc+"murnaghan.";
        break;
      default:
        msg = "Nonexistent EOS method was passed to " + function;
        throw aurostd::xerror(_AFLOW_FILE_NAME_, QHA_ARUN_MODE, msg, _INPUT_UNKNOWN_);
        break;
    }
    filename += DEFAULT_QHA_THERMO_FILE;

    file.precision(10);

    file << AFLOWIN_SEPARATION_LINE << std::endl;
    file << "[SCQHA_THERMAL_PROPERTIES]START" << std::endl;

    // print header
    file << setw(5)  << "# T[K]"               << setw(SW) << ' ' <<
      setw(TW) << "Veq[ev/atom]"         << setw(SW) << ' ' <<
      setw(TW) << "F(Veq)[eV/atom]"      << setw(SW) << ' ' <<
      setw(TW) << "B[GPa]"               << setw(SW) << ' ' <<
      setw(TW) << "beta[10^-6/K]"        << setw(SW) << ' ' <<
      setw(TW) << "Cv[kB/atoms]"         << setw(SW) << ' ' <<
      setw(TW) << "Cp[kB/atoms]"         << setw(SW) << ' ' <<
      setw(TW) << "gamma"                << setw(SW) << ' ' <<
      std::endl;

    double T = Temperatures[0];

    // self-consistent loop for the initial temperature (defined by user)
    int iter=1;
    while (iter <= max_scqha_iteration){
      Pe   = dEOSpoly(V, fit_params);
      VPg  = VPgamma(T, V);
      Vnew = VPg/Pe;
      if (abs(V - Vnew)/V > Vtol) V += (Vnew - V) * dV; else break;
      iter++;
    }
    if (iter == max_scqha_iteration){
      msg="Maximum number of iterations in self consistent loop is reached";
      msg += " at T="+aurostd::utype2string<double>(T)+"K.";
      pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *p_FileMESSAGE, 
          *p_oss, _LOGGER_MESSAGE_);
    }

    double dT = (Temperatures[Ntemperatures-1]-Temperatures[0])/(Ntemperatures-1);

    // a set of variables predefined here
    double Cvi     = 0.0;  // mode-dependent specific heat at V=const
    double w       = 0.0;  // frequency for a given volume
    double w0K     = 0.0;  // frequency for a volume at T=0K
    double expx    = 0.0;  // temperature-dependent exponential factor
    double expx0   = 0.0;  // temperature-dependent exponential factor
    double ui      = 0.0;  // mode-dependent internal energy
    double Bgamma  = 0.0;  // contribution to bulk modulus of 2nd order gamma component
    double Bdgamma = 0.0;  // check https://doi.org/10.1016/j.commatsci.2016.04.012
    double Belec   = 0.0;  // "electronic" bulk modulus
    double B       = 0.0;  // total bulk modulus
    double Pgamma  = 0.0;  // "phononic" pressure
    double gamma   = 0.0;  // Grueneisen parameter
    double fi      = 0.0;  // mode-dependet free energy
    double Feq     = 0.0;  // total free energy for equilibrium volume at given T
    double CP      = 0.0;  // isobaric specific heat
    double CV      = 0.0;  // isochoric specific heat
    double GP      = 0.0;  // average Grueneisen parameter
    double beta    = 0.0;  // coefficient of thermal expansion
    double d2wdV2  = 0.0;  // second derivative of frequency w.r.t volume
    double betaT   = 1.0/KBOLTZEV/T;

    uint NIrQpoints = omegaV_mesh.size(); // number of irreducible q-points
    int NQpoints = 0; // the total number of q-points (to be determined in the loop)

    for (int Tid=0; Tid<Ntemperatures; Tid++){
      T = Temperatures[Tid];
      betaT = 1.0/KBOLTZEV/T;

      // calculate the next equilibrium volume using a self-consistent loop if beta=0 or
      // if the user wants so
      if (all_iterations_self_consistent || !(abs(beta)>0)){
        iter=1;
        while (iter <= max_scqha_iteration){
          Pe   = dEOSpoly(V, fit_params);
          VPg  = VPgamma(T, V);
          Vnew = VPg/Pe;
          if (abs(V - Vnew)/V > Vtol) V += (Vnew - V) * dV; else break;
          iter++;
        }
        if (iter == max_scqha_iteration){
          string msg="Maximum number of iterations in self consistent loop is reached";
          msg += " at T="+aurostd::utype2string<double>(T)+"K";
          pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *p_FileMESSAGE, 
              *p_oss, _LOGGER_MESSAGE_);
        }
      }
      else{
        V *= (1 + beta*dT);
      }

      // calculate thermodynamic properties
      Cvi = 0.0; CP = 0.0; CV = 0.0; 
      GP = 0.0; gamma = 0.0;
      w = 0.0; w0K = 0.0; expx = 0.0; expx0 = 0.0;
      ui = 0.0;
      Bgamma = 0.0; Bdgamma = 0.0; Belec = 0.0; Pgamma = 0.0; B = 0.0;
      fi = 0.0; Feq = 0.0;

      for (uint q=0; q<NIrQpoints; q++){
        for (uint branch=0; branch<omegaV_mesh[q].size(); branch++){
          xvector<double> xomega = aurostd::vector2xvector(omegaV_mesh[q][branch]);
          w = extrapolateFrequency(V, xomega);
          w *= THz2Hz*PLANCKSCONSTANTEV_h; // [THz] -> [eV]

          w0K = extrapolateFrequency(V0K, xomega);
          w0K *= THz2Hz*PLANCKSCONSTANTEV_h; // [THz] -> [eV]

          ui = 0.5*w;
          fi = 0.5*w;
          if (w > 1e-12 && T > 0){
            expx = exp(w*betaT);

            // use volume at T=0K for calculation of isochoric specific heat
            expx0 = exp(w0K*betaT);
            Cvi = pow(w,2)*expx0/pow(expx0-1.0,2) * qpWeights[q];

            ui  += w/(expx - 1.0);
            ui  *= qpWeights[q];

            fi += KBOLTZEV*T*log(1-exp(-w*betaT));
            fi *= qpWeights[q];

            d2wdV2 = (xomega[1]+xomega[3]-2.0*xomega[2])/
              pow(0.5*(GPvolumes[0]-GPvolumes[2]),2);
            d2wdV2 *= THz2Hz*PLANCKSCONSTANTEV_h;

            gamma = extrapolateGamma(V, xomega);

            Bgamma  += (ui - T*Cvi*pow(betaT,2)*KBOLTZEV)*pow(gamma, 2);
            Bdgamma -= ui*((1+gamma)*gamma - pow(V,2)/w*d2wdV2);

            GP += extrapolateGamma(V, xomega) * Cvi;
            CV += Cvi;
            Feq  += fi;
          }
        }
        NQpoints += qpWeights[q];
      }

      Feq /= NQpoints;
      Feq /= NatomsOrigCell;
      Feq += evalEOSmodel(V, fit_params, method);

      Bgamma /= NQpoints; Bdgamma /= NQpoints;
      Bgamma /= NatomsOrigCell; Bdgamma /= NatomsOrigCell;
      Bgamma /= V; Bdgamma /= V;

      Belec  = BulkModulus(V, fit_params);
      Pgamma = VPgamma(T, V)/V;

      B = Belec + (Bgamma + Bdgamma + Pgamma)*eV2GPa;

      GP /= CV;
      CV /= NQpoints; CV /= Nbranches;
      CV *= 3*pow(betaT,2); // [kB/atom]

      beta = KBOLTZEV*CV*GP/V/(B/eV2GPa); // [K^-1]

      CP = CV + V*T*B*pow(beta,2)/eV2GPa/KBOLTZEV; // [kB/atom]
      file << setw(5)  << T                   << setw(SW) << ' ' <<
        setw(TW) << V                   << setw(SW) << ' ' <<
        setw(TW) << Feq                 << setw(SW) << ' ' <<
        setw(TW) << B                   << setw(SW) << ' ' <<
        setw(TW) << beta * 1e6          << setw(SW) << ' ' << //[10^-6/K]
        setw(TW) << CV                  << setw(SW) << ' ' <<
        setw(TW) << CP                  << setw(SW) << ' ' <<
        setw(TW) << GP                  << setw(SW) << ' ' <<
        std::endl;
    }
    file << "[SCQHA_THERMAL_PROPERTIES]STOP" << std::endl;
    file << AFLOWIN_SEPARATION_LINE << std::endl;

    if (!aurostd::stringstream2file(file, filename)){
      msg = "Error writing to " + filename + "file.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function,msg,_FILE_ERROR_);
    }
  }

  ////////////////////////////////////////////////////////////////////////////////
  // output-related functions

  /// Writes temperature-dependent properties to file.
  ///
  /// List of properties: temperature, volume, free energy, bulk modulus,
  /// thermal expansion coefficient, specific heat at const V, specific heat at
  /// const P, average Grueneisen parameter calculated from mode-dependent
  /// Grueneisen parameters and average Grueneisen calculated from the thermal
  /// expansion coefficient.
  ///
  void QHAN::writeThermalProperties(EOSmethod eos_method, QHAmethod qha_method)
  {
    string function = "QHAN::writeThermalProperties():";
    string msg = "";

    // type of qha calculation
    string qha = "";
    switch(qha_method){
      case (QHA3P_CALC):
        qha = "qha3p";
        break;
      case(QHA_CALC):
        qha = "qha";
        break;
      default:
        msg = "Nonexistent QHA method was passed to " + function;
        throw aurostd::xerror(_AFLOW_FILE_NAME_, QHA_ARUN_MODE, msg, _INPUT_UNKNOWN_);
        break;
    }

    // the name of the output file depends on the EOS fit method and on the type of
    // QHA calculation
    stringstream file;
    string filename = DEFAULT_QHA_FILE_PREFIX;
    switch(eos_method){
      case(EOS_POLYNOMIAL):
        filename += qha + ".polynomial.";
        break;
      case(EOS_BIRCH_MURNAGHAN):
        filename += qha + ".birch-murnaghan.";
        break;
      case(EOS_MURNAGHAN):
        filename += qha + ".murnaghan.";
        break;
      default:
        msg = "Nonexistent EOS method was passed to " + function;
        throw aurostd::xerror(_AFLOW_FILE_NAME_, QHA_ARUN_MODE, msg, _INPUT_UNKNOWN_);
        break;
    }
    filename += DEFAULT_QHA_THERMO_FILE;

    msg = "Writing T-dependent properties to "+filename;
    pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *p_FileMESSAGE, *p_oss,
        _LOGGER_MESSAGE_);

    file.precision(10);

    file << AFLOWIN_SEPARATION_LINE << std::endl;
    file << "[QHA_THERMAL_PROPERTIES]START" << std::endl;
    // write header
    file << setw(5)  << "# T[K]"               << setw(SW) << ' ' <<
      setw(TW) << "Veq[ev/atom]"         << setw(SW) << ' ' <<
      setw(TW) << "F(V0)[eV/atom]"       << setw(SW) << ' ' <<
      setw(TW) << "B[GPa]"               << setw(SW) << ' ' <<
      setw(TW) << "beta[10^-6/K]"        << setw(SW) << ' ' <<
      setw(TW) << "Cv[kB/atoms]"         << setw(SW) << ' ' <<
      setw(TW) << "Cp[kB/atoms]"         << setw(SW) << ' ' <<
      setw(TW) << "gamma(beta,B,Cv)"     << setw(SW) << ' ' <<
      setw(TW) << "Bprime";
    // the following properties are calculated only with a regular QHA calculation
    if (qha_method==QHA_CALC){
      file  << setw(SW) << ' ' <<
        setw(TW) << "beta_mesh[10^-6/K]"   << setw(SW) << ' ' <<
        setw(TW) << "Cv_mesh[kB/atoms]"    << setw(SW) << ' ' <<
        setw(TW) << "Cp_mesh[kB/atoms]"    << setw(SW) << ' ' <<
        setw(TW) << "gamma_mesh";
    }
    file << std::endl;

    xvector<double> F(N_EOSvolumes); // free energy
    xvector<double> xvolumes = aurostd::vector2xvector(EOSvolumes);

    switch(qha_method){
      case(QHA3P_CALC):
        for (int Vid=0; Vid<N_EOSvolumes; Vid++) F[Vid+1] = FreeEnergyTaylorExpansion(0, Vid);
        break;
      case(QHA_CALC):
        for (int Vid=0; Vid<N_EOSvolumes; Vid++) F[Vid+1] = FreeEnergy(0, Vid);
        break;
      default:
        msg = "Nonexistent QHA method was passed to " + function;
        throw aurostd::xerror(_AFLOW_FILE_NAME_, QHA_ARUN_MODE, msg, _INPUT_UNKNOWN_);
        break;
    }

    fitToEOSmodel(F, eos_method);
    double V0K = EOS_volume_at_equilibrium; // equilibrium volume at 0K

    double T = 0.0, Veq = 0.0, Feq = 0.0, B = 0.0, Bp = 0.0, beta = 0.0, CV = 0.0, CP = 0.0, GP = 0.0;
    double CV_mesh = 0.0, GP_mesh = 0.0, CP_mesh = 0.0, beta_mesh = 0.0; // these properties are calculated 
    // by weighted sum over q-points mesh
    for (int Tid=0; Tid<Ntemperatures; Tid++){
      T = Temperatures[Tid];

      switch(qha_method){
        case(QHA3P_CALC):
          for (int Vid=0; Vid<N_EOSvolumes; Vid++) F[Vid+1] = FreeEnergyTaylorExpansion(T, Vid);
          break;
        case(QHA_CALC):
          for (int Vid=0; Vid<N_EOSvolumes; Vid++) F[Vid+1] = FreeEnergy(T, Vid);
          break;
        default:
          msg = "Nonexistent QHA method was passed to " + function;
          throw aurostd::xerror(_AFLOW_FILE_NAME_, QHA_ARUN_MODE, msg, _INPUT_UNKNOWN_);
          break;
      }

      if (includeElectronicContribution) F += electronicFreeEnergySommerfeld(T);

      // stop if energy minimum is no longer within a given set of volumes
      if (!isMinimumWithinBounds(F)){
        msg = "Calculation is stopped at T=" + aurostd::utype2string<double>(T) + " [K]";
        msg+= "since there is no free energy minimum within a given volume range.";
        pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *p_FileMESSAGE,
            *p_oss, _LOGGER_MESSAGE_);
        break;
      }

      fitToEOSmodel(F, eos_method);
      // note that the state of the EOS fit is changed by the ThermalExpansion and/or
      // IsochoricSpecificHeat functions, so save Veq, Feq and B for future use
      Veq = EOS_volume_at_equilibrium;
      Feq = EOS_energy_at_equilibrium;
      B   = EOS_bulk_modulus_at_equilibrium;  // [GPa]
      Bp  = EOS_Bprime_at_equilibrium;
      beta = ThermalExpansion(T, eos_method, qha_method); // [K^-1]
      CV   = IsochoricSpecificHeat(T, V0K, eos_method, qha_method)/KBOLTZEV; // [kB/atom]
      CP   = CV + Veq*T*B*pow(beta,2)/eV2GPa/KBOLTZEV; // [kB/atom]
      GP   = (beta/CV)*B*Veq/eV2GPa/KBOLTZEV;

      // the following properties are calculated only with a regular QHA calculation
      if (qha_method==QHA_CALC){
        calcCVandGPfit(T, Veq, CV_mesh, GP_mesh);
        beta_mesh = KBOLTZEV*CV_mesh*GP_mesh/Veq/(B/eV2GPa); // [K^-1]
        CP_mesh   = CV_mesh + Veq*T*B*pow(beta_mesh,2)/eV2GPa/KBOLTZEV; // [kB/atom]
      }

      // write values to file
      file << setw(5)  << T                   << setw(SW) << ' ' <<
        setw(TW) << Veq                 << setw(SW) << ' ' <<
        setw(TW) << Feq                 << setw(SW) << ' ' <<
        setw(TW) << B                   << setw(SW) << ' ' <<
        setw(TW) << beta * 1e6          << setw(SW) << ' ' << //[10^-6/K]
        setw(TW) << CV                  << setw(SW) << ' ' <<
        setw(TW) << CP                  << setw(SW) << ' ' <<
        setw(TW) << GP                  << setw(SW) << ' ' <<
        setw(TW) << Bp;
      // the following properties are calculated only with a regular QHA calculation
      if (qha_method==QHA_CALC){
        file << setw(SW) << ' ' <<
          setw(TW) << beta_mesh * 1e6     << setw(SW) << ' ' << //[10^-6/K]
          setw(TW) << CV_mesh             << setw(SW) << ' ' <<
          setw(TW) << CP_mesh             << setw(SW) << ' ' <<
          setw(TW) << GP_mesh             << setw(SW);
      }
      file << std::endl;
    }

    file << "[QHA_THERMAL_PROPERTIES]STOP" << std::endl;
    file << AFLOWIN_SEPARATION_LINE << std::endl;

    if (!aurostd::stringstream2file(file, filename)){
      msg = "Error writing to " + filename + "file.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function,msg,_FILE_ERROR_);
    }
  }

  /// Writes the F(V,T) data to aflow.qha.FVT.out file
  ///
  void QHAN::writeFVT()
  {
    string function = "QHA::writeFVT():";
    string msg = "Writing F(V,T) relations to file.";
    pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *p_FileMESSAGE, *p_oss,
        _LOGGER_MESSAGE_);

    stringstream file;
    string filename = DEFAULT_QHA_FILE_PREFIX + DEFAULT_QHA_FVT_FILE;
    file.precision(10);

    file << AFLOWIN_SEPARATION_LINE << std::endl;
    file << "[QHA_FVT]START" << std::endl;

    xvector<double> Felec(N_EOSvolumes);
    double T = 0.0;

    for (int Tid = 0; Tid < Ntemperatures; Tid++){
      for (int Vid = 0; Vid < N_EOSvolumes; Vid++){
        T = Temperatures[Tid];
        if (includeElectronicContribution) Felec = electronicFreeEnergySommerfeld(T);
        file << setw(TW) << EOSvolumes[Vid]    << setw(SW) << ' '
          << setw(TW) << FreeEnergy(T, Vid) << setw(SW) << ' '
          << setw(TW) << Felec[Vid+1]       << setw(SW) << ' '
          << setw(TW) << E0_V[Vid] <<
          std::endl;
      }
      file << std::endl << std::endl;
    }
    file << "[QHA_FVT]STOP" << std::endl;
    file << AFLOWIN_SEPARATION_LINE << std::endl;

    if (!aurostd::stringstream2file(file, filename)){
      msg = "Error writing to " + filename + "file.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function,msg,_FILE_ERROR_);
    }
  }

  /// Writes mode-dependent Grueneisen parameter along a given path in k-space.
  ///
  void QHAN::writeGPpath(double V, const string &directory)
  {
    string function = "QHAN::writeGPpath():";
    string msg = "Calculating Grueneisen parameters along the path in";
    msg += " that was used to for the phonon dispersion.";
    pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *p_FileMESSAGE, *p_oss,
        _LOGGER_MESSAGE_);
    // we will save bands-projected Grueneisen parameter in xEIGENVAL
    xEIGENVAL GPpath(gp_ph_dispersions.front());
    GPpath.Vol = V;
    GPpath.temperature = 0;

    xstructure struc = origStructure;
    struc.SetVolume(V);
    xvector<double> lattice(3);
    lattice[1] = struc.a * 1E-10;
    lattice[2] = struc.b * 1E-10;
    lattice[3] = struc.c * 1E-10;
    GPpath.lattice = lattice;
    GPpath.carstring = "GRUEN_PATH";

    xvector<double> xomega(N_GPvolumes);
    //venergy.at(kpoint number).at(band number).at(spin number)
    for (uint q=0; q<GPpath.venergy.size(); q++){
      for (int branch=0; branch<Nbranches; branch++){
        for (int Vid=0; Vid<N_GPvolumes; Vid++){
          xomega[Vid+1] = gp_ph_dispersions[Vid].venergy.at(q).at(branch).at(0);
        }

        //GPpath.venergy.at(q).at(branch).at(0) = calcGrueneisen(V, xomega);
        GPpath.venergy.at(q).at(branch).at(0) = calcGrueneisenFD(xomega);
      }
    }

    stringstream eigenval;
    eigenval << GPpath;

    string filename = directory+'/'+DEFAULT_QHA_FILE_PREFIX+DEFAULT_QHA_GP_PATH_FILE;
    aurostd::stringstream2file(eigenval, filename);
    if (!aurostd::FileExist(filename)){
      msg = "Cannot open "+filename+" file.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function,msg,_FILE_ERROR_);
    }
  }


  /// Writes the average Grueneisen parameter, which is calculated using the finite
  /// difference method and isochoric specific heat.
  /// Output file: "aflow.qha.gp.avg.out"
  ///
  void QHAN::writeAverageGPfiniteDifferences()
  {
    string function = "QHAN::writeAverageGPfiniteDifferences(): ";
    string msg = "Writing T-dependence of the average Grueneisen parameter.";
    pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *p_FileMESSAGE, *p_oss,
        _LOGGER_MESSAGE_);

    stringstream file;
    string filename = DEFAULT_QHA_FILE_PREFIX + DEFAULT_QHA_GP_AVG_FILE;
    file.precision(10);

    file << AFLOWIN_SEPARATION_LINE << std::endl;
    file << "[QHA_AVG_GP]START" << std::endl;

    file << setw(5)  << "# T[K]"           << setw(SW) << ' ' <<
      setw(TW) << "gamma"            << setw(SW) << ' ' <<
      setw(TW) << "CV"               << setw(SW) << ' ' <<
      std::endl;

    double T = 0.0, CV =0.0, GP = 0.0;
    for (int Tid=0; Tid<Ntemperatures; Tid++){
      T = Temperatures[Tid];

      calcCVandGP(T, CV, GP);

      file << setw(5)  << T                   << setw(SW) << ' ' <<
        setw(TW) << GP                  << setw(SW) << ' ' <<
        setw(TW) << CV                  << setw(SW) << ' ' <<
        std::endl;
    }

    file << "[QHA_AVG_GP]STOP" << std::endl;
    file << AFLOWIN_SEPARATION_LINE << std::endl;

    if (!aurostd::stringstream2file(file, filename)){
      msg = "Error writing to " + filename + "file.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function,msg,_FILE_ERROR_);
    }
  }

  /// Writes the mode-dependent Grueneisen parameters calculated at each q-point in the IBZ.
  ///
  void QHAN::writeGPmeshFD()
  {
    string function = "QHAN::writeGPmeshFD():";
    string msg = "Writing Grueneisen parameters calculated on the mesh of q-points.";
    pflow::logger(QHA_ARUN_MODE, function, msg, currentDirectory, *p_FileMESSAGE, *p_oss,
        _LOGGER_MESSAGE_);

    stringstream file;
    string filename = DEFAULT_QHA_FILE_PREFIX + DEFAULT_QHA_GP_MESH_FILE;
    file.precision(6);

    xvector<double> xomega;
    uint NIrQpoints = omegaV_mesh.size();

    xEIGENVAL xeigen;
    xeigen.number_atoms = NatomsOrigCell;
    xeigen.number_loops = 0;
    xeigen.spin = 1;
    xeigen.Vol = origStructure.GetVolume()/NatomsOrigCell;
    xvector<double> lattice(3);
    lattice[1] = origStructure.a * 1e-10;
    lattice[2] = origStructure.b * 1e-10;
    lattice[3] = origStructure.c * 1e-10;
    xeigen.lattice = lattice;
    xeigen.POTIM = 0;
    xeigen.temperature = 0;
    xeigen.carstring = "GRUEN_MESH";
    xeigen.title = system_title;
    xeigen.number_electrons = 0;
    xeigen.number_kpoints = NIrQpoints;
    xeigen.number_bands = Nbranches;

    xeigen.vweight.resize(NIrQpoints);
    xeigen.vkpoint.resize(NIrQpoints, xvector<double>(3));
    xeigen.venergy.resize(NIrQpoints, deque<deque<double> >(xeigen.number_bands,
          deque<double>(2)));

    for (uint q=0; q<NIrQpoints; q++){
      xeigen.vkpoint[q] = qPoints[q];
      xeigen.vweight[q] = qpWeights[q];
      for (int branch=0; branch<Nbranches; branch++){
        xomega = aurostd::vector2xvector(omegaV_mesh[q][branch]);
        xeigen.venergy[q][branch][0] = xomega[2];
        xeigen.venergy[q][branch][1] = calcGrueneisenFD(xomega);
      }
    }

    file << xeigen;
    if (!aurostd::stringstream2file(file, filename)){
      msg = "Error writing to " + filename + "file.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function,msg,_FILE_ERROR_);
    }
  }

  /// Writes the volume-dependent phonon frequencies obtained from the series of EOS APL
  /// calculations.
  ///
  void QHAN::writeFrequencies()
  {
    string function = "writeFrequencies():", msg = "";
    stringstream file;
    string filename = DEFAULT_QHA_FILE_PREFIX + DEFAULT_QHA_FREQS_FILE;

    file << AFLOWIN_SEPARATION_LINE << std::endl;
    file << "[QHA_FREQUENCIES]START" << std::endl;
    for (uint q=0; q<omegaV_mesh_EOS.size(); q++){
      for (int branch=0; branch<Nbranches; branch++){
        for (int Vid=0; Vid<N_EOSvolumes; Vid++){
          file << EOSvolumes[Vid] << " " << omegaV_mesh_EOS[q][branch][Vid] << std::endl;
        }
        file << std::endl << std::endl;
      }
    }
    file << "[QHA_FREQUENCIES]STOP" << std::endl;
    file << AFLOWIN_SEPARATION_LINE << std::endl;

    if (!aurostd::stringstream2file(file, filename)){
      msg = "Error writing to " + filename + "file.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function,msg,_FILE_ERROR_);
    }
  }
}
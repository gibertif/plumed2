/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "ActionRegister.h"
#include "Function.h"
#include "tools/Exception.h"
#include "tools/Communicator.h"
#include "tools/BiasRepresentation.h"
#include "tools/File.h"
#include "tools/Tools.h"
#include "tools/Stopwatch.h"
#include "tools/Matrix.h"
#include "tools/Grid.h"
#include <iostream>
#include <memory>
#include <numeric>

#define DP2CUTOFF 6.25

using namespace std;

namespace PLMD {
namespace function {


//+PLUMEDOC FUNCTION FUNCSUMHILLS
/*
Insert description of the Iterative c(t) function
*/
//+ENDPLUMEDOC
/*
class FilesHandlerCT {
  vector <string> filenames;
  Action *action;
  Log *log;
  bool parallelread;
  unsigned beingread;
  bool isopen;
};
*/

double  mylogCT( double v1 ) {
  return log(v1);
}

class FuncIterativeCT :
  public Function
{
  vector<string> hillsFiles,histoFiles;
  vector <std::unique_ptr<IFile>>  ifiles;
  vector<string> proj;
  vector<string> names;
  bool iscltool,integratehills,parallelread;
  double beta;
  string outct,fmt;
  vector<std::vector<double> > c_t;
  int n_iteration;
  unsigned beingread=0;
  bool isopen=false;
public:
  explicit FuncIterativeCT(const ActionOptions&);
  void calculate(); // this probably is not needed
  bool checkFilesAreExisting(const vector<string> & hills_files );
  static void registerKeywords(Keywords& keys);
  int n_eval;
  int stride;
  vector<double> inst_pot;
  struct Gaussian {
    vector<double> center;
    vector<double> sigma;
    double height;
    bool   multivariate; // this is required to discriminate the one dimensional case
    vector<double> invsigma;
    Gaussian(const vector<double> & center,const vector<double> & sigma,double height, bool multivariate ):
      center(center),sigma(sigma),height(height),multivariate(multivariate),invsigma(sigma) {
      // to avoid troubles from zero element in flexible hills
      for(unsigned i=0; i<invsigma.size(); ++i) abs(invsigma[i])>1.e-20?invsigma[i]=1.0/invsigma[i]:0.;
    }
  };
  vector<Gaussian> hills;
private:
  bool readBunch(vector<Gaussian> &hills, vector<string> &names);
  bool scanOneHill(IFile *ifile,  vector<Value> &v, vector<double> &center, vector<double>  &sigma, double &height, bool &multivariate);
  void writeToFile(OFile& ofile);
  double calculateOverlap(const vector<double>& cv, const Gaussian& g2);
};

PLUMED_REGISTER_ACTION(FuncIterativeCT,"FUNCITERATIVECT")

void FuncIterativeCT::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys);
  keys.use("ARG");
  keys.add("optional","HILLSFILES"," source file for hills creation(may be the same as HILLS)"); // this can be a vector!
  keys.add("optional","KT"," the kt factor");
  keys.add("optional","OUTCT"," output file for histogram ");
  keys.add("optional","STRIDE"," stride when you do it on the fly ");
  keys.add("optional","PREFACTOR"," prefactor for the initial bias (e.g. the height of the first hill in Well-Tempered Metadynamics) ");
  keys.addFlag("ISCLTOOL",true,"use via plumed command line: calculate at read phase and then go");
  keys.addFlag("PARALLELREAD",false,"read parallel HILLS file");
  keys.add("optional","FMT","the format that should be used to output real numbers");
}

FuncIterativeCT::FuncIterativeCT(const ActionOptions&ao):
  Action(ao),
  Function(ao),
  iscltool(false),
  integratehills(true),
  parallelread(false),
  beta(-1.),
  fmt("%14.9f")
{
  // format
  parse("FMT",fmt);
  log<<"  Output format is "<<fmt<<"\n";
  // hills file:
  parseVector("HILLSFILES",hillsFiles);
  if(hillsFiles.size()==0) {
    integratehills=false; // default behaviour
  } else {
    integratehills=true;
    for(unsigned i=0; i<hillsFiles.size(); i++) log<<"  hillsfile  : "<<hillsFiles[i]<<"\n";
  }
  // is a cltool: then you start and then die
  parse("KT",beta);
  plumed_massert(beta>0.,"KT cannot be 0!");
  beta=1./beta;
  log<<"  beta is "<<beta<<"\n";
  parseFlag("ISCLTOOL",iscltool);
  //
  parseFlag("PARALLELREAD",parallelread);
  if(parallelread) {
    plumed_merror("parallelread is not yet implemented !!!");
  }

  if(iscltool) {
    std::vector<Value*> tmphillsvalues, tmphistovalues;
      for(unsigned i=0; i<getNumberOfArguments(); i++) {
        // allocate a new value from the old one: no deriv here
        // if we are summing hills then all the arguments are needed
        tmphillsvalues.push_back( getPntrToArgument(i) );
      }
    // check if the files exists
      for (int i=0;i<tmphillsvalues.size();i++){
        names.push_back(tmphillsvalues[i]->getName());
      }
      checkFilesAreExisting(hillsFiles);
    }

    string outct = "c_t.dat";
    parse("OUTCT",outct);
    log<<"  output file for c(t) is :  "<<outct<<"\n";
    stride = 100;
    parse("STRIDE",stride);
    double prefactor = 1.0;
    parse("PREFACTOR",prefactor);
    // Check all the keyword
    checkRead();

// Stopwatch is logged when it goes out of scope
    Stopwatch sw(log);
// Stopwatch is stopped when swh goes out of scope
    auto swh=sw.startStop("0 Summing hills");
    // read a number of hills
    int nfiles=0;
    bool ibias=integratehills;

    for(unsigned i=0; i<hillsFiles.size(); i++) {
      std::unique_ptr<IFile> ifile(new IFile());
      ifile->link(*this);
      plumed_massert((ifile->FileExist(hillsFiles[i])), "the file "+hillsFiles[i]+" does not exist " );
      ifiles.emplace_back(std::move(ifile));
    }


    while(true) {
      if(  integratehills  && ibias  ) {
        log<<"  reading hills: \n";
        ibias=readBunch(hills,names) ; log<<"\n";
      }

      if(!ibias)integratehills=false;// once you get to the final bunch just give up

      if ( !ibias) break; //when both are over then just quit

      nfiles++;
    }

    log<<"  using a stride to evaluate c(t) of :  "<<stride<<"\n";

    n_eval = int(hills.size()/stride);

    log<< "  corresponding to  " << n_eval << " number of evaluation\n";

    // dump: need to project?
    std::ostringstream ostr; ostr<<nfiles;
    log<<"\n";
    log<<"  Now calculating...\n";
    log<<"\n";

    log<<"\n";
    log<<"  Evaluating the Potential Matrix V(s(t),t').\n";
    log<<"\n";

    Matrix<double> pot_matrix(n_eval,n_eval);
    inst_pot.resize(n_eval);

    // this should go in a function for the matrix multiplication
    prefactor /= hills[0].height;

    for (int index_1=0 ; index_1<n_eval ; index_1++){
      int time_1 = index_1*stride;
      if (index_1 % 100 == 0) {
        log << "  Done " << index_1 << " iterations over " <<  n_eval << "... \n";
      }
      double alpha = 0.0;
      for (int index_2=0 ; index_2<=index_1 ; index_2++){
        int time_2 =index_2*stride;

        double V_st_tp = 0.0;
        for (int time=0 ; time<time_1 ; time++){
          V_st_tp += calculateOverlap(hills[time_2].center,hills[time]);
        }
        alpha = exp(-beta*V_st_tp*prefactor) ;
        pot_matrix[index_1][index_2] = alpha ;
      }
      inst_pot[index_1] = 1.0/alpha ;
    }

    log << "  Calculation of V(s(t),t') is finished! \n";
    log << "\n" ;

    n_iteration = 20 ;
    c_t.resize(n_iteration);
    for (size_t it = 0; it < n_iteration; it++) {
      c_t[it].resize(n_eval);
      std::fill(c_t[it].begin(),c_t[it].end(),1.0);
    }

    vector<double> denominator(n_eval);
    vector<double> numerator(n_eval);
    vector<double> nom(n_eval);

    log << "  Iteratively calculating the weight! \n";

    for (size_t it = 1; it < n_iteration; it++) {
      log << "  Doing iteration number: " << it << " out of " << n_iteration << "\n";
      denominator[0] = c_t[it-1][0] * inst_pot[0] ;
      nom[0] = denominator[0];
      for (size_t index_1 = 1; index_1 < n_eval; index_1++) {
        denominator[index_1] = c_t[it-1][index_1] * inst_pot[index_1] ;
        nom[index_1] = nom[index_1-1] + denominator[index_1] ;
      }

      for (int index_1=0 ; index_1<n_eval ; index_1++){
        numerator[index_1] = 0.0 ;
        for (int index_2=0 ; index_2<n_eval ; index_2++){
          numerator[index_1] += pot_matrix[index_1][index_2]*denominator[index_2] ;
        }
        c_t[it][index_1] =  numerator[index_1]/nom[index_1] ;
      }
    }

    log << "  Done. Printing the result on file!" << "\n";

    // print out the last CT
    OFile myout;
    myout.link(*this);
    myout.open(outct);
    writeToFile(myout);

    return;
  }


  double FuncIterativeCT::calculateOverlap(const vector<double>& cv, const Gaussian& g2)
  {
    double dp2=0.0;
    double bias=0.0;
    // I use a pointer here because cv is const (and should be const)
    // but when using doInt it is easier to locally replace cv[0] with
    // the upper/lower limit in case it is out of range
    const double *pcv=NULL; // pointer to cv
    double tmpcv[1]; // tmp array with cv (to be used with doInt_)
    if(cv.size()>0) pcv=&cv[0];

    if(g2.multivariate) {
      unsigned k=0;
      unsigned ncv=cv.size();
      // recompose the full sigma from the upper diag cholesky
      Matrix<double> mymatrix(ncv,ncv);
      for(unsigned i=0; i<ncv; i++) {
        for(unsigned j=i; j<ncv; j++) {
          mymatrix(i,j)=mymatrix(j,i)=g2.sigma[k]; // recompose the full inverse matrix
          k++;
        }
      }

      for(unsigned i=0; i<cv.size(); ++i) {
        double dp_i=difference(i,g2.center[i],pcv[i]);
        for(unsigned j=i; j<cv.size(); ++j) {
          if(i==j) {
            dp2+=dp_i*dp_i*mymatrix(i,j)*0.5;
          } else {
            double dp_j=difference(j,g2.center[j],pcv[j]);
            dp2+=dp_i*dp_j*mymatrix(i,j);
          }
        }
      }
      if(dp2<DP2CUTOFF) {
        bias=g2.height*exp(-dp2);
      }
    } else {
      for(unsigned i=0; i<cv.size(); ++i) {
        double dp=difference(i,g2.center[i],pcv[i])*g2.invsigma[i];
        dp2+=dp*dp;
      }
      dp2*=0.5;
      if(dp2<DP2CUTOFF) {
        bias=g2.height*exp(-dp2);
      }
    }
    return bias;
  }

  bool FuncIterativeCT::readBunch(vector<Gaussian> &hills, vector<string> &names) {
    bool morefiles; morefiles=true;
    // read one by one hills
    IFile *ff;
    ff=ifiles[beingread].get();

    if(!isopen) {
      log<<"  opening file "<< hillsFiles[beingread]<<"\n";
      ff->open(hillsFiles[beingread]); isopen=true;
    }

    int n=0;
    while(true) {
      bool fileisover=true;
      unsigned ncv=getNumberOfArguments();
      vector<double> center(ncv);
      vector<double> sigma(ncv);
      double height;
      int nhills=0;
      bool multivariate=false;

      std::vector<Value> tmpvalues;
      for(unsigned j=0; j<getNumberOfArguments(); ++j){
        tmpvalues.push_back( Value( this, getPntrToArgument(j)->getName(), false ) );
      }

      while(scanOneHill(ff,tmpvalues,center,sigma,height,multivariate)) {
        hills.emplace_back(Gaussian(center,sigma,height,multivariate));
        n=hills.size();
      }
      if(fileisover) {
        log<<"  closing file "<<hillsFiles[beingread]<<"\n";
        ff->close();
        isopen=false;
        log<<"  now total "<<hills.size()<<" kernels \n";
        beingread++;
        if(beingread<ifiles.size()) {
          ff=ifiles[beingread].get(); ff->open(hillsFiles[beingread]);
          log<<"  opening file "<<hillsFiles[beingread]<<"\n";
          isopen=true;
        } else {
          morefiles=false;
          log<<"  final chunk: now with "<<n<<" kernels  \n";
          break;
        }
      }
      // if there are no more files to read and this file is over then quit
      if(fileisover && !morefiles) {break;}
      // if you are in the middle of a file and you are here
      // then means that you read what you need to read
      if(!fileisover ) {break;}
    }
    return morefiles;
  }


  bool FuncIterativeCT::scanOneHill(IFile *ifile,  vector<Value> &tmpvalues, vector<double> &center, vector<double>  &sigma, double &height, bool &multivariate)
  {
    double dummy;
    multivariate=false;
    if(ifile->scanField("time",dummy)) {
      unsigned ncv; ncv=tmpvalues.size();
      for(unsigned i=0; i<ncv; ++i) {
        ifile->scanField( &tmpvalues[i] );
        if( tmpvalues[i].isPeriodic() && ! getPntrToArgument(i)->isPeriodic() ) {
          error("in hills file periodicity for variable " + tmpvalues[i].getName() + " does not match periodicity in input");
        } else if( tmpvalues[i].isPeriodic() ) {
          std::string imin, imax; tmpvalues[i].getDomain( imin, imax );
          std::string rmin, rmax; getPntrToArgument(i)->getDomain( rmin, rmax );
          if( imin!=rmin || imax!=rmax ) {
            error("in hills file periodicity for variable " + tmpvalues[i].getName() + " does not match periodicity in input");
          }
        }
        center[i]=tmpvalues[i].get();
      }
      // scan for kerneltype
      std::string ktype="gaussian";
      if( ifile->FieldExist("kerneltype") ) ifile->scanField("kerneltype",ktype);
      // scan for multivariate label: record the actual file position so to eventually rewind
      std::string sss;
      ifile->scanField("multivariate",sss);
      if(sss=="true") multivariate=true;
      else if(sss=="false") multivariate=false;
      else plumed_merror("cannot parse multivariate = "+ sss);
      if(multivariate) {
        sigma.resize(ncv*(ncv+1)/2);
        Matrix<double> upper(ncv,ncv);
        Matrix<double> lower(ncv,ncv);
        for(unsigned i=0; i<ncv; i++) {
          for(unsigned j=0; j<ncv-i; j++) {
            ifile->scanField("sigma_"+getPntrToArgument(j+i)->getName()+"_"+getPntrToArgument(j)->getName(),lower(j+i,j));
            upper(j,j+i)=lower(j+i,j);
          }
        }
        Matrix<double> mymult(ncv,ncv);
        Matrix<double> invmatrix(ncv,ncv);
        mult(lower,upper,mymult);
        // now invert and get the sigmas
        Invert(mymult,invmatrix);
        // put the sigmas in the usual order: upper diagonal (this time in normal form and not in band form)
        unsigned k=0;
        for(unsigned i=0; i<ncv; i++) {
          for(unsigned j=i; j<ncv; j++) {
            sigma[k]=invmatrix(i,j);
            k++;
          }
        }
      } else {
        for(unsigned i=0; i<ncv; ++i) {
          ifile->scanField("sigma_"+getPntrToArgument(i)->getName(),sigma[i]);
        }
      }

      ifile->scanField("height",height);
      ifile->scanField("biasf",dummy);
      if(ifile->FieldExist("clock")) ifile->scanField("clock",dummy);
      if(ifile->FieldExist("lower_int")) ifile->scanField("lower_int",dummy);
      if(ifile->FieldExist("upper_int")) ifile->scanField("upper_int",dummy);
      ifile->scanField();
      return true;
    } else {
      return false;
    }
  }

void FuncIterativeCT::calculate() {
  // this should be connected only with a grid representation to metadynamics
  // at regular time just dump it
  plumed_merror("You should have never got here: this stuff is not yet implemented!");
}

void FuncIterativeCT::writeToFile(OFile& ofile){
  for (unsigned int i=0; i<n_eval ; i++){
    ofile.printField("neval",int(i*stride));
    for (size_t l = 0; l < n_iteration; l++) {
      ofile.printField("c_t_"+std::to_string(l),-std::log(c_t[l][i])/beta);
    }
    ofile.printField("inst_pot",std::log(inst_pot[i])/beta);
    ofile.printField();
  }
  ofile.flush();
}

bool FuncIterativeCT::checkFilesAreExisting(const vector<string> & hills_files ) {
  plumed_massert(hills_files.size()!=0,"the number of  files provided should be at least one" );
  std::unique_ptr<IFile> ifile(new IFile());
  ifile->link(*this);
  for(unsigned i=0; i< hills_files.size(); i++) {
    plumed_massert(ifile->FileExist(hills_files[i]),"missing file "+hills_files[i]);
  }
  return true;

}

}
}

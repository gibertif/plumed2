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
#include "CLTool.h"
#include "CLToolRegister.h"
#include "tools/Tools.h"
#include "core/Action.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Communicator.h"
#include "tools/Random.h"
#include <cstdio>
#include <string>
#include <vector>
#include <iostream>
#include "tools/File.h"
#include "core/Value.h"
#include "tools/Matrix.h"

using namespace std;

namespace PLMD {
namespace cltools {

//+PLUMEDOC TOOLS sum_hills
/*
Write doc for Iterative c(t)
*/
//+ENDPLUMEDOC

class CLToolIterativeCT : public CLTool {
public:
  static void registerKeywords( Keywords& keys );
  explicit CLToolIterativeCT(const CLToolOptions& co );
  int main(FILE* in,FILE*out,Communicator& pc);
  string description()const;
/// find a list of variables present, if they are periodic and which is the period
/// return false if the file does not exist
  static bool findCvsAndPeriodic(std::string filename, std::vector< std::vector <std::string> > &cvs,std::vector<std::string> &pmin,std::vector<std::string> &pmax, bool &multivariate, string &lowI_, string &uppI_);
};

void CLToolIterativeCT::registerKeywords( Keywords& keys ) {
  CLTool::registerKeywords( keys );
  keys.addFlag("--help-debug",false,"print special options that can be used to create regtests");
  keys.add("optional","--hills","specify the name of the hills file");
  keys.add("optional","--stride","specify the stride for integrating hills file (default 100)");
  keys.add("optional","--outfile","specify the output file for iterative_c_t");
  keys.add("optional","--kt","specify temperature in energy units for integrating out variables");
  keys.add("optional","--prefactor","specify the height of the hill at the beginning of the calculation");
  keys.add("optional","--fmt","specify the output format");
}

CLToolIterativeCT::CLToolIterativeCT(const CLToolOptions& co ):
  CLTool(co)
{
  inputdata=commandline;
}

string CLToolIterativeCT::description()const { return " calculate the iterative c(t) estimator."; }

int CLToolIterativeCT::main(FILE* in,FILE*out,Communicator& pc) {

// Read the hills input file name
  vector<string> hillsFiles;
  bool dohills;
  dohills=parseVector("--hills",hillsFiles);

  plumed_massert(dohills,"you should use --hills command");

  vector< vector<string> > vcvs;
  vector<string> vpmin;
  vector<string> vpmax;
  string lowI_, uppI_;
// parse it as it was a restart
  bool vmultivariate;
  findCvsAndPeriodic(hillsFiles[0], vcvs, vpmin, vpmax, vmultivariate, lowI_, uppI_);

  vector<std::string> sigma;

  // now put into a neutral vector

  vector< vector<string> > cvs;
  vector<string> pmin;
  vector<string> pmax;

  cvs=vcvs;
  pmin=vpmin;
  pmax=vpmax;

  PlumedMain plumed;
  std::string ss;
  unsigned nn=1;
  ss="setNatoms";
  plumed.cmd(ss,&nn);
  if(Communicator::initialized())  plumed.cmd("setMPIComm",&pc.Get_comm());
  plumed.cmd("init",&nn);
  vector <bool> isdone(cvs.size(),false);
  for(unsigned i=0; i<cvs.size(); i++) {
    if(!isdone[i]) {
      isdone[i]=true;
      std::vector<std::string> actioninput;
      std::vector <unsigned> inds;
      actioninput.push_back("FAKE");
      actioninput.push_back("ATOMS=1");
      actioninput.push_back("LABEL="+cvs[i][0]);
      std::vector<std::string> comps, periods;
      if(cvs[i].size()>1) {comps.push_back(cvs[i][1]); inds.push_back(i);}
      periods.push_back(pmin[i]); periods.push_back(pmax[i]);
      for(unsigned j=i+1; j<cvs.size(); j++) {
        if(cvs[i][0]==cvs[j][0] && !isdone[j]) {
          if(cvs[i].size()==1 || cvs[j].size()==1  )plumed_merror("you cannot have twice the same label and no components ");
          if(cvs[j].size()>1) {
            comps.push_back(cvs[j][1]);
            periods.push_back(pmin[j]); periods.push_back(pmax[j]);
            isdone[j]=true; inds.push_back(j);
          }
        }

      }
      // drain all the components
      std::string addme;
      if(comps.size()>0) {
        addme="COMPONENTS=";
        for(unsigned i=0; i<comps.size()-1; i++)addme+=comps[i]+",";
        addme+=comps.back();
        actioninput.push_back(addme);
      }
      // periodicity (always explicit here)
      addme="PERIODIC=";
      for(unsigned j=0; j<periods.size()-1; j++) {
        addme+=periods[j]+",";
      }
      addme+=periods.back();
      actioninput.push_back(addme);
      plumed.readInputWords(actioninput);
    }
  }

  unsigned ncv=cvs.size();
  std::vector<std::string> actioninput;

  std::string kt; kt=std::string("1.");// assign an arbitrary value just in case that idw.size()==cvs.size()

  plumed_massert(parse("--kt",kt)," You need to define --kt ");

  std::string addme;

  actioninput.push_back("FUNCITERATIVECT");
  actioninput.push_back("ISCLTOOL");

  // set names
  std::string outfile;
  if(parse("--outfile",outfile)) {
    actioninput.push_back("OUTCT="+outfile);
  }else{
    actioninput.push_back("OUTCT=c_t.dat");
  }
  std::string prefactor;
  if(parse("--prefactor",prefactor)) {
    actioninput.push_back("PREFACTOR="+prefactor);
  }else{
    actioninput.push_back("PREFACTOR=1.0");
  }

  addme="ARG=";
  for(unsigned i=0; i<(ncv-1); i++) {
    if(cvs[i].size()==1) {
      addme+=std::string(cvs[i][0])+",";
    } else {
      addme+=std::string(cvs[i][0])+"."+std::string(cvs[i][1])+",";
    }
  }
  if(cvs[ncv-1].size()==1) {
    addme+=std::string(cvs[ncv-1][0]);
  } else {
    addme+=std::string(cvs[ncv-1][0])+"."+std::string(cvs[ncv-1][1]);
  }
  actioninput.push_back(addme);
  //for(unsigned i=0;i< actioninput.size();i++){
  //  cerr<<"AA "<<actioninput[i]<<endl;
  //}
  addme="HILLSFILES="; for(unsigned i=0; i<hillsFiles.size()-1; i++)addme+=hillsFiles[i]+","; addme+=hillsFiles[hillsFiles.size()-1];
  actioninput.push_back(addme);
  std::string  stride; stride="";

  if(parse("--stride",stride)) {
    actioninput.push_back("STRIDE="+stride);
  }else{
    actioninput.push_back("STRIDE=100");
  }

  actioninput.push_back("KT="+kt);

  std::string fmt; fmt="";
  parse("--fmt",fmt);
  if(fmt!="")actioninput.push_back("FMT="+fmt);


//  for(unsigned i=0;i< actioninput.size();i++){
//   cerr<<"AA "<<actioninput[i]<<endl;
//  }
  plumed.readInputWords(actioninput);
  // if not a grid, then set it up automatically
  return 0;
}

bool CLToolIterativeCT::findCvsAndPeriodic(std::string filename, std::vector< std::vector<std::string>  > &cvs, std::vector<std::string> &pmin,std::vector<std::string> &pmax, bool &multivariate, string &lowI_, string &uppI_) {
  IFile ifile;
  ifile.allowIgnoredFields();
  std::vector<std::string> fields;
  if(ifile.FileExist(filename)) {
    cvs.clear(); pmin.clear(); pmax.clear();
    ifile.open(filename);
    ifile.scanFieldList(fields);
    bool before_sigma=true;
    for(unsigned i=0; i<fields.size(); i++) {
      size_t pos = 0;
      size_t founds,foundm,foundp;
      //found=(fields[i].find("sigma_", pos) || fields[i].find("min_", pos) || fields[i].find("max_", pos) ) ;
      founds=fields[i].find("sigma_", pos)  ;
      foundm=fields[i].find("min_", pos)  ;
      foundp=fields[i].find("max_", pos)  ;
      if (founds!=std::string::npos || foundm!=std::string::npos ||  foundp!=std::string::npos )before_sigma=false;
      // cvs are after time and before sigmas
      size_t  found;
      found=fields[i].find("time", pos);
      if( found==std::string::npos && before_sigma) {
        // separate the components
        size_t dot=fields[i].find_first_of('.');
        std::vector<std::string> ss;
        // this loop does not take into account repetitions
        if(dot!=std::string::npos) {
          std::string a=fields[i].substr(0,dot);
          std::string name=fields[i].substr(dot+1);
          ss.push_back(a);
          ss.push_back(name);
          cvs.push_back(ss);
        } else {
          std::vector<std::string> ss;
          ss.push_back(fields[i]);
          cvs.push_back(ss);
        }
        //std::cerr<<"found variable number  "<<cvs.size()<<" :  "<<cvs.back()[0]<<std::endl;
        //if((cvs.back()).size()!=1){
        //	std::cerr<<"component    "<<(cvs.back()).back()<<std::endl;
        //}
        // get periodicity
        pmin.push_back("none");
        pmax.push_back("none");
        std::string mm; if((cvs.back()).size()>1) {mm=cvs.back()[0]+"."+cvs.back()[1];} else {mm=cvs.back()[0];}
        if(ifile.FieldExist("min_"+mm)) {
          std::string val;
          ifile.scanField("min_"+mm,val);
          pmin[pmin.size()-1]=val;
          // std::cerr<<"found min   :  "<<pmin.back()<<std::endl;
        }
        //std::cerr<<"found min   :  "<<pmin.back()<<std::endl;
        if(ifile.FieldExist("max_"+mm)) {
          std::string val;
          ifile.scanField("max_"+mm,val);
          pmax[pmax.size()-1]=val;
          // std::cerr<<"found max   :  "<<pmax.back()<<std::endl;
        }
        //std::cerr<<"found max   :  "<<pmax.back()<<std::endl;
      }
    }
    // is multivariate ???
    std::string sss;
    multivariate=false;
    if(ifile.FieldExist("multivariate")) {
      ;
      ifile.scanField("multivariate",sss);
      if(sss=="true") { multivariate=true;}
      else if(sss=="false") { multivariate=false;}
    }
    // do interval?
    if(ifile.FieldExist("lower_int")) {
      ifile.scanField("lower_int",lowI_);
      ifile.scanField("upper_int",uppI_);
    } else {
      lowI_="-1.";
      uppI_="-1.";
    }
    ifile.scanField();
    return true;
  } else {
    return false;
  }
}


PLUMED_REGISTER_CLTOOL(CLToolIterativeCT,"iterative_c_t")
}
}

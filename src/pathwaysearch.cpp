#include<boost/property_tree/ptree.hpp>
#include<boost/property_tree/xml_parser.hpp>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include "pathwaysearch.hpp"
#include "db_manager.hpp"
#include "mysqlrepository.hpp"
#include "ramrepository.hpp"
//#include "folate_dbclient.hpp"
//#include "mygo_dbclient.hpp"
#include "utility.hpp"

#include<cstring>
//short unsigned settings->max_order = 2;


//int MAX_ADDABLE=1000;

struct byDouble{
    bool operator()(const double & i1, const double  & i2) const{
        return (i1<i2);
    }
};

inline bool getbool(const string & str){
  return !str.compare("true");
}

inline bool filehandle(ifstream & ifs, const char * filename){
  ifs.open(filename);
  if (!ifs.is_open()){
    cerr<<"File "<<filename<<" not found. ";
    return false;
  }
  return true;
}


template <class Type>
void loadvector(ifstream & ifs, Type * vec, uint rows){
  string line;
  for(uint i=0;i<rows;++i){
    getline(ifs,line);
    istringstream iss(line);
    iss>>vec[i];
  }
}

void loadmatrix(ifstream & ifs, double ** matrix, uint rows, uint cols,bool transpose){
  string line;
  for(uint i=0;i<rows;++i){
    getline(ifs,line);
    istringstream iss(line);
    for(uint j=0;j<cols;++j){
      if (transpose){
        iss>>matrix[j][i];
      }else{
        iss>>matrix[i][j];
      }
    }
  }
}

int colcount(const char * filename){
  int c=0;
  ifstream ifs(filename);
  if (!ifs.is_open()){
     cerr<<filename<<" not found.\n";
  }else{
    string line,token;
    getline(ifs,line);
    istringstream iss(line);
    while(iss>>token) ++c;
    cerr<<"Found "<<c<<" cols in "<<filename<<endl;
  }
  return c;
}

int rowcount(const char * filename){
  ifstream ifs(filename);
  if (!ifs.is_open()){
     cerr<<filename<<" not found.\n";
     throw "File input error."; 
  }
  int c=0;
  string line;
  while(getline(ifs,line)){
    ++c;
  }
  ifs.close();
  return c;
}


MCMCSampler::MCMCSampler(pathway_settings_t * settings){
  this->settings = settings;
  cerr<<"Allocating memory\n";
  iTotalStudyPersons = rowcount(settings->file_trait.data());
  iTotalSnps = rowcount(settings->file_snplist.data());
  iTotalEnvCov = colcount(settings->file_study_env_var.data());
  iTotalExposures = iTotalSnps+iTotalEnvCov;
  expGeno_obs = new double[iTotalSnps];
  expGeno_prior = new double[iTotalSnps];
  dPhenoVector = new double[iTotalStudyPersons];
  cStudyGenoMatrix = new char*[iTotalSnps];
  fStudyEnvMatrix = new double *[iTotalEnvCov];
  for (unsigned int i=0;i<iTotalSnps;++i){
     cStudyGenoMatrix[i] = new char[iTotalStudyPersons];
  }
  fStudyEnvMatrix = new double *[iTotalEnvCov];
  for (unsigned int i=0;i<iTotalEnvCov;++i){
    fStudyEnvMatrix[i] = new double[iTotalStudyPersons];
  }
  if (settings->usedb && !settings->marginal_prior){
    data = new MySqlRepository(this,settings->dbhost.data(),settings->dbuser.data(),
      settings->dbpw.data(),settings->dbname.data());
  }else{
    data = new RamRepository(settings->marginal_prior);
  }
  math = new MathUtils(time(NULL));
  if (settings->dynamic_prior){
    iTotalColInA = colcount(settings->file_prior_endopheno.data());
    iTotalColInZ = rowcount(settings->file_selected_endo.data())+1;
    iTotalPriorPersons = rowcount(settings->file_prior_endopheno.data());
    cPriorGenoMatrix = new char*[iTotalSnps];
    for (unsigned int i=0;i<iTotalSnps;++i){
      cPriorGenoMatrix[i] = new char[iTotalPriorPersons];
    }
    fPriorEnvMatrix = new double *[iTotalEnvCov];
    for (unsigned int i=0;i<iTotalEnvCov;++i){
      fPriorEnvMatrix[i] = new double[iTotalPriorPersons];
    }
    fBioMarkerMatrix = new double*[iTotalColInA];
    for(uint i=0;i<iTotalColInA;++i){
      fBioMarkerMatrix[i] = new double[iTotalPriorPersons];
    }
    selectedBiomarkersInZ = new uint[iTotalColInZ];
    //if (settings->usedb) data = new MySqlFolateClientRepository(this,math,settings->dbname.data());
  }else{
    iTotalColInA = colcount(settings->file_a_matrix.data())-1;
    iTotalColInZ = colcount(settings->file_z_matrix.data());
    dStaticZ = new double * [iTotalExposures];
    for(uint i=0;i<iTotalExposures;++i){
      dStaticZ[i] = new double[iTotalColInZ];
    }
    dStaticA = new double * [iTotalExposures];
    for(uint i=0;i<iTotalExposures;++i){
      dStaticA[i] = new double[iTotalColInA];
    }
    
    //if (settings->usedb) data = new MySqlMyGoClientRepository(this,math,dStaticA);
  }
  termWeights = new double[iTotalColInA];
  //if (settings->usedb) data = new MySqlFolateClientRepository
  cerr<<"Totalcol in Z: "<<iTotalColInZ<<endl;
  cerr<<"Totalcol in A: "<<iTotalColInA<<endl;
//  if (data==NULL) data = new RamRepository(math,false);
  currentModel = new ModelParams();
  prevModel = new ModelParams();
  cerr<<"Successfully allocated memory\n";
  out_models.open("models.out");
  out_models<<"MODEL\tMAINEFFS\tSIZE\tLOG-LIKELIHOOD\n";
  out_prior_mean.open("prior_mean.out");
  out_prior_mean<<"MODEL\tINTERCEPT\tSE";
  for(uint j=0;j<iTotalColInZ-1;++j){
    out_prior_mean<<"\tCOL"<<j<<"\tSE";
  }
  out_prior_mean<<endl;
  out_prior_var.open("prior_var.out");
  out_prior_var<<"MODEL\tTAU^2\tSIGMA^2\n";
  out_beta.open("betas.out");
  out_beta<<"MODEL\tVAR_ID\tPOST_BETA\tPOST_SE\tMAINEFFS\n";
  out_prior_beta.open("prior_betas.out");
  out_prior_beta<<"MODEL\tMAINEFFS\tMLE_BETA\tSE\tZ_MEAN\tA_MEAN\tPRIOR_BETA\tSE\tADD_PROB\n";

}

MCMCSampler::~MCMCSampler(){
  cerr<<"Deleting MCMC sampler\n"; 
  cerr<<"Closing file handles\n";
  out_beta.close();
  out_prior_beta.close();
  out_prior_mean.close();
  out_prior_var.close();
  out_models.close();
  cerr<<"Cleaning up data \n";
  delete[] expGeno_obs;
  delete[] expGeno_prior;
  delete[] dPhenoVector;
  for (unsigned int i=0;i<iTotalEnvCov;++i){
    delete[] fStudyEnvMatrix[i];
  }
  delete[] fStudyEnvMatrix;
  for (unsigned int i=0;i<iTotalSnps;++i){
     delete[] cStudyGenoMatrix[i];
  }
  delete[] cStudyGenoMatrix;
  delete data;
  delete math;
  if (settings->dynamic_prior){
    for (unsigned int i=0;i<iTotalEnvCov;++i){
      delete[] fPriorEnvMatrix[i];
    }
    delete[] fPriorEnvMatrix;
    for (unsigned int i=0;i<iTotalSnps;++i){
      delete[] cPriorGenoMatrix[i];
    }
    delete[] cPriorGenoMatrix;    
    for(uint i=0;i<iTotalColInA;++i){
      delete[] fBioMarkerMatrix[i];
    }
    delete[] fBioMarkerMatrix;
    delete[] selectedBiomarkersInZ;
    //if (settings->usedb) data = new MySqlFolateClientRepository(this,math,settings->dbname.data());
  }else{
    for(uint i=0;i<iTotalExposures;++i){
      delete[] dStaticZ[i];
    }
    delete[]dStaticZ;
    for(uint i=0;i<iTotalExposures;++i){
      delete[] dStaticA[i];
    }
    delete[]dStaticA;
  }
  delete[] termWeights;
  delete currentModel;
  delete prevModel;
  cerr<<"Cleaned up data \n";
}


void MCMCSampler::init(){

  srand(time(NULL));
  dBICPenalty = sqrt(log(iTotalStudyPersons)); // standard BIC

  ifstream ifs;
  if (!filehandle(ifs,settings->file_trait.data())) throw "trait file required";
  cerr<<"Reading in "<<iTotalStudyPersons<<" phenotypes from "<<settings->file_trait<<endl;
  uint i=0;
  string line;
  while(getline(ifs,line)){
    istringstream iss(line);
    iss>>dPhenoVector[i];
    ++i;
  }
  ifs.close();
  cerr<<"Reading in genotypes for study database...\n";
  readGenoFile(settings->file_study_geno.data(),cStudyGenoMatrix,iTotalStudyPersons,iTotalSnps);
  // impute missing genotypes simply by their expected value (no LD)
  computeExpGeno(cStudyGenoMatrix,iTotalStudyPersons,expGeno_obs);
  cerr<<"Reading in "<<iTotalEnvCov<<" environmental covariates in study...\n";
  if(filehandle(ifs, settings->file_study_env_var.data())){
    loadmatrix(ifs,fStudyEnvMatrix,iTotalStudyPersons,iTotalEnvCov,true);
    ifs.close();
  }
  if (settings->dynamic_prior){
    cerr<<"Reading in genotypes for prior database...\n";
    readGenoFile(settings->file_prior_geno.data(),cPriorGenoMatrix,iTotalPriorPersons,iTotalSnps);
    computeExpGeno(cPriorGenoMatrix,iTotalPriorPersons,expGeno_prior);
    cerr<<"Reading in "<<iTotalEnvCov<<" environmental covariates in prior...\n";
    if(filehandle(ifs, settings->file_prior_env_var.data())){
      loadmatrix(ifs,fPriorEnvMatrix,iTotalPriorPersons,iTotalEnvCov,true);
      ifs.close();
    }
    cerr<<"Reading in "<<iTotalColInA<<" biomarkers...\n";
    if(filehandle(ifs, settings->file_prior_endopheno.data())){
      loadmatrix(ifs,fBioMarkerMatrix,iTotalPriorPersons,iTotalColInA,true);
      ifs.close();
    }
    cerr<<"Reading in "<<(iTotalColInZ-1)<<" selected biomarkers to be used in Z matrix...\n";
    if(filehandle(ifs, settings->file_selected_endo.data())){
      loadvector<uint>(ifs,selectedBiomarkersInZ,iTotalColInZ-1);
      ifs.close();
    }
  } else{ // END FOLATE CASE
    cerr<<"Static prior:\n";
    cerr<<"Reading in "<<iTotalColInZ<<" cols in Z...\n";
    if(filehandle(ifs, settings->file_z_matrix.data())){
      loadmatrix(ifs,dStaticZ,iTotalExposures,iTotalColInZ,false);
      ifs.close();
    }
    cerr<<"Reading in "<<iTotalColInA<<" cols in A...\n";
    if(filehandle(ifs, settings->file_a_matrix.data())){
      loadmatrix(ifs,dStaticA,iTotalExposures,iTotalColInA,false);
      ifs.close();
    }
    // set the term weihts
    for(uint i=0;i<iTotalColInA;++i){
      termWeights[i] = 0.;
    }
    cerr<<"Initializing term weights for snp:\n";
    for(uint snp=0;snp<iTotalExposures;++snp){
      //if (snp % 1000==0) cerr<<" "<<snp;
      double ind[iTotalColInA];
      for(uint i=0;i<iTotalColInA;++i){
        ind[i] = dStaticA[snp][i];
        termWeights[i] +=ind[i];
      }
    }
    termWeightsTotal = 0.;
    for(uint i=0;i<iTotalColInA;++i){
     termWeights[i] =(iTotalExposures-termWeights[i]+1.)/iTotalExposures;
     termWeightsTotal+=termWeights[i];
     //cerr<<"GO Term weight for index "<<i<<": "<<termWeights[i]<<endl;
    }
  } // end cases
  uint initSnps[settings->max_modelsize];
  cerr<<"Reading in selected effects to initialize model...\n";
  if(!filehandle(ifs, settings->file_initmodel.data())){
    throw "Must enter an initial model file.";
  }
  i=0;
  while(getline(ifs,line)){
    istringstream iss(line);
    uint snp;
    iss>>snp;
    initSnps[i++]=snp;
    if (i>=settings->max_modelsize){
      cerr<<"Max model size is "<<settings->max_modelsize<<endl;
      throw "Number of variables exceeds the maximum model size";
    }
  }
  ifs.close();
  int modeldim = i;

  currentModel->iModelSize=0;
  int mainEff[1];
  EffectPtr newEff;
  for(int i=-1;i< modeldim;++i){
    if (i==-1){
      mainEff[0] = i;
    }else{
      mainEff[0] = initSnps[i];
    }
    createNewEffect(mainEff,1,false,newEff);
    newEff->indexInModel = currentModel->iModelSize;
    currentModel->effectsInModel[currentModel->iModelSize++]=newEff;
    newEff->debugOut();
    currentModel->effectsInModelSet.insert(newEff);
    EffectPtr & intercept = currentModel->effectsInModel[0];
    if (i>=0){
      newEff->child1Id = intercept->getId();
      ++intercept->referenceCount;
    }else{
      intercept = newEff;
    }
  }
  checkDataIntegrity();
  cerr<<"Initialized the model with "<<currentModel->iModelSize<<" variables.\n";
}

void MCMCSampler::readGenoFile(const char * filename, char ** mat,uint iMaxRows,uint iMaxCols){
  string line;
  ifstream genoFile(filename,ios::in|ios::binary);
  cerr<<"Reading in "<<filename<<endl;
  cerr<<"max cols: "<<iMaxCols<<", max persons: "<<iMaxRows<<endl;
  uint linesread=0,matrixrow=0;
  if (!genoFile)
    throw "Geno file not found\n";
  if (genoFile.is_open()) {
    // check that number of iTotalSnps matches length of first row
    string testfirst;
    getline (genoFile,testfirst);
    if (testfirst.length()!=iMaxCols){
      cerr<<"Cols specified: "<<iMaxCols<<", found: "<<testfirst.length()<<endl;
      throw "Number of iTotalSnps does not match number of cols in inputfile";
    }
    genoFile.seekg(0,ios::beg);
    cerr<<"Reading geno file...\n";
    while (linesread<iMaxRows){
        char line[iMaxCols+1];
        genoFile.read(line,iMaxCols+1);
        for (uint matrixcol=0;matrixcol<iMaxCols;++matrixcol){
          char genoChar;
          switch(line[matrixcol]){
                case '0':
                    genoChar = UNDEFINED_GENO;
                    break;
                case '1':
                    genoChar = 0;
                    break;
                case '2':
                    genoChar = 1;
                    break;
                case '3':
                    genoChar = 2;
                    break;
                default:
                    cerr<<"Found char "<<genoChar<<endl;
                    throw "Genotype char not recognized";
          }  // end switch
          mat[matrixcol][matrixrow] = genoChar;
        } // end for loop
        ++matrixrow;
      linesread++;
    }
    genoFile.close();
    cerr<<"Geno file closed\n";
    cerr<<"Total rows filled: "<<matrixrow<<endl;
  }
}




double lratio(double ln1,double ln0){
  //cerr<<"ln1:"<<ln1<<"ln0:"<<ln0<<endl;
  if (isinf(ln1)||isnan(ln1)) ln1 = -1e50;
  if (isinf(ln0)||isnan(ln0)) ln0 = -1e50;
  double lratio=exp(ln1-ln0);
  if (isinf(lratio)){
    if (ln1>ln0) return 1e50;
    else return 1e-50;
  }
  return lratio;
}
// simple imputation, can be replaced by more sophisticated method later.

void MCMCSampler::computeExpGeno(char ** genomatrix, uint totalpersons,double * expGeno){
  for(uint i=0;i<iTotalSnps;++i){
    uint sum = 0;
    uint rows = 0;
    for(uint j=0;j<totalpersons;++j){
      if (cStudyGenoMatrix[i][j]!=UNDEFINED_GENO){
        ++rows;
        sum+=cStudyGenoMatrix[i][j];
      }
    }
    expGeno[i] = 1.*sum/rows;
  }
}

void outputAdjacency(double ** mat, uint dim){
  //cerr<<"Adjacency:\n";
  for (uint i=0;i<dim;++i){
    for (uint j=0;j<dim;++j){
      if (j) {
         //cerr<<" ";
      }
      //cerr<<mat[i][j];
    }
    //cerr<<endl;
  }
}

void MCMCSampler::updateCurrentModelString(string & currentModelString){
    ostringstream ossmodel;
    int modelcount=0;
    EffectSet::iterator it;
    for(it=currentModel->effectsInModelSet.begin();
    it!=currentModel->effectsInModelSet.end();it++){
      EffectPtr eff = *it;
      suint dim;
      const int * maineffs=eff->getMainEffects(dim);
      if (modelcount++) ossmodel<<";";
      for(suint i=0;i<dim;++i){
        if (i) {
          ossmodel<<" ";
        }
        ossmodel<<maineffs[i];
      }
    }
    currentModelString = ossmodel.str();
}

void MCMCSampler::doAnalysis(){
  // initialize current model
  init();
  uint & ms = currentModel->iModelSize;
  if (settings->use_a) updateAdjacencyMatrix(-1);
  if (drawFromPosterior(true)){
    cerr<<"Successfully initialized model\n";
  }else{
    throw "Could not fit initial model. Exiting.";
  }
  //cerr<<"Initial model has adjacency structure:\n";
  outputAdjacency(currentModel->adjacencyMatrix,ms);
  for(iIteration=1;iIteration<=settings->iterations;++iIteration){
    if (iIteration % 1000 == 0){
      cerr<<endl<<"Sample iteration: "<<iIteration<<endl;
    }else{
      cerr<<".";
    }
    // backup all parameters in case anything fails
    if (drawFromPosterior(false)){
      prevModel->copy(currentModel);
      checkDataIntegrity();
      if (!updateModelForm()){  // change the model by adding/deleting variables
      //cerr<<"Rolling back to previous model parameters"<<endl;
        currentModel->copy(prevModel);
      }else{
        //cerr<<"Successfully changed model parameters"<<endl;
      }
    }else{
     throw "Base model could not be fit\n";
    }
    checkDataIntegrity();
    // dump to output the current effect estimates
    //cerr<<"Model size is now: "<<currentModel->iModelSize<<endl;
    ostringstream ossmodel;
    int modelcount=0;
    EffectSet::iterator it;
    //cerr<<"First stage effects:\n";
    for(it=currentModel->effectsInModelSet.begin();
    it!=currentModel->effectsInModelSet.end();it++){
      EffectPtr eff = *it;
      uint effId = eff->getId();
      uint i = eff->indexInModel;
      //cerr<<"Beta/SE/Zpi effect:\t"<<effId<<"\t"<<currentModel->betas[i]<<"\t"<<sqrt(currentModel->invInfoMatrix[i*currentModel->iModelSize+i])<<"\t"<<currentModel->Zpi[i]<<endl;
      //eff->debugOut();
      suint dim;
      const int * maineffs=eff->getMainEffects(dim);
      out_beta<<iIteration<<"\t"<<effId<<"\t"<<currentModel->post_betas[i]<<"\t"<<sqrt(currentModel->post_betas_var[i])<<"\t";
      if (modelcount++) ossmodel<<";";
      for(suint i=0;i<dim;++i){
        if (i) {
          ossmodel<<" ";
          out_beta<<" ";
        }
        ossmodel<<maineffs[i];
        out_beta<<maineffs[i];
      }
      out_beta<<endl;
    }
    for (uint i=0;i<ms;++i){
      out_prior_beta<<iIteration<<"\t";
      //const double * Z = currentModel->effectsInModel[i]->getPriorInZ();
      suint dim;
      EffectPtr eff1 = currentModel->effectsInModel[i];
      const int * maineffs=eff1->getMainEffects(dim);
      for(suint j=0;j<dim;++j){
        if (j) {
          out_prior_beta<<" ";
        }
        out_prior_beta<<maineffs[j];
      }
      out_prior_beta<<"\t"<<currentModel->betas[i]<<"\t"<<sqrt(currentModel->invInfoMatrix[i*ms+i])<<"\t"<<currentModel->hm_mean[i]<<"\t"<<currentModel->car_mean[i]<<"\t"<<currentModel->prior_betas[i]<<"\t"<<sqrt(currentModel->prior_betas_var[i])<<"\t"<<currentModel->effectsInModel[i]->alpha<<endl;
    }
    string currmodel = ossmodel.str();
    //cerr<<"Model string:\t"<<currmodel<<endl;
    out_models<<iIteration<<"\t"<<currmodel<<"\t"<<ms<<"\t"<<currentModel->posteriorLikelihood<<endl;
    out_prior_var<<iIteration<<"\t"<<currentModel->tau2<<"\t"<<currentModel->sigma2<<endl;
    if (settings->use_z){
      out_prior_mean<<iIteration;
      for (uint z1=0; z1<iTotalColInZ; ++z1){
        if (currentModel->zColDefined[z1]){
          out_prior_mean<<"\t"<<currentModel->pi[z1];
          out_prior_mean<<"\t"<<currentModel->pi_se[z1];
        }else{
          out_prior_mean<<"\t\\N\t\\N";
        }
      }
      out_prior_mean<<"\n";
    }
  }
}

// this procedure should draw from the prior and then from the
// likelihood

bool MCMCSampler::drawFromPosterior(bool maximize){
  //cerr<<"Drawing from posterior\n";
  //bool & includePrior = currentModel->bIncludePrior;
  double & logPost = currentModel->posteriorLikelihood;
  uint ms = currentModel->iModelSize;
  // compute the beta hat
  double logLike;
  bool fitPi = true;
  if (settings->logistic){
    fitPi = maximize?(fitLogisticModel(false,true)||fitLogisticModel(true,true)):true;
  }else{
    fitPi = maximize?fitLinearModel():true;
  }
  if (fitPi && updatePi()){
    if (settings->marginal_prior){
      logLike = 0;
    }else{
      logLike = 2.*currentModel->logLikelihood;
    }
  }else{
    return false;
  }
  double & logPrior = currentModel->priorLikelihood;
  logPost = logPrior = 0.;
  if (!settings->marginal_prior){  // we don't care about shrinkage for marginal prior
    // now draw from sigma and tau conditional on beta hat
    updateSigmaTau(logPrior);
    //cerr<<"LogLike: "<<logLike<<" Log prior: "<<logPrior<<endl;
    //if (ms<=1) logPrior =0.;
    // generate a posterior draw for beta
    for(uint i=1;i<ms;++i){
      currentModel->post_betas_var[i] = 1./(1.0/currentModel->prior_betas_var[i] + 1.0/currentModel->invInfoMatrix[i*ms+i]);
      //double pbeta = currentModel->betas[i]<0?-currentModel->prior_betas[i]:
      //currentModel->prior_betas[i];
      double pbeta = currentModel->prior_betas[i];
      //double posterior_mean = posterior_variance * (currentModel->prior_betas[i]/
      //currentModel->prior_betas_var[i] + currentModel->betas[i]/currentModel->invInfoMatrix[i*ms+i]);
      double posterior_mean = currentModel->post_betas_var[i] * (pbeta/currentModel->prior_betas_var[i] + currentModel->betas[i]/currentModel->invInfoMatrix[i*ms+i]);
      // Sample from the posterior now.
      currentModel->post_betas[i] = math->stdNormal(currentModel->post_betas_var[i])+posterior_mean;
    }
  }
  logPost = logLike ;
  //logPost = logLike + logPrior;
  //cerr<<"Log Posterior: "<<logPost<<endl;
  return true;
}

bool MCMCSampler::fitLinearModel(){
  double & logL = currentModel->logLikelihood;
  if (settings->marginal_prior){
    //getRandomBetas(ms,currentModel->betas);
    logL = 0.;
    return true;
  }
  int rank = currentModel->iModelSize;
  int samplesize = iTotalStudyPersons;
  double aff[samplesize];
  double fitted[samplesize];
  for(int i=0;i<samplesize;++i){
    aff[i] = dPhenoVector[i];
  }
  double designMat[samplesize * rank];
  for(int j=0;j<rank;++j){
    EffectPtr & currEff = currentModel->effectsInModel[j];
    const double * effSizes = currEff->getResiduals();
    for(int i=0;i<samplesize;++i){
      designMat[i*rank+j] = effSizes[i];
    }
  }
  return math->linReg(aff,designMat,currentModel->betas,currentModel->invInfoMatrix,samplesize,rank,logL,fitted);
}

// computes the MLE of Beta

bool MCMCSampler::fitLogisticModel(bool bInitBetas,bool print){
  double & logL = currentModel->logLikelihood;
  //uint ms = currentModel->iModelSize;
  if (settings->marginal_prior){
    //getRandomBetas(ms,currentModel->betas);
    logL = 0.;
    return true;
  }
  string modelStr;
  updateCurrentModelString(modelStr);
  try{
    //cerr<<"Searching for model "<<modelStr<<" of size "<<currentModel->iModelSize<<endl;
    data->retrieveModelFit(modelStr,logL);
    for(uint i = 0;i<currentModel->iModelSize;++i){
      unsigned long int effid = currentModel->effectsInModel[i]->getId();
      //cerr<<" id "<<effid<<endl;
      data->retrieveModelParam(modelStr,effid,currentModel->betas[i],currentModel->invInfoMatrix[i*currentModel->iModelSize+i]);
    }
  }catch(uint error){
    if (error==DataManager::ERROR_CODE_DATA_NOT_FOUND){
      //cerr<<"Entering fitLogisticModel\n";
      bool bConverged = false;
      if (bInitBetas){
        //cerr<<"Resetting betas to zero.\n";
        reInit(currentModel->betas,currentModel->iModelSize);
      }
      //double curtime = clock();
      uint currentIter=0;
      double chisq;
      bool nan = false;
      chisq = 0.;
      double designMat[iTotalStudyPersons * currentModel->iModelSize];
      for(uint i=0;i<currentModel->iModelSize;++i){
        EffectPtr & currEff = currentModel->effectsInModel[i];
        const double * effSizes = currEff->getResiduals();
        for(uint j=0;j<iTotalStudyPersons;++j){
          designMat[j*currentModel->iModelSize+i] = effSizes[j];
        }
      }
      do{
        if (!math->scoreTest(dPhenoVector,designMat,currentModel->betas,currentModel->invInfoMatrix,iTotalStudyPersons,currentModel->iModelSize,currentModel->logLikelihood,chisq) || isnan(chisq)) bConverged = false;
        //cerr<<"Chisq: "<<chisq<<endl;
        bConverged = chisq<CONVERGED;
        ++currentIter;
      }while(!bConverged && !nan&& currentIter<MAX_NEWTON_RAPHSON);
      if (currentIter==MAX_NEWTON_RAPHSON||isnan(chisq)){
          //cerr<<"Current model could not converge.\n";
        //bInitBetas = true;
      }
      if (bConverged){
        data->storeModelFit(modelStr,logL);
        for(uint i = 0;i<currentModel->iModelSize;++i){
          unsigned long int effid = currentModel->effectsInModel[i]->getId();
          data->storeModelParam(modelStr,effid,currentModel->betas[i],currentModel->invInfoMatrix[i*currentModel->iModelSize+i]);
        }
      }
      return bConverged;
    }
  }
  //cerr<<"Found cached logL of "<<logL<<endl;
  return true;
}


void MCMCSampler::updateAdjacencyMatrix(int newindex){
  if (settings->marginal_prior) return;
  int ms=currentModel->iModelSize;
  double nullValue=0.;
  double rhoDiag = nullValue;
  const double * corr1;
  const double * corr2;
  for(int i=0;i<ms;++i){
    corr1 = currentModel->effectsInModel[i]->getPriorInA();
    for(int j=i;j<ms;++j){
      corr2 = currentModel->effectsInModel[j]->getPriorInA();
      if (newindex==-1 || i==newindex || j==newindex){
        double rho;
        if (i==j){
          rho = rhoDiag;
        }else{
          if (!i || !j){ // this is the first row or the first column
            rho = nullValue;
            if (j>i) currentModel->adjacencyMatrix[j][i]
             = rho;
          }else{ // this is every other row,column
            try{
              rho = data->retrieveRho(currentModel->effectsInModel[i]->getId(), currentModel->effectsInModel[j]->getId());
            }catch(uint error){
              if (error==DataManager::ERROR_CODE_DATA_NOT_FOUND){
                if (settings->dynamic_prior){
                  rho = math->corr(corr1,corr2,iTotalColInA);
                  rho = rho>.8?1:0;
                }else{
                  rho = math->corrIndicator(corr1,corr2,
                  iTotalColInA,termWeights,termWeightsTotal);
                  rho = rho>.0?1:0;
                }
                data->storeRho(currentModel->effectsInModel[i]->getId(),
                currentModel->effectsInModel[j]->getId(),rho);
              }
            }
            rho=(rho<=nullValue)?nullValue:rho;
            if (j>i) currentModel->adjacencyMatrix[j][i]
             = rho;
          }
        }
        currentModel->adjacencyMatrix[i][j] = rho;
      }
    } // inner for loop
  } //outer for loop
  outputAdjacency(currentModel->adjacencyMatrix,ms);
}

void MCMCSampler::updateBetaBar(uint i){
    uint & ms = currentModel->iModelSize;
    EffectPtr & eff = currentModel->effectsInModel[i];
    suint iTotalMainEff;
    eff->getMainEffects(iTotalMainEff);
    double * weights = currentModel->adjacencyMatrix[i];
    double rowsum=0.,betasum=0.;
    for(uint j=1;j<ms;++j){
      if (i!=j){
        if (weights[j]>0){
          betasum+=currentModel->betas[j]*weights[j];
          rowsum+=weights[j];
        }
      }
    }
    currentModel->totaldistance[i] = (rowsum==0)?1:rowsum+1;
    currentModel->beta_bars[i] = (rowsum==0)?0:betasum/rowsum;
}

void MCMCSampler::updateBetaBars(){
    uint & ms = currentModel->iModelSize;
    // if converged, compute the mean betas for neighbors
    for(uint i=0;i<ms;++i){
      updateBetaBar(i);
    }
}

void MCMCSampler::checkDataIntegrity(){
    if (currentModel->iModelSize!=currentModel->effectsInModelSet.size()){
        //cerr<<"Set has "<<currentModel->effectsInModelSet.size()<<" elements\n";
        throw "Array and set for model not in sync";
    }
    for (int i=0;i<(int)currentModel->iModelSize;++i){
      EffectPtr & eff = currentModel->effectsInModel[i];
      if (eff->indexInModel!=i) {
        //cerr<<"Saw "<<eff->indexInModel<<" expecting "<<i<<endl;
        throw "Mismatch in model index\n";
      }
      //eff->debugOut();
      EffectSet::iterator it= currentModel->effectsInModelSet.find(eff);
      if (it==currentModel->effectsInModelSet.end()){
          eff->debugOut();
          throw "This effect is supposed to be in the model set";
      }
    }
}

void MCMCSampler::getRandomBetas(uint size,double * betas){
  double sigmamat[size*size];
  double means[size];
  for(uint i=0;i<size;++i){
    for(uint j=0;j<size;++j){
      if (i==j) sigmamat[i*size+j] = 0.1;
      else sigmamat[i*size+j] = 0.;
    }
    means[i] = 0.;
  }
  math->ranMVNormal(sigmamat,means,betas,size);
}

bool MCMCSampler::updateModelForm(){
  //cerr<<"Entering updateModelForm()\n";
  uint & ms = currentModel->iModelSize;
  bool atIntercept = ms<=1?true:false;
  bool atMax = ms==settings->max_modelsize?true:false;
  bool returnCode1 = true;//,returnCode2 = true;
  // try to add
  double dDelProb,dSwapProb;
  dDelProb = 1.*(ms-1)/(settings->max_modelsize-1);
  dSwapProb = .5;
  if (currentModel->iModelSize>1 && math->RandomUniform()>dSwapProb){
    returnCode1 = swapEffect();
  }else{
    if (math->RandomUniform()>dDelProb){
      if (currentModel->iModelSize<settings->max_modelsize)
        returnCode1 = addEffect(atIntercept);
      else returnCode1 = deleteEffect(atMax);
    }else{
      if (!atIntercept)
      returnCode1 = deleteEffect(atMax);
      else returnCode1 = addEffect(atIntercept);
    }
  }
  return returnCode1;
}

bool MCMCSampler::swapEffect(){
  double oldLogPost = currentModel->posteriorLikelihood;
  EffectPtr add,del;
  bool found = false;
  int counter=10;
  uint randIndex;

  double probOld,probNew;
  while(!found && --counter>0){
    randIndex = 1+static_cast<uint>((currentModel->iModelSize-1)*math->RandomUniform());
    del = currentModel->effectsInModel[randIndex];
    suint effsize;
    del->getMainEffects(effsize);
    if (del->referenceCount==0){
      probOld =  probAddEffect(randIndex);
      found = true;
    }
      probOld = probAddEffect(randIndex);
  }
  if (!found) {
    //cerr<<"No main effects found to swap out\n";
    rollbackAdjacency(randIndex);
    return false;
  }else{
    //cerr<<counter<<": Found effect to swap out: "; del->debugOut();
  }
  int counter2=10;
  bool found2 = false;
  int indices[1];
  //double start=clock();
  while(!found2 && --counter2>0){
    indices[0] = static_cast<uint>(iTotalExposures*
    math->RandomUniform());
    if (createNewEffect(indices,1,false,add)==
    Effect::RETURNCODE_ADDABLE){
      currentModel->effectsInModel[randIndex] = add;
      probNew = probAddEffect(randIndex);
      found2 = true;
    }
  }
  if (!found2){ // roll back
    //cerr<<"No effects found to swap in\n";
    currentModel->effectsInModel[randIndex] = del;
    rollbackAdjacency(randIndex);
    return false;
  }else{
    //cerr<<"Swapped in new main effect:\n"; add->debugOut();
    add->newestMainEff = randIndex;
    add->indexInModel = randIndex;
    add->child1Id = 0;
    add->child2Id = -1; // this means it is undefined.
  }
  currentModel->effectsInModel[randIndex] = add;
  add->indexInModel = randIndex;
  currentModel->effectsInModelSet.erase(del);
  currentModel->effectsInModelSet.insert(add);
  if (!drawFromPosterior(true)) {
    rollbackAdjacency(randIndex);
    return false;
  }
  double newLogPost = currentModel->posteriorLikelihood;
  double proposalRatio = 1.0;
  proposalRatio = probNew/probOld;
  //}
  double LR,HastingsRatio;
  LR = lratio(newLogPost,oldLogPost);
  HastingsRatio = 1.0*LR*proposalRatio;
  double sa_prob = HastingsRatio>0?HastingsRatio:0;

  if (HastingsRatio>1 || sa_prob > math->RandomUniform()){
    //cerr<<"swap Committing:\n";
    decrementDescendantCounts(del);
    incrementDescendantCounts(add);
  }else{
    //cerr<<"swap Rolling back"<<endl;
    currentModel->effectsInModel[randIndex] = del;
    rollbackAdjacency(randIndex);
    currentModel->effectsInModelSet.erase(add);
    currentModel->effectsInModelSet.insert(del);
    return false;
  }
  return true;
}

bool MCMCSampler::randomAdd(EffectPtr & add,bool interaction,bool randInter){
  //cerr<<"Entering randomAdd\n";
  bool found = false;
  int counter=5;
  while(!found&&--counter>0){
    int indices[1];
    if (!interaction){
      uint randIndex = static_cast<uint>(iTotalExposures*math->RandomUniform());
//       check if it can be added as a main effect
      indices[0] = randIndex;
      // take a 50-50 prob of main effect or higher order effect
      if (createNewEffect(indices,1,false,add)==Effect::RETURNCODE_ADDABLE){
        add->newestMainEff = randIndex;
        add->child1Id = 0;
        add->child2Id = -1; // this means it is undefined.
        found = true;
      }
    }else{ // propose a higher order effect
      int counter2 = 10;
      do{
        // select an effect from the model at random
        uint modelindex = currentModel->iModelSize-1;
        if (randInter){
          modelindex = static_cast<uint>(1+(currentModel->iModelSize-1)*
          math->RandomUniform());
        }
        EffectPtr & currentEffect = currentModel->effectsInModel[modelindex];
  //       check if it can be added as a main effect
        suint dim;
        const int * mainEff = currentEffect->getMainEffects(dim);
        if (dim<settings->max_order){
          EffectSet::iterator it1;
          uint randIndex;
          int counter3 = 10;
          do{
            // select a random main effect to append.
            randIndex = static_cast<uint>(iTotalExposures*math->RandomUniform());
            indices[0] = randIndex;
            EffectPtr simpleQuery = EffectPtr(new Effect(indices,1));
            it1 = currentModel->effectsInModelSet.find(simpleQuery);

            if(it1!=currentModel->effectsInModelSet.end() &&
            !DataManager::inIntArr(mainEff,dim,randIndex)){
              // The effect is eligible since it is hierarchical
              int expandEffs[dim+1];
              for(unsigned short int k=0;k<dim;++k) expandEffs[k] = mainEff[k];
              expandEffs[dim] = randIndex;
              // Check if the effect is redundant
              if (createNewEffect(expandEffs,dim+1,false,add)==Effect::RETURNCODE_ADDABLE){
                add->newestMainEff = randIndex;
                add->child1Id = currentEffect->getId();
                add->child2Id = (*it1)->getId();
                    found = true;
              }
            }
          }while(!found && --counter3>0);
        }
        ++modelindex;
      }while(!found&&--counter2>0);
    }
  }
  return found;
}

bool MCMCSampler::randomDelete(EffectPtr & del,uint & randIndex){
  // handle the random deletion
  bool found = false;
  int counter=10;
  while(!found && --counter>0){
    randIndex = 1+static_cast<uint>((currentModel->iModelSize-1)*math->RandomUniform());
    del = currentModel->effectsInModel[randIndex];
    if (del->referenceCount==0){
      suint effsize;
      del->getMainEffects(effsize);
      found = true;
    }
  }
  if (!found) return false;
  return true;
}





double MCMCSampler::probAddEffect(uint index){
  double & prob = currentModel->effectsInModel[index]->alpha = 0.;
  if (settings->marginal_prior){
    prob = .5;
    return prob;
  }

  uint ms = currentModel->iModelSize;
  double rhoArr[ms];
  rhoArr[0] = 0.;
  double sum=0.;
  double priorvar = currentModel->sigma2+currentModel->tau2;
  if (settings->use_z && currentModel->totaldefinedcolinz>1){
    const double * Z = currentModel->effectsInModel[index]->getPriorInZ();
    uint definedcount=0;
    for(uint i=0;i<iTotalColInZ;++i){
      if (currentModel->zColDefined[i]){
        ++definedcount;
        sum+=Z[i] * currentModel->pi[i];
      }
    }
    if (ms<=definedcount) {
      sum = 0.;
    }
  }
  if (settings->use_a){
    updateAdjacencyMatrix(index);
    updateBetaBar(index);
    sum+=currentModel->beta_bars[index];
    //cerr<<"With CAR, sum is "<<sum<<endl;
  }

  // ignore sign
  sum = math->stdNormal(priorvar)+abs(sum);
  // shift everything by one std
  sum-=dBICPenalty*sqrt(priorvar);
  prob = math->cdfNormal(sum,priorvar,0.);
  return prob;
}

void MCMCSampler::rollbackAdjacency(uint index1){
    updateAdjacencyMatrix(index1);
    updateBetaBar(index1);
}

void MCMCSampler::rollbackAdjacency(uint index1,uint index2){
    updateAdjacencyMatrix(index1);
    updateAdjacencyMatrix(index2);
    updateBetaBar(index1);
    updateBetaBar(index2);
}
bool MCMCSampler::addEffect(bool force){
  // Establish the prior and likelihood for the base model
  double oldLogPost=currentModel->posteriorLikelihood;
  EffectPtr addEff1,addEff2;
  uint timeout=1000;
  double proposalRatio =  1.;
  double prob1 = .5,prob2 = .5;
  int added = 0;
  double mainEffProb = settings->max_order==1?1.0:math->beta(settings->main_effect_pref,1.);
  bool interaction = (mainEffProb>math->RandomUniform()||currentModel->iModelSize<=2)?false:true;
  if (interaction){
    if (math->RandomUniform()>.5){ // random interaction
      if (!randomAdd(addEff1,true,true)) return false;
      currentModel->effectsInModel[currentModel->iModelSize++] = addEff1;
      addEff1->indexInModel = currentModel->iModelSize-1;
      ++added;
      prob1 = probAddEffect(currentModel->iModelSize-1);
    }else{ // a new main effect followed by its interaction
      //cerr<<"Adding main + inter\n";
      if (!randomAdd(addEff1,false,false)) return false;
      currentModel->effectsInModel[currentModel->iModelSize++] = addEff1;
      addEff1->indexInModel = currentModel->iModelSize-1;
      prob1 = probAddEffect(currentModel->iModelSize-1);
      if (!randomAdd(addEff2,true,false)) return false;
      currentModel->effectsInModel[currentModel->iModelSize++] = addEff2;
      addEff2->indexInModel = currentModel->iModelSize-1;
      prob2 = probAddEffect(currentModel->iModelSize-1);
      added+=2;
    }
  }else{ //main effect
    //cerr<<"Adding main\n";
    if (!randomAdd(addEff1,false,false)) return false;
    currentModel->effectsInModel[currentModel->iModelSize++] = addEff1;
    addEff1->indexInModel = currentModel->iModelSize-1;
    ++added;
    prob1 = probAddEffect(currentModel->iModelSize-1);
  }
  if (!timeout) return false;

  proposalRatio*=(prob1*prob2)/((1.-prob1)*(1.-prob2));
  if (isinf(proposalRatio)) proposalRatio = 1e50;
  currentModel->effectsInModelSet.insert(addEff1);
  if (added==2) currentModel->effectsInModelSet.insert(addEff2);
  if (!drawFromPosterior(true)) return false;
  double newLogPost = currentModel->posteriorLikelihood;
  double HastingsRatio;
  double LR = lratio(newLogPost,oldLogPost);
  HastingsRatio = 1.0*LR*proposalRatio;
  double sa_prob = HastingsRatio>0?HastingsRatio:0;
  if (HastingsRatio>1 || sa_prob > math->RandomUniform()){
    incrementDescendantCounts(addEff1);
    if (added==2){
      incrementDescendantCounts(addEff2);
    }
  }else{
    currentModel->effectsInModelSet.erase(addEff1);
    if (added==2){
      currentModel->effectsInModelSet.erase(addEff2);
    }
    return false;
  }
  return true;
}


bool MCMCSampler::deleteEffect(bool force){
  double oldLogPost = currentModel->posteriorLikelihood;
  EffectPtr delEff;
  uint delIndex;
  double proposalRatio =  1.;
  bool found = false;
  double prob;
  while(!found){
    if (!randomDelete(delEff,delIndex)) {
      rollbackAdjacency(delIndex);
      return false;
    }
    prob = probAddEffect(delIndex);
    found = prob<math->RandomUniform()?true:false;
  }
  proposalRatio*= (1.-prob)/(prob);
  if (isinf(proposalRatio)) proposalRatio = 1e50;

  bool insertion = false;
  if (delIndex<(currentModel->iModelSize-1)){ // re-arrange for performance
    insertion = true;
    // put the last element into the removed index
    currentModel->effectsInModel[delIndex] =
    currentModel->effectsInModel[currentModel->iModelSize-1];
    currentModel->effectsInModel[delIndex]->indexInModel = delIndex;
  }
  --currentModel->iModelSize;
  // make sure the adjacency structure is correct
  //updateAdjacencyMatrix(delIndex);
  currentModel->effectsInModelSet.erase(delEff);
  if (!drawFromPosterior(true)) {
    if (insertion){
      rollbackAdjacency(delIndex,currentModel->iModelSize-1);
    }
    return false;
  }
  double newLogPost = currentModel->posteriorLikelihood;
  double LR,HastingsRatio;
  LR = lratio(newLogPost,oldLogPost);
  HastingsRatio = 1.0*LR*proposalRatio;
  double sa_prob = HastingsRatio>0?HastingsRatio:0;
  suint mainEffDim;
  delEff->getMainEffects(mainEffDim);
  if (HastingsRatio>1 || sa_prob > math->RandomUniform()){
    decrementDescendantCounts(delEff);
  }else{
    ++currentModel->iModelSize;
    if (insertion){
      currentModel->effectsInModel[currentModel->iModelSize-1] = currentModel->effectsInModel[delIndex];
      currentModel->effectsInModel[currentModel->iModelSize-1]->indexInModel = currentModel->iModelSize-1;
      currentModel->effectsInModel[delIndex] = delEff;
      rollbackAdjacency(delIndex,currentModel->iModelSize-1);
    }
    currentModel->effectsInModelSet.insert(delEff);
    return false;
  }
  return true;
}

bool MCMCSampler::updatePi(){
  
  if (!settings->use_z||settings->marginal_prior) return true;
  // remove the intercept of the first level
  uint rows = currentModel->iModelSize - 1;
  double zmatbycol[iTotalColInZ * rows];
  // populate Z matrix by col major
  //cerr<<"Zmatrix:\n";
  //bool univariate_debug = false;
  for(uint i=0;i<rows;++i){
    const double * Z = currentModel->effectsInModel[i+1]->getPriorInZ();
    //cerr<<currentModel->betas[i+1];
    for (uint j=0;j<iTotalColInZ;++j){
      // Any very small value of Z is likely to be spurious.
      //zmatbycol[j*rows+i] = abs(Z[j])<.1?0:Z[j];
      zmatbycol[j*rows+i] =  Z[j];
      //zmatbycol[j*rows+i] = (!i && j)?0.:Z[j];
      //cerr<<" "<<zmatbycol[j*rows+i];
    }
    //cerr<<endl;
  }
  // intercept of Z is always defined
  currentModel->zColDefined[0] = true;
  uint & totaldefinedcol= currentModel->totaldefinedcolinz = 1;
  // check that there is variation in all cols
  for (uint j=1;j<iTotalColInZ;++j){
    //if (currentModel->zColDefined[j]){
    //cerr<<"Checking col "<<j<<endl;
    uint k=0;
    bool indep = true;
    while(indep && k<j){
      double * col1 = zmatbycol+k*rows;
      double * col2 = zmatbycol+j*rows;
      bool diff = false;
      uint row=0;
      while(row<rows && !diff){
        diff=(col1[row]!=col2[row]);
        ++row;
      }
      if (!diff) indep = false;
      ++k;
    }
    currentModel->zColDefined[j] = indep;
    totaldefinedcol+=indep;
  }
  //cerr<<"Total defined cols "<<totaldefinedcol<<endl;
  // populate the trimmed Z matrix
  //double  pi[totaldefinedcol];
  // transpose the Z matrix
  double zmatbyrow[rows * totaldefinedcol];
  double invinfomatrix[totaldefinedcol*totaldefinedcol];
  for(uint i=0;i<rows;++i){
    int k=0;
    for (uint j=0;j<iTotalColInZ;++j){
      if (currentModel->zColDefined[j]){
        zmatbyrow[i*totaldefinedcol+k] = zmatbycol[j*rows+i];
        k++;
      }
    }
  }
  double absbeta[rows];
  for(uint i=0;i<rows;++i){
    // testing
    //absbeta[i] = fabs(currentModel->betas[i]);
    // original
    absbeta[i] = currentModel->betas[i+1];
    //absbeta[i] = !i?0.:(currentModel->betas[i]);
    //cerr<<absbeta[i]<<" "<<zmatbyrow[i*2+1]<<endl;
  }
  bool linRegOK = false;
  if (rows>totaldefinedcol){
    double pi_max[totaldefinedcol];
    //if (print)
    //cerr<<"Running linear regression total defined cols: "<<totaldefinedcol<<"\n";
    double fitted[rows];
    double logL=0;
    linRegOK =  math->linReg(absbeta, zmatbyrow,pi_max, invinfomatrix, rows, totaldefinedcol,logL,fitted);
    //cerr<<linRegOK;
    int z2=0;
    for (uint z1=0; z1<iTotalColInZ; ++z1){
      if ( !currentModel->zColDefined[z1] || !linRegOK ||  isnan(pi_max[z2]) || abs(pi_max[z2])>10.){
        currentModel->zColDefined[z1] = false;
      }else{
       //cerr<<z1<<" "<<currentModel->zColDefined[z1]<<" "<<linRegOK<<" "<<pi_max[z2]<<endl;
        currentModel->pi[z1] = math->stdNormal(invinfomatrix[z2*totaldefinedcol+z2])+pi_max[z2];
        currentModel->pi_se[z1] = sqrtf(invinfomatrix[z2*totaldefinedcol+z2]);
        ++z2;
        //currentModel->pi[z1] = 0.;
      }
    }
  }else{
    for (uint j=1;j<iTotalColInZ;++j){
      currentModel->zColDefined[j] = false;
    }
  }
  totaldefinedcol = 0;
  for (uint j=0;j<iTotalColInZ;++j){
     if (currentModel->zColDefined[j]) ++totaldefinedcol;
  }
  //cerr<<"Updated pi: "<<currentModel->totaldefinedcolinz<<endl;
  return true;
}


// OK follow these steps
// 1. draw from the conditional distributions of sigma and tau based
//  on the MLE of beta
// 2. draw betas from the conditional distributions of CAR and a univariate normal
// based on tau and sigma respectively
// update sigma2 and tau2 using the CAR
void MCMCSampler::updateSigmaTau(double & LL){
  uint ms = currentModel->iModelSize;
  //outputAdjacency(currentModel->adjacencyMatrix,ms);
  updateBetaBars();
  //double numer_sigma=0.,numer_tau=0.,denom=0.;
  // skip the intercept for computing tau sigma
  double denom=0.;
  double s_numer=0,t_numer=0;
  //double hm_mean[ms];
  //double car_mean[ms];

  for (uint i=0;i<ms;++i){
    denom+=pow(math->stdNormal(1),2);
    // compute zpi
    const double * Z=currentModel->effectsInModel[i]->getPriorInZ();
    currentModel->Zpi[i]=0;
    for(uint j=0;j<iTotalColInZ;++j){
      if (currentModel->zColDefined[j]){
        currentModel->Zpi[i]+=Z[j] * currentModel->pi[j];
        //cerr<<"Z[j]"<<Z[j]<<" pi "<<currentModel->pi[j]<<endl;
      }
    }
    currentModel->hm_mean[i] = settings->use_z?currentModel->Zpi[i]:0;
    s_numer+=pow(currentModel->hm_mean[i]-currentModel->betas[i],2);
    //cerr<<"hm_mean"<<i<<":"<<currentModel->hm_mean[i]<<endl;
    currentModel->car_mean[i] = settings->use_a?currentModel->beta_bars[i]:0;
    t_numer+=pow(currentModel->car_mean[i]-currentModel->betas[i],2);
  }
  //cerr<<"s_numer"<<s_numer<<"t_numer"<<t_numer<<"denom"<<endl;
  currentModel->sigma2 = s_numer/denom;
  currentModel->tau2 = t_numer/denom;
  //currentModel->tau2 = settings->use_a?t_numer/denom:0;
  LL=0.;
//  if (!settings->marginal_prior) //cerr<<"tag,index,mle,se,hm,car,priorbeta,priorvar\n";
  for (uint i=0;i<ms;++i){
    double priormean = currentModel->hm_mean[i]+currentModel->car_mean[i];
    currentModel->prior_betas_var[i] = currentModel->sigma2+currentModel->tau2/currentModel->totaldistance[i];
    currentModel->prior_betas[i]=math->stdNormal(currentModel->prior_betas_var[i])+priormean;
    //if (!settings->marginal_prior) //cerr<<"hyper,"<<currentModel->effectsInModel[i]->getId()<<","<<currentModel->betas[i]<<","<<currentModel->invInfoMatrix[i*ms+i]<<","<<currentModel->hm_mean[i]<<","<<currentModel->car_mean[i]<<","<<currentModel->prior_betas[i]<<","<<currentModel->prior_betas_var[i]<<endl;
  }
//  //cerr<<"updatesigmatau final: tau2: "<<currentModel->tau2<<" new sigma2: "<<currentModel->sigma2<<endl;
}

suint MCMCSampler::createNewEffect(const int * mainEff,int iTotalMainEff,
bool bTestAgainstModel,EffectPtr & newEff){

//    EffectPtr nullEff;
    EffectPtr testEff = EffectPtr(new Effect(mainEff,iTotalMainEff));
    // if this effect is in the pool of considered effects return;
    EffectSet::iterator it= currentModel->effectsInModelSet.find(testEff);
    if(it!=currentModel->effectsInModelSet.end()){
        newEff = *it;
        return Effect::RETURNCODE_IN_MODEL_ALREADY;
    }
    // if we explored this effect in the past, re-use it

    newEff = data->loadCachedEffect(testEff);
    if (newEff==NULL) {
      //cerr<<"Generating a new effect\n";
        // initialize everything in preparation for caching this effect
        double * obsEffects = new double[iTotalStudyPersons];
        // pre-compute the multiplicative effect size for the interaction
        // in the case control data
        double tempeffects[settings->max_order][iTotalStudyPersons];
        for(unsigned short int j=0;j<iTotalMainEff;++j){
          if(mainEff[j]==-1){
            for(uint i=0;i<iTotalStudyPersons;++i){
              tempeffects[j][i] = 1.0;
            }
          }else{
            if (static_cast<uint>(mainEff[j])>=iTotalSnps){
              cerr<<"Loading from environmental covariates\n";
              for(uint i=0;i<iTotalStudyPersons;++i){
                tempeffects[j][i] = 1.0*fStudyEnvMatrix[mainEff[j]-iTotalSnps][i];
                //cerr<<tempeffects[j][i];
              }
              //cerr<<endl;
            }else{
              for(uint i=0;i<iTotalStudyPersons;++i){
                double genotype = cStudyGenoMatrix[mainEff[j]][i]==
                UNDEFINED_GENO?expGeno_obs[mainEff[j]]:
                cStudyGenoMatrix[mainEff[j]][i];
                tempeffects[j][i] = 1.0*genotype;
              }
            }
            //if (!bMyGoEnabled){
              math->standardize(tempeffects[j],iTotalStudyPersons);
            //}
          }
        }

        for(uint i=0;i<iTotalStudyPersons;++i){
          double effsize = 1.0;
          for(unsigned short int j=0;j<iTotalMainEff;++j){
            effsize*=tempeffects[j][i];
          }
          obsEffects[i] = effsize;
        }

        //double * corrTotal = NULL;
        //uint goterms=0,totalbiomarkers=0;
        //uint iTotalColInA = data->getTotalColInA();
        double * corrInZ = NULL;
        double * corrInA = NULL;
        if (!settings->marginal_prior){
          corrInZ = new double[iTotalColInZ];
          corrInA = new double[iTotalColInA];
          corrInZ[0] = 1;
          double colsum;
          uint validrows;
          if(!settings->dynamic_prior){
            double temp2d[settings->max_order][iTotalColInA];
            validrows=0;
            if (mainEff[0] > -1){
              for(suint i=0;i<iTotalMainEff;++i){
                // annotations relevant only for SNPs
                //if (mainEff[i]<static_cast<int>(iTotalSnps)){
                  for(uint j=0;j<iTotalColInA;++j){
                    temp2d[validrows][j] = dStaticA[mainEff[i]][j];
                  }
                  ++validrows;
                //}
              }
            }

            colsum = 0.;
            for(uint j=0;j<iTotalColInA;++j){
              double min=validrows?1:0;
              for(suint i=0;i<validrows;++i){
                if (temp2d[i][j]<min) min = temp2d[i][j];
                colsum+=temp2d[i][j];
              }
              corrInA[j] = min;

            }
            for(uint col=0;col<iTotalColInZ-1;++col){
              double score=mainEff[0]==-1?-1:1;
  //            int rows=0;
              for(suint i=0;i<iTotalMainEff;++i){
                if ( mainEff[i] >=0 ){
                  if (score>dStaticZ[mainEff[i]][col]) score = dStaticZ[mainEff[i]][col];
                }
              }
              //if (rows) score/=rows;
              corrInZ[1+col] = score;
            }
          }else if (settings->dynamic_prior){
            //totalbiomarkers = iTotalColInA;
            // do the same for the prior observations but we don't need to save it.
            double priorEffects[iTotalPriorPersons];
            for(uint i=0;i<iTotalPriorPersons;++i){
                double effsize = 1.0;
                for(unsigned short int j=0;j<iTotalMainEff;++j){
                    if (mainEff[j]==-1){
                        effsize*=1.;
                    }else if (static_cast<uint>(mainEff[j])>=iTotalSnps){
                        effsize*=1.*fPriorEnvMatrix[mainEff[j]-iTotalSnps][i];
                    }else{
                        double genotype = cPriorGenoMatrix[mainEff[j]][i]==
                        UNDEFINED_GENO?
                        expGeno_prior[mainEff[j]]:cPriorGenoMatrix[mainEff[j]][i];
                        effsize*=1.*genotype;
                    }
                }
                priorEffects[i] = effsize;
            }
            for (uint col=0;col<iTotalColInA;++col){
              double * metCol = fBioMarkerMatrix[col];

              corrInA[col] = math->corr(priorEffects,metCol,iTotalPriorPersons);
            }
            for(uint col=0;col<iTotalColInZ-1;++col){
              corrInZ[col+1] = corrInA[selectedBiomarkersInZ[col]];
            }
          }
        }

        //cerr<<"Maineff: "<<mainEff[0]<<endl;
        unsigned long int lastID = data->getLastEffId();
        newEff = EffectPtr(new Effect(lastID,mainEff,iTotalMainEff,
        obsEffects,corrInA,corrInZ));
        data->saveCachedEffect(newEff);
        //cerr<<"Maineff saved: "<<endl;
    }
    return Effect::RETURNCODE_ADDABLE;

}


bool MCMCSampler::findEffInArray(unsigned long int iid,EffectPtr & eff){
    uint index=0;
    bool found = false;
//    EffectPtr eff;
    while(!found && index<currentModel->iModelSize){
        if (currentModel->effectsInModel[index]->getId()==iid){
            found = true;
            eff = currentModel->effectsInModel[index];
            //eff->debugOut();
        }
        ++index;
    }
    if (!found){
       //cerr<<"Not found: iid "<<iid<<endl;
       //eff->debugOut();
       //throw "Could not find effect in model array";
    }
    return found;

//    return eff;
}

bool MCMCSampler::incrementDescendantCounts(EffectPtr eff){
    //cerr<<"Incrementing for eff: ";
    //eff->debugOut();
    bool found = false;
    eff->referenceCount = 0;
    if (eff->iTotalMainEff>1){
        EffectPtr child2;
        found = findEffInArray(eff->child2Id,child2);
        ++child2->referenceCount;
    }
    EffectPtr child1;
    found = findEffInArray(eff->child1Id,child1);
    ++child1->referenceCount;
    return found;
}



bool MCMCSampler::decrementDescendantCounts(EffectPtr eff){
    //cerr<<"Decrementing for eff: ";
    //eff->debugOut();
    bool found = false;
    if (eff->iTotalMainEff>1){
        EffectPtr child2;
        found = findEffInArray(eff->child2Id,child2);
        --child2->referenceCount;
    }
    EffectPtr child1;
    found = findEffInArray(eff->child1Id,child1);
    --child1->referenceCount;
    return found;
}

void MCMCSampler::reInit(double * vec,uint dim){
//    //cerr<<"Initializing\n";
    for(uint i=0;i<dim;++i) vec[i] = 0.;
}



Effect::Effect(const int*  mainEffects,unsigned short int iTotalMainEff){
    init(0,mainEffects,iTotalMainEff,NULL,NULL,NULL);
}


struct byInt{
    bool operator()(const int & i1, const int  & i2) const{
        return (i1<i2);
    }
};

Effect::Effect(unsigned long int id,const int * mainEffects, unsigned short int iTotalMainEff, const double *  residuals,const double * prior_in_a,const double * prior_in_z){
    init(id,mainEffects,iTotalMainEff,residuals,prior_in_a,prior_in_z);
}

void Effect::init(unsigned long int id,const int * mainEffects, unsigned short int iTotalMainEff,const double *  residuals,const double * prior_in_a, const double * prior_in_z){

    this->id = id;
    // we want to make sure we have our own copy
    // before sorting
    list<int> toSort;
    for(uint i=0;i<iTotalMainEff;++i) toSort.push_back(mainEffects[i]);
    toSort.sort(byInt());
    int * sortedArr = new int[iTotalMainEff];
    this->iTotalMainEff = iTotalMainEff;
    uint count=0;
    for(list<int>::iterator it=toSort.begin();
    it!=toSort.end();it++){
        sortedArr[count++] = *it;
    }
    this->mainEff = sortedArr;
    this->residuals = residuals;
    this->prior_in_a = prior_in_a;
    this->prior_in_z = prior_in_z;
    this->referenceCount = 0;
    this->indexInModel = -1;
    this->alpha = 0.;
    child1Id=0;
    child2Id=0;
}


void Effect::debugOut(){
    cerr<<"Effect ID: "<<id<<", maineff:";
    for (uint i=0;i<iTotalMainEff;++i) cerr<<" "<<mainEff[i];
    cerr<<",modindex:"<<indexInModel<<
    ",refcount:"<<referenceCount<<endl;
}



const double * Effect::getPriorInA(){
    return prior_in_a;
}

const double * Effect::getPriorInZ(){
    return prior_in_z;
}

const double * Effect::getResiduals(){
    return residuals;
}

Effect::~Effect(){

    if (mainEff!=NULL){
//        cerr<<"Destroying maineff:\n";
//        debugOut();
        delete[] mainEff;
    }
    if (residuals!=NULL){
//      cerr<<"Destroying obseff:\n";
       delete[] residuals;
    }
    if (prior_in_a!=NULL){
//      cerr<<"Destroying corpriorall:\n";
       delete[] prior_in_a;
    }
    if (prior_in_z!=NULL){
//      cerr<<"Destroying corrinpi:\n";
       delete[] prior_in_z;
    }
}

unsigned long int Effect::getId(){
    return this->id;
}


//intlist Effect::getMainEffects(){
//    return this->mainEffects;
//}

const int * Effect::getMainEffects(unsigned short int & iTotalMainEff){
    iTotalMainEff = this->iTotalMainEff;
    return this->mainEff;
}


ModelParams::ModelParams(){
  adjacencyMatrix = new double * [MAX_MODELSIZE];
  for(uint i=0;i<MAX_MODELSIZE;++i){
    adjacencyMatrix[i] = new double [MAX_MODELSIZE];
  }
}

ModelParams::~ModelParams(){
//  if (Zmatrix!=NULL) delete[] Zmatrix;
//  if (pi!=NULL) delete[] pi;
  for(uint i=0;i<MAX_MODELSIZE;++i){
    delete[] adjacencyMatrix[i];
  }
  delete[] adjacencyMatrix;
}

//uint ModelParams::getTotalColInZ(){
//    return iTotalColInZ;
//}

//void ModelParams::setTotalColInZ(uint cols){
//    if (cols<1) throw "Must have at least 1 col in Z";
//    iTotalColInZ = cols;
//    Zmatrix = new double[iTotalColInZ * MAX_MODELSIZE];
//    pi = new double[iTotalColInZ];
//    dimTallies = new uint[Effect::MAX_DIM+1];
//    for(uint i=1;i<=Effect::MAX_DIM;++i){
//      dimTallies[i] = 0;
//    }
//    for(uint i=0;i<MAX_MODELSIZE;++i){
//      betas[i] = 0.;
//    }
//}


void ModelParams::copy(ModelParams * & old){
    //sigma2 = old->sigma2;
    //tau2 = old->tau2;
    //phi2 = old->phi2;
    iModelSize = old->iModelSize;
    // copy just the necessary components from the first level
    logLikelihood = old->logLikelihood;
    //bIncludePrior = old->bIncludePrior;
    //posteriorLikelihood = old->posteriorLikelihood;
    sigma2 = old->sigma2;
    tau2 = old->tau2;
    effectsInModelSet.clear();
    totaldefinedcolinz = old->totaldefinedcolinz;
    for(uint i=0;i<old->iModelSize;++i){
      effectsInModel[i] = old->effectsInModel[i];
      effectsInModel[i]->indexInModel = i;
      effectsInModelSet.insert(old->effectsInModel[i]);
      betas[i] = old->betas[i];
      alpha[i] = old->alpha[i];
      hm_mean[i] = old->hm_mean[i];
      car_mean[i] = old->car_mean[i];
      Zpi[i] = old->Zpi[i];
      beta_bars[i] = old->beta_bars[i];
      totaldistance[i] = old->totaldistance[i];
      prior_betas[i] = old->prior_betas[i];
      prior_betas_var[i] = old->prior_betas_var[i];
      post_betas[i] = old->post_betas[i];
      post_betas_var[i] = old->post_betas_var[i];
      zColDefined[i] = old->zColDefined[i];
      pi[i] = old->pi[i];
      pi_se[i] = old->pi_se[i];
      for(uint j=0;j<old->iModelSize;++j){
        uint index = i*old->iModelSize+j;
        adjacencyMatrix[i][j] = old->adjacencyMatrix[i][j];
        invInfoMatrix[index] = old->invInfoMatrix[index];
      }
    }
}


void convertgeno(char * g,uint len,double * vec){
  for(uint i=0;i<len;++i){
    vec[i]  = g[i];
  }
}

int main(int argc,char * argv[]){
  try{
    if (argc<2){
      cout<<"Usage: <XML settings file>\n";
      exit(1);
    }
    const char * filename=argv[1];
    ifstream ifs(filename);
    if (!ifs.is_open()){
      cout<<"Configuration file "<<filename<<" not found.\n";
      exit(1);
    }
    boost::property_tree::ptree pt;
    read_xml(filename, pt);
    ifs.close();
    pathway_settings_t * settings = new pathway_settings_t;
    settings->file_initmodel = pt.get<string>("settings.inputfiles.initmodel");
    settings->file_snplist = pt.get<string>("settings.inputfiles.snplist");
    settings->file_trait = pt.get<string>("settings.inputfiles.trait");
    settings->file_study_geno = pt.get<string>("settings.inputfiles.study_geno");
    settings->file_study_env_var = pt.get<string>("settings.inputfiles.study_env_var");
    settings->file_prior_geno = pt.get<string>("settings.inputfiles.endo_files.prior_geno");
    settings->file_prior_env_var = pt.get<string>("settings.inputfiles.endo_files.prior_env_var");
    settings->file_prior_endopheno = pt.get<string>("settings.inputfiles.endo_files.prior_endopheno");
    settings->file_selected_endo = pt.get<string>("settings.inputfiles.endo_files.selected_endo");
    settings->file_z_matrix = pt.get<string>("settings.inputfiles.annotation_files.z_matrix");
    settings->file_a_matrix = pt.get<string>("settings.inputfiles.annotation_files.a_matrix");
    settings->usedb = getbool(pt.get<string>("settings.database.enable"));
    settings->dbhost = pt.get<string>("settings.database.host");
    settings->dbuser = pt.get<string>("settings.database.user");
    settings->dbpw = pt.get<string>("settings.database.pw");
    settings->dbname = pt.get<string>("settings.database.name");
    settings->dynamic_prior = getbool(pt.get<string>("settings.sampling.use_endoprior"));
    settings->use_a = getbool(pt.get<string>("settings.sampling.use_a"));
    settings->use_z = getbool(pt.get<string>("settings.sampling.use_z"));
    settings->logistic = getbool(pt.get<string>("settings.sampling.logistic"));
    settings->marginal_prior = getbool(pt.get<string>("settings.sampling.marginal_prior"));
    settings->max_order = pt.get<int>("settings.sampling.max_order");
    settings->max_modelsize = pt.get<int>("settings.sampling.max_modelsize");
    settings->iterations = pt.get<int>("settings.sampling.iterations");
    settings->main_effect_pref = pt.get<double>("settings.sampling.main_effect_pref");
     
    MCMCSampler * pAnalyzer = new MCMCSampler(settings);
    // This begins the analysis
    double dtime = clock();
    pAnalyzer->doAnalysis();
    cerr<<"Total time for analysis: "<<(clock()-dtime)/CLOCKS_PER_SEC<<" seconds.\n";
    delete pAnalyzer;
    delete settings;
  }catch(const char* message){
      cerr<<"Exception at main:"<<endl<<message<<endl;
  }
  cerr<<"Program is complete.\n";
  exit(EXIT_SUCCESS);
}


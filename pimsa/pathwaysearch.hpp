#include<math.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<list>
#include<queue>
#include<map>
#include<set>
#include<algorithm>
#include<iostream>
#include<fstream>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_matrix.h>
#include<boost/shared_ptr.hpp>

using namespace std;

class MathUtils;
class DataManager;
class ModelParams;


typedef short unsigned int suint;

static const uint MIN_MODELSIZE = 1;
static const uint MAX_MODELSIZE = 10000;
//static const uint MIN_MODELSIZE = 1;
//static const uint MAX_MODELSIZE = 10000;

class pathway_settings_t{
public:
  bool dynamic_prior;
  bool usedb;
  bool use_a;
  bool use_z;
  bool logistic;
  bool marginal_prior;
  double main_effect_pref;
  uint iterations;  
  uint max_order;
  uint max_modelsize;  
  string dbhost;
  string dbuser;
  string dbpw;
  string dbname;
  string file_initmodel;
  string file_snplist;
  string file_trait;
  string file_study_geno;
  string file_study_env_var;
  string file_prior_geno;
  string file_prior_env_var;
  string file_prior_endopheno;
  string file_selected_endo;
  string file_z_matrix;
  string file_a_matrix;
};

//class Effect;


class Effect{
public:
    Effect(unsigned long int id,const int * mainEff,unsigned short int iTotalMainEff, const double * residuals,const double * prior_in_a, const double * prior_in_z);
    Effect(const int * mainEff,unsigned short int iTotalMainEff);
    ~Effect();
    unsigned long int getId();
    const int * getMainEffects(unsigned short int & iTotalMainEff);
    void debugOut();

    //static short unsigned MAX_DIM;
    static const unsigned short int STATE_MODEL_IN = 1;
    static const unsigned short int STATE_MODEL_OUT = 2;
    static const unsigned short int STATE_MODEL_CACHED= 3;
    static const suint RETURNCODE_ADDABLE = 0;
    static const suint RETURNCODE_NOT_ADDABLE = 1;
    static const suint RETURNCODE_IN_MODEL_ALREADY = 2;
//    static const uint GO_TERMS = 5831;
    int indexInModel;
    unsigned long int child1Id,child2Id;
    uint referenceCount;
    int weight;
    unsigned short int iTotalMainEff;
    int newestMainEff;
    //double zpi;
    double alpha;
    const double * getResiduals();
    const double * getPriorInA();
    const double * getPriorInZ();
private:
    void init(unsigned long int id,const int * mainEffects,
    unsigned short int iTotalMainEff, const double *  residuals,
    const double * prior_in_a,const double * prior_in_z);
    unsigned long int id;
    const int * mainEff;
    // dimension should be iobssamplesize
    const double * residuals;
    // dimension should be cols in Z + 1 (intercept)
    const double * prior_in_z;
    // dimension should be cols in TOTAL_METABOLITES + 1 (intercept)
    const double * prior_in_a;
};

typedef boost::shared_ptr<Effect> EffectPtr;

//struct byEffect;

struct byEffect{
    bool operator()(const EffectPtr & eff1, const EffectPtr & eff2) const{
        unsigned short int dim1,dim2;
        const int * effdata1 = eff1->getMainEffects(dim1);
        const int * effdata2 = eff2->getMainEffects(dim2);
        if (dim1!=dim2){
            return (dim1<dim2);
        }else{
            for (unsigned short int i=0;i<dim1;++i){
                const int & i1 = effdata1[i];
                const int & i2 = effdata2[i];
                // looop from most to least sig digit
                if (i1<i2){
                    return true;
                }else if (i2<i1){
                    return false;
                }else{
//                    cerr<<"Tie breaker\n";
                }

                // this digit must otherwise be equal, iterate
                // to the next digit.
//                it1++;
//                it2++;
            }
//            cerr<<"Looks like they are the same\n";
            return false;
        }
    }
};



typedef set<EffectPtr,byEffect> EffectSet;
typedef list<EffectPtr> EffectList;




typedef vector<char> charvector;
typedef list<float> inputrowtype;
typedef list<inputrowtype> inputmatrixtype;

// proposed parameters
class ModelParams{
public:
    ModelParams();
    ~ModelParams();
    void copy(ModelParams * &);
    uint totaldefinedcolinz;
    bool zColDefined[MAX_MODELSIZE];
    double pi[MAX_MODELSIZE];
    double pi_se[MAX_MODELSIZE];
    double priorLikelihood,logLikelihood,posteriorLikelihood;
    double alpha[MAX_MODELSIZE];
    double hm_mean[MAX_MODELSIZE];
    double car_mean[MAX_MODELSIZE];
    double Zpi[MAX_MODELSIZE];
    double beta_bars[MAX_MODELSIZE];
  //these are the total connections for the cluster of betas
    double totaldistance[MAX_MODELSIZE];
    double lassoOmegas[MAX_MODELSIZE];
    double invInfoMatrix[MAX_MODELSIZE*MAX_MODELSIZE];
    double invinfoMatDecomp[MAX_MODELSIZE*MAX_MODELSIZE];
    double ** adjacencyMatrix;
    // dimension of betas is dependent on size of model
    double post_betas[MAX_MODELSIZE];
    double post_betas_var[MAX_MODELSIZE];
    double betas[MAX_MODELSIZE];
    double prior_betas[MAX_MODELSIZE];
    double prior_betas_var[MAX_MODELSIZE];
    double deviance[MAX_MODELSIZE];
    double tau2,sigma2;

    uint iModelSize;
    EffectSet effectsInModelSet;
    EffectPtr effectsInModel[MAX_MODELSIZE];
private:
    uint iTotalColInZ;
};




class MCMCSampler:public Analyzer{
public:
    void init(const ptree & pt);
    MCMCSampler();
    ~MCMCSampler();
    void run();
    friend class DataManager;
    friend class MySqlRepository;
private:
    pathway_settings_t * settings;
    //MCMCSampler(GWA * pGWA,CommandArguments args);
    //void doAnalysis();
    //void doAnalysis2();
    // accessed by data managers so publicized
    ModelParams *currentModel,*prevModel,*proposed1,*proposed2;

    // DIMENSIONS

    uint iTotalStudyPersons;
    uint iTotalPriorPersons;
    uint iTotalSnps;
    uint iTotalEnvCov;
    uint iTotalColInA;
    uint iTotalColInZ;

    // DATA STRUCTURES
    
    double * dPhenoVector;
    char ** cStudyGenoMatrix;
    char ** cPriorGenoMatrix;
    // imputed genotypes
    double * expGeno_obs,*expGeno_prior;
    double ** fBioMarkerMatrix,**fStudyEnvMatrix,**fPriorEnvMatrix;
    double ** dStaticZ;
    double ** dStaticA;
    double * termWeights;
    double termWeightsTotal;
    uint * selectedBiomarkersInZ;

    static const double CONVERGED = .00001;
    static const uint MAX_NEWTON_RAPHSON = 100;
    static const double UNDEFINED_PHENO=-99;
    static const int UNDEFINED_GENO = -1;

    // AIC/BIC penalty
    double dBICPenalty;
    // handle to the data repository
    DataManager * data;
    // math utils
    MathUtils * math;
    // outfiles used for summary statistics
    ofstream out_beta;
    ofstream out_prior_mean;
    ofstream out_prior_var;
    ofstream out_models;
    ofstream out_prior_beta;

    uint iTotalExposures;
    uint iIteration;

    void init();
    bool drawFromPosterior(bool maximize);
    bool fitLinearModel();
    bool fitLogisticModel(bool resetBeta,bool print);
    bool updatePi();
    void updateAdjacencyMatrix(int newindex);
    void updateBetaBar(uint modelindex);
    void updateBetaBars();
    void updateSigmaTau(double & LL);
    bool updateModelForm();
    bool addEffect(bool force);
    bool randomAdd(EffectPtr & add,bool interaction,bool randomInteraction);
    bool swapEffect();
    bool deleteEffect(bool force);
    bool randomDelete(EffectPtr & del,uint & delIndex);
    double probAddEffect(uint indexInModel);

    // utility functions
    suint createNewEffect(const int * mainEff,int iTotalMainEff,
    bool bTestAgainstModel,EffectPtr & newEff);
    void reInit(double * vector,uint dim);
    void copyEffects(const EffectPtr* & source,EffectPtr *& dest);
    bool incrementDescendantCounts(EffectPtr eff);
    bool decrementDescendantCounts(EffectPtr eff);
    bool findEffInArray(unsigned long int iid,EffectPtr & eff);
    void getRandomBetas(uint size,double * betas);
    void computeExpGeno(char ** genomatrix, uint totalpersons,double * expGeno);
    void checkDataIntegrity();
    void rollbackAdjacency(uint index);
    void rollbackAdjacency(uint index1,uint index2);
    void updateCurrentModelString(string & modelStr);
    void readGenoFile(const char * filename, char ** genoMatrix,uint rows,uint cols);


};

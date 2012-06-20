#include<vector>
#include<set>
#include<map>

typedef vector<const double * > matrixDouble;
typedef set<EffectPtr,byEffect> EffectSet;
typedef vector<EffectPtr> EffectVector;

struct RhoContainer{
    int i1,i2;
    double rho;
};

struct byIntPair{
    bool operator()(const RhoContainer& pair1, const RhoContainer& pair2) const{
        if(pair1.i1<pair2.i1) return true;
        else if(pair1.i1>pair2.i1) return false;
        else{
            // tie breaker
            if (pair1.i2<pair2.i2) return true;
        }
        return false;
    }
};
typedef set<RhoContainer,byIntPair> RhoContainerSet;
typedef map<string,double> ModelPropertyMap;


class RamRepository:public DataManager{
public:
    RamRepository(bool isnull);
    ~RamRepository();
private:
    virtual void saveCachedEffect(const EffectPtr & eff);
    virtual EffectPtr loadCachedEffect(const EffectPtr & eff);
    virtual void listCachedEffects();
    virtual unsigned long int getLastEffId();

    virtual void storeRho(unsigned long int effectID1,unsigned long int effectID2,double rho);
    // returns -1 if rho does not exist.
    virtual double retrieveRho(unsigned long int effectID1,unsigned long int effectID2);
    virtual void retrieveModelParam(string & modelstr, unsigned long int effid, double & beta, double & beta_var);
    virtual void storeModelParam(string &  modelstr, unsigned long int effid, double & beta, double & beta_var);
    virtual void retrieveModelFit(string & modelstr,  double & logL);
    virtual void storeModelFit(string & modelstr,  double & logL);

    RhoContainerSet rhoSet;
    ModelPropertyMap modelLogLikelihoodMap;
    ModelPropertyMap modelBetaMap;
    ModelPropertyMap modelBetaVarMap;
    matrixDouble corrVX;
    EffectSet allEffectSet;
    EffectSet delEffectSet;
    EffectVector delEffectVector;
    unsigned long int iDeletedEffects;
    unsigned long int iLastEffectId;
    bool isnull;


};


#include<mysql.h>

class MCMCSampler;

class MySqlRepository:public DataManager{
public:
   MySqlRepository(MCMCSampler * mcmc, const char * hostname, const char * username, const char * password, const char * dbname);
   ~MySqlRepository();
private:
   MYSQL * mysql;
   MCMCSampler * mcmc;
   int total;
   void sortKey(unsigned long int & effID1,unsigned long int & effID2);
   virtual void saveCachedEffect(const EffectPtr & eff);
   virtual EffectPtr loadCachedEffect(const EffectPtr & eff);
   virtual void listCachedEffects();
   //virtual void addHit(unsigned long int effID);
   //virtual void addIteration();
   virtual unsigned long int getLastEffId();
   virtual void storeRho(unsigned long int effectID1,unsigned long int effectID2,double rho);
   // throw if rho does not exist.
   virtual double retrieveRho(unsigned long int effectID1,unsigned long int effectID2);
   virtual void retrieveModelParam(string & modelstr, unsigned long int effid, double & beta, double & beta_var);
   virtual void storeModelParam(string & modelstr, unsigned long int effid, double & beta, double & beta_var);
   virtual void retrieveModelFit(string & modelstr,  double & logL);
   virtual void storeModelFit(string & modelstr,  double & logL);
};

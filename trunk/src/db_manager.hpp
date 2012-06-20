class DataManager{
public:
    DataManager();
    virtual void saveCachedEffect(const EffectPtr & eff)=0;
    virtual EffectPtr loadCachedEffect(const EffectPtr & eff)=0;
    virtual unsigned long int getLastEffId()=0;
    virtual void storeRho(unsigned long int effectID1,unsigned long int effectID2,double rho)=0;
    virtual double retrieveRho(unsigned long int effectID1,unsigned long int effectID2)=0;
    virtual void retrieveModelParam(string & modelstr, unsigned long int effid, double & beta, double & beta_var)=0;
    virtual void storeModelParam(string & modelstr, unsigned long int effid, double & beta, double & beta_var)=0;
    virtual void retrieveModelFit(string & modelstr,  double & logL)=0;
    virtual void storeModelFit(string & modelstr,  double & logL)=0;
    static bool inIntArr(const int * arr,uint dim,int key);
    virtual ~DataManager();
    static const uint ERROR_CODE_DATA_NOT_FOUND  =1;
protected:
  void getStringForm(const double * arr,uint size,string & output);
  void getArrForm(const string & input,double * arr,uint size);
};


#include"../common/main.hpp"
#include"../common/analyzer.hpp"
#include"../common/io.hpp"
#include"../common/utility.hpp"
#include "pathwaysearch.hpp"
#include "db_manager.hpp"
#include "ramrepository.hpp"

RamRepository::RamRepository(bool isnull):DataManager(){
    this->isnull = isnull;
    cerr<<"Creating a RamRepository\n";
    iDeletedEffects = 0;
    iLastEffectId = -1;
}

RamRepository::~RamRepository(){
    cerr<<"Cleaning up RamRepository\n";
}

RhoContainer makeKey(unsigned long int eff1Id,unsigned long int eff2Id){
    RhoContainer rhoContainer;
    rhoContainer.rho=0.;
    if (eff1Id>eff2Id){
        rhoContainer.i1 = eff2Id;
        rhoContainer.i2 = eff1Id;
    }else{
        rhoContainer.i1 = eff1Id;
        rhoContainer.i2 = eff2Id;
    }
    return rhoContainer;
}

void RamRepository::storeRho(unsigned long int eff1Id,unsigned long int eff2Id,double rho){
    if (isnull) return ;
    RhoContainer rhoContainer = makeKey(eff1Id,eff2Id);
    rhoContainer.rho = rho;
//    cerr<<"Inserting rho "<<rho <<" container with key "<<eff1Id<<" "<<eff2Id<<"\n";
    rhoSet.insert(rhoContainer);
}

double RamRepository::retrieveRho(unsigned long int eff1Id,unsigned long int eff2Id){
    if (isnull) throw ERROR_CODE_DATA_NOT_FOUND;
    RhoContainer rhoContainer = makeKey(eff1Id,eff2Id);
    RhoContainerSet::iterator it = rhoSet.find(rhoContainer);
    if (it==rhoSet.end()){
        throw ERROR_CODE_DATA_NOT_FOUND;
//        cerr<<"Rho could not be found";
//        exit(0);
    }
    return it->rho;
}

unsigned long int RamRepository::getLastEffId(){
  //cerr<<"last effect id is "<<iLastEffectId<<endl;
  return ++iLastEffectId;
}

void RamRepository::saveCachedEffect(const EffectPtr & eff){
  allEffectSet.insert(eff);
}

//void RamRepository::addHit(unsigned long effid){
//}

void RamRepository::listCachedEffects(){
    for(EffectSet::iterator it =allEffectSet.begin();
    it!=allEffectSet.end();it++){
        EffectPtr eff = *it;
        eff->debugOut();
    }
}



EffectPtr RamRepository::loadCachedEffect(const EffectPtr & eff){
    EffectPtr nullEff;
    //if (isnull) return nullEff;
    EffectSet::iterator it = allEffectSet.find(eff);
    if(it!=allEffectSet.end()){
        //cerr<<"Found eff in cached repository\n";
        return *it;
    }
    //cerr<<"Did not find eff in cached repository\n";
    return nullEff;
}

    void RamRepository::retrieveModelParam(string & modelstr, unsigned long int effid, double & beta, double & beta_var){
      if (isnull) throw ERROR_CODE_DATA_NOT_FOUND;
      ostringstream oss;
      oss<<modelstr<<" "<<effid<<endl;
      ModelPropertyMap::iterator it1 = modelBetaMap.find(oss.str());
      ModelPropertyMap::iterator it2 = modelBetaVarMap.find(oss.str());
      if (it1!=modelBetaMap.end()){
        //cerr<<"Retrieving beta hat\n";
        beta = it1->second;
      }else{
        throw ERROR_CODE_DATA_NOT_FOUND;
      }
      if (it2!=modelBetaVarMap.end()){
        //cerr<<"Retrieving se of beta hat\n";
        beta_var = it2->second;
      }else{
        throw ERROR_CODE_DATA_NOT_FOUND;
      }
    }

    void RamRepository::storeModelParam(string & modelstr, unsigned long int effid, double & beta, double & beta_var){
      if (isnull) return;
      ostringstream oss;
      oss<<modelstr<<" "<<effid<<endl;
      modelBetaMap[oss.str()]=beta; 
      modelBetaVarMap[oss.str()]=beta_var; 
    }

    void RamRepository::retrieveModelFit(string & modelstr,  double & logL){
      if (isnull) throw ERROR_CODE_DATA_NOT_FOUND;
      ModelPropertyMap::iterator it1 = modelLogLikelihoodMap.find(modelstr);
      if (it1!=modelLogLikelihoodMap.end()){
        //cerr<<"Loading logL for "<<modelstr<<endl;
        logL = it1->second;
      }else{
        throw ERROR_CODE_DATA_NOT_FOUND;
      }
    }


    void RamRepository::storeModelFit(string & modelstr,  double & logL){
      if (isnull) return;
       modelLogLikelihoodMap[modelstr]=logL; 
      //cerr<<"Storing logL for "<<modelstr<<endl;
    }


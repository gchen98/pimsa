#include <mysql.h>
//#include "gwa.hpp"
#include "pathwaysearch.hpp"
#include "db_manager.hpp"
#include "mysqlrepository.hpp"
//#include "utility.hpp"
//#include "folate_dbclient.hpp"

//select a.rs,b.term,c.term from (snplist as a,term as b)
//left join (snp_terms as c) on (c.rs=a.rs and c.term=b.term)
//where a.rs='rs981'


MySqlRepository::MySqlRepository(MCMCSampler * mcmc,
  const char * hostname, const char * username, const char * password, const char * dbname){
  this->mcmc = mcmc;
  mysql=NULL;
  if (mysql_library_init(0, NULL, NULL)) {
    throw "Could not initialize MySQL library";
  }
  mysql = mysql_init(mysql);
  if (mysql==NULL){
    throw "MySQL cannot init";
  }
  mysql = mysql_real_connect(mysql,hostname,username,password,dbname,3306,NULL,0);
  if (mysql==NULL){
    throw "MySQL cannot connect";
  }
  const char * effstr="CREATE TEMPORARY TABLE effect (effid int unsigned, \
  main1 mediumint(8) default NULL, \
  main2 mediumint(8) default NULL, \
  main3 mediumint(8) default NULL,\
  main4 mediumint(8) default NULL,\
  main5 mediumint(8) default NULL,\
  totalmaineff tinyint(3) unsigned default NULL,\
  residuals mediumtext,\
  prior_in_a mediumtext,\
  prior_in_z mediumtext,\
  PRIMARY KEY  (effid),\
  KEY index_maineff (main1,main2,main3,main4,main5)) ENGINE=MyISAM DEFAULT CHARSET=latin1";
  const char * rhostr="  CREATE TABLE `rho` (\
  `eff1` int(10) unsigned NOT NULL default '0',\
  `eff2` int(10) unsigned NOT NULL default '0',\
  `rho` float default NULL,\
  PRIMARY KEY  (`eff1`,`eff2`)\
) ENGINE=MEMORY DEFAULT CHARSET=latin1";
  mysql_query(mysql,effstr);
  mysql_query(mysql,rhostr);

    
  cerr<<"Initialized MySQL\n";
}

MySqlRepository::~MySqlRepository(){
  mysql_query(mysql,"drop table effect");
  mysql_query(mysql,"drop table rho");
  mysql_close(mysql);
  mysql_library_end();
  cerr<<"MySQL obj destructed\n";
}


void MySqlRepository::sortKey(unsigned long int & effID1,unsigned long int & effID2){
  if (effID2<effID1){
    unsigned long int temp = effID1;
    effID1 = effID2;
    effID2 = temp;
  }
}

void MySqlRepository::storeRho(unsigned long int eff1Id,unsigned long int eff2Id,double rho){
    sortKey(eff1Id,eff2Id);
    ostringstream insert;
    insert<<"insert into rho (eff1,eff2,rho) values ("<<eff1Id<<","<<eff2Id<<","<<rho<<")";
    mysql_query(mysql,insert.str().data());
}

double MySqlRepository::retrieveRho(unsigned long int eff1Id,unsigned long int eff2Id){
    sortKey(eff1Id,eff2Id);
    ostringstream select;
    select<<"select rho from rho where eff1="<<eff1Id<<" and eff2="<<eff2Id;
    //cerr<<select.str()<<endl;
    mysql_query(mysql,select.str().data());
    MYSQL_RES * res = mysql_store_result(mysql);
    double rho=-99;
    if (res){
      MYSQL_ROW row = mysql_fetch_row(res);
      if (row){
        rho = atof(row[0]);
      }
    }
    mysql_free_result(res);
    if (rho==-99) throw ERROR_CODE_DATA_NOT_FOUND;
    return rho;

}

void MySqlRepository::saveCachedEffect(const EffectPtr & eff){
    ostringstream oss;
    int maineffs[]={0,0,0,0,0};
    suint itotalmaineff;
    const int * mainEff = eff->getMainEffects(itotalmaineff);
    for(uint i=0;i<itotalmaineff;++i){
        maineffs[i] = mainEff[i];
    }
    string str_residuals,str_prior_in_a="",str_prior_in_z;
    getStringForm(eff->getResiduals(),mcmc->iTotalStudyPersons,
    str_residuals);
    //cerr<<"totalcol in z"<<mcmc->iTotalColInZ<<endl;
    getStringForm(eff->getPriorInZ(),mcmc->iTotalColInZ,
    str_prior_in_z);
    getStringForm(eff->getPriorInA(),mcmc->iTotalColInA,
    str_prior_in_a);
    oss<<"insert into effect(effid,main1,main2,main3,main4,main5,totalmaineff,"<<
    "residuals,prior_in_a,prior_in_z) values ("
    <<eff->getId();
    for(uint i=0;i<5;++i){
        oss<<","<<maineffs[i];
    }
    oss<<","<<itotalmaineff<<",'"<<str_residuals<<"','"<<str_prior_in_a<<"','"<<str_prior_in_z<<"')";
    //cerr<<oss.str()<<endl;

    mysql_query(mysql,oss.str().data());
}


EffectPtr MySqlRepository::loadCachedEffect(const EffectPtr & eff){
    EffectPtr newEff;
    suint itotalmaineff;
    const int * mainEff = eff->getMainEffects(itotalmaineff);
    int main5[]={0,0,0,0,0};
    for(uint i=0;i<itotalmaineff;++i){
      main5[i] = mainEff[i];
    }
    ostringstream select;
    select<<"select * from effect where ";
    for(uint i=0;i<5;++i){
        if (i) select<<" and ";
        select<<"main"<<i+1<<"="<<main5[i];
//        maineffs[i] = mainEff[i];
    }
    mysql_query(mysql,select.str().data());
//   cerr<<select.str().data()<<endl;
    MYSQL_RES * res = mysql_store_result(mysql);

    if (!res){
    }else{
      MYSQL_ROW row = mysql_fetch_row(res);
      int count=0;
      if (row){
        long int effId=atoi(row[count++]);
        int maineff[mcmc->settings->max_order];
        for(uint i=0;i<5;++i){
            maineff[i] = atoi(row[count++]);
        }
        int totaleff=atoi(row[count++]);
        double * effsize = new double[mcmc->iTotalStudyPersons];
        getArrForm(string(row[count++]),effsize,mcmc->iTotalStudyPersons); //residuals

        double * corrpriorall = new double[mcmc->iTotalColInA];
        string str_prior_in_a = string(row[count++]); // prior in a
        getArrForm(str_prior_in_a,corrpriorall,mcmc->iTotalColInA);

        double * corrpriorinz = new double[mcmc->iTotalColInZ];
        getArrForm(string(row[count++]),corrpriorinz,mcmc->iTotalColInZ);

        newEff = EffectPtr(new Effect(effId,maineff,totaleff,effsize,
        corrpriorall,corrpriorinz));

      }else{
//        cerr<<"Did not find eff in db\n";
      }
    }
    mysql_free_result(res);
    return newEff;
}

unsigned long int MySqlRepository::getLastEffId(){
    mysql_query(mysql,"select count(effid) from effect");
    MYSQL_RES * res = mysql_store_result(mysql);
    MYSQL_ROW row = mysql_fetch_row(res);
    unsigned long int count;
    if (row){
        count = atoi(row[0]);
    }
    mysql_free_result(res);
    return count;
}
 
void MySqlRepository::retrieveModelParam(string &  modelstr, unsigned long int effid, double & beta, double & beta_var){
  throw ERROR_CODE_DATA_NOT_FOUND;
}
void MySqlRepository::storeModelParam(string & modelstr, unsigned long int effid, double & beta, double & beta_var){}
void MySqlRepository::retrieveModelFit(string & modelstr,  double & logL){
  throw ERROR_CODE_DATA_NOT_FOUND;
}
void MySqlRepository::storeModelFit(string & modelstr,  double & logL){}
void MySqlRepository::listCachedEffects(){}

#include"../common/main.hpp"
#include"../common/analyzer.hpp"
#include"../common/io.hpp"
#include"../common/utility.hpp"
#include "pathwaysearch.hpp"
#include "db_manager.hpp"

DataManager::DataManager(){
}

DataManager::~DataManager(){
}

bool DataManager::inIntArr(const int * arr,uint dim,int key){
    bool found = false;
    unsigned short int j=0;
    while(!found && j<dim){
        if (arr[j++]==key) found=true;
    }
    return found;
}

void DataManager::getStringForm(const double * arr,uint size,string & output){
  ostringstream oss;
  for(uint i=0;i<size;++i){
    if (i) oss<<" ";
    oss<<arr[i];
  }
  output = oss.str();
}

void DataManager::getArrForm(const string & input,double * arr,uint size){

  istringstream iss(input);

  for(uint i=0;i<size;++i){
    iss>>arr[i];
  }
}


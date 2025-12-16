#include <iostream>
#include <string>
#include <map>
#include "OptionInfo.hpp"

class TestOptions {
public:
    TestOptions();
    void setOpt(const std::string &name, const double value);
    double getOpt(const std::string &name);
    
    std::map<std::string, OptionInfo> info_map;
    
    double doubleOpt;
    int intOpt;
    bool boolOpt;
};

TestOptions::TestOptions() :
    info_map()
{
    info_map["doubleOpt"]= OptionInfo(&doubleOpt);
    info_map["intOpt"]= OptionInfo(&intOpt);
    info_map["boolOpt"]= OptionInfo(&boolOpt);
}

void TestOptions::setOpt(const std::string &name, const double value) {
    auto it = info_map.find(name);
    if(it==info_map.end())    
        throw std::runtime_error("Unknown parameter \"" + name + "\".");
    
    it->second.setOpt(value);
}

double TestOptions::getOpt(const std::string &name) {
    auto it = info_map.find(name);
    if(it==info_map.end())    
        throw std::runtime_error("Unknown parameter \"" + name + "\".");
    
    return it->second.getOpt();    
}




int main()
{
  TestOptions o;

  
  o.setOpt("doubleOpt", 2.3);
  o.setOpt("intOpt", 42);
  o.setOpt("boolOpt", 1);
  
  std::cout << "doubleOpt: " << o.doubleOpt << ", " << o.getOpt("doubleOpt") << std::endl;
  std::cout << "intOpt: " << o.intOpt << ", " << o.getOpt("intOpt") << std::endl;
  std::cout << "boolOpt: " << o.boolOpt << ", " << o.getOpt("boolOpt") << std::endl;
  
}

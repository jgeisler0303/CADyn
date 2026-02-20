#ifndef OPTIONINFO_HPP_
#define OPTIONINFO_HPP_

#include <limits>
#include <stdexcept>
#include <string>
#include <map>
#include <fstream>
#include <sstream>

class OptionException: public std::exception {
public:
    OptionException(const std::string& msg= "General option exception") :
        m_msg(msg)
    {  }
    
    virtual const char* what() const throw () {
        return m_msg.c_str();
    }
    
    const std::string m_msg;
};

enum data_type_id_t {
    type_id_unknown,
    type_id_double,
    type_id_int,
    type_id_bool
};

class OptionInfo {
public:
    OptionInfo() :
        data_type_id(type_id_unknown),
        data(nullptr)
        {}
    
    OptionInfo(double *data_, double min_val= -std::numeric_limits<double>::infinity(), double max_val= std::numeric_limits<double>::infinity()) :
        min_val(min_val),
        max_val(max_val),
        data_type_id(type_id_double),
        data((void*)data_)
        {}
    OptionInfo(int *data_, double min_val= -std::numeric_limits<double>::infinity(), double max_val= std::numeric_limits<double>::infinity()) :
        min_val(min_val),
        max_val(max_val),
        data_type_id(type_id_int),
        data((void*)data_)
        {}
    OptionInfo(bool *data_, double min_val= -std::numeric_limits<double>::infinity(), double max_val= std::numeric_limits<double>::infinity()) :
        min_val(min_val),
        max_val(max_val),
        data_type_id(type_id_bool),
        data((void*)data_)
        {}

    void setOpt(const double value) const;
    void setOpt(const void* value) const;
    double getOpt() const;
    void getOpt(void *p) const;

    double min_val, max_val;
    
private:
   data_type_id_t data_type_id;
   void *data;
};

void OptionInfo::setOpt(const double value) const {
    if(value<min_val || value>max_val)
        throw OptionException("Option value must not be less than " + std::to_string(min_val) + " and not greater than " + std::to_string(max_val));
    
    switch(data_type_id) {
        case type_id_unknown:
            throw OptionException("Unknown data_type_id in OptionInfo.setOpt");
        case type_id_double:
            ((double*)data)[0]= value;
            break;
        case type_id_int:
            ((int*)data)[0]= value;
            break;
        case type_id_bool:
            ((bool*)data)[0]= value;
            break;
                    default:
                        throw OptionException("Unimplemented typeid in OptionInfo.setOpt");
    }
}

void OptionInfo::setOpt(const void* value) const {
    switch(data_type_id) {
        case type_id_unknown:
            throw OptionException("Unknown data_type_id in OptionInfo.setOpt");
        case type_id_double:
            if(((double*)value)[0]<min_val || ((double*)value)[0]>max_val)
                throw OptionException("Option value must not be less than " + std::to_string(min_val) + " and not greater than " + std::to_string(max_val));
            ((double*)data)[0]= ((double*)value)[0];
            break;
        case type_id_int:
            if(((int*)value)[0]<min_val || ((int*)value)[0]>max_val)
                throw OptionException("Option value must not be less than " + std::to_string(min_val) + " and not greater than " + std::to_string(max_val));
            ((int*)data)[0]= ((int*)value)[0];
            break;
        case type_id_bool:
            ((bool*)data)[0]= ((bool*)value)[0];
            break;
        default:
            throw OptionException("Unimplemented typeid in OptionInfo.setOpt");
    }
}


double OptionInfo::getOpt() const {
    switch(data_type_id) {
        case type_id_unknown:
            throw OptionException("Unknown data_type_id in OptionInfo.getOpt");
        case type_id_double:
            return ((double*)data)[0];
        case type_id_int:
            return ((int*)data)[0];
        case type_id_bool:
            return ((bool*)data)[0];
        default:
            throw OptionException("Unimplemented typeid in OptionInfo.getOpt");
    }
}


void OptionInfo::getOpt(void *value) const {
    switch(data_type_id) {
        case type_id_unknown:
            throw OptionException("Unknown data_type_id in OptionInfo.getOpt");
        case type_id_double:
            ((double*)value)[0]= ((double*)data)[0];
            break;
        case type_id_int:
            ((int*)value)[0]= ((int*)data)[0];
            break;
        case type_id_bool:
            ((bool*)value)[0]= ((bool*)data)[0];
            break;
        default:
            throw OptionException("Unimplemented typeid in OptionInfo.getOpt");
    }
}

class OptionsAccessor {
public:
    OptionsAccessor() : options_map() {}
    
    void setOpt(const std::string &name, const double value) {
        auto it = options_map.find(name);
        if(it==options_map.end())    
            throw OptionException("Unknown option \"" + name + "\".");
        try {
            it->second.setOpt(value);
        } catch(const OptionException& e) {
            throw OptionException("Option \"" + name + "\" must not be less than " + std::to_string(it->second.min_val) + " and not greater than " + std::to_string(it->second.max_val) + ".");
        }
    }
    
    double getOpt(const std::string &name) {
        auto it = options_map.find(name);
        if(it==options_map.end())    
            throw OptionException("Unknown option \"" + name + "\".");
        
        return it->second.getOpt();    
    }

    void setOptionsFromFile(const std::string &fileName)  {
        std::ifstream infile(fileName);
        double value;
        
        if(infile.is_open()) {
            int i= 0;
            for(std::string line; std::getline(infile, line);) {
                ++i;
                std::istringstream iss(line);
                
                std::string optionName;
                iss >> optionName;
                
                if(optionName.empty() || (optionName[0]=='*' || optionName[0]=='#')) // comment or empty line
                    continue;
                
                iss >> value;
                if(iss.fail()) {
                    fprintf(stderr, "Could not read value for option \"%s\" in line %d.\n", optionName.c_str(), i);
                }                
                
                setOpt(optionName, value);
            }            
        } else
            throw OptionException("Could not open options file \"" + fileName + "\". Valid options are AbsTol, RelTol, StepTol, HMin, JacRecalc, MaxSteps");
    }
    
    std::map<std::string, OptionInfo> options_map;
};

#endif /* OPTIONINFO_HPP_ */

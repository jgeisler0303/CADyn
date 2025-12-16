#ifndef OPTIONINFO_HPP_
#define OPTIONINFO_HPP_

#include <limits>
#include <stdexcept>

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
        data_type_id(type_id_double),
        data((void*)data_),
        min_val(min_val),
        max_val(max_val)
        {}
    OptionInfo(int *data_, double min_val= -std::numeric_limits<double>::infinity(), double max_val= std::numeric_limits<double>::infinity()) :
        data_type_id(type_id_int),
        data((void*)data_),
        min_val(min_val),
        max_val(max_val)
        {}
    OptionInfo(bool *data_, double min_val= -std::numeric_limits<double>::infinity(), double max_val= std::numeric_limits<double>::infinity()) :
        data_type_id(type_id_bool),
        data((void*)data_),
        min_val(min_val),
        max_val(max_val)
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
        throw std::range_error("Option value must not be less than " + std::to_string(min_val) + " and not greater than " + std::to_string(max_val));
    
    switch(data_type_id) {
        case type_id_unknown:
            throw std::domain_error("unknown data_type_id in OptionInfo.setOpt");
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
            throw std::domain_error("unimplemented typeid in OptionInfo.setOpt");
    }
}

void OptionInfo::setOpt(const void* value) const {
    switch(data_type_id) {
        case type_id_unknown:
            throw std::domain_error("unknown data_type_id in OptionInfo.setOpt");
        case type_id_double:
            if(((double*)value)[0]<min_val || ((double*)value)[0]>max_val)
                throw std::range_error("Option value must not be less than " + std::to_string(min_val) + " and not greater than " + std::to_string(max_val));
            ((double*)data)[0]= ((double*)value)[0];
            break;
        case type_id_int:
            if(((int*)value)[0]<min_val || ((int*)value)[0]>max_val)
                throw std::range_error("Option value must not be less than " + std::to_string(min_val) + " and not greater than " + std::to_string(max_val));
            ((int*)data)[0]= ((int*)value)[0];
            break;
        case type_id_bool:
            ((bool*)data)[0]= ((bool*)value)[0];
            break;
        default:
            throw std::domain_error("unimplemented typeid in OptionInfo.setOpt");
    }
}


double OptionInfo::getOpt() const {
    switch(data_type_id) {
        case type_id_unknown:
            throw std::domain_error("unknown data_type_id in OptionInfo.getOpt");
        case type_id_double:
            return ((double*)data)[0];
        case type_id_int:
            return ((int*)data)[0];
        case type_id_bool:
            return ((bool*)data)[0];
        default:
            throw std::domain_error("unimplemented typeid in OptionInfo.getOpt");
    }
}


void OptionInfo::getOpt(void *value) const {
    switch(data_type_id) {
        case type_id_unknown:
            throw std::domain_error("unknown data_type_id in OptionInfo.getOpt");
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
            throw std::domain_error("unimplemented typeid in OptionInfo.setOpt");
    }
}

#endif /* OPTIONINFO_HPP_ */

#ifndef PARAMETER_INFO_HPP_
#define PARAMETER_INFO_HPP_


class ParameterInfo {
public:
    ParameterInfo() :
        isSet(false),
        nrows(0),
        ncols(0),
        offset(0),
        data(nullptr)
        {}
        
    ParameterInfo(real_type *data_, int nrows_, int ncols_, int offset_) :
        isSet(false),
        nrows(nrows_),
        ncols(ncols_),
        offset(offset_),
        data(data_)
        {}

    bool setParam(const real_type *value);
    real_type getParam() const;
    void getParam(real_type *p) const;
    bool isSet;
    int nrows, ncols;
    int offset;
    
private:
    real_type *data;
};

bool ParameterInfo::setParam(const real_type *value) {
    bool wasSet= isSet;
    
    memcpy(data, value, sizeof(real_type)*nrows*ncols);
    isSet= true;
    return wasSet;
}

real_type ParameterInfo::getParam() const {
    return data[0];
}

void ParameterInfo::getParam(real_type *value) const {
    memcpy(value, data, sizeof(real_type)*nrows*ncols);
}

#endif /* PARAMETER_INFO_HPP_ */

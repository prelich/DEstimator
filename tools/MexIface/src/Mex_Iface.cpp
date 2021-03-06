/** @file Mex_Iface.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 04-01-2014
 * @brief The class definition for Mex_Iface.
 *
 * @copyright Mark J. Olah and The Regents of the University of New Mexico (2014).
 *            This code is free for non-commercial use and modification, provided
 *            this copyright notice remains unmodified and attached to the code
 */

#include "Mex_Iface.h"
#include <algorithm>
#include <sstream>
#include <cstdint>

/**
 * @brief Checks that a matlab mxArray object has the correct 2D dimensions
 * @param m A pointer to the mxArray to check
 * @param rows the expected number of rows should be >0
 * @param cols the expected number of colss should be >0
 *
 * Throws an exception if the number of rows or cols do not match
 */
void Mex_Iface::checkDim(const mxArray *m, int rows, int cols) const
{
    std::ostringstream msg;
    int nrows=mxGetDimensions(m)[0];
    int ncols=mxGetDimensions(m)[1];
    if ((rows>0) && (nrows!=rows)) {
        msg<<"Expected M="<<rows<<" rows, Got: "<<nrows;
        error("InputShape",msg.str());
    }
    if ((cols>0) && (ncols!=cols)) {
        msg<<"Expected "<<cols<<" cols, Got: "<<ncols;
        error("InputShape",msg.str());
    }
}

/**
 * @brief Checks that two matlab mxArray objects have the same sized last dimension.
 * @param m1 A pointer to the first mxArray to check
 * @param m2 A pointer to the second mxArray to check
 *
 * Throws an exception if the last dimensions do not match.
 */
void Mex_Iface::checkSameLastDim(const mxArray *m1, const mxArray *m2) const
{
    std::ostringstream msg;
    int nd1=mxGetNumberOfDimensions(m1);
    int nd2=mxGetNumberOfDimensions(m2);
    int last1=mxGetDimensions(m1)[nd1-1];
    int last2=mxGetDimensions(m2)[nd2-1];
    if (last1 != last2){
        msg<<"Got last dim1:"<<last1<<" last dim2:"<<last2;
        error("InputShape",msg.str());
    }
}

/**
 * @brief Checks that the Iface mex function was called with a minimum
 *        number of input and output arguments.
 * @param min_nlhs The minimum number of Left-hand-side (output) arguments to require.
 * @param min_nrhs The minimum number of Right-hand-side (input) arguments to require.
 *
 * If min_nlhs or min_nrhs arguments are negative the respective side is not checked.
 *
 * Throws an exception if there are not enough lhs or rhs arguments
 */
void Mex_Iface::checkMinNumArgs(int min_nlhs, int min_nrhs) const
{
    std::ostringstream msg;
    if ((min_nlhs>=0) && (nlhs< (uint32_t)min_nlhs)) {
        msg<<"Got:"<<nlhs<<" LHS Expected: >="<<min_nlhs;
        error("NumOutputArgs",msg.str());
    }
    if ((min_nrhs>=0) && (nrhs< (uint32_t)min_nrhs)) {
        msg<<"Got:"<<nrhs<<" RHS Expected: >="<<min_nrhs;
        error("NumInputArgs",msg.str());
    }
}


/**
 * @brief Checks that the Iface mex function was called with a maximum
 *        number of input and output arguments.
 * @param max_nlhs The maximum number of Left-hand-side (output) arguments to allow.
 * @param max_nrhs The maximum number of Right-hand-side (input) arguments to allow.
 *
 * If max_nlhs or max_nrhs arguments are negative the respective side is not checked.
 *
 * Throws an exception if there are too many lhs or rhs arguments
 */
void Mex_Iface::checkMaxNumArgs(int max_nlhs,int max_nrhs) const
{
    std::ostringstream msg;
    if ((max_nlhs>=0) && (nlhs> (uint32_t)max_nlhs)) {
        msg<<"Got:"<<nlhs<<" LHS Expected: <="<<max_nlhs;
        error("NumOutputArgs",msg.str());
    }
    if ((max_nrhs>=0) && (nrhs> (uint32_t)max_nrhs)) {
        msg<<"Got:"<<nrhs<<" RHS Expected: <="<<max_nrhs;
        error("NumInputArgs",msg.str());
    }
}


/**
 * @brief Checks that the Iface mex function was called with exactly the expected
 *        number of input and output arguments.
 * @param expected_nlhs The expected number of Left-hand-side (output) arguments to require.
 * @param expected_nrhs The expected number of Right-hand-side (input) arguments to require.
 *
 * If expected_nlhs or expected_nrhs arguments are negative the respective side is not checked.
 *
 * Throws an exception if there are an incorrect number of lhs or rhs arguments.
 */
void Mex_Iface::checkNumArgs(int expected_nlhs, int expected_nrhs) const
{
    std::ostringstream msg;
    if ((expected_nlhs>=0) && (nlhs!= (uint32_t) expected_nlhs)) {
        msg<<"Got:"<<nlhs<<" LHS Expected:"<<expected_nlhs;
        error("NumOutputArgs",msg.str());
    }
    if ((expected_nrhs>=0) && (nrhs!= (uint32_t)expected_nrhs)) {
        msg<<"Got:"<<nrhs<<" RHS Expected:"<<expected_nrhs;
        error("NumInputArgs",msg.str());
    }
}

/**
 * @brief Calls a named member function on the instance of the wrapped class.
 *
 * @param name The name of the method to call, as given to the mexFunction call.
 *
 * Throws an error if the name is not in the methodmap std::map data structure.
 */
void Mex_Iface::callMethod(std::string name)
{
    auto it=methodmap.find(name);
    if (it==methodmap.end()){ 
        mexPrintf("[IFACE] Got unknown method: %s\n", name.c_str());
        exploreMexArgs(nrhs, rhs);
        component_error("callMethod","UnknownMethod","Got bad method name.");
    } else {
        try {
            it->second();
        } catch (std::exception &e) {
            component_error(name.c_str(),"C++Exception",e.what());
        }
    }
}

/**
 * @brief Calls a named static member function of the wrapped class.
 *
 * @param name The name of the static method to call, as given to the mexFunction call.
 *
 * Throws an error if the name is not in the staticmethodmap std::map data structure.
 */
void Mex_Iface::callStaticMethod(std::string name)
{
    auto it=staticmethodmap.find(name);
    if (it==staticmethodmap.end()){ 
        mexPrintf("[IFACE] Got unknown static method: %s\n", name.c_str());
        exploreMexArgs(nrhs, rhs);
        component_error("callMethod","UnknownStaticMethod","Got bad method name.");
    } else {
        it->second();
    }
}

/**
 * @brief Reads a mxArray as a scalar C++ int32_t type.
 *
 * @param mxdata The pointer to the mxArray to interpret.
 *
 * The mxArray must be a signed 32 or 64 bit integer type.
 *
 * Throws an error if the conversion cannot be made.
 */
int32_t Mex_Iface::getInt(const mxArray *mxdata)
{
    if(mxdata==nullptr) mxdata=rhs[rhs_idx++];
    switch (mxGetClassID(mxdata)) {
        case mxINT32_CLASS: return *static_cast<uint32_t *>(mxGetData(mxdata));
        case mxINT64_CLASS: return *static_cast<uint64_t *>(mxGetData(mxdata));
        default: break;
    }
    std::ostringstream msg;
    msg<<"Expected class int32 or int64.  Got:"<<get_mx_class_name(mxdata);
    error("InputType",msg.str());
    return 0; //never get here
}

/**
 * @brief Reads a mxArray as a scalar C++ uint32_t type.
 *
 * @param mxdata The pointer to the mxArray to interpret.
 *
 * The mxArray must be an unsigned 32 or 64 bit integer type.
 *
 * Throws an error if the conversion cannot be made.
 */
uint32_t Mex_Iface::getUnsigned(const mxArray *mxdata)
{
    if(mxdata==nullptr) mxdata=rhs[rhs_idx++];
    switch (mxGetClassID(mxdata)) {
        case mxUINT32_CLASS: return *static_cast<uint32_t *>(mxGetData(mxdata));
        case mxUINT64_CLASS: return *static_cast<uint64_t *>(mxGetData(mxdata));
        default: break;
    }
    std::ostringstream msg;
    msg<<"Expected class uint32 or uint64.  Got:"<<get_mx_class_name(mxdata);
    error("InputType",msg.str());
    return 0; //never get here
}

/**
 * @brief Reads a mxArray as a scalar C++ float type.
 *
 * @param mxdata The pointer to the mxArray to interpret.
 *
 * The mxArray must be a floating point type.
 *
 * Throws an error if the conversion cannot be made.
 */
float Mex_Iface::getFloat(const mxArray *mxdata) 
{
    if(mxdata==nullptr) mxdata=rhs[rhs_idx++];
    switch (mxGetClassID(mxdata)) {
        case mxSINGLE_CLASS: return *static_cast<float *>(mxGetData(mxdata));
        case mxDOUBLE_CLASS: return *static_cast<double *>(mxGetData(mxdata));
        default: break;
    }
    std::ostringstream msg;
    msg<<"Expected floating point.  Got:"<<get_mx_class_name(mxdata);
    error("InputType",msg.str());
    return 0; //never get here
}


/**
 * @brief Reads a mxArray as a scalar C++ double type.
 *
 * @param mxdata The pointer to the mxArray to interpret.
 *
 * The mxArray must be a floating point type.
 *
 * Throws an error if the conversion cannot be made.
 */
double Mex_Iface::getDouble(const mxArray *mxdata)
{
    if(mxdata==nullptr) mxdata=rhs[rhs_idx++];
    switch (mxGetClassID(mxdata)) {
        case mxSINGLE_CLASS: return *static_cast<float *>(mxGetData(mxdata));
        case mxDOUBLE_CLASS: return *static_cast<double *>(mxGetData(mxdata));
        default: break;
    }
    std::ostringstream msg;
    msg<<"Expected floating point.  Got:"<<get_mx_class_name(mxdata);
    error("InputType",msg.str());
    return 0; //never get here
}


/**
 * @brief Reads a mxArray as a string.
 *
 * @param mxdata The pointer to the mxArray to interpret.
 *
 * Throws an error if the conversion cannot be made.
 */
std::string Mex_Iface::getString(const mxArray *mxdata)
{
    if(mxdata==nullptr) mxdata=rhs[rhs_idx++];
    char str[MAX_STR_LEN];
    if (mxGetString(mxdata, str, MAX_STR_LEN)) {
        std::ostringstream msg;
        msg<<"Expected string.  Got:"<<get_mx_class_name(mxdata);
        error("InputType",msg.str());
    }
    return std::string(str);
}


/**
 * @brief Reports an error condition to Matlab using the mexErrMsgIdAndTxt function
 *
 * @param condition A string describing the error condition encountered.
 * @param message An informative message to accompany the error.
 */
void Mex_Iface::error(std::string condition, std::string message) const
{
    std::ostringstream msgid;
    msgid<<mex_name<<":"<<condition;
    mexErrMsgIdAndTxt(msgid.str().c_str(), message.c_str());
}

/**
 * @brief Reports an error condition in a specified component to Matlab using the mexErrMsgIdAndTxt function
 *
 * @param condition A string describing the component in which the error was encountered.
 * @param condition A string describing the error condition encountered.
 * @param message An informative message to accompany the error.
 */
void Mex_Iface::component_error(std::string component, std::string condition, std::string message) const
{
    std::ostringstream msgid;
    msgid<<mex_name<<":"<<component<<":"<<condition;
    mexErrMsgIdAndTxt(msgid.str().c_str(), message.c_str());
}

/**
 * @brief A helper function for processing strings into matlab structures.
 *
 * @param p is the parameter name in the StatsT map.
 */
bool isSubStruct(Mex_Iface::StatsT::value_type &p)
{
    return p.first.find_first_of('.')>0;
}

/**
 * @brief A helper function for processing strings into matlab structures.
 *
 * @param stats A StatsT type maping from strings to scalar doubles.
 */
void Mex_Iface::outputStatsToStruct(const Mex_Iface::StatsT &stats)
{
    int nfields=stats.size();
    const char **fnames=new const char*[nfields];
    int i=0;
    for(auto &stat: stats) 
        fnames[i++]=stat.first.c_str();
    
    mxArray *st=mxCreateStructMatrix(1,1,nfields,fnames);
    delete[] fnames;
    for(auto &stat: stats)
        mxSetField(st, 0, stat.first.c_str(), makeDouble(stat.second));
    lhs[lhs_idx++]=st;
}

/**
 * @brief A helper function for processing strings into matlab structures.
 *
 * @param stats A StatsT type maping from strings to scalar doubles.
 */
Mex_Iface::StatsT 
Mex_Iface::getDoubleStruct(const mxArray *mxdata)
{
    if(mxdata==nullptr) mxdata=rhs[rhs_idx++];
    if(!mxIsStruct(mxdata)) {
        std::ostringstream msg;
        msg<<"Expected Struct.  Got:"<<get_mx_class_name(mxdata);
        error("InputType",msg.str());
    }
    StatsT stats;
    for(int i=0; i<mxGetNumberOfFields(mxdata); i++)
        stats[mxGetFieldNameByNumber(mxdata,i)] = getDouble(mxGetFieldByNumber(mxdata,0,i));
    return stats;
}


/**
 * @brief The mexFunction that will be exposed as the entry point for the .mex file
 *
 * @param[in] _nlhs The number of left-hand-side (input) arguments passed from the Matlab side of the Iface.
 * @param[in] _lhs The input arguments passed from the Matlab side of the Iface.
 * @param[in] _nrhs The number of right-hand-side (output) arguments requested from the Matlab side of the Iface.
 * @param[in,out] _lhs The output arguments requested from the Matlab side of the Iface to be filled in.
 *
 * This command is the main entry point for the .mex file, and allows the mexFunction to act like a class interface.
 * Special \@new, \@delete, \@static strings allow objects to be created and destroyed and static functions to be called
 * otherwise the command is interpreted as a member function to be called on the given object handle which is expected
 * to be the second argument.
 *
 */
void Mex_Iface::mexFunction(uint32_t _nlhs, mxArray *_lhs[], uint32_t _nrhs, const mxArray *_rhs[])
{
    setArguments(_nlhs,_lhs,_nrhs,_rhs);
    checkMinNumArgs(0,1);
    auto command=getString(rhs[0]);
    popRhs();//remove command from RHS
    if (command=="@new") { objConstruct(); }
    else if (command=="@delete") { objDestroy(); }
    else if (command=="@static") {
        checkMinNumArgs(0,1);
        auto command=getString(rhs[0]);
        popRhs();//remove handle from RHS
        callStaticMethod(command);
    } else {
        checkMinNumArgs(0,1);
        getObjectFromHandle(rhs[0]);
        popRhs();//remove handle from RHS
        callMethod(command);
    }
}


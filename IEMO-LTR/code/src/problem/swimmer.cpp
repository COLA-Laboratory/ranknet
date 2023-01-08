//
// Created by gylai on 22-9-28.
// Contact: Lai Guiyu <guiyulai.chn@gmail.com>
// COLA-Lab@UESTC
//

#include "problem/swimmer.h"

#include <cmath>

#include "core/global.h"

namespace emoc {

    SWIMMER::SWIMMER(int dec_num, int obj_num) :Problem(dec_num, obj_num)
    {
        std::cout<<"<info> SWIMMER::SWIMMER"<<std::endl;
        for (int i = 0; i < dec_num; ++i)
        {
            lower_bound_[i] = 0.0;
            upper_bound_[i] = 1.0;
        }
    }

    SWIMMER::~SWIMMER()
    {

    }

    void SWIMMER::CalObj(Individual *ind)
    {
        std::cout<<"<info> [SWIMMER::CalObj]"<<std::endl;
        PyRun_SimpleString("import sys");
        PyRun_SimpleString("sys.path.append('./')");

        wchar_t * s2[] = { L" " }; // 宽字符，长度为2字节
        PySys_SetArgv(1, s2);   // TODO: 加入argv参数 否则出错.PyAPI_FUNC(void) PySys_SetArgv(int, wchar_t **);
        //导入python文件
        PyObject *pModule = PyImport_ImportModule("MuJoCo");
        if (!pModule)
        {
            std::cerr << "[ERROR] Can not open python file MuJoCo.py!"<<std::endl;
            return;
        }
        PyObject *pSwimmer = NULL;
        if (obj_num_ == 2) {
            pSwimmer = PyObject_GetAttrString(pModule, "Swimmer");
        } else if (obj_num_ == 3) {
            pSwimmer = PyObject_GetAttrString(pModule, "SwimmerV3");
        } else {
            exit(-1);
            std::cerr<<"<error> [SWIMMER::CalObj] wrong obj number!"<<std::endl;
        }
        PyObject *x = PyList_New(dec_num_);
        for (int i = 0; i < dec_num_; ++i)
        {
            PyList_SetItem(x, i, PyFloat_FromDouble(ind->dec_[i]));
        }

        PyObject *pRet = PyObject_CallFunctionObjArgs(pSwimmer, x, NULL);
        if (pRet != NULL)  // 验证是否调用成功
        {
            int listSize = PyList_Size(pRet);
            PyObject *float_py = NULL;
            for (int i = 0; i < listSize; ++i)
            {
                float_py = PyList_GetItem(pRet, i);
                ind->obj_[i] = PyFloat_AsDouble(float_py);
            }
        }
        else
        {
            std::cerr<<"[ERROR] Function excuted unsuccessfully."<<std::endl;
            return;
        }

        /**
         * 当函数的返回值是New reference时，
         * 需要对PyObject * 变量使用Py_DECREF(），
         * 返回值是 Borrowed reference时，无需使用Py_DECREF(）
         */
        Py_DECREF(x);
        Py_DECREF(pSwimmer);
        Py_DECREF(pRet);
        Py_DECREF(pModule);

        return;
    }

}